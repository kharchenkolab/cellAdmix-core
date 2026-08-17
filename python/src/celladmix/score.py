"""Score result objects."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ._score_utils import rules_from_annotation, score_annotation


def _rule_support(rules: pd.DataFrame, pair_sets: dict) -> list:
    """Fraction of ensemble members keeping a rule for each pair."""
    n = max(len(pair_sets), 1)
    keys = zip(rules["source_cell_type"].astype(str), rules["target_cell_type"].astype(str))
    return [sum(key in kept for kept in pair_sets.values()) / n for key in keys]


class CellAdmixScore:
    """Membrane or bridge scoring result."""

    def __init__(self, fit, method: str, result: dict, params: dict | None = None):
        self.fit = fit
        self.method = method
        self.result = result
        # Scoring parameters are kept so ensemble members can be scored with
        # the exact same settings.
        self.params = dict(params or {})
        # Per-member kept-rule cache for the ensemble correction.
        self._member_rules_cache: dict = {}
        # Native bindings return row-wise dicts; convert once for plotting and
        # rule helpers while keeping the raw result available.
        self.pairs = pd.DataFrame(result.get("scores", []))
        self.summary = pd.DataFrame(result.get("summary", []))

    def __repr__(self) -> str:
        return f"CellAdmixScore(method={self.method!r}, n_summary={len(self.summary)})"

    def annotation(self, *, p_thresh: float = 0.1, adjust_p: bool = False) -> dict:
        """Return source and target calls by factor."""
        return score_annotation(self.summary, p_thresh=p_thresh, adjust_p=adjust_p)

    def rules(
        self,
        *,
        p_thresh: float = 0.1,
        adjust_p: bool = False,
        targets=None,
        max_rules: int | None = None,
        native_check: bool = True,
        native_median_thresh: float = 0.1,
        native_expr_thresh: float = 0.05,
        native_outlier_min_frac: float = 0.1,
        neighbor_k: int = 15,
    ) -> pd.DataFrame:
        """Return factor/target correction rules using the strongest source per factor.

        With ``native_check`` (the default), each rule is tested against
        target cells that have no source-type cells among their nearest
        neighbors; rules whose factor persists in those source-distant cells
        are flagged ``keep=False`` (with the reason in ``native_check``) and
        skipped by :meth:`correct`.
        """
        if self.summary.empty:
            return pd.DataFrame(
                columns=["factor", "source_cell_type", "target_cell_type", "p_value", "neg_log10_p", "rule_id"]
            )
        # Rules are derived from summary-level calls, not individual pair rows.
        rules = rules_from_annotation(
            self.annotation(p_thresh=p_thresh, adjust_p=adjust_p),
            target_cell_types=targets,
        )
        rules = rules.sort_values(["p_value", "factor", "target_cell_type"])
        if max_rules is not None:
            rules = rules.head(int(max_rules))
        if native_check and not rules.empty:
            from ._score_utils import apply_native_check

            rules = apply_native_check(
                rules,
                self.fit,
                median_thresh=native_median_thresh,
                expr_thresh=native_expr_thresh,
                outlier_min_frac=native_outlier_min_frac,
                neighbor_k=neighbor_k,
            )
        return rules

    def correct(
        self,
        rules: pd.DataFrame | None = None,
        *,
        p_thresh: float = 0.1,
        name: str | None = None,
        ensemble: int | None = None,
        vote: float = 0.3,
    ):
        """Apply score-derived correction rules to the parent fit.

        By default the correction is a molecule-vote ensemble: every NMF
        restart retained by the fit is scored and vetted independently, and
        a molecule is removed when at least ``vote`` of the members remove
        it. Pass ``ensemble=1`` for a single-fit correction using the
        selected restart only. Rules flagged ``keep=False`` by the
        native-factor check are skipped; when ``rules`` is passed
        explicitly, its source→target pairs restrict every member.
        """
        import math

        from . import _core

        user_rules = rules is not None
        rules = self.rules(p_thresh=p_thresh) if rules is None else rules
        if "keep" in getattr(rules, "columns", ()):
            rules = rules[rules["keep"].fillna(True).astype(bool)]
        name = name or f"{self.method}_clean"

        ensemble = 10 if ensemble is None else max(1, int(ensemble))
        members: list[int] = []
        member_rules = None
        if ensemble > 1:
            available = _core.ensemble_prepare(
                str(self.fit.run_path), int(self.fit.dataset.num_threads)
            )
            if available < 2:
                print(
                    "Run carries no ensemble member pool (single-restart fit, or a "
                    "run cached before member pools were stored; refit with "
                    "overwrite=True to enable it); applying a single-fit correction."
                )
            else:
                selected = int(self.fit.manifest["nmf_diagnostics"]["selected_run"]) - 1
                members = [selected] + [m for m in range(available) if m != selected]
                members = members[: min(ensemble, len(members))]
                restrict = (
                    rules[["source_cell_type", "target_cell_type"]].drop_duplicates()
                    if user_rules and len(rules)
                    else None
                )
                member_rules = self._ensemble_member_rules(
                    members, selected, rules, p_thresh=p_thresh, restrict_pairs=restrict
                )

        if member_rules is None:
            return self.fit.correct(rules, name=name)

        min_votes = max(1, math.ceil(vote * len(members)))
        correction = self.fit.correct(
            member_rules["rules"],
            name=name,
            rule_members=member_rules["member"],
            min_votes=min_votes,
        )
        rules = rules.copy()
        rules["support"] = _rule_support(rules, member_rules["pair_sets"])
        correction.rules = rules
        print(
            f"Ensemble correction over {len(members)} members: molecules removed "
            f"by >= {min_votes} members (vote >= {vote:.2f}) are dropped "
            f"({correction.manifest.get('n_removed', 0):,} molecules)."
        )
        return correction

    def _score_member(self, member: int) -> "CellAdmixScore":
        """Score one ensemble member with this score's method and parameters."""
        fn = self.fit.score_membrane if self.method == "membrane" else self.fit.score_bridge
        return fn(ensemble_member=int(member), **self.params)

    def _ensemble_member_rules(
        self,
        members,
        selected,
        primary_rules: pd.DataFrame,
        *,
        p_thresh: float,
        restrict_pairs: pd.DataFrame | None = None,
    ) -> dict:
        """Per-member kept rules, mirroring the primary rule derivation."""
        from . import _core
        from ._score_utils import apply_native_check

        def pair_keys(df):
            return set(zip(df["source_cell_type"].astype(str), df["target_cell_type"].astype(str)))

        frames = []
        frame_members: list[int] = []
        pair_sets: dict[int, set] = {}
        for member in members:
            if member == selected:
                member_df = primary_rules
            else:
                cache_key = (member, p_thresh)
                member_df = self._member_rules_cache.get(cache_key)
                if member_df is None:
                    member_score = self._score_member(member)
                    member_df = member_score.rules(p_thresh=p_thresh, native_check=False)
                    if not member_df.empty:
                        raw = _core.ensemble_member_fractions(str(self.fit.run_path), int(member))
                        fr = raw["fractions"]
                        cells = self.fit.cell_factors()
                        frame = pd.DataFrame(
                            np.asarray(fr["data"], dtype=float).reshape(fr["rows"], fr["cols"]),
                            columns=[f"factor_{k + 1}_fraction" for k in range(fr["cols"])],
                        )
                        frame["cell_id"] = [str(c) for c in raw["cell_id"]]
                        merged = cells.drop(
                            columns=[c for c in cells.columns if c.endswith("_fraction")]
                        ).merge(frame, on="cell_id", how="left")
                        member_df = apply_native_check(member_df, self.fit, cell_factors=merged)
                    if "keep" in member_df.columns:
                        member_df = member_df[member_df["keep"].fillna(True).astype(bool)]
                    self._member_rules_cache[cache_key] = member_df
            if restrict_pairs is not None and len(member_df):
                allowed = pair_keys(restrict_pairs)
                keys = list(
                    zip(
                        member_df["source_cell_type"].astype(str),
                        member_df["target_cell_type"].astype(str),
                    )
                )
                member_df = member_df[[k in allowed for k in keys]]
            pair_sets[member] = pair_keys(member_df) if len(member_df) else set()
            if len(member_df):
                frames.append(member_df[["factor", "target_cell_type"]])
                frame_members.extend([int(member)] * len(member_df))
        rules = (
            pd.concat(frames, ignore_index=True)
            if frames
            else pd.DataFrame(columns=["factor", "target_cell_type"])
        )
        return {"rules": rules, "member": frame_members, "pair_sets": pair_sets}

    def plot_heatmap(self, *, p_thresh: float = 0.1, adjust_p: bool = False, **kwargs):
        from .plotting import plot_score_heatmap

        return plot_score_heatmap(self.summary, p_thresh=p_thresh, adjust_p=adjust_p, **kwargs)

    def plot_pairs(self, *, factors=None, factor=None, **kwargs):
        from .plotting import plot_score_pairs

        if factors is None and factor is not None:
            factors = factor
        if factors is None:
            factors = sorted(self.summary["factor"].dropna().unique().astype(int).tolist())
        if np.isscalar(factors):
            factors = [int(factors)]
        return plot_score_pairs(self.summary, factors=[int(x) for x in factors], **kwargs)

    def examples(self, **kwargs) -> pd.DataFrame:
        """Select target-cell examples with strong non-native factor evidence."""
        from .examples import select_example_cells

        return select_example_cells(self, **kwargs)

    def plot_example(self, example, *, p_thresh: float = 0.1, adjust_p: bool = False, **kwargs):
        """Plot one score-selected example cell, or a cell given by its ID."""
        from .examples import plot_cell_example, prepare_cell_example, select_example_cells

        score_annotation = kwargs.pop("score_annotation", None)
        score_annotation = score_annotation or self.annotation(p_thresh=p_thresh, adjust_p=adjust_p)
        if isinstance(example, str):
            example = select_example_cells(
                self, score_annotation=score_annotation, p_thresh=p_thresh,
                adjust_p=adjust_p, cells=example)
        prepare_keys = {
            "cell_data", "boundaries", "stains", "cell_types", "markers",
            "padding", "min_side", "max_pixels",
        }
        prepare_kwargs = {k: kwargs.pop(k) for k in list(kwargs) if k in prepare_keys}
        prepared = prepare_cell_example(
            self.fit, example, score_annotation=score_annotation, **prepare_kwargs)
        return plot_cell_example(prepared, **kwargs)

    def plot_examples(self, examples=None, *, p_thresh: float = 0.1, adjust_p: bool = False, cells=None, **kwargs):
        """Plot a grid of score-selected example cells.

        Pass ``cells`` (or a list of cell IDs as ``examples``) to plot specific
        cells instead of the automatically selected ones.
        """
        from .examples import plot_examples

        score_annotation = kwargs.pop("score_annotation", None)
        score_annotation = score_annotation or self.annotation(p_thresh=p_thresh, adjust_p=adjust_p)
        if cells is None and examples is not None and not hasattr(examples, "columns"):
            cells = examples
            examples = None
        if examples is None:
            examples = self.examples(
                score_annotation=score_annotation, p_thresh=p_thresh,
                adjust_p=adjust_p, cells=cells)
        return plot_examples(examples, fit=self.fit, score_annotation=score_annotation, **kwargs)
