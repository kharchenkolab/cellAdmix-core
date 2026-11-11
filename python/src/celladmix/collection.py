"""Multi-sample collections: consensus factorization across samples.

A collection groups per-sample datasets that share one cell-type annotation
scheme so that factorization, scoring, and correction treat every sample by
the same procedure: one shared factor space (a consensus of the per-sample
factorizations), one set of removal rules derived from the evidence pooled
across samples, and sample-specific molecule labels and removal amounts.
"""

from __future__ import annotations

import hashlib
from typing import Optional

import numpy as np
import pandas as pd

from . import _core


class CellAdmixCollection:
    """Per-sample datasets of a multi-sample analysis and the fits made on them.

    The workflow is ``fit()`` (one ordinary fit per sample, cached like any
    fit), ``fit_consensus()`` (match factors across samples, average them into
    shared profiles, and label every sample against those profiles), then
    scoring, rules, and correction on the returned :class:`CellAdmixConsensusFit`.
    """

    def __init__(self, datasets: dict, *, sample_allocation: str = "proportional"):
        if not isinstance(datasets, dict) or not datasets:
            raise ValueError("datasets must be a non-empty dict of sample name -> CellAdmix")
        for name, dataset in datasets.items():
            if not str(name):
                raise ValueError("sample names must be non-empty")
            if getattr(dataset, "annotation", None) is None:
                raise ValueError(f"dataset {name!r} has no annotation")
        if sample_allocation not in ("proportional", "equal"):
            raise ValueError("sample_allocation must be 'proportional' or 'equal'")
        self.datasets = {str(k): v for k, v in datasets.items()}
        self.sample_allocation = sample_allocation
        self.fits: dict = {}

    def __repr__(self) -> str:
        return (f"CellAdmixCollection(samples={list(self.datasets)}, "
                f"allocation={self.sample_allocation!r}, fitted={list(self.fits)})")

    def samples(self) -> pd.DataFrame:
        rows = []
        for name, dataset in self.datasets.items():
            rows.append({"sample": name, "annotated_cells": int(len(dataset.annotation)),
                         "cell_types": int(dataset.annotation.nunique()),
                         "output_dir": str(dataset.output_dir)})
        return pd.DataFrame(rows)

    def fit(self, **kwargs):
        """One ordinary fit per sample; cached runs are reused like any fit."""
        for name, dataset in self.datasets.items():
            self.fits[name] = dataset.fit(**kwargs)
        return self

    def fit_consensus(self, *, match_min: float = 0.3, min_restarts: int = 2, min_stability: float = 0.0,
                      name: str = "consensus", num_threads: Optional[int] = None, **fit_kwargs):
        """Cluster the pooled restart factors of all samples into shared
        profiles and project every sample onto them."""
        state = self._prepare_profiles(**fit_kwargs)
        min_restarts = max(1, int(min_restarts))
        min_stability = min(1.0, max(0.0, float(min_stability)))
        # Every restart of every sample enters the pool; each sample's primary
        # restart is marked so the matching table can describe its own factors.
        pools = {s: [h[state["genes"]] for h in _fit_member_profiles(state["fits"][s])]
                 for s in state["samples"]}
        selected = {s: int(state["fits"][s].manifest["nmf_diagnostics"]["selected_run"])
                    for s in state["samples"]}
        consensus = _consensus_profiles(pools, state["weights"], selected, match_min=match_min,
                                        min_restarts=min_restarts, min_stability=min_stability)
        if consensus["profiles"].empty:
            raise ValueError(f"no factor recurs across restarts at match_min = {match_min}; "
                             "lower match_min, min_restarts, or min_stability")
        # The matching table describes each sample's primary factors: the
        # consensus factor whose cluster they joined, their ownership correlation
        # with that consensus profile, their within-cluster agreement, and their
        # restart stability within the sample.
        rows = []
        components = consensus["components"]
        for s in state["samples"]:
            comp = components[(components["sample"] == s) & components["selected"]].sort_values("factor")
            corr = _profile_correlation(pools[s][selected[s] - 1], consensus["profiles"])
            for _, row in comp.iterrows():
                k = row["consensus_factor"]
                matched = not pd.isna(k)
                rows.append({"sample": s, "sample_factor": int(row["factor"]),
                             "consensus_factor": int(k) if matched else pd.NA,
                             "correlation": float(corr[int(row["factor"]) - 1, int(k) - 1]) if matched else np.nan,
                             "matched": matched, "agreement": row["agreement"]})
        matching = pd.DataFrame(rows)
        matching["own_stability"] = _own_factor_stability(state["fits"], matching["sample"], matching["sample_factor"])
        # One member pool shared by every sample: the restart profiles of all
        # samples aligned to the consensus factor order and interleaved, so an
        # ensemble of any size draws from each sample alike and every member
        # votes with the same decisions in every sample. The consensus itself
        # leads the pool as member 0.
        aligned = [_align_members(pools[s], consensus["profiles"], match_min=match_min) for s in state["samples"]]
        shared = _interleave(aligned)
        members = {s: shared for s in state["samples"]}
        run_id = f"{name}_{_profiles_hash(consensus['profiles'], shared)}"
        return CellAdmixConsensusFit(
            self, self._project(consensus["profiles"], members, run_id, num_threads),
            state["fits"], consensus["profiles"], matching, consensus["factors"],
            consensus["components"], consensus["unreproduced"],
            self._loss_table(state, consensus["profiles"]),
            params={"match_min": match_min, "min_restarts": min_restarts, "min_stability": min_stability,
                    "weights": state["weights"], "restarts_total": consensus["restarts_total"],
                    "name": name, "run_id": run_id, "source": "consensus"})

    def fit_profiles(self, profiles, *, members=None, name: str = "profiles",
                     match_min: float = 0.3, num_threads: Optional[int] = None, **fit_kwargs):
        """Project every sample onto given profiles (a factor x gene frame with
        gene columns): profiles of a joint factorization, or one sample's own,
        so alternatives to the consensus run through the same scoring, rules,
        and diagnostics. ``members``, when given, is a list of frames shaped
        like ``profiles`` that becomes every sample's member pool after
        alignment to the profiles' factor order."""
        state = self._prepare_profiles(**fit_kwargs)
        profiles = pd.DataFrame(profiles).copy()
        profiles.columns = [str(g) for g in profiles.columns]
        genes = _common_genes({"profiles": profiles, **state["profiles"]})
        profiles = profiles[genes]
        profiles.index = [f"F{i + 1}" for i in range(len(profiles))]
        profiles_norm = _normalize_profiles(profiles)
        rows = []
        for s in state["samples"]:
            own = _normalize_profiles(state["profiles"][s])
            match = _match_factors(own, profiles_norm)
            for i in range(own.shape[0]):
                j, r = match["match"][i], match["correlation"][i]
                ok = j >= 0 and np.isfinite(r) and r >= match_min
                rows.append({"sample": s, "sample_factor": i + 1,
                             "consensus_factor": j + 1 if ok else pd.NA,
                             "correlation": r, "matched": ok})
        matching = pd.DataFrame(rows)
        matching["own_stability"] = _own_factor_stability(state["fits"], matching["sample"], matching["sample_factor"])
        aligned = []
        if members:
            frames = []
            for member in members:
                member = pd.DataFrame(member).copy()
                member.columns = [str(g) for g in member.columns]
                frames.append(member)
            aligned = _align_members(frames, profiles, match_min=match_min)
        member_pools = {s: aligned for s in state["samples"]}
        run_id = f"{name}_{_profiles_hash(profiles, aligned)}"
        # With given profiles the per-factor table counts the samples whose
        # primary factors matched each profile.
        counts = matching[matching["matched"]]["consensus_factor"].astype(int).value_counts()
        factors = pd.DataFrame({
            "consensus_factor": range(1, len(profiles) + 1),
            "cluster": pd.NA,
            "samples": [int(counts.get(k, 0)) for k in range(1, len(profiles) + 1)]})
        factors["recurrence"] = factors["samples"] / len(state["samples"])
        for column in ("restarts", "components", "stability", "agreement"):
            factors[column] = pd.NA
        components = pd.DataFrame({
            "sample": matching["sample"], "restart": pd.NA, "factor": matching["sample_factor"],
            "group": pd.NA, "selected": True, "cluster": pd.NA, "agreement": np.nan,
            "consensus_factor": matching["consensus_factor"]})
        return CellAdmixConsensusFit(
            self, self._project(profiles, member_pools, run_id, num_threads),
            state["fits"], profiles, matching, factors, components, components.iloc[0:0],
            self._loss_table(state, profiles),
            params={"match_min": match_min, "name": name, "run_id": run_id, "source": "profiles"})

    def _prepare_profiles(self, **fit_kwargs) -> dict:
        """Per-sample fits with their raw-space profiles over the common gene
        panel and the averaging weights."""
        if not self.fits:
            self.fit(**fit_kwargs)
        fits = self.fits
        samples = list(fits)
        if len(samples) < 2:
            raise ValueError("a collection fit needs at least two samples")
        molecules = {s: float(fits[s].manifest.get("n_transcripts", 0)) for s in samples}
        total = sum(molecules.values()) or 1.0
        if self.sample_allocation == "equal":
            weights = {s: 1.0 for s in samples}
        else:
            weights = {s: molecules[s] / total for s in samples}
        profiles = {s: _fit_profiles_raw(fits[s]) for s in samples}
        genes = _common_genes(profiles)
        profiles = {s: h[genes] for s, h in profiles.items()}
        return {"fits": fits, "samples": samples, "molecules": molecules,
                "weights": weights, "profiles": profiles, "genes": genes}

    def _loss_table(self, state: dict, shared: pd.DataFrame) -> pd.DataFrame:
        """How well shared profiles represent each sample's own training rows."""
        rows = []
        for s in state["samples"]:
            run_path = str(state["fits"][s].run_path)
            own = _training_loss(run_path, state["profiles"][s])
            under_shared = _training_loss(run_path, shared)
            objectives = list(state["fits"][s].manifest.get("nmf_diagnostics", {}).get("candidate_final_objectives", []) or [])
            spread = ((max(objectives) - min(objectives)) / min(objectives)
                      if len(objectives) > 1 and min(objectives) > 0 else np.nan)
            rows.append({"sample": s, "own_loss": own, "consensus_loss": under_shared,
                         "loss_gap": under_shared / own - 1 if own > 0 else np.nan,
                         "restart_spread": spread})
        return pd.DataFrame(rows)

    def _project(self, shared: pd.DataFrame, member_pools: dict, run_id: str,
                 num_threads: Optional[int]) -> dict:
        """Fixed-profile projection of every sample under one run id."""
        projected = {}
        for s, own in self.fits.items():
            dataset = self.datasets[s]
            variant = own.manifest.get("pipeline_options", {}).get("nmf_variant", "ls_nmf")
            pool = member_pools.get(s) or []
            projected[s] = dataset.fit(
                nmf_fixed_h=shared,
                nmf_fixed_h_members=pool if pool else None,
                run_id=run_id,
                nmf_variant=variant,
                nmf_n_runs=1,
                num_threads=num_threads or dataset.num_threads,
            )
        return projected


class CellAdmixConsensusFit:
    """Shared gene profiles across samples, the factor matching that produced
    them, one fixed-profile fit per sample (molecule labels and ensemble
    members in the consensus factor index), and representation diagnostics."""

    def __init__(self, collection, fits, own_fits, profiles, matching, factors, components,
                 unreproduced, loss, params=None):
        self.collection = collection
        self.fits = fits
        self.own_fits = own_fits
        self.profiles = profiles
        self.matching = matching
        self.factors = factors
        self.components = components
        self.unreproduced = unreproduced
        self.loss = loss
        self.params = dict(params or {})

    def __repr__(self) -> str:
        single = int((self.factors["samples"] == 1).sum()) if len(self.factors) else 0
        return (f"CellAdmixConsensusFit(samples={list(self.fits)}, n_factors={len(self.profiles)}, "
                f"single_sample_factors={single}, unreproduced_components={len(self.unreproduced)})")

    def samples(self) -> list:
        return list(self.fits)

    def factor_stability(self) -> pd.DataFrame:
        """Per consensus factor and sample: the number of that sample's restarts
        with a component in the factor's cluster, and that number as a share of
        the sample's restarts (the program's stability within the sample). Zero
        rows mark samples that do not carry the factor at all."""
        comp = self.components[self.components["consensus_factor"].notna()]
        totals = self.params.get("restarts_total") or {}
        rows = []
        for k in range(1, len(self.profiles) + 1):
            for s in self.fits:
                sub = comp[(comp["consensus_factor"] == k) & (comp["sample"] == s)]
                restarts = int(sub["restart"].nunique()) if len(sub) else 0
                total = totals.get(s)
                rows.append({"consensus_factor": k, "sample": s, "restarts": restarts,
                             "restarts_total": total if total is not None else pd.NA,
                             "stability": restarts / total if total else np.nan})
        return pd.DataFrame(rows)

    def loadings(self) -> pd.DataFrame:
        """Consensus profiles, one row per factor and one column per gene."""
        return self.profiles

    def stability(self) -> pd.DataFrame:
        """Within-sample restart stability of each own fit next to the
        agreement of that sample's aligned restarts with the consensus."""
        rows = []
        for s, fit in self.fits.items():
            own = self.own_fits[s].manifest.get("nmf_diagnostics", {})
            shared = fit.manifest.get("nmf_diagnostics", {})
            rows.append({"sample": s,
                         "within_sample_stability": own.get("candidate_best_match_correlation_mean", np.nan),
                         "consensus_agreement": shared.get("candidate_best_match_correlation_mean", np.nan),
                         "members": int(fit.manifest.get("fixed_profile_members", 0))})
        return pd.DataFrame(rows)

    def sample_usage(self) -> pd.DataFrame:
        """Share of each cell type's molecules labeled with each factor, per sample."""
        frames = []
        for s, fit in self.fits.items():
            cells = fit.cell_factors()
            annotation = fit.dataset.annotation
            cells = cells.assign(cell_type=[annotation.get(str(c)) for c in cells["cell_id"].astype(str)])
            cells = cells[cells["cell_type"].notna()]
            fraction_cols = sorted(
                [c for c in cells.columns if c.startswith("factor_") and c.endswith("_fraction")],
                key=lambda c: int(c.split("_")[1]))
            molecules = cells[fraction_cols].to_numpy(dtype=float) * cells["transcript_count"].to_numpy(dtype=float)[:, None]
            by_type = pd.DataFrame(molecules, columns=range(1, len(fraction_cols) + 1)).groupby(cells["cell_type"].to_numpy()).sum()
            share = by_type.div(by_type.sum(axis=1).clip(lower=1), axis=0)
            long = by_type.stack().rename("molecules").reset_index()
            long.columns = ["cell_type", "factor", "molecules"]
            long["share"] = share.stack().to_numpy()
            long.insert(0, "sample", s)
            frames.append(long)
        return pd.concat(frames, ignore_index=True)

    def usage_divergence(self, *, by_sample: bool = False) -> pd.DataFrame:
        """Per cell type: half the L1 distance between a sample's factor usage
        and the pooled usage, maximized over samples (or, with ``by_sample``,
        one row per cell type and sample)."""
        usage = self.sample_usage()
        pooled = usage.groupby(["cell_type", "factor"], as_index=False)["molecules"].sum()
        totals = pooled.groupby("cell_type")["molecules"].sum()
        pooled["pooled_share"] = pooled["molecules"] / pooled["cell_type"].map(totals).clip(lower=1)
        merged = usage.merge(pooled[["cell_type", "factor", "pooled_share"]], on=["cell_type", "factor"])
        merged["deviation"] = (merged["share"] - merged["pooled_share"]).abs()
        per_sample = merged.groupby(["cell_type", "sample"], as_index=False)["deviation"].sum()
        per_sample["divergence"] = per_sample.pop("deviation") / 2.0
        per_sample["molecules"] = per_sample["cell_type"].map(totals).to_numpy()
        if by_sample:
            return per_sample.sort_values(["sample", "divergence"], ascending=[True, False]).reset_index(drop=True)
        out = per_sample.groupby("cell_type", as_index=False)["divergence"].max()
        out["molecules"] = out["cell_type"].map(totals).to_numpy()
        return out.sort_values("divergence", ascending=False).reset_index(drop=True)

    def diagnostics(self) -> pd.DataFrame:
        """One row per sample: how well the shared factor space represents it
        (the training loss under the shared profiles relative to its own fit,
        with the sample's restart spread as the yardstick), how many of its
        primary factors entered the shared space, its restart stability next
        to its agreement with the consensus, and the largest factor-usage
        divergence of any of its cell types. Flags are reports:
        ``representation_flag`` when the loss gap exceeds the restart spread,
        ``shared_factors_flag`` when fewer than half of the sample's primary
        factors are in the shared space."""
        loss = self.loss.set_index("sample") if len(self.loss) else pd.DataFrame()
        stability = self.stability().set_index("sample")
        deviation = self.usage_divergence(by_sample=True)
        rows = []
        for s in self.fits:
            m = self.matching[self.matching["sample"] == s]
            gap = float(loss.loc[s, "loss_gap"]) if s in loss.index else np.nan
            spread = float(loss.loc[s, "restart_spread"]) if s in loss.index else np.nan
            d = deviation[deviation["sample"] == s]
            worst = d.loc[d["divergence"].idxmax()] if len(d) else None
            shared = int(m["matched"].sum())
            rows.append({
                "sample": s, "loss_gap": gap, "restart_spread": spread,
                "representation_flag": bool(np.isfinite(gap) and np.isfinite(spread) and gap > spread),
                "primary_factors": int(len(m)), "primary_factors_shared": shared,
                "shared_fraction": shared / len(m) if len(m) else np.nan,
                "shared_factors_flag": bool(len(m) and shared / len(m) < 0.5),
                "within_sample_stability": float(stability.loc[s, "within_sample_stability"]) if s in stability.index else np.nan,
                "consensus_agreement": float(stability.loc[s, "consensus_agreement"]) if s in stability.index else np.nan,
                "max_usage_divergence": float(worst["divergence"]) if worst is not None else np.nan,
                "divergent_cell_type": str(worst["cell_type"]) if worst is not None else None})
        return pd.DataFrame(rows)

    def cells(self) -> pd.DataFrame:
        frames = []
        for s, fit in self.fits.items():
            cells = fit.cell_factors().copy()
            cells["sample"] = s
            cells["orig_cell_id"] = cells["cell_id"].astype(str)
            cells["cell_id"] = s + "_" + cells["orig_cell_id"]
            frames.append(cells)
        return pd.concat(frames, ignore_index=True)

    def counts(self, sample: Optional[str] = None):
        if sample is not None:
            return self.fits[sample].counts()
        return _bind_sample_counts({s: fit.counts() for s, fit in self.fits.items()})

    def audit_admixture(self, **kwargs):
        return CellAdmixCollectionAudit(self, {s: fit.audit_admixture(**kwargs) for s, fit in self.fits.items()})

    def score(self, method: str = "bridge", *, name: Optional[str] = None, **kwargs):
        if method not in ("bridge", "membrane"):
            raise ValueError("method must be 'bridge' or 'membrane'")
        return self._score_samples(method, name or method, **kwargs)

    def score_bridge(self, *, name: str = "bridge", **kwargs):
        return self._score_samples("bridge", name, **kwargs)

    def score_membrane(self, *, name: str = "membrane", **kwargs):
        return self._score_samples("membrane", name, **kwargs)

    def correct(self, score, rules=None, **kwargs):
        if not isinstance(score, CellAdmixCollectionScore):
            raise TypeError("score must be a CellAdmixCollectionScore")
        return score.correct(rules=rules, **kwargs)

    def _allocation_shares(self) -> dict:
        samples = list(self.fits)
        if self.collection.sample_allocation == "equal":
            return {s: 1.0 for s in samples}
        molecules = {s: float(self.fits[s].manifest.get("n_transcripts", 0)) for s in samples}
        total = sum(molecules.values()) or 1.0
        return {s: len(samples) * molecules[s] / total for s in samples}

    def _score_samples(self, method: str, name: str, **kwargs):
        defaults = ({"candidate_pairs_per_type_pair": 1000, "max_cells_per_type_pair": 400}
                    if method == "bridge" else
                    {"candidate_pairs_per_type_pair": 400, "max_cells_per_type_pair": 400})
        shares = self._allocation_shares()
        scores = {}
        for s, fit in self.fits.items():
            args = dict(kwargs)
            # The candidate budget of k single-sample analyses is split across
            # the samples by the collection's allocation rule.
            for cap, base in defaults.items():
                args[cap] = max(1, int(round(kwargs.get(cap, base) * shares[s])))
            fn = fit.score_membrane if method == "membrane" else fit.score_bridge
            scores[s] = fn(**args)
        return CellAdmixCollectionScore(self, method, name, scores, params=kwargs)


class CellAdmixCollectionScore:
    """Per-sample scores of a consensus fit together with the pooled summary.

    Per target/source/factor, pair counts are summed, mean scores are
    pair-weighted, and the one-sided p-values are combined by Stouffer's
    method with weights proportional to the square root of each sample's pair
    count. Rules derived from the pooled summary carry the fraction of samples
    whose own evidence supports them (``sample_support``), as a report, never
    a gate.
    """

    def __init__(self, consensus, method, name, sample_scores, params=None):
        self.consensus = consensus
        self.method = method
        self.name = name
        self.sample_scores = sample_scores
        self.params = dict(params or {})
        self._pooled = None

    def __repr__(self) -> str:
        return f"CellAdmixCollectionScore(method={self.method!r}, samples={list(self.sample_scores)})"

    def samples(self) -> list:
        return list(self.sample_scores)

    def sample_summaries(self) -> pd.DataFrame:
        frames = [score.summary.assign(sample=s) for s, score in self.sample_scores.items() if len(score.summary)]
        if not frames:
            return pd.DataFrame(columns=["sample", "target_cell_type", "source_cell_type", "factor",
                                         "n_pairs", "mean_score", "p_value", "neg_log10_p"])
        out = pd.concat(frames, ignore_index=True)
        return out[["sample"] + [c for c in out.columns if c != "sample"]]

    @property
    def pairs(self) -> pd.DataFrame:
        frames = [score.pairs.assign(sample=s) for s, score in self.sample_scores.items() if len(score.pairs)]
        if not frames:
            return pd.DataFrame()
        out = pd.concat(frames, ignore_index=True)
        return out[["sample"] + [c for c in out.columns if c != "sample"]]

    @property
    def summary(self) -> pd.DataFrame:
        if self._pooled is None:
            self._pooled = _pool_score_summaries(self.sample_summaries())
        return self._pooled

    def annotation(self, *, p_thresh: float = 0.1, adjust_p: bool = False) -> dict:
        from ._score_utils import score_annotation

        return score_annotation(self.summary, p_thresh=p_thresh, adjust_p=adjust_p)

    def rules(self, *, p_thresh: float = 0.1, adjust_p: bool = False, targets=None,
              native_check: bool = True, native_median_thresh: float = 0.1,
              native_expr_thresh: float = 0.05, native_outlier_min_frac: float = 0.1,
              neighbor_k: int = 15) -> pd.DataFrame:
        from ._score_utils import rules_from_annotation

        if self.summary.empty:
            return pd.DataFrame(columns=["factor", "source_cell_type", "target_cell_type",
                                         "p_value", "neg_log10_p", "rule_id", "sample_support"])
        rules = rules_from_annotation(self.annotation(p_thresh=p_thresh, adjust_p=adjust_p),
                                      target_cell_types=targets)
        rules = rules.sort_values(["p_value", "factor", "target_cell_type"]).reset_index(drop=True)
        rules = self._add_sample_support(rules, p_thresh)
        if native_check:
            if rules.empty:
                from ._score_utils import apply_native_check

                rules = apply_native_check(rules, None)
            else:
                rules = self._pooled_native_check(rules, native_median_thresh, native_expr_thresh,
                                                  native_outlier_min_frac, neighbor_k)
        return rules

    def plot_heatmap(self, *, p_thresh: float = 0.1, adjust_p: bool = False,
                     native_check: bool = True, **kwargs):
        """The pooled evidence drawn like a single sample's scoring matrix:
        rows cell types, columns shared factors; removal calls the
        native-factor check rejects are marked ``(*)``."""
        from ._score_utils import vetoed_calls_from_rules
        from .plotting import plot_score_heatmap

        vetoed = None
        if native_check:
            try:
                flagged = self.rules(p_thresh=p_thresh, adjust_p=adjust_p)
            except Exception:
                flagged = None
            vetoed = vetoed_calls_from_rules(flagged)
        return plot_score_heatmap(self.summary, p_thresh=p_thresh, adjust_p=adjust_p,
                                  vetoed_calls=vetoed, **kwargs)

    def correct(self, rules=None, *, name: Optional[str] = None, p_thresh: float = 0.1,
                adjust_p: bool = False, targets=None, ensemble: Optional[int] = None,
                vote: float = 0.3, apply_rules: str = "everywhere",
                num_threads: Optional[int] = None):
        """Apply the pooled rules to every sample.

        The molecule-vote ensemble is shared across samples: the member pool
        is the same in every sample's run, each member is scored in every
        sample, its evidence pooled like the primary evidence, and one member
        rule set votes in every sample. By default every rule is applied in
        every sample (``apply_rules="everywhere"``); with
        ``apply_rules="where_present"`` a rule is applied in a sample only
        when at least one of that sample's own NMF restarts carries the
        rule's program, so a sample is not cleaned of a program its own data
        never finds (the molecules labeled with it there are nearest-factor
        artifacts). The returned rules list the samples whose restarts do not
        carry each rule's program in ``samples_without_program``. The
        per-sample steps (labeling each
        sample's molecules with every member, scoring every member in every
        sample, and the vote in each sample) are independent jobs that run
        concurrently in forked worker processes, as many at a time as fit a
        thread budget given the threads each job uses; ``num_threads`` sets
        the budget (default: the ``CELLADMIX_NUM_THREADS`` environment
        variable, or else the number of available cores).
        """
        import math

        from . import _core

        if apply_rules not in ("everywhere", "where_present"):
            raise ValueError("apply_rules must be 'everywhere' or 'where_present'")
        if rules is None:
            rules = self.rules(p_thresh=p_thresh, adjust_p=adjust_p, targets=targets)
        rules = rules.copy()
        kept = rules
        if "keep" in getattr(kept, "columns", ()):
            kept = kept[kept["keep"].fillna(True).astype(bool)]
        if kept.empty:
            raise ValueError("no rules to apply")
        name = name or f"{self.name}_clean"
        fits = self.consensus.fits
        samples = list(fits)

        # The programs each sample's own restarts carry, from the consensus
        # components (None for a fit whose profiles came from elsewhere).
        carried = self._carried_factors()
        if apply_rules == "where_present" and carried is None:
            raise ValueError("apply_rules='where_present' needs a consensus fit, whose restart "
                             "components tell which programs each sample carries")
        if carried is not None:
            rules["samples_without_program"] = [
                ",".join(s for s in samples if int(f) not in carried[s]) for f in rules["factor"]]

        ensemble = 10 if ensemble is None else max(1, int(ensemble))
        member_rules = None
        members: list[int] = []
        if ensemble > 1:
            available = min(self._prepare_members(num_threads))
            if available < 2:
                print("Runs carry no shared ensemble member pool; applying single-fit corrections.")
            else:
                selected = int(fits[samples[0]].manifest["nmf_diagnostics"]["selected_run"]) - 1
                members = [selected] + [m for m in range(available) if m != selected]
                members = members[: min(ensemble, len(members))]
                member_rules = self._pooled_member_rules(members, selected, kept, p_thresh=p_thresh,
                                                         adjust_p=adjust_p, targets=targets,
                                                         num_threads=num_threads)
        min_votes = max(1, math.ceil(vote * len(members))) if members else 1
        # The rule rows applied in a sample: all of them, or only those whose
        # program the sample's own restarts carry.
        def in_sample(s):
            table = kept if member_rules is None else member_rules["rules"]
            member = None if member_rules is None else member_rules["member"]
            if apply_rules == "where_present":
                return _rules_in_sample(table, member, carried[s])
            return {"rules": table, "member": member}

        def vote_in_sample(s):
            applied = in_sample(s)
            if member_rules is None:
                return fits[s]._correct_run(applied["rules"], name=name)
            return fits[s]._correct_run(applied["rules"], name=name, rule_members=applied["member"],
                                        min_votes=min_votes)

        # The vote in each sample is an independent single-threaded job; the
        # samples are processed concurrently within the thread budget.
        manifests = _run_jobs(vote_in_sample, samples,
                              workers=_parallel_workers(num_threads, 1, len(samples)))
        corrections = {s: fits[s]._make_correction(name, manifest, in_sample(s)["rules"])
                       for s, manifest in zip(samples, manifests)}
        if member_rules is not None:
            from .score import _rule_support

            rules["support"] = _rule_support(rules, member_rules["pair_sets"])
            removed = ", ".join(f"{s}: {c.manifest.get('n_removed', 0):,}" for s, c in corrections.items())
            skipped = ("; rules on programs a sample's restarts do not carry were skipped there"
                       if apply_rules == "where_present" else "")
            print(f"Pooled ensemble correction over {len(members)} members shared by {len(samples)} "
                  f"samples: molecules removed by >= {min_votes} members (vote >= {vote:.2f}) are "
                  f"dropped ({removed}){skipped}.")
        for correction in corrections.values():
            correction.rules = rules
        return CellAdmixCollectionCorrection(self, name, corrections, rules=rules,
                                             params={"apply_rules": apply_rules})

    def _carried_factors(self):
        """Per sample, the consensus factors at least one of the sample's own
        restarts carries; None when the fit's profiles did not come from a
        consensus of the samples' restarts."""
        consensus = self.consensus
        if consensus.params.get("source") != "consensus":
            return None
        stability = consensus.factor_stability()
        return {s: set(stability[(stability["sample"] == s) & (stability["restarts"] > 0)]
                       ["consensus_factor"].astype(int)) for s in consensus.fits}

    def _member_summary(self, sample: str, member: int) -> pd.DataFrame:
        """One ensemble member's score summary in one sample, with the same
        method and parameters as the primary score; persisted in the sample's
        run so re-corrections reuse it."""
        fit = self.consensus.fits[sample]
        score = self.sample_scores[sample]
        from .dataset import _annotation_hash

        key = "|".join([self.method, str(member), _annotation_hash(fit.dataset.annotation),
                        repr(sorted(score.params.items()))])
        cache_dir = fit.run_path / "scores"
        cache_dir.mkdir(exist_ok=True)
        path = cache_dir / f"collection_member_summary_{self.method}_m{member}.pkl"
        if path.exists():
            try:
                stored = pd.read_pickle(path)
            except Exception:
                stored = None
            if stored is not None and stored.get("key") == key:
                return stored["summary"]
        summary = score._score_member(member).summary
        pd.to_pickle({"key": key, "summary": summary}, path)
        return summary

    def _member_cell_factors(self, fit, member: int) -> pd.DataFrame:
        """Cell factor fractions of a run under one member's labeling."""
        from . import _core

        raw = _core.ensemble_member_fractions(str(fit.run_path), int(member))
        fr = raw["fractions"]
        cells = fit.cell_factors()
        frame = pd.DataFrame(np.asarray(fr["data"], dtype=float).reshape(fr["rows"], fr["cols"]),
                             columns=[f"factor_{k + 1}_fraction" for k in range(fr["cols"])])
        frame["cell_id"] = [str(c) for c in raw["cell_id"]]
        return cells.drop(columns=[c for c in cells.columns if c.endswith("_fraction")]).merge(
            frame, on="cell_id", how="left")

    def _threads_per_sample(self) -> int:
        """Threads one job on a sample uses: the largest dataset thread count."""
        return max(int(getattr(fit.dataset, "num_threads", 1) or 1) for fit in self.consensus.fits.values())

    def _prepare_members(self, num_threads: Optional[int] = None) -> list:
        """Make sure every sample's run carries one labeling of its molecules
        per member of the shared pool (written once, then reused) and return
        the number of members available per sample. The samples are prepared
        concurrently, each with its dataset's thread count."""
        from . import _core

        fits = self.consensus.fits

        def prepare(sample):
            fit = fits[sample]
            return int(_core.ensemble_prepare(str(fit.run_path), int(fit.dataset.num_threads)))

        return _run_jobs(prepare, list(fits),
                         workers=_parallel_workers(num_threads, self._threads_per_sample(), len(fits)))

    def _score_members(self, members, num_threads: Optional[int] = None) -> None:
        """Score the given ensemble members in every sample. The (sample,
        member) jobs are independent and run concurrently, as many as the
        thread budget allows given each dataset's own thread count; each
        worker writes its summary to the run's cache, which the parent then
        reads."""
        jobs = [(s, int(m)) for m in members for s in self.sample_scores]
        if not jobs:
            return
        _run_jobs(lambda job: self._member_summary(*job), jobs,
                  workers=_parallel_workers(num_threads, self._threads_per_sample(), len(jobs)))

    def _pooled_member_rules(self, members, selected, primary_rules, *, p_thresh, adjust_p, targets,
                             num_threads: Optional[int] = None) -> dict:
        """Per-member rule sets for the shared molecule-vote ensemble: each
        member's summaries from every sample are pooled like the primary
        summaries, thresholded, vetted by the pooled native-factor check on
        the member's own cell factor fractions, and restricted to the primary
        rules' source -> target pairs."""
        from ._score_utils import rules_from_annotation, score_annotation

        self._score_members([m for m in members if m != selected], num_threads)

        def pair_keys(df):
            return set(zip(df["source_cell_type"].astype(str), df["target_cell_type"].astype(str)))

        restrict = pair_keys(primary_rules)
        frames, frame_members, pair_sets = [], [], {}
        for member in members:
            if member == selected:
                member_df = primary_rules
            else:
                summaries = [self._member_summary(s, member).assign(sample=s) for s in self.sample_scores]
                summaries = [f for f in summaries if len(f)]
                pooled = _pool_score_summaries(pd.concat(summaries, ignore_index=True)) if summaries else pd.DataFrame()
                if pooled.empty:
                    member_df = primary_rules.iloc[0:0]
                else:
                    member_df = rules_from_annotation(score_annotation(pooled, p_thresh=p_thresh, adjust_p=adjust_p),
                                                      target_cell_types=targets)
                    if not member_df.empty:
                        member_df = self._pooled_native_check(
                            member_df, 0.1, 0.05, 0.1, 15,
                            cell_factors_fn=lambda fit, m=member: self._member_cell_factors(fit, m))
                        member_df = member_df[member_df["keep"].fillna(True).astype(bool)]
            if len(member_df):
                keys = list(zip(member_df["source_cell_type"].astype(str), member_df["target_cell_type"].astype(str)))
                member_df = member_df[[k in restrict for k in keys]]
            pair_sets[member] = pair_keys(member_df) if len(member_df) else set()
            if len(member_df):
                frames.append(member_df[["factor", "target_cell_type"]])
                frame_members.extend([int(member)] * len(member_df))
        rules = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["factor", "target_cell_type"])
        return {"rules": rules, "member": frame_members, "pair_sets": pair_sets}

    def _add_sample_support(self, rules: pd.DataFrame, p_thresh: float) -> pd.DataFrame:
        """Fraction of samples whose own summary reaches p_thresh for the
        rule's factor and target (any source), with their names. Report only."""
        summaries = self.sample_summaries()
        n_samples = len(self.sample_scores)
        support, names = [], []
        for _, rule in rules.iterrows():
            rows = summaries[(summaries["factor"].astype(int) == int(rule["factor"]))
                             & (summaries["target_cell_type"].astype(str) == str(rule["target_cell_type"]))]
            if rows.empty:
                support.append(0.0)
                names.append("")
                continue
            best = rows.groupby("sample")["p_value"].min()
            hits = sorted(best.index[best < p_thresh].tolist())
            support.append(len(hits) / n_samples)
            names.append(",".join(hits))
        rules = rules.copy()
        rules["sample_support"] = support
        rules["samples_significant"] = names
        return rules

    def _pooled_native_check(self, rules, median_thresh, expr_thresh, outlier_min_frac, neighbor_k,
                             cell_factors_fn=None):
        """Source exposure computed within each sample and stacked, so the
        check's statistics pool the samples without neighborhoods crossing.
        ``cell_factors_fn`` supplies the per-sample cell factor fractions: the
        run's own (default) or one ensemble member's."""
        from ._score_utils import apply_native_check, source_exposure_counts

        if cell_factors_fn is None:
            cell_factors_fn = lambda fit: fit.cell_factors()  # noqa: E731
        fits = self.consensus.fits
        annotations = {s: fit.dataset.annotation for s, fit in fits.items()}
        types = sorted({str(t) for ann in annotations.values() for t in ann.dropna().unique()})
        cell_frames, count_blocks, ann_parts = [], [], []
        for s, fit in fits.items():
            cells = cell_factors_fn(fit).copy()
            counts, _ = source_exposure_counts(cells, annotations[s], neighbor_k=neighbor_k, types=types)
            cells["cell_id"] = s + "_" + cells["cell_id"].astype(str)
            ann = annotations[s].copy()
            ann.index = s + "_" + ann.index.astype(str)
            cell_frames.append(cells)
            count_blocks.append(counts)
            ann_parts.append(ann)
        cells = pd.concat(cell_frames, ignore_index=True)
        counts = np.vstack(count_blocks)
        annotation = pd.concat(ann_parts)
        return apply_native_check(rules, None, median_thresh=median_thresh, expr_thresh=expr_thresh,
                                  outlier_min_frac=outlier_min_frac, neighbor_k=neighbor_k,
                                  cell_factors=cells, annotation=annotation, exposure=(counts, types))


class CellAdmixCollectionCorrection:
    """The corrections of every sample under one pooled rule set;
    ``params["apply_rules"]`` records how the rules were applied."""

    def __init__(self, score, name, corrections, rules=None, params=None):
        self.score = score
        self.name = name
        self.corrections = corrections
        self.rules = rules
        self.params = dict(params or {})

    def __repr__(self) -> str:
        removed = {s: c.manifest.get("n_removed") for s, c in self.corrections.items()}
        return f"CellAdmixCollectionCorrection(name={self.name!r}, n_removed={removed})"

    def samples(self) -> list:
        return list(self.corrections)

    def diagnostics(self, report=None) -> pd.DataFrame:
        """One row per sample: the number of rules applied in the sample;
        among them, the number that lack local evidence (the sample is absent
        from the rule's ``samples_significant``) and the fraction of the
        sample's removed molecules attributable to them; the number of rules
        whose program none of the sample's own restarts carries (applied
        there under ``apply_rules="everywhere"``, skipped under
        ``"where_present"``) and the fraction of its removed molecules
        attributable to those (a removed molecule is attributed to the rule
        of its cell's type and its primary factor label); and, when a cleanup
        report is given, the audit's leakage removed, median pair
        sensitivity, and own-marker false removal, with
        ``false_removal_flag`` when a sample's false removal exceeds twice
        the median across samples. Flags are reports."""
        rules = self.rules
        if rules is not None and "keep" in rules.columns:
            rules = rules[rules["keep"].fillna(True).astype(bool)]
        consensus = self.score.consensus
        fits = consensus.fits
        stability = consensus.factor_stability() if consensus.params.get("source") == "consensus" else None
        where_present = self.params.get("apply_rules") == "where_present"
        rows = []
        for s, correction in self.corrections.items():
            removed = _removed_by_type_factor(fits[s], correction)
            total = float(removed["removed"].sum()) if len(removed) else 0.0
            if rules is None or rules.empty:
                rows.append({"sample": s, "rules_applied": 0, "rules_without_local_evidence": 0,
                             "removed_without_local_evidence": np.nan, "rules_on_absent_programs": pd.NA,
                             "removed_via_absent_programs": np.nan, "molecules_removed": total})
                continue
            lookup = {(str(t), int(f)): float(r) for t, f, r in zip(removed["cell_type"], removed["factor"], removed["removed"])}

            def attributed(subset):
                return sum(lookup.get((str(t), int(f)), 0.0) for t, f in zip(subset["target_cell_type"], subset["factor"])) / max(total, 1.0)

            absent = None
            applied = rules
            if stability is not None:
                carried = set(stability[(stability["sample"] == s) & (stability["restarts"] > 0)]["consensus_factor"].astype(int))
                on_carried = rules["factor"].astype(int).isin(carried).to_numpy(dtype=bool)
                absent = rules[~on_carried]
                if where_present:
                    applied = rules[on_carried]
            significant = applied["samples_significant"].fillna("").astype(str) if "samples_significant" in applied.columns \
                else pd.Series([""] * len(applied), index=applied.index)
            has_local = significant.map(lambda x: s in x.split(",") if x else False).to_numpy(dtype=bool)
            without = applied[~has_local]
            rows.append({"sample": s, "rules_applied": int(len(applied)),
                         "rules_without_local_evidence": int(len(without)),
                         "removed_without_local_evidence": attributed(without),
                         "rules_on_absent_programs": int(len(absent)) if absent is not None else pd.NA,
                         "removed_via_absent_programs": attributed(absent) if absent is not None else np.nan,
                         "molecules_removed": total})
        out = pd.DataFrame(rows)
        if report is not None:
            summary = report.summary()
            if "sample" in summary.columns:
                summary = summary.set_index("sample")
                for column in ("leakage_removed_overall", "median_pair_sensitivity", "own_marker_false_removal"):
                    if column in summary.columns:
                        out[column] = out["sample"].map(summary[column]).to_numpy()
                if "own_marker_false_removal" in out.columns:
                    reference = out["own_marker_false_removal"].median()
                    out["false_removal_flag"] = (len(out) > 1) & (out["own_marker_false_removal"] > 2 * reference)
        return out

    def summary(self, *, by_sample: bool = True) -> pd.DataFrame:
        frames = [c.summary().assign(sample=s) for s, c in self.corrections.items()]
        per_sample = pd.concat(frames, ignore_index=True)
        per_sample = per_sample[["sample"] + [c for c in per_sample.columns if c != "sample"]]
        if by_sample:
            return per_sample
        pooled = per_sample.groupby("cell_type", as_index=False)[
            ["n_cells", "n_modified_cells", "molecules_before", "molecules_after", "molecules_removed"]].sum()
        pooled["fraction_cells_modified"] = pooled["n_modified_cells"] / pooled["n_cells"].clip(lower=1)
        pooled["fraction_molecules_removed"] = pooled["molecules_removed"] / pooled["molecules_before"].clip(lower=1)
        return pooled

    def ensemble(self) -> dict:
        return {s: c.ensemble() for s, c in self.corrections.items()}

    def cells(self) -> pd.DataFrame:
        frames = []
        for s, correction in self.corrections.items():
            cells = pd.read_parquet(correction.manifest["paths"]["cells_parquet"]).copy()
            cells["sample"] = s
            cells["orig_cell_id"] = cells["cell_id"].astype(str)
            cells["cell_id"] = s + "_" + cells["orig_cell_id"]
            frames.append(cells)
        return pd.concat(frames, ignore_index=True)

    def counts(self, sample: Optional[str] = None):
        if sample is not None:
            return self.corrections[sample].counts()
        return _bind_sample_counts({s: c.counts() for s, c in self.corrections.items()})


class CellAdmixCollectionAudit:
    """One audit per sample, with pooled per-pair estimates and per-sample
    cleanup verification. ``plot_map()`` draws each sample's admixture map
    before any correction, one row per sample, on one color scale; the map
    after a correction is ``evaluate(correction).plot_map()``."""

    def __init__(self, consensus, audits):
        self.consensus = consensus
        self.audits = audits

    def samples(self) -> list:
        return list(self.audits)

    def pairs(self, *, by_sample: bool = False, detected_only: bool = False) -> pd.DataFrame:
        frames = [audit.pairs(detected_only=detected_only).assign(sample=s) for s, audit in self.audits.items()]
        per_sample = pd.concat(frames, ignore_index=True)
        per_sample = per_sample[["sample"] + [c for c in per_sample.columns if c != "sample"]]
        pooled_columns = ["source", "target", "admixed_molecules", "rate", "detected",
                          "n_samples", "n_samples_detected"]
        if by_sample:
            return per_sample
        if per_sample.empty or "rate" not in per_sample.columns:
            return pd.DataFrame(columns=pooled_columns)
        # Pooled estimate per pair: admixed molecules summed, the rate taken
        # over the samples' summed target molecule counts, detection when any
        # sample detects the pair.
        rate = per_sample["rate"].astype(float)
        per_sample["target_molecules"] = np.where(rate > 0, per_sample["admixed_molecules"] / rate.where(rate > 0, 1), np.nan)
        rows = []
        for (source, target), group in per_sample.groupby(["source", "target"], sort=False):
            target_molecules = float(np.nansum(group["target_molecules"]))
            admixed = float(group["admixed_molecules"].sum())
            rows.append({"source": source, "target": target, "admixed_molecules": admixed,
                         "rate": admixed / target_molecules if target_molecules > 0 else 0.0,
                         "detected": bool(group["detected"].any()), "n_samples": int(len(group)),
                         "n_samples_detected": int(group["detected"].sum())})
        pooled = pd.DataFrame(rows)
        if detected_only and not pooled.empty:
            pooled = pooled[pooled["detected"]]
        return pooled.sort_values("admixed_molecules", ascending=False).reset_index(drop=True)

    def plot_map(self, value: str = "rate", detected_only: bool = True, samples=None, axes=None):
        """Each sample's admixture map before any correction, one row per
        sample (``samples`` selects a subset), on one color scale; returns the
        axes array."""
        from .audit import _audit_pairs, _plot_pair_maps

        chosen = list(self.audits) if samples is None else [s for s in samples if s in self.audits]
        return _plot_pair_maps([[(sample, _audit_pairs(self.audits[sample], detected_only))]
                                for sample in chosen], value, axes)

    def evaluate(self, correction, **kwargs):
        corrections = correction.corrections if isinstance(correction, CellAdmixCollectionCorrection) else correction
        reports = {s: audit.evaluate(corrections[s], **kwargs) for s, audit in self.audits.items()}
        return CellAdmixCollectionCleanupReport(reports)


class CellAdmixCollectionCleanupReport:
    """The cleanup report of every sample, with stacked summaries.
    ``plot_map()`` draws, one row per sample, the admixture map before the
    correction next to the admixture that remains after it, all panels on one
    color scale."""

    def __init__(self, reports):
        self.reports = reports

    def samples(self) -> list:
        return list(self.reports)

    def summary(self) -> pd.DataFrame:
        frames = []
        for s, report in self.reports.items():
            summary = report.summary()
            frame = summary if isinstance(summary, pd.DataFrame) else pd.DataFrame([summary])
            frames.append(frame.assign(sample=s))
        out = pd.concat(frames, ignore_index=True)
        return out[["sample"] + [c for c in out.columns if c != "sample"]]

    def pairs(self) -> pd.DataFrame:
        frames = [report.pairs().assign(sample=s) for s, report in self.reports.items()]
        out = pd.concat(frames, ignore_index=True)
        return out[["sample"] + [c for c in out.columns if c != "sample"]]

    def _selected(self, samples) -> list:
        return list(self.reports) if samples is None else [s for s in samples if s in self.reports]

    def plot_map(self, value: str = "rate", samples=None, axes=None):
        """Per sample, the admixture map before the correction next to what
        remains after it, all panels on one color scale; returns the axes array."""
        from .audit import _audit_pairs, _plot_pair_maps, _remaining_pairs

        return _plot_pair_maps([[(f"{s}: before correction", _audit_pairs(self.reports[s].audit)),
                                 (f"{s}: after correction", _remaining_pairs(self.reports[s]))]
                                for s in self._selected(samples)], value, axes)

    def plot_remaining(self, samples=None, ax=None):
        """The estimated admixture of each sample's molecules before and after
        the correction, as bars grouped by sample."""
        import matplotlib.pyplot as plt

        from .audit import _draw_remaining

        frames = []
        for s in self._selected(samples):
            report = self.reports[s]
            frames.append(report.audit.remaining({"corrected": report.correction}).assign(sample=s))
        if ax is None:
            _, ax = plt.subplots(figsize=(max(4.4, 0.9 * len(frames) + 2.5), 3.6))
        return _draw_remaining(ax, pd.concat(frames, ignore_index=True))

    def plot_exposure(self, source=None, target=None, strict: bool = True, samples=None,
                      ncols: int = 4, axes=None):
        """Each sample's exposure profile before and after the correction, one
        panel per sample in a grid of ``ncols`` columns; samples without a
        detected pair are left out. Returns the axes array."""
        import matplotlib.pyplot as plt

        from .audit import _draw_exposure

        frames = {}
        for s in self._selected(samples):
            report = self.reports[s]
            try:
                frames[s] = report.audit.exposure(source, target, report.correction, strict)
            except ValueError:
                continue
        if not frames:
            raise ValueError("No detected admixture pairs")
        n = len(frames)
        ncols = max(1, min(ncols, n))
        nrows = -(-n // ncols)
        if axes is None:
            fig, axes = plt.subplots(nrows, ncols, figsize=(4.8 * ncols, 3.8 * nrows),
                                     squeeze=False, layout="constrained")
        else:
            axes = np.asarray(axes).reshape(nrows, ncols)
        flat = axes.ravel()
        for ax, (s, df) in zip(flat, frames.items()):
            _draw_exposure(ax, df, pooled=source is None, source=source, target=target)
            ax.set_title(s if source is None else f"{s}: {source} \u2192 {target}", fontsize=10)
        for ax in flat[n:]:
            ax.set_axis_off()
        return axes


# ---- consensus helpers -------------------------------------------------------

def _removed_by_type_factor(fit, correction) -> pd.DataFrame:
    """Molecules removed from a sample per annotated cell type and primary
    factor label: each cell's per-factor molecule counts before the
    correction (factor fractions times molecule count in the fit) minus those
    after it."""
    before = fit.cell_factors()
    after = pd.read_parquet(correction.manifest["paths"]["cells_parquet"])
    annotation = fit.dataset.annotation
    fraction_cols = sorted([c for c in before.columns if c.startswith("factor_") and c.endswith("_fraction")],
                           key=lambda c: int(c.split("_")[1]))
    b = before.set_index(before["cell_id"].astype(str))
    a = after.set_index(after["cell_id"].astype(str)).reindex(b.index)
    mass_before = b[fraction_cols].to_numpy(dtype=float) * b["transcript_count"].to_numpy(dtype=float)[:, None]
    mass_after = a[fraction_cols].fillna(0.0).to_numpy(dtype=float) * a["transcript_count"].fillna(0.0).to_numpy(dtype=float)[:, None]
    removed = np.clip(mass_before - mass_after, 0.0, None)
    types = pd.Series([annotation.get(c) for c in b.index], index=b.index)
    keep = types.notna().to_numpy()
    if not keep.any():
        return pd.DataFrame(columns=["cell_type", "factor", "removed"])
    frame = pd.DataFrame(removed[keep], columns=range(1, len(fraction_cols) + 1))
    frame["cell_type"] = types[keep].to_numpy()
    by_type = frame.groupby("cell_type").sum()
    by_type.columns.name = "factor"
    return by_type.stack().rename("removed").reset_index()


def _fit_profiles_raw(fit, h=None) -> pd.DataFrame:
    """A fit's gene profiles (factor x gene frame) in raw gene space.

    Weighted least-squares fits store raw-space profiles. The inverse-square-
    root KL variant fits in a feature space where every gene is scaled by a
    weight recorded in the run; dividing by that weight returns the profiles
    to raw space, where profiles of different samples can be compared and
    averaged.
    """
    if h is None:
        frame = fit.factor_loadings().T.copy()
    else:
        frame = pd.DataFrame(h)
    frame.columns = [str(g) for g in frame.columns]
    variant = fit.manifest.get("pipeline_options", {}).get("nmf_variant", "ls_nmf")
    weights = list(fit.manifest.get("nmf_gene_weights", []) or [])
    if variant == "invsqrt_kl" and len(weights) == frame.shape[1]:
        scale = np.array([1.0 / w if w > 0 else 0.0 for w in weights])
        frame = frame * scale[None, :]
    frame.index = [f"F{i + 1}" for i in range(len(frame))]
    return frame


def _own_factor_stability(fits: dict, samples, factors) -> list:
    """Restart stability of given factors within their own sample's fit (the
    fit's matched ownership stability per factor), NaN when not recorded."""
    out = []
    for sample, factor in zip(list(samples), list(factors)):
        stability = fits[sample].manifest.get("nmf_diagnostics", {}).get("selected_factor_stability", []) or []
        k = int(factor) if not pd.isna(factor) else 0
        out.append(float(stability[k - 1]) if 1 <= k <= len(stability) else np.nan)
    return out


def _fit_member_profiles(fit) -> list:
    """Every persisted restart profile matrix of a fit, in raw gene space."""
    raw = _core.ensemble_h(str(fit.run_path))
    genes = [str(g) for g in raw["genes"]]
    out = []
    for member in raw["members"]:
        matrix = np.asarray(member["data"], dtype=float).reshape(int(member["rows"]), int(member["cols"]))
        out.append(_fit_profiles_raw(fit, pd.DataFrame(matrix, columns=genes)))
    return out


def _common_genes(profiles: dict) -> list:
    """Gene panel shared by all profile frames; the panels must be identical."""
    first = next(iter(profiles.values()))
    genes = list(first.columns)
    for frame in profiles.values():
        if set(frame.columns) != set(genes):
            differ = sorted(set(frame.columns) ^ set(genes))
            raise ValueError(f"samples do not share one gene panel ({len(differ)} genes differ, e.g. {differ[:3]})")
    return genes


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Column-wise weighted median of a (components x genes) array: the
    midpoint of the lower weighted median (the smallest value whose
    cumulative weight reaches half the total) and the upper one (the largest
    value whose cumulative weight from above reaches half), which reduces to
    the ordinary median for equal weights. The R binding uses the same
    definition."""
    order = np.argsort(values, axis=0, kind="stable")
    sorted_values = np.take_along_axis(values, order, axis=0)
    sorted_weights = weights[order]
    half = weights.sum() / 2.0
    columns = np.arange(values.shape[1])
    lower = sorted_values[(np.cumsum(sorted_weights, axis=0) >= half).argmax(axis=0), columns]
    from_above = np.cumsum(sorted_weights[::-1], axis=0)[::-1]
    upper_index = values.shape[0] - 1 - (from_above[::-1] >= half).argmax(axis=0)
    upper = sorted_values[upper_index, columns]
    return (lower + upper) / 2.0


def _normalize_profiles(frame: pd.DataFrame) -> pd.DataFrame:
    sums = frame.sum(axis=1).replace(0, 1)
    return frame.div(sums, axis=0)


def _match_factors(profiles: pd.DataFrame, reference: pd.DataFrame) -> dict:
    """One-to-one matching of the rows of ``profiles`` to those of
    ``reference`` by gene-ownership correlation (0-based, -1 for none)."""
    raw = _core.match_factors(
        profiles.to_numpy(dtype=float).ravel(order="C").tolist(), int(profiles.shape[0]),
        reference.to_numpy(dtype=float).ravel(order="C").tolist(), int(reference.shape[0]))
    return {"match": [int(m) for m in raw["match"]], "correlation": [float(c) for c in raw["correlation"]]}


def _profile_correlation(profiles: pd.DataFrame, reference: pd.DataFrame) -> np.ndarray:
    """Gene-ownership correlation between every row of ``profiles`` and every
    row of ``reference`` (each set normalized within itself)."""
    raw = _core.profile_correlation(
        profiles.to_numpy(dtype=float).ravel(order="C").tolist(), int(profiles.shape[0]),
        reference.to_numpy(dtype=float).ravel(order="C").tolist(), int(reference.shape[0]))
    return np.asarray(raw["data"], dtype=float).reshape(int(raw["rows"]), int(raw["cols"]))


def _consensus_profiles(pools: dict, weights: dict, selected: dict, *, match_min: float,
                        min_restarts: int, min_stability: float = 0.0) -> dict:
    """Shared profiles by pooled clustering of every restart of every sample.

    ``pools`` maps each sample to the list of its restart profile frames (raw
    gene space, common gene columns); ``selected`` gives each sample's
    primary restart (1-based). The pooled components are clustered by average
    linkage on gene-ownership correlation, cut at ``match_min``; a cluster
    whose components come from at least ``min_restarts`` distinct restarts, of
    any samples, and from at least a share ``min_stability`` of the
    restarts of some sample, becomes a consensus factor, and the remaining
    components are reported as unreproduced. The consensus profile of a cluster is the
    weighted per-gene median of its unit-sum component profiles (every
    component weighted by its sample's weight, so a program found by more of
    a sample's restarts carries more of that sample's evidence), renormalized
    to unit sum and rescaled to the weighted mean of the original row sums so
    that labeling sees the usual profile scale; the median keeps the profile
    of a loosely agreeing cluster from blurring.
    Factors are ordered by the number of samples covered, then restarts,
    then profile mass."""
    genes = list(next(iter(pools.values()))[0].columns)
    blocks, meta = [], []
    group = 0
    for s, pool in pools.items():
        for r, h in enumerate(pool, start=1):
            group += 1
            blocks.append(h.to_numpy(dtype=float))
            meta.append(pd.DataFrame({"sample": s, "restart": r, "factor": range(1, h.shape[0] + 1),
                                      "group": group, "selected": r == selected[s]}))
    stacked = np.vstack(blocks)
    components = pd.concat(meta, ignore_index=True)
    clustering = _core.cluster_factors(stacked.ravel(order="C").tolist(), int(stacked.shape[0]),
                                       [int(g) for g in components["group"]], float(match_min))
    components["cluster"] = [int(c) + 1 for c in clustering["cluster"]]
    components["agreement"] = [float(a) for a in clustering["within_correlation"]]
    row_sums = stacked.sum(axis=1)
    restarts_total = {s: len(pool) for s, pool in pools.items()}
    rows = []
    for k, idx in components.groupby("cluster").indices.items():
        per_sample = components.iloc[idx].groupby("sample")["group"].nunique()
        rows.append({"cluster": int(k), "samples": int(components["sample"].iloc[idx].nunique()),
                     "restarts": int(components["group"].iloc[idx].nunique()), "components": int(len(idx)),
                     "stability": float(max(per_sample[s] / restarts_total[s] for s in per_sample.index)),
                     "agreement": float(np.nanmean(components["agreement"].iloc[idx])) if len(idx) > 1 else np.nan,
                     "mass": float(row_sums[idx].mean())})
    summary = pd.DataFrame(rows).sort_values(["samples", "restarts", "mass"], ascending=False, kind="stable")
    reproduced = summary[(summary["restarts"] >= min_restarts)
                         & (summary["stability"] >= min_stability)].reset_index(drop=True)
    sums = stacked.sum(axis=1, keepdims=True)
    sums[sums <= 0] = 1.0
    normalized = stacked / sums
    consensus = np.zeros((len(reproduced), len(genes)))
    for i, cluster in enumerate(reproduced["cluster"]):
        idx = np.flatnonzero(components["cluster"].to_numpy() == cluster)
        w = np.array([float(weights[s]) for s in components["sample"].iloc[idx]])
        profile = _weighted_median(normalized[idx], w)
        profile = profile / max(float(profile.sum()), np.finfo(float).eps)
        consensus[i] = profile * (float((w * row_sums[idx]).sum()) / float(w.sum()))
    lookup = {int(c): i + 1 for i, c in enumerate(reproduced["cluster"])}
    components["consensus_factor"] = [lookup.get(int(c), pd.NA) for c in components["cluster"]]
    factors = pd.DataFrame({
        "consensus_factor": range(1, len(reproduced) + 1),
        "cluster": reproduced["cluster"].astype(int),
        "samples": reproduced["samples"].astype(int),
        "recurrence": reproduced["samples"] / len(pools),
        "restarts": reproduced["restarts"].astype(int),
        "components": reproduced["components"].astype(int),
        "stability": reproduced["stability"],
        "agreement": reproduced["agreement"]})
    unreproduced = components[components["consensus_factor"].isna()].reset_index(drop=True)
    frame = pd.DataFrame(consensus, index=[f"F{i + 1}" for i in range(len(reproduced))], columns=genes)
    return {"profiles": frame, "components": components, "factors": factors,
            "unreproduced": unreproduced, "restarts_total": restarts_total}


def _align_members(pool: list, consensus: pd.DataFrame, *, match_min: float) -> list:
    """Restart profiles of one sample reordered to the consensus factor index;
    consensus factors without a counterpart in a restart take the consensus
    profile, so every member is a complete profile set in the shared index."""
    consensus_norm = _normalize_profiles(consensus)
    aligned_members = []
    for h in pool:
        h = h[list(consensus.columns)]
        match = _match_factors(_normalize_profiles(h), consensus_norm)
        aligned = consensus.copy()
        for i in range(h.shape[0]):
            j, r = match["match"][i], match["correlation"][i]
            if j >= 0 and np.isfinite(r) and r >= match_min:
                aligned.iloc[j] = h.iloc[i].to_numpy(dtype=float)
        aligned_members.append(aligned)
    return aligned_members


def _training_loss(run_path: str, profiles: pd.DataFrame, num_threads: int = 1) -> float:
    return float(_core.training_loss(
        run_path, profiles.to_numpy(dtype=float).ravel(order="C").tolist(),
        int(profiles.shape[0]), [str(g) for g in profiles.columns], int(num_threads)))


_FORKED_JOB = None  # the function a forked worker inherits and applies


def _forked_job(item):
    return _FORKED_JOB(item)


def _run_jobs(fn, items, workers: int) -> list:
    """Apply ``fn`` to every item in forked worker processes (``workers`` at a
    time; in this process when ``workers`` is 1) and return the results in
    order. An error in any job is raised here."""
    global _FORKED_JOB

    items = list(items)
    if not items:
        return []
    if workers <= 1:
        return [fn(item) for item in items]
    import multiprocessing

    context = multiprocessing.get_context("fork")
    _FORKED_JOB = fn
    try:
        with context.Pool(workers, maxtasksperchild=1) as pool:
            return pool.map(_forked_job, items, chunksize=1)
    finally:
        _FORKED_JOB = None


def _parallel_workers(num_threads: Optional[int], per_job: int, n_jobs: int) -> int:
    """Number of concurrent worker processes for independent jobs that each
    use ``per_job`` threads, within a total thread budget: ``num_threads``, or
    by default the ``CELLADMIX_NUM_THREADS`` environment variable, or else
    every available core. Jobs run sequentially on Windows, which cannot
    fork."""
    import os
    import sys

    if sys.platform.startswith("win"):
        return 1
    if num_threads is not None:
        budget = int(num_threads)
    else:
        env = os.environ.get("CELLADMIX_NUM_THREADS", "")
        budget = int(env) if env.strip().isdigit() and int(env) > 0 else 0
        if budget < 1:
            try:
                budget = len(os.sched_getaffinity(0))
            except (AttributeError, OSError):
                budget = os.cpu_count() or 1
    budget = max(1, budget)
    return max(1, min(n_jobs, budget // max(1, per_job)))


def _rules_in_sample(rules: pd.DataFrame, member, carried) -> dict:
    """The rule rows applied in one sample under ``apply_rules="where_present"``:
    those whose program is among ``carried``, the factors the sample's own
    restarts carry; ``member`` (one entry per row, or None) is subset
    alongside."""
    rows = rules["factor"].astype(int).isin({int(f) for f in carried}).to_numpy(dtype=bool)
    subset = rules[rows]
    if member is None:
        return {"rules": subset, "member": None}
    return {"rules": subset, "member": [m for m, keep in zip(list(member), rows) if keep]}


def _interleave(lists: list) -> list:
    """Round-robin interleaving of several lists (one per sample) into one."""
    lists = [items for items in lists if items]
    out = []
    for i in range(max((len(items) for items in lists), default=0)):
        for items in lists:
            if i < len(items):
                out.append(items[i])
    return out


def _profiles_hash(profiles: pd.DataFrame, members: list) -> str:
    """Content hash of shared profiles together with their member pool, so a
    cached fixed-profile run is reused only for the same profiles and members."""
    stacked = pd.concat([profiles] + [m[list(profiles.columns)] for m in members]) if members else profiles
    return _matrix_hash(stacked)


def _matrix_hash(frame: pd.DataFrame) -> str:
    # Rounded to five decimals first, so that floating-point noise between two
    # computations of the same consensus does not change the id.
    values = np.round(frame.to_numpy(dtype=float), 5)
    return hashlib.sha1(values.tobytes()).hexdigest()[:8]


def _pool_score_summaries(summaries: pd.DataFrame) -> pd.DataFrame:
    """Pool per-sample summaries into one per target/source/factor: pair
    counts summed, mean scores pair-weighted, one-sided p-values combined by
    Stouffer's method with weights proportional to sqrt(pairs)."""
    from scipy.stats import norm

    columns = ["target_cell_type", "source_cell_type", "factor", "n_pairs", "mean_score",
               "p_value", "neg_log10_p", "n_samples"]
    if summaries is None or summaries.empty:
        return pd.DataFrame(columns=columns)
    has_null = "mean_null_score" in summaries.columns
    rows = []
    for (target, source, factor), group in summaries.groupby(["target_cell_type", "source_cell_type", "factor"], sort=True):
        n = np.array(group["n_pairs"].astype(float).to_numpy(), dtype=float, copy=True)
        n[~np.isfinite(n) | (n <= 0)] = 1.0
        p = np.clip(group["p_value"].astype(float).to_numpy(), 1e-300, 1 - 1e-16)
        z = norm.isf(p)
        z_pooled = float(np.sum(np.sqrt(n) * z) / np.sqrt(np.sum(n)))
        p_pooled = float(norm.sf(z_pooled))
        row = {"target_cell_type": target, "source_cell_type": source, "factor": int(factor),
               "n_pairs": float(group["n_pairs"].astype(float).sum()),
               "mean_score": float(np.sum(n * group["mean_score"].astype(float).to_numpy()) / np.sum(n))}
        if has_null:
            row["mean_null_score"] = float(np.sum(n * group["mean_null_score"].astype(float).to_numpy()) / np.sum(n))
        row["p_value"] = p_pooled
        row["neg_log10_p"] = float(-np.log10(max(p_pooled, 1e-300)))
        row["n_samples"] = int(len(group))
        rows.append(row)
    return pd.DataFrame(rows)


def _bind_sample_counts(counts_by_sample: dict):
    """Column-bind per-sample (matrix, genes, cells) count triples with
    sample-prefixed cell ids, aligning gene rows by name."""
    from scipy import sparse

    first_genes = None
    blocks, cells = [], []
    for s, (matrix, genes, cell_ids) in counts_by_sample.items():
        genes = pd.Index([str(g) for g in genes])
        if first_genes is None:
            first_genes = genes
        if set(genes) != set(first_genes):
            raise ValueError("samples do not share one gene panel")
        order = genes.get_indexer(first_genes)
        blocks.append(sparse.csc_matrix(matrix)[order, :])
        cells.extend(f"{s}_{c}" for c in cell_ids)
    return sparse.hstack(blocks).tocsc(), pd.Index(first_genes, name="gene"), pd.Index(cells, name="cell_id")
