"""Score result objects."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ._score_utils import rules_from_annotation, score_annotation


class CellAdmixScore:
    """Membrane or bridge scoring result."""

    def __init__(self, fit, method: str, result: dict):
        self.fit = fit
        self.method = method
        self.result = result
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
    ) -> pd.DataFrame:
        """Return factor/target correction rules using the strongest source per factor."""
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
        return rules

    def correct(self, rules: pd.DataFrame | None = None, *, p_thresh: float = 0.1, name: str | None = None):
        """Apply score-derived correction rules to the parent fit."""
        rules = self.rules(p_thresh=p_thresh) if rules is None else rules
        return self.fit.correct(rules, name=name or f"{self.method}_clean")

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
