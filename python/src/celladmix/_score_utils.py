"""Shared score annotation logic used by rules and plotting."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd


def score_annotation(summary: pd.DataFrame, *, p_thresh: float = 0.1, adjust_p: bool = False) -> dict:
    """Summarize score evidence into target matrix, source calls, and removal calls."""
    required = {"target_cell_type", "source_cell_type", "factor", "p_value", "neg_log10_p"}
    missing = required.difference(summary.columns)
    if missing:
        raise ValueError(f"score summary is missing required columns: {', '.join(sorted(missing))}")
    data = summary.dropna(subset=["target_cell_type", "source_cell_type", "factor"]).copy()
    if data.empty:
        return {
            "score_matrix": pd.DataFrame(),
            "source_score": pd.DataFrame(),
            "source_calls": {},
            "remove_calls": [],
            "threshold": -math.log10(p_thresh),
            "summary": data,
        }
    data["factor"] = data["factor"].astype(int)
    data["plot_p_value"] = data["p_value"].astype(float)
    if adjust_p:
        # Adjustment is applied over the provided summary rows, matching the
        # plotting/rule context rather than a global experiment-wide universe.
        data["plot_p_value"] = _benjamini_hochberg(data["plot_p_value"])
        data["plot_score"] = -np.log10(data["plot_p_value"])
    else:
        data["plot_score"] = data["neg_log10_p"].astype(float)
    data.loc[np.isnan(data["plot_score"]), "plot_score"] = np.nan

    cell_types = sorted(set(data["target_cell_type"]) | set(data["source_cell_type"]))
    factors = sorted(data["factor"].unique())
    columns = [f"F{factor}" for factor in factors]
    score_matrix = pd.DataFrame(np.nan, index=cell_types, columns=columns)
    source_score = pd.DataFrame(np.nan, index=cell_types, columns=columns)

    for cell_type in cell_types:
        for factor in factors:
            column = f"F{factor}"
            target_values = data.loc[
                (data["target_cell_type"] == cell_type) & (data["factor"] == factor),
                "plot_score",
            ]
            target_values = target_values[np.isfinite(target_values)]
            if len(target_values):
                score_matrix.loc[cell_type, column] = float(target_values.max())

            source_values = data.loc[
                (data["source_cell_type"] == cell_type) & (data["factor"] == factor),
                "plot_score",
            ]
            source_values = source_values[np.isfinite(source_values)]
            if len(source_values):
                source_score.loc[cell_type, column] = float(source_values.mean())

    threshold = -math.log10(p_thresh)
    source_calls: dict[int, str] = {}
    remove_calls: list[tuple[int, str]] = []
    for factor in factors:
        column = f"F{factor}"
        # Source calls come from source evidence; removal calls come from
        # significant target evidence after excluding the inferred source type.
        src = source_score[column].dropna()
        if not src.empty:
            source_calls[int(factor)] = str(src.idxmax())
        targets = score_matrix.index[
            np.isfinite(score_matrix[column].to_numpy())
            & (score_matrix[column].to_numpy() > threshold)
        ].tolist()
        source = source_calls.get(int(factor))
        for target in targets:
            if source is None or target != source:
                remove_calls.append((int(factor), str(target)))

    return {
        "score_matrix": score_matrix,
        "source_score": source_score,
        "source_calls": source_calls,
        "remove_calls": remove_calls,
        "threshold": threshold,
        "summary": data,
    }


def rules_from_annotation(annotation: dict, *, target_cell_types=None) -> pd.DataFrame:
    """Convert score annotation output into correction rules."""
    summary = annotation["summary"]
    rows = []
    target_set = None if target_cell_types is None else set(map(str, target_cell_types))
    for factor, target in annotation["remove_calls"]:
        if target_set is not None and target not in target_set:
            continue
        candidates = summary[
            (summary["factor"].astype(int) == int(factor))
            & (summary["target_cell_type"].astype(str) == str(target))
        ]
        if candidates.empty:
            p_value = np.nan
            neg_log10_p = np.nan
        elif "plot_p_value" in candidates.columns:
            best = candidates["plot_p_value"].astype(float).idxmin()
            p_value = float(candidates.loc[best, "p_value"])
            neg_log10_p = float(candidates.loc[best, "neg_log10_p"])
        else:
            best = candidates["p_value"].astype(float).idxmin()
            p_value = float(candidates.loc[best, "p_value"])
            neg_log10_p = float(candidates.loc[best, "neg_log10_p"])
        rows.append(
            {
                "factor": int(factor),
                "source_cell_type": annotation["source_calls"].get(int(factor), pd.NA),
                "target_cell_type": target,
                "p_value": p_value,
                "neg_log10_p": neg_log10_p,
                "rule_id": f"{int(factor)}_{target}",
            }
        )
    return pd.DataFrame(
        rows,
        columns=["factor", "source_cell_type", "target_cell_type", "p_value", "neg_log10_p", "rule_id"],
    )


def _benjamini_hochberg(values: pd.Series) -> np.ndarray:
    p = values.astype(float).to_numpy()
    out = np.full_like(p, np.nan, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return out
    idx = np.where(finite)[0]
    order = idx[np.argsort(p[finite])]
    ranked = p[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out[order] = np.minimum(ranked, 1.0)
    return out

def source_exposure_counts(cells: pd.DataFrame, annotation, *, neighbor_k: int = 15):
    """Count annotated cell types among each cell's k nearest neighbor cells.

    Returns ``(counts, types)`` where counts is an ``len(cells) x len(types)``
    array aligned to the rows of ``cells``. Unannotated cells occupy space as
    neighbors but contribute no counts.
    """
    from . import _core

    types = sorted(map(str, pd.Series(annotation).dropna().unique()))
    type_index = {t: i for i, t in enumerate(types)}
    codes = [
        type_index.get(str(annotation.get(str(c), None)), -1)
        for c in cells["cell_id"].astype(str)
    ]
    raw = _core.cell_neighbor_type_counts(
        cells["x"].astype(float).tolist(),
        cells["y"].astype(float).tolist(),
        codes,
        len(types),
        int(neighbor_k),
    )
    counts = np.asarray(raw["data"], dtype=float).reshape(raw["rows"], raw["cols"])
    return counts, types


def apply_native_check(
    rules: pd.DataFrame,
    fit,
    *,
    median_thresh: float = 0.1,
    expr_thresh: float = 0.05,
    outlier_min_frac: float = 0.1,
    neighbor_k: int = 15,
) -> pd.DataFrame:
    """Flag likely-native factor/target rules using source-distant cells.

    For each rule (factor f, target T, source S), target cells with zero
    S-type cells among their nearest neighbors form the source-distant group.
    The rule is flagged ``keep=False`` when the factor persists in that group
    (median fraction above ``median_thresh``), when the distant expression is
    a cross-type outlier, when there is no positive exposure gradient, or when
    the check cannot be evaluated safely.
    """
    rules = rules.copy()
    checks = {
        "keep": [], "native_check": [], "native_distant_n": [],
        "native_exposed_n": [], "native_distant_median": [],
        "native_exposure_gradient": [], "native_distant_expr_frac": [],
    }
    if rules.empty:
        for key, values in checks.items():
            rules[key] = values
        return rules

    annotation = getattr(getattr(fit, "dataset", None), "annotation", None)
    cells = fit.cell_factors()
    if annotation is None or cells.empty:
        raise ValueError("native_check requires an annotation and cell factors")
    counts, types = source_exposure_counts(cells, annotation, neighbor_k=neighbor_k)
    cell_types = np.array(
        [str(annotation.get(str(c), None)) for c in cells["cell_id"].astype(str)]
    )

    def factor_column(f):
        return cells[f"factor_{int(f)}_fraction"].to_numpy(dtype=float)

    for _, rule in rules.iterrows():
        target = str(rule["target_cell_type"])
        source = rule["source_cell_type"]
        record = dict(keep=True, native_check="pass", native_distant_n=np.nan,
                      native_exposed_n=np.nan, native_distant_median=np.nan,
                      native_exposure_gradient=np.nan, native_distant_expr_frac=np.nan)

        def flag(reason):
            record["keep"] = False
            record["native_check"] = reason

        if pd.isna(source) or str(source) not in types:
            flag("source_not_in_annotation")
        elif target not in types:
            flag("target_not_in_annotation")
        else:
            source_col = types.index(str(source))
            fractions = factor_column(rule["factor"])
            distant_mask = counts[:, source_col] == 0
            in_target = cell_types == target
            distant = in_target & distant_mask
            exposed = in_target & ~distant_mask
            record["native_distant_n"] = int(distant.sum())
            record["native_exposed_n"] = int(exposed.sum())
            if distant.sum() <= 1:
                flag("no_distant_cells")
            elif exposed.sum() <= 1:
                flag("no_exposed_cells")
            else:
                fr_d = fractions[distant]
                fr_e = fractions[exposed]
                record["native_distant_median"] = float(np.median(fr_d))
                record["native_exposure_gradient"] = float(fr_e.mean() - fr_d.mean())
                record["native_distant_expr_frac"] = float((fr_d > expr_thresh).mean())
                if np.all(fr_d == 0) and np.all(fr_e == 0):
                    flag("no_expression")
                else:
                    # Cross-type outlier test on source-distant expression rates.
                    rates = []
                    for cell_type in types:
                        if cell_type == str(source):
                            continue
                        members = (cell_types == cell_type) & distant_mask
                        rates.append(
                            float((fractions[members] > expr_thresh).mean())
                            if members.sum() > 1 else 0.0
                        )
                    rates = np.asarray(rates)
                    q3 = np.quantile(rates, 0.75)
                    upper = q3 + 1.5 * (q3 - np.quantile(rates, 0.25))
                    target_rate = record["native_distant_expr_frac"]
                    is_outlier = target_rate > upper and target_rate > outlier_min_frac
                    if record["native_distant_median"] > median_thresh:
                        flag("native_median")
                    elif is_outlier:
                        flag("native_outlier")
                    elif record["native_exposure_gradient"] < 0:
                        flag("no_gradient")
        for key in checks:
            checks[key].append(record[key])

    for key, values in checks.items():
        rules[key] = values
    return rules
