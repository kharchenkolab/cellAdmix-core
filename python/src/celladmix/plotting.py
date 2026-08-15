"""Lightweight plotting helpers for Python notebooks."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from ._score_utils import score_annotation
from .factor_sources import factor_source_annotation


def plot_loadings(loadings: pd.DataFrame, *, n_genes: int = 8, ncol: int = 3, ax=None):
    import matplotlib.pyplot as plt

    factors = loadings.columns.tolist()
    ncol = max(1, int(ncol))
    nrow = math.ceil(len(factors) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.4, max(2.0, 1.75 * nrow)), squeeze=False)
    for ax_i, factor in zip(axes.ravel(), factors):
        top = loadings[factor].sort_values(ascending=False).head(n_genes).iloc[::-1]
        ax_i.barh(top.index, top.values, color="black")
        ax_i.set_title(str(factor), fontsize=9)
        ax_i.tick_params(axis="y", labelsize=8)
        ax_i.tick_params(axis="x", labelsize=7)
    for ax_i in axes.ravel()[len(factors) :]:
        ax_i.axis("off")
    fig.tight_layout()
    return fig


def plot_stability(fit, *, min_stability: float = 0.3, label: bool = True, ax=None):
    """Plot NMF restart stability versus factor molecule importance."""
    import matplotlib.pyplot as plt

    diagnostics = fit.manifest.get("nmf_diagnostics", {}) or {}
    n_runs = int((fit.manifest.get("pipeline_options", {}) or {}).get("nmf_n_runs", 1))
    stability = np.asarray(diagnostics.get("selected_factor_stability", []), dtype=float)
    n_factors = int(fit.manifest.get("n_factors") or len(stability))
    if len(stability) < n_factors:
        stability = np.concatenate([stability, np.full(n_factors - len(stability), np.nan)])
    stability = stability[:n_factors]
    cells = fit.cell_factors()
    factor_cols = [f"factor_{i}_fraction" for i in range(1, n_factors + 1)]
    if all(col in cells.columns for col in factor_cols) and "transcript_count" in cells.columns:
        fractions = cells[factor_cols].to_numpy(dtype=float)
        weights = cells["transcript_count"].to_numpy(dtype=float)
        # Weight cell-level fractions by transcript counts to approximate the
        # fraction of assigned molecules represented by each factor.
        importance = np.nansum(fractions * weights[:, None], axis=0)
        importance = importance / max(np.nansum(importance), np.finfo(float).eps)
        x_label = "Factor importance (% of assigned molecules)"
    else:
        loading = fit.factor_loadings().sum(axis=0).to_numpy(dtype=float)
        importance = loading / max(np.nansum(loading), np.finfo(float).eps)
        x_label = "Factor importance (% of loading mass)"
    if ax is None:
        _, ax = plt.subplots(figsize=(6.2, 4.2))
    if n_runs <= 1 or not np.isfinite(stability).any():
        ax.text(0.5, 0.5, "NMF restart stability is unavailable", ha="center", va="center")
        ax.axis("off")
        return ax
    x = 100 * np.maximum(importance, 0)
    y = np.clip(stability, -1, 1)
    colors = np.where(y >= 0.6, "#1B9E77", np.where(y >= min_stability, "#7570B3", "#D95F02"))
    max_importance = max(np.nanmax(importance), np.finfo(float).eps)
    sizes = 60 + 240 * np.sqrt(np.maximum(importance, 0) / max_importance)
    ax.axhline(min_stability, linestyle="--", color="#8c8c8c", linewidth=0.8)
    ax.axhline(0.6, linestyle=":", color="#b3b3b3", linewidth=0.8)
    scatter = ax.scatter(x, y, s=sizes, c=colors, alpha=0.9)
    # Dot size encodes factor importance; invert the size formula for the key.
    handles, size_labels = scatter.legend_elements(
        prop="sizes",
        num=3,
        func=lambda s: 100 * max_importance * np.square(np.maximum(s - 60, 0) / 240),
        fmt="{x:.1f}%",
    )
    for handle in handles:
        handle.set_markerfacecolor("#8c8c8c")
        handle.set_markeredgecolor("none")
    size_legend = ax.legend(
        handles, size_labels, loc="lower right", frameon=True,
        fontsize=6.5, title="Importance", title_fontsize=7,
        labelspacing=1.1, borderpad=0.9)
    size_legend.get_frame().set_facecolor("white")
    size_legend.get_frame().set_alpha(0.55)
    size_legend.get_frame().set_linewidth(0)
    if label:
        for i, (xi, yi) in enumerate(zip(x, y), start=1):
            if np.isfinite(xi) and np.isfinite(yi):
                ax.text(xi, yi + 0.015, f"F{i}", fontsize=8, ha="center")
    ax.set_title("NMF factor stability vs. importance")
    ax.set_xlabel(x_label)
    ax.set_ylabel("Matched ownership correlation across restarts")
    ax.set_ylim(max(-0.1, np.nanmin(y) - 0.05), 1.02)
    ax.figure.tight_layout()
    return ax


def plot_spatial(data: pd.DataFrame, *, color_by: str, x: str = "x", y: str = "y", title=None, s: float = 1.0, ax=None):
    """Plot cell-level spatial coordinates colored by a column."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(6.0, 4.8))
    frame = data.dropna(subset=[x, y, color_by]).copy()
    values = frame[color_by]
    if pd.api.types.is_numeric_dtype(values) and values.nunique(dropna=True) > 12:
        sc = ax.scatter(frame[x], frame[y], c=values, s=s, cmap="viridis", linewidths=0, alpha=0.85)
        ax.figure.colorbar(sc, ax=ax, fraction=0.046, pad=0.03, label=color_by)
    else:
        labels = sorted(values.astype(str).unique())
        cmap = plt.get_cmap("tab20")
        for i, label_i in enumerate(labels):
            sub = frame[values.astype(str) == label_i]
            ax.scatter(sub[x], sub[y], s=s, color=cmap(i % 20), linewidths=0, alpha=0.85, label=label_i)
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=7, markerscale=4)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title or color_by)
    ax.figure.tight_layout()
    return ax


def plot_molecule_removal(summary: pd.DataFrame):
    """Plot per-cell molecule removal counts and fractions."""
    import matplotlib.pyplot as plt

    if summary.empty:
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, "No correction summary available", ha="center", va="center")
        ax.axis("off")
        return fig
    data = summary.copy()
    data["cell_type"] = data["cell_type"].fillna("unknown").astype(str)
    order = (
        data.groupby("cell_type")["n_removed"].median().sort_values(ascending=False).index.tolist()
        if "n_removed" in data.columns
        else sorted(data["cell_type"].unique())
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), squeeze=False)
    for ax, col, ylabel in [
        (axes[0, 0], "n_removed", "Molecules removed per cell"),
        (axes[0, 1], "fraction_removed", "Fraction removed per cell"),
    ]:
        series = [data.loc[data["cell_type"] == label, col].astype(float).to_numpy() for label in order]
        ax.boxplot(series, labels=order, showfliers=False)
        ax.set_ylabel(ylabel)
        ax.tick_params(axis="x", rotation=45, labelsize=8)
    fig.tight_layout()
    return fig


def plot_score_heatmap(
    summary: pd.DataFrame,
    *,
    p_thresh: float = 0.1,
    adjust_p: bool = False,
    max_neg_log10: float = 6,
    source_prior=None,
    source_prior_min_score: float = 0.15,
    source_prior_min_margin: float = 0.03,
    ax=None,
):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    if ax is None:
        _, ax = plt.subplots(figsize=(6.8, 4.4))
    ann = score_annotation(summary, p_thresh=p_thresh, adjust_p=adjust_p)
    prior = factor_source_annotation(
        source_prior,
        min_score=source_prior_min_score,
        min_margin=source_prior_min_margin,
    )
    gene_sources = {}
    if prior is not None and not prior.empty:
        called = prior[
            prior["called"]
            & prior["source_cell_type"].notna()
            & (prior["source_cell_type"].astype(str) != "")
        ]
        gene_sources = dict(zip(called["factor"].astype(int), called["source_cell_type"].astype(str)))
    matrix = ann["score_matrix"]
    if matrix.empty:
        ax.text(0.5, 0.5, "No score summaries to plot", ha="center", va="center")
        ax.axis("off")
        return ax

    values = matrix.to_numpy(dtype=float)
    plot_values = np.array(values, copy=True)
    plot_values[~np.isfinite(plot_values)] = np.nan
    plot_values = np.clip(plot_values, 0, max_neg_log10)
    # Match the R display: evidence below threshold stays white.
    plot_values[plot_values <= ann["threshold"]] = 0
    threshold = max(ann["threshold"], np.finfo(float).eps)
    denom = max(max_neg_log10 - threshold, np.finfo(float).eps)
    scaled = np.zeros_like(plot_values)
    finite = np.isfinite(plot_values) & (plot_values > threshold)
    scaled[finite] = (plot_values[finite] - threshold) / denom
    cmap = LinearSegmentedColormap.from_list("celladmix_white_red", ["white", "#B2182B"])
    image = ax.imshow(np.nan_to_num(scaled, nan=0.0), aspect="auto", cmap=cmap, vmin=0, vmax=1)

    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=0)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index)
    ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="#BDBDBD", linestyle="-", linewidth=0.7)
    ax.tick_params(which="minor", bottom=False, left=False)

    remove_calls = set(ann["remove_calls"])
    for i, cell_type in enumerate(matrix.index):
        for j, column in enumerate(matrix.columns):
            factor = int(str(column).lstrip("F"))
            spatial_source = ann["source_calls"].get(factor)
            gene_source = gene_sources.get(factor)
            label = ""
            if spatial_source == cell_type and gene_source == cell_type:
                label = "SG"
            elif spatial_source == cell_type:
                label = "S"
            elif gene_source == cell_type:
                label = "G"
            if (factor, cell_type) in remove_calls:
                label = f"{label}*"
            if label:
                ax.text(
                    j,
                    i,
                    label,
                    ha="center",
                    va="center",
                    fontsize=8.5 if len(label) > 1 else 10.5,
                    fontweight="bold",
                )

    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("-log10(p)")
    ticks = [0, 1]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{ann['threshold']:.1f}", f"{max_neg_log10:.1f}"])
    ax.set_xlabel("Factor")
    ax.set_ylabel("")
    ax.set_title(f"Factor annotation, p < {p_thresh:g}")
    if prior is not None:
        ax.text(
            0.0,
            -0.16,
            "S spatial source; G gene-content source; SG agreement; * target",
            transform=ax.transAxes,
            fontsize=7,
            ha="left",
            va="top",
        )
    ax.figure.tight_layout()
    return ax


def plot_factor_source_heatmap(
    source_score,
    *,
    min_score: float = 0.15,
    min_margin: float = 0.03,
    max_score: float | None = None,
    ax=None,
):
    """Plot marker-weighted gene-content source evidence."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    if ax is None:
        _, ax = plt.subplots(figsize=(6.8, 4.4))
    matrix = source_score.score_matrix.T.copy()
    if matrix.empty:
        ax.text(0.5, 0.5, "No factor source scores to plot", ha="center", va="center")
        ax.axis("off")
        return ax
    values = matrix.to_numpy(dtype=float)
    max_score = float(max_score) if max_score is not None else float(np.nanmax(values))
    max_score = max(max_score, np.finfo(float).eps)
    scaled = np.clip(values / max_score, 0, 1)
    cmap = LinearSegmentedColormap.from_list("celladmix_source_blue", ["white", "#2166AC"])
    image = ax.imshow(scaled, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=0)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index)
    ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="#BDBDBD", linestyle="-", linewidth=0.7)
    ax.tick_params(which="minor", bottom=False, left=False)

    calls = source_score.annotation(min_score=min_score, min_margin=min_margin)
    for _, row in calls.iterrows():
        if not bool(row["called"]) or pd.isna(row["source_cell_type"]):
            continue
        if row["source_cell_type"] not in matrix.index or row["factor_label"] not in matrix.columns:
            continue
        ax.text(
            matrix.columns.get_loc(row["factor_label"]),
            matrix.index.get_loc(row["source_cell_type"]),
            "S",
            ha="center",
            va="center",
            fontsize=10.5,
            fontweight="bold",
        )
    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("source score")
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(["0", f"{max_score:.2f}"])
    ax.set_xlabel("Factor")
    ax.set_ylabel("")
    ax.set_title("Factor source gene-content score")
    ax.figure.tight_layout()
    return ax


def plot_factor_source_top_genes(source_score, *, n_genes: int = 6, ncol: int = 3):
    """Plot marker-weighted genes supporting each factor source call."""
    import matplotlib.pyplot as plt

    data = source_score.top_factor_genes(n_genes=n_genes)
    if data.empty:
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, "No marker-weighted factor genes available", ha="center", va="center")
        ax.axis("off")
        return fig
    panels = data[["factor_label", "source_cell_type"]].drop_duplicates().reset_index(drop=True)
    ncol = max(1, int(ncol))
    nrow = math.ceil(len(panels) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.4, max(2.0, 1.9 * nrow)), squeeze=False)
    for ax_i, row in zip(axes.ravel(), panels.itertuples(index=False)):
        sub = data[
            (data["factor_label"] == row.factor_label)
            & (data["source_cell_type"] == row.source_cell_type)
        ].sort_values("rank", ascending=False)
        ax_i.barh(sub["gene"], sub["contribution"], color="black")
        ax_i.set_title(f"{row.factor_label}\\n{row.source_cell_type}", fontsize=8.5)
        ax_i.tick_params(axis="y", labelsize=7.5)
        ax_i.tick_params(axis="x", labelsize=7)
    for ax_i in axes.ravel()[len(panels) :]:
        ax_i.axis("off")
    fig.tight_layout()
    return fig


def plot_score_pairs(
    summary: pd.DataFrame,
    *,
    factors,
    p_thresh: float = 0.1,
    adjust_p: bool = False,
    max_abs_significance: float = 6,
    ncol: int = 2,
):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

    ann = score_annotation(summary, p_thresh=p_thresh, adjust_p=adjust_p)
    data = ann["summary"]
    cell_types = list(ann["score_matrix"].index)
    ncol = max(1, int(ncol))
    nrow = math.ceil(len(factors) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(9.0, max(3.8, 3.6 * nrow)), squeeze=False)
    sig_cmap = LinearSegmentedColormap.from_list("celladmix_signed_sig", ["#2166AC", "white", "#B2182B"])
    src_cmap = LinearSegmentedColormap.from_list("celladmix_source", ["white", "#525252"])
    for ax, factor in zip(axes.ravel(), factors):
        factor = int(factor)
        sub = data[data["factor"].astype(int) == factor].copy()
        score_col = "mean_score" if "mean_score" in sub.columns else "neg_log10_p"
        sub["score_value"] = sub[score_col].astype(float)
        sub["significance_value"] = sub["plot_score"].astype(float)
        # Pair plots encode direction with sign and confidence with magnitude.
        sub["signed_significance"] = np.sign(sub["score_value"]) * sub["significance_value"]
        pivot = sub.pivot_table(
            index="target_cell_type",
            columns="source_cell_type",
            values="signed_significance",
            aggfunc="mean",
        ).reindex(index=list(reversed(cell_types)), columns=cell_types)

        values = pivot.to_numpy(dtype=float)
        source = ann["source_score"].get(f"F{factor}", pd.Series(index=cell_types, dtype=float)).reindex(cell_types)
        combined = np.vstack([source.to_numpy(dtype=float), values])
        combined = np.clip(combined, -max_abs_significance, max_abs_significance)

        norm = TwoSlopeNorm(vcenter=0, vmin=-max_abs_significance, vmax=max_abs_significance)
        im = ax.imshow(combined, cmap=sig_cmap, norm=norm, aspect="equal")
        source_vals = source.to_numpy(dtype=float)
        finite_src = source_vals[np.isfinite(source_vals)]
        if len(finite_src):
            lo, hi = finite_src.min(), finite_src.max()
            scaled_src = (source_vals - lo) / max(hi - lo, np.finfo(float).eps)
            for j, val in enumerate(scaled_src):
                color = src_cmap(0 if not np.isfinite(val) else val)
                ax.add_patch(plt.Rectangle((j - 0.5, -0.5), 1, 1, facecolor=color, edgecolor="white", linewidth=0.6))
            best = np.nanargmax(source_vals) if np.isfinite(source_vals).any() else None
            if best is not None:
                ax.add_patch(plt.Rectangle((best - 0.5, -0.5), 1, 1, fill=False, edgecolor="black", linewidth=1.1))
        for j, val in enumerate(source_vals):
            if np.isfinite(val):
                ax.text(j, 0, f"{val:.1f}", ha="center", va="center", fontsize=6.5)

        ax.set_xticks(np.arange(len(cell_types)))
        ax.set_xticklabels(cell_types, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(np.arange(len(cell_types) + 1))
        ax.set_yticklabels(["source evidence"] + list(reversed(cell_types)), fontsize=7)
        ax.set_xticks(np.arange(-0.5, len(cell_types), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(cell_types) + 1, 1), minor=True)
        ax.grid(which="minor", color="white", linestyle="-", linewidth=0.8)
        ax.tick_params(which="minor", bottom=False, left=False)
        ax.set_title(f"F{factor}", fontsize=10)
        ax.set_xlabel("Source cell type", fontsize=8)
        ax.set_ylabel("Target cell type", fontsize=8)
        ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.03, label="Signed -log10(p)")
    for ax in axes.ravel()[len(factors) :]:
        ax.axis("off")
    fig.tight_layout()
    return fig
