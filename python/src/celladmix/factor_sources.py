"""Gene-content source priors for fitted cellAdmix factors."""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import _core


class CellAdmixFactorSourceScore:
    """Marker-weighted factor-to-cell-type source evidence."""

    def __init__(self, fit, result: dict):
        self.fit = fit
        self.result = result
        self.summary = result["scores"]
        self.score_matrix = result["score_matrix"]
        self.marker_weights = result["marker_weights"]
        self.top_markers = result["top_markers"]

    def __repr__(self) -> str:
        return (
            "CellAdmixFactorSourceScore("
            f"n_factors={self.summary['factor'].nunique()}, "
            f"n_cell_types={self.summary['cell_type'].nunique()})"
        )

    def annotation(self, *, min_score: float = 0.15, min_margin: float = 0.03) -> pd.DataFrame:
        """Return the best gene-content source call for each factor."""
        best = self.summary[self.summary["is_best"]].copy()
        best["called"] = (best["score"] >= float(min_score)) & (best["margin"] >= float(min_margin))
        best["source_cell_type"] = np.where(best["called"], best["cell_type"], pd.NA)
        return best[
            ["factor", "factor_label", "source_cell_type", "score", "margin", "called"]
        ].reset_index(drop=True)

    def top_factor_genes(self, *, n_genes: int = 8, sources=None) -> pd.DataFrame:
        """Return marker-weighted genes supporting each factor source call."""
        return factor_source_top_genes(self.result, n_genes=n_genes, sources=sources)

    def plot_heatmap(self, **kwargs):
        from .plotting import plot_factor_source_heatmap

        return plot_factor_source_heatmap(self, **kwargs)

    def plot_top_genes(self, **kwargs):
        from .plotting import plot_factor_source_top_genes

        return plot_factor_source_top_genes(self, **kwargs)


def score_factor_sources(
    fit,
    *,
    annotation=None,
    counts=None,
    normalize: str = "log1p_cpm",
    scale_factor: float = 10000.0,
    specificity_power: float = 1.0,
    factor_power: float = 1.0,
    min_marker_logfc: float = 0.0,
    top_markers_per_type: int | None = None,
    top_genes: int = 25,
) -> CellAdmixFactorSourceScore:
    """Score each NMF factor against cell-type-specific marker signatures."""
    from scipy import sparse

    labels = fit.dataset.annotation if annotation is None else annotation
    if labels is None:
        raise ValueError("score_factor_sources() requires a cell annotation")
    labels = pd.Series(labels).dropna().astype(str)
    labels.index = labels.index.astype(str)
    labels = labels[labels != ""]
    if labels.empty:
        raise ValueError("annotation does not contain usable labels")

    if counts is None:
        counts, genes, cells = fit.counts()
    else:
        counts, genes, cells = counts
    genes = pd.Index(genes.astype(str) if hasattr(genes, "astype") else [str(x) for x in genes], name="gene")
    cells = pd.Index(cells.astype(str) if hasattr(cells, "astype") else [str(x) for x in cells], name="cell_id")
    counts = counts.tocsc() if sparse.issparse(counts) else sparse.csc_matrix(counts)

    loadings = fit.factor_loadings()
    loadings.index = loadings.index.astype(str)
    loadings.columns = [str(col) if str(col).startswith("F") else f"F{int(col)}" for col in loadings.columns]

    gene_pos = pd.Series(np.arange(len(genes)), index=genes)
    common_genes = [gene for gene in loadings.index if gene in gene_pos.index]
    if not common_genes:
        raise ValueError("No overlap between factor-loading genes and count-matrix genes")
    cell_pos = pd.Series(np.arange(len(cells)), index=cells)
    common_cells = labels.index.intersection(cell_pos.index)
    if common_cells.empty:
        raise ValueError("No overlap between count-matrix cells and annotation cells")

    gene_idx = gene_pos.loc[common_genes].to_numpy(dtype=int)
    cell_idx = cell_pos.loc[common_cells].to_numpy(dtype=int)
    labels = labels.loc[common_cells].astype(str)
    matrix = counts[gene_idx, :][:, cell_idx]

    cell_types = sorted(labels.unique().tolist())
    if len(cell_types) < 2:
        raise ValueError("Factor source scoring requires at least two annotated cell types")
    type_pos = {label: i for i, label in enumerate(cell_types)}
    type_idx = np.array([type_pos[x] for x in labels], dtype=int)
    design = sparse.csr_matrix(
        (np.ones(len(type_idx)), (np.arange(len(type_idx)), type_idx)),
        shape=(len(type_idx), len(cell_types)),
    )
    bulk = matrix @ design
    if normalize not in {"log1p_cpm", "cpm", "none"}:
        raise ValueError("normalize must be one of 'log1p_cpm', 'cpm', or 'none'")
    if normalize != "none":
        lib = np.asarray(bulk.sum(axis=0)).ravel()
        scale = float(scale_factor) / np.maximum(lib, 1.0)
        bulk = bulk @ sparse.diags(scale)
    if normalize == "log1p_cpm":
        bulk = bulk.tocsc(copy=True)
        bulk.data = np.log1p(bulk.data)
    profile = np.asarray(bulk.toarray(), dtype=float)

    weights = np.zeros_like(profile, dtype=float)
    for j in range(len(cell_types)):
        if len(cell_types) == 2:
            other = profile[:, 1 - j]
        else:
            other = np.max(np.delete(profile, j, axis=1), axis=1)
        margin = np.maximum(profile[:, j] - other, 0.0)
        margin[margin < float(min_marker_logfc)] = 0.0
        if top_markers_per_type is not None:
            keep_n = max(1, int(top_markers_per_type))
            keep = np.argsort(-margin)[:keep_n]
            mask = np.zeros_like(margin, dtype=bool)
            mask[keep] = True
            margin[~mask] = 0.0
        weights[:, j] = margin
    weights = weights ** float(specificity_power)

    factor_values = np.maximum(loadings.loc[common_genes].to_numpy(dtype=float).T, 0.0)
    factor_values = factor_values ** float(factor_power)
    numerator = factor_values @ weights
    denom = np.sqrt((factor_values**2).sum(axis=1))[:, None] * np.sqrt((weights**2).sum(axis=0))[None, :]
    score_matrix_np = np.divide(numerator, denom, out=np.zeros_like(numerator), where=denom > 0)
    factors = loadings.columns.tolist()
    score_matrix = pd.DataFrame(score_matrix_np, index=factors, columns=cell_types)

    rows = []
    for i, factor_label in enumerate(factors, start=1):
        scores = score_matrix.loc[factor_label]
        order = scores.sort_values(ascending=False)
        second = float(order.iloc[1]) if len(order) > 1 else 0.0
        for rank, (cell_type, score) in enumerate(order.items(), start=1):
            rows.append(
                {
                    "factor": i,
                    "factor_label": factor_label,
                    "cell_type": cell_type,
                    "score": float(score),
                    "rank": rank,
                    "best_score": float(order.iloc[0]),
                    "second_score": second,
                    "margin": float(score) - second,
                    "is_best": rank == 1,
                }
            )
    scores = pd.DataFrame(rows)

    marker_rows = []
    for j, cell_type in enumerate(cell_types):
        values = weights[:, j]
        hit = np.where(values > 0)[0]
        for idx in hit:
            marker_rows.append({"gene": common_genes[idx], "cell_type": cell_type, "weight": float(values[idx])})
    marker_weights = pd.DataFrame(marker_rows, columns=["gene", "cell_type", "weight"])

    top_marker_rows = []
    for j, cell_type in enumerate(cell_types):
        values = weights[:, j]
        order = np.argsort(-values)[: int(top_genes)]
        rank = 1
        for idx in order:
            if values[idx] <= 0:
                continue
            top_marker_rows.append(
                {
                    "cell_type": cell_type,
                    "rank": rank,
                    "gene": common_genes[idx],
                    "marker_weight": float(values[idx]),
                    "expression": float(profile[idx, j]),
                }
            )
            rank += 1
    top_markers = pd.DataFrame(
        top_marker_rows,
        columns=["cell_type", "rank", "gene", "marker_weight", "expression"],
    )

    result = {
        "method": "factor_sources",
        "scores": scores,
        "score_matrix": score_matrix,
        "marker_weights": marker_weights,
        "top_markers": top_markers,
        "profile": pd.DataFrame(profile, index=common_genes, columns=cell_types),
        "factor_loadings": pd.DataFrame(factor_values, index=factors, columns=common_genes),
        "params": {
            "normalize": normalize,
            "scale_factor": float(scale_factor),
            "specificity_power": float(specificity_power),
            "factor_power": float(factor_power),
            "min_marker_logfc": float(min_marker_logfc),
            "top_markers_per_type": top_markers_per_type,
            "n_genes": len(common_genes),
            "n_cells": len(common_cells),
            "n_cell_types": len(cell_types),
        },
    }
    return CellAdmixFactorSourceScore(fit, result)


def factor_source_annotation(source_prior=None, *, min_score: float = 0.15, min_margin: float = 0.03):
    """Normalize a source-prior object into a factor-source annotation table."""
    if source_prior is None:
        return None
    if isinstance(source_prior, CellAdmixFactorSourceScore):
        return source_prior.annotation(min_score=min_score, min_margin=min_margin)
    if isinstance(source_prior, dict) and source_prior.get("method") == "factor_sources":
        return CellAdmixFactorSourceScore(None, source_prior).annotation(
            min_score=min_score, min_margin=min_margin
        )
    if isinstance(source_prior, pd.DataFrame):
        out = source_prior.copy()
        if "factor" not in out.columns and "factor_label" in out.columns:
            out["factor"] = out["factor_label"].astype(str).str.replace("^F", "", regex=True).astype(int)
        if "factor_label" not in out.columns and "factor" in out.columns:
            out["factor_label"] = "F" + out["factor"].astype(int).astype(str)
        if "source_cell_type" not in out.columns and "cell_type" in out.columns:
            out["source_cell_type"] = out["cell_type"]
        if "called" not in out.columns:
            out["called"] = out["source_cell_type"].notna() & (out["source_cell_type"].astype(str) != "")
        for col in ("score", "margin"):
            if col not in out.columns:
                out[col] = np.nan
        return out[["factor", "factor_label", "source_cell_type", "score", "margin", "called"]]
    if isinstance(source_prior, (pd.Series, dict)):
        series = pd.Series(source_prior)
        factor_label = series.index.astype(str)
        factor = pd.Index(factor_label).str.replace("^F", "", regex=True).astype(int)
        return pd.DataFrame(
            {
                "factor": factor,
                "factor_label": [f"F{x}" for x in factor],
                "source_cell_type": series.astype(str).to_numpy(),
                "score": np.nan,
                "margin": np.nan,
                "called": series.notna().to_numpy(),
            }
        )
    raise TypeError("source_prior must be a CellAdmixFactorSourceScore, DataFrame, Series, dict, or None")


def factor_source_top_genes(result: dict, *, n_genes: int = 8, sources=None) -> pd.DataFrame:
    """Return marker-weighted genes supporting each selected factor source."""
    annotation = factor_source_annotation(result)
    if sources is not None:
        source_map = pd.Series(sources)
        annotation["source_cell_type"] = annotation["factor_label"].map(source_map)
    loadings = result["factor_loadings"]
    marker_weights = result["marker_weights"]
    rows = []
    for _, call in annotation.iterrows():
        source = call["source_cell_type"]
        if pd.isna(source) or str(source) == "":
            continue
        weights = marker_weights[marker_weights["cell_type"] == source]
        if weights.empty:
            continue
        factor_label = call["factor_label"]
        values = loadings.loc[factor_label, weights["gene"].to_numpy()].to_numpy(dtype=float)
        contribution = values * weights["weight"].to_numpy(dtype=float)
        order = np.argsort(-contribution)[: int(n_genes)]
        for rank, pos in enumerate(order, start=1):
            rows.append(
                {
                    "factor": int(call["factor"]),
                    "factor_label": factor_label,
                    "source_cell_type": source,
                    "rank": rank,
                    "gene": weights.iloc[pos]["gene"],
                    "contribution": float(contribution[pos]),
                    "loading": float(values[pos]),
                    "marker_weight": float(weights.iloc[pos]["weight"]),
                }
            )
    return pd.DataFrame(rows)
