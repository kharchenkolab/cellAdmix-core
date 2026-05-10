"""SpatialData integration helpers.

The core still reads high-throughput Xenium molecules directly from the bundle.
SpatialData is used as the first-class source for user-facing cell annotations
and as a destination for fitted/corrected results.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from .dataset import CellAdmix


def _default_table(sdata, table: Optional[str]):
    if table is not None:
        return sdata.tables[table]
    if len(sdata.tables) != 1:
        raise ValueError("Specify table= because the SpatialData object has multiple tables")
    return next(iter(sdata.tables.values()))


def from_spatialdata(
    sdata,
    *,
    xenium_dir,
    output_dir="out",
    table: Optional[str] = None,
    annotation: Optional[str | pd.Series] = None,
    cell_id_key: Optional[str] = None,
    num_threads: Optional[int] = None,
):
    """Construct a ``CellAdmix`` dataset from a SpatialData Xenium object.

    The Xenium bundle path is still explicit so cellAdmix can stream molecule
    tables and image files from disk. If ``annotation`` is a string, it is read
    from the selected AnnData table's ``obs``.
    """
    adata = _default_table(sdata, table)
    cell_ids = (
        adata.obs[cell_id_key].astype(str)
        if cell_id_key is not None
        else pd.Series(adata.obs_names.astype(str), index=adata.obs_names)
    )
    if annotation is None:
        for candidate in ("cell_type", "annotation", "cluster", "cluster_id"):
            if candidate in adata.obs.columns:
                annotation = candidate
                break
    if isinstance(annotation, str):
        if annotation not in adata.obs.columns:
            raise ValueError(f"annotation column {annotation!r} is absent from SpatialData table")
        # Re-key annotations by the Xenium/cellAdmix cell ids, which may be
        # stored separately from AnnData obs_names in some SpatialData tables.
        values = pd.Series(adata.obs[annotation].to_numpy(), index=cell_ids.to_numpy(), name=annotation)
        ann = values.dropna()
    elif annotation is None:
        ann = None
    else:
        ann = annotation
    return CellAdmix(Path(xenium_dir), output_dir=output_dir, annotation=ann, num_threads=num_threads)


def add_fit_to_spatialdata(sdata, fit, *, table: Optional[str] = None, prefix: str = "celladmix_"):
    """Add per-cell factor summaries to a SpatialData table's obs."""
    adata = _default_table(sdata, table)
    cells = fit.cell_factors().set_index("cell_id")
    aligned = cells.reindex(adata.obs_names.astype(str))
    for column in aligned.columns:
        if column.startswith("factor_") or column in ("dominant_factor", "dominant_fraction"):
            # AnnData obs is cell-major; native fit output is keyed by cell_id.
            adata.obs[prefix + column] = aligned[column].to_numpy()
    return sdata


def add_corrected_counts_to_spatialdata(
    sdata,
    correction,
    *,
    table: Optional[str] = None,
    layer: str = "celladmix_corrected_counts",
):
    """Add corrected counts to a SpatialData table as an AnnData layer."""
    from scipy import sparse

    adata = _default_table(sdata, table)
    matrix, genes, cells = correction.counts()
    row_order = pd.Index(adata.obs_names.astype(str))
    col_order = pd.Index(adata.var_names.astype(str))
    cell_pos = cells.get_indexer(row_order)
    gene_pos = genes.get_indexer(col_order)
    valid_cells = cell_pos >= 0
    valid_genes = gene_pos >= 0
    if not valid_cells.any() or not valid_genes.any():
        raise ValueError("Corrected counts do not overlap the SpatialData table")
    # correction.counts() returns genes x cells; AnnData layers are cells x genes.
    sub = matrix[gene_pos[valid_genes], :][:, cell_pos[valid_cells]].T.tocoo()
    row_index = pd.RangeIndex(len(row_order))[valid_cells].to_numpy()[sub.row]
    col_index = pd.RangeIndex(len(col_order))[valid_genes].to_numpy()[sub.col]
    adata.layers[layer] = sparse.csr_matrix(
        (sub.data, (row_index, col_index)),
        shape=adata.shape,
    )
    return sdata
