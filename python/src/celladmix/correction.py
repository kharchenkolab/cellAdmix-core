"""Corrected run objects."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import _core
from .state import clustering_result_to_frame


class CellAdmixCorrection:
    """A corrected cellAdmix run."""

    def __init__(self, run_path: str, manifest: dict, rules=None):
        self.run_path = Path(run_path)
        self.manifest = manifest
        self.rules = rules

    def __repr__(self) -> str:
        removed = self.manifest.get("n_removed", "unknown")
        return f"CellAdmixCorrection(run_path={str(self.run_path)!r}, n_removed={removed})"

    def summary(self) -> pd.DataFrame:
        """Return a compact per-cell-type correction summary."""
        return _compact_correction_summary(self.cell_summary())

    def ensemble(self) -> dict:
        """Return molecule-vote ensemble diagnostics for this correction."""
        histogram = self.manifest.get("vote_histogram")
        return {
            "members": self.manifest.get("n_members", 1),
            "min_votes": self.manifest.get("min_votes", 1),
            "n_removed": self.manifest.get("n_removed"),
            "vote_histogram": histogram,
        }

    def cell_summary(self) -> pd.DataFrame:
        """Return per-cell molecule removal statistics."""
        path = self.manifest.get("correction_summary_parquet")
        if path is None:
            path = self.run_path / "correction_summary.parquet"
        return pd.read_parquet(path)

    def counts(self):
        """Return corrected counts as a scipy CSC matrix with genes x cells orientation."""
        from scipy import sparse

        # The native layer returns CSC arrays directly; keep orientation explicit
        # here because AnnData expects the transpose.
        raw = _core.collect_counts(str(self.run_path))
        data = np.asarray(raw["data"], dtype=float)
        indices = np.asarray(raw["indices"], dtype=np.int32)
        indptr = np.asarray(raw["indptr"], dtype=np.int64)
        genes = list(raw["genes"])
        cells = list(raw["cells"])
        matrix = sparse.csc_matrix((data, indices, indptr), shape=(len(genes), len(cells)))
        return matrix, pd.Index(genes, name="gene"), pd.Index(cells, name="cell_id")

    def to_anndata(self):
        """Return corrected counts as an AnnData object with cells x genes orientation."""
        import anndata as ad

        matrix, genes, cells = self.counts()
        return ad.AnnData(X=matrix.T.tocsr(), obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))

    def cell_state_umap(
        self,
        *,
        annotation: pd.Series | None = None,
        cells_max: int = 5000,
        min_molecules: int = 10,
        min_genes: int = 5,
        n_variable_genes: int = 1000,
        pca_dims: int = 30,
        graph_k: int = 15,
        umap_neighbors: int = 15,
        umap_epochs: int = 200,
        normalization_scale: float = 5000.0,
        seed: int = 1,
        num_threads: int = 1,
    ) -> pd.DataFrame:
        """Compute a cell-state UMAP from corrected counts."""
        # Corrected runs are already materialized on disk, so clustering can
        # reuse the same native count-collection path as uncorrected runs.
        result = _core.cluster_run_counts(
            str(self.run_path),
            min_molecules=min_molecules,
            min_genes=min_genes,
            cells_max=int(cells_max),
            n_variable_genes=n_variable_genes,
            pca_dims=pca_dims,
            graph_k=graph_k,
            compute_umap=True,
            umap_neighbors=umap_neighbors,
            umap_epochs=umap_epochs,
            num_threads=int(num_threads),
            normalization_scale=normalization_scale,
            seed=int(seed),
        )
        return clustering_result_to_frame(result, annotation=annotation)

    def plot_removed_molecules(self):
        """Plot per-cell molecule removal counts and fractions."""
        from .plotting import plot_molecule_removal

        return plot_molecule_removal(self.cell_summary())


def _compact_correction_summary(summary: pd.DataFrame) -> pd.DataFrame:
    required = {
        "cell_type",
        "n_molecules_before",
        "n_molecules_after",
        "n_removed",
        "fraction_removed",
    }
    if summary.empty:
        return pd.DataFrame()
    missing = required.difference(summary.columns)
    if missing:
        raise ValueError(f"correction summary is missing columns: {', '.join(sorted(missing))}")
    data = summary.copy()
    data["cell_type"] = data["cell_type"].fillna("unknown").replace("", "unknown").astype(str)
    data["modified"] = data["n_removed"].astype(float) > 0

    def summarize(group: pd.DataFrame, label: str) -> dict:
        modified = group[group["modified"]]
        molecules_before = float(group["n_molecules_before"].sum())
        molecules_removed = float(group["n_removed"].sum())
        return {
            "cell_type": label,
            "n_cells": int(len(group)),
            "n_modified_cells": int(group["modified"].sum()),
            "fraction_cells_modified": float(group["modified"].sum() / max(1, len(group))),
            "molecules_before": molecules_before,
            "molecules_after": float(group["n_molecules_after"].sum()),
            "molecules_removed": molecules_removed,
            "fraction_molecules_removed": float(molecules_removed / max(1.0, molecules_before)),
            "median_removed_per_modified_cell": float(modified["n_removed"].median()) if len(modified) else 0.0,
            "median_fraction_removed_per_modified_cell": (
                float(modified["fraction_removed"].median()) if len(modified) else 0.0
            ),
        }

    rows = [summarize(data, "all")]
    for cell_type, group in data.groupby("cell_type", sort=True):
        rows.append(summarize(group, str(cell_type)))
    return pd.DataFrame(rows)
