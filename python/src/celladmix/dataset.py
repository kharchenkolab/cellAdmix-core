"""Dataset-level Python API."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

from . import _core
from .fit import CellAdmixFit
from .io import annotation_vectors, ensure_list, read_annotation, recommended_rank
from .state import clustering_result_to_frame


def _default_threads() -> int:
    return max(1, min(10, os.cpu_count() or 1))


class CellAdmix:
    """A cellAdmix dataset backed by a Xenium bundle and an on-disk cache."""

    def __init__(
        self,
        source,
        *,
        output_dir: str | os.PathLike = "out",
        format: str = "xenium",
        annotation=None,
        annotation_col: Optional[str] = None,
        cell_id_col: str = "cell_id",
        num_threads: Optional[int] = None,
        cell_filter=None,
        gene_filter=None,
        min_qv: float = -1.0,
        keep_unassigned: bool = False,
        keep_non_gene: bool = False,
    ):
        self.source = Path(source)
        self.output_dir = Path(output_dir)
        self.format = format
        self.num_threads = int(num_threads or _default_threads())
        self.cell_filter = ensure_list(cell_filter)
        self.gene_filter = ensure_list(gene_filter)
        self.min_qv = float(min_qv)
        self.keep_unassigned = bool(keep_unassigned)
        self.keep_non_gene = bool(keep_non_gene)
        self.annotation = (
            read_annotation(annotation, annotation_col=annotation_col, cell_id_col=cell_id_col)
            if annotation is not None
            else None
        )
        self.store_dir = self.output_dir / "input_store"
        self.runs_dir = self.output_dir / "runs"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.runs_dir.mkdir(parents=True, exist_ok=True)

    def __repr__(self) -> str:
        ann = "none" if self.annotation is None else f"{self.annotation.nunique()} labels"
        return (
            f"CellAdmix(format={self.format!r}, source={str(self.source)!r}, "
            f"output_dir={str(self.output_dir)!r}, annotation={ann}, "
            f"num_threads={self.num_threads})"
        )

    def set_annotation(self, annotation, *, annotation_col: Optional[str] = None, cell_id_col: str = "cell_id"):
        """Set or replace the active cell annotation."""
        self.annotation = read_annotation(annotation, annotation_col=annotation_col, cell_id_col=cell_id_col)
        return self

    def ensure_store(self, *, force: bool = False, materialize_molecules: bool = True, verbose: bool = True):
        """Build or reuse the on-disk input store."""
        if self.format != "xenium":
            raise NotImplementedError("Python v1 currently implements Xenium stores")
        # The C++ store builder owns parsing and parquet materialization; the
        # Python layer only normalizes user options and paths.
        return _core.build_xenium_store(
            str(self.source),
            str(self.store_dir),
            cell_filter=self.cell_filter,
            gene_filter=self.gene_filter,
            min_qv=self.min_qv,
            keep_unassigned=self.keep_unassigned,
            keep_non_gene=self.keep_non_gene,
            materialize_molecules=materialize_molecules,
            force=force,
            num_threads=self.num_threads,
        )

    def fit(
        self,
        *,
        rank: Optional[int] = None,
        rank_multiplier: float = 1.2,
        rank_cap: int = 30,
        nmf_variant: str = "invsqrt_kl",
        nmf_init: str = "auto",
        nmf_n_runs: Optional[int] = None,
        molecule_scoring: str = "gene_loadings",
        run_id: Optional[str] = None,
        overwrite: bool = False,
        num_threads: Optional[int] = None,
        verbose: bool = True,
        **kwargs,
    ) -> CellAdmixFit:
        """Fit NMF factors and molecule labels.

        By default, rank is inferred from the active annotation and independent
        NMF restarts use the dataset thread count.
        """
        self.ensure_store(force=False)
        threads = int(num_threads or self.num_threads)
        rank = int(rank or recommended_rank(self.annotation, multiplier=rank_multiplier, cap=rank_cap))
        nmf_n_runs = int(nmf_n_runs or threads)
        run_id = run_id or f"fit_rank{rank}_{nmf_variant}"
        run_dir = self.runs_dir / run_id
        if (run_dir / "run.json").exists() and not overwrite:
            # Rehydrate existing native runs instead of refitting by default.
            manifest = _core.read_run_manifest(str(run_dir))
            return CellAdmixFit(self, str(run_dir), manifest)

        cell_ids, labels = annotation_vectors(self.annotation)
        # Training strata are aligned to the store cell order inside C++.
        training_strata = _core.training_strata_for_cells(str(self.store_dir), cell_ids, labels)
        result = _core.fit_store(
            str(self.store_dir),
            str(run_dir),
            rank,
            nmf_variant=nmf_variant,
            nmf_init=nmf_init,
            molecule_scoring=molecule_scoring,
            nmf_n_runs=nmf_n_runs,
            num_threads=threads,
            training_cell_strata=training_strata,
            verbose=verbose,
            **kwargs,
        )
        return CellAdmixFit(self, str(run_dir), result)

    def cell_state_umap(
        self,
        *,
        cells_max: int | None = 5000,
        min_molecules: int = 10,
        min_genes: int = 5,
        n_variable_genes: int = 1000,
        pca_dims: int = 30,
        graph_k: int = 15,
        cluster_resolution: float = 1.0,
        compute_umap: bool = True,
        umap_neighbors: int = 15,
        umap_epochs: int = 200,
        umap_parallel_optimization: bool = True,
        normalization_scale: float = 5000.0,
        seed: int = 1,
        num_threads: Optional[int] = None,
    ):
        """Compute a cell-state UMAP from original input-store counts."""
        self.ensure_store(force=False, materialize_molecules=False)
        result = _core.cluster_store_counts(
            str(self.store_dir),
            min_molecules=min_molecules,
            min_genes=min_genes,
            cells_max=-1 if cells_max is None else int(cells_max),
            n_variable_genes=n_variable_genes,
            pca_dims=pca_dims,
            graph_k=graph_k,
            cluster_resolution=cluster_resolution,
            compute_umap=compute_umap,
            umap_neighbors=umap_neighbors,
            umap_epochs=umap_epochs,
            num_threads=int(num_threads or self.num_threads),
            umap_parallel_optimization=umap_parallel_optimization,
            normalization_scale=normalization_scale,
            seed=int(seed),
        )
        return clustering_result_to_frame(result, annotation=self.annotation)
