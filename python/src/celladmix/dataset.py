"""Dataset-level Python API."""

from __future__ import annotations

import hashlib
import os
import shutil
from pathlib import Path
from typing import Optional

from . import _core
from .fit import CellAdmixFit
from .io import annotation_vectors, ensure_list, read_annotation, recommended_rank
from .state import clustering_result_to_frame


def _default_threads() -> int:
    return max(1, min(10, os.cpu_count() or 1))


def _annotation_hash(annotation) -> str:
    if annotation is None:
        return ""
    digest = hashlib.sha1()
    for cell, label in sorted(zip(annotation.index.astype(str), annotation.astype(str))):
        digest.update(cell.encode())
        digest.update(b"=")
        digest.update(label.encode())
        digest.update(b";")
    return digest.hexdigest()


_COMPARABLE_FIT_PARAMS = (
    "rank", "nmf_variant", "nmf_init", "molecule_scoring", "ncv_k", "graph_k",
    "same_label_ratio", "nmf_iterations", "nmf_n_runs", "nmf_train_max_rows",
    "nmf_min_molecules", "seed", "training_scope_cell_types",
)


def _fit_param_diff(manifest: dict, requested: dict, annotation_hash: str) -> list[str]:
    """Describe requested fit parameters that differ from a cached run.

    Only caller-specified parameters participate; automatically resolved
    values (e.g. an omitted ncv_k) match whatever the cached run recorded.
    """

    def normalize(value):
        if isinstance(value, (list, tuple, set)):
            return ",".join(sorted(map(str, value)))
        return str(value)

    recorded_options = manifest.get("pipeline_options", {}) or {}
    diff = []
    for name in _COMPARABLE_FIT_PARAMS:
        if name not in requested or name not in recorded_options:
            continue
        recorded = normalize(recorded_options[name])
        wanted = normalize(requested[name])
        if recorded != wanted:
            diff.append(f"{name}: {recorded} -> {wanted}")
    recorded_hash = manifest.get("annotation_hash") or ""
    if recorded_hash and annotation_hash and recorded_hash != annotation_hash:
        diff.append("annotation content")
    return diff


class CellAdmix:
    """A cellAdmix dataset backed by a Xenium bundle and an on-disk cache."""

    def __init__(
        self,
        source=None,
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
        self.source = None if source is None else Path(source)
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
        source = None if self.source is None else str(self.source)
        return (
            f"CellAdmix(format={self.format!r}, source={source!r}, "
            f"output_dir={str(self.output_dir)!r}, annotation={ann}, "
            f"num_threads={self.num_threads})"
        )

    @classmethod
    def attach_existing(
        cls,
        output_dir: str | os.PathLike,
        *,
        source=None,
        format: Optional[str] = None,
        annotation=None,
        annotation_col: Optional[str] = None,
        cell_id_col: str = "cell_id",
        num_threads: Optional[int] = None,
    ):
        """Attach to an existing cellAdmix output directory.

        This does not build or scan the input store. If ``source`` is omitted,
        the method uses the source path recorded in ``input_store/store.json``
        when available. Bundle-dependent methods such as membrane image
        discovery need a usable source path.
        """
        output_dir = Path(output_dir)
        store_dir = output_dir / "input_store"
        store_manifest = None
        if store_dir.exists():
            try:
                store_manifest = _core.read_input_store_manifest(str(store_dir))
            except Exception:
                store_manifest = None
        if source is None and store_manifest is not None:
            source = store_manifest.get("source_path") or None
        if format is None and store_manifest is not None:
            format = store_manifest.get("source_type") or "xenium"
        return cls(
            source,
            output_dir=output_dir,
            format=format or "xenium",
            annotation=annotation,
            annotation_col=annotation_col,
            cell_id_col=cell_id_col,
            num_threads=num_threads,
        )

    def load_fit(self, run_id_or_path) -> CellAdmixFit:
        """Load a persisted fit by run id or path."""
        run_path = Path(run_id_or_path)
        if not run_path.exists():
            run_path = self.runs_dir / str(run_id_or_path)
        manifest = _core.read_run_manifest(str(run_path))
        run_dir = Path(manifest["paths"]["root_dir"])
        return CellAdmixFit(self, str(run_dir), manifest)

    def set_annotation(self, annotation, *, annotation_col: Optional[str] = None, cell_id_col: str = "cell_id"):
        """Set or replace the active cell annotation."""
        self.annotation = read_annotation(annotation, annotation_col=annotation_col, cell_id_col=cell_id_col)
        return self

    def ensure_store(self, *, force: bool = False, materialize_molecules: bool = True, verbose: bool = True):
        """Build or reuse the on-disk input store."""
        if self.format != "xenium":
            raise NotImplementedError("Python v1 currently implements Xenium stores")
        if self.source is None:
            raise ValueError(
                "Cannot build an input store without a source path. "
                "Pass source=... or attach to an existing completed run."
            )
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
        nmf_variant: str = "ls_nmf",
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
        explicit_n_runs = nmf_n_runs
        rank = int(rank or recommended_rank(self.annotation, multiplier=rank_multiplier, cap=rank_cap))
        nmf_n_runs = int(nmf_n_runs or threads)
        run_id = run_id or f"fit_rank{rank}_{nmf_variant}"
        run_dir = self.runs_dir / run_id
        annotation_hash = _annotation_hash(self.annotation)
        if (run_dir / "run.json").exists():
            if overwrite:
                shutil.rmtree(run_dir)
            else:
                manifest = _core.read_run_manifest(str(run_dir))
                requested = dict(kwargs)
                requested.update(
                    rank=rank, nmf_variant=nmf_variant, nmf_init=nmf_init,
                    molecule_scoring=molecule_scoring)
                if explicit_n_runs is not None:
                    requested["nmf_n_runs"] = nmf_n_runs
                diff = _fit_param_diff(manifest, requested, annotation_hash)
                if not diff:
                    # Cached runs are reused only when the request matches them.
                    if verbose:
                        recorded = manifest.get("package_version") or ""
                        note = (
                            f"; fitted with package {recorded}, pass overwrite=True to refit"
                            if recorded and recorded != _core.core_version()
                            else ""
                        )
                        print(f"Reusing cached run {run_id} (parameters match{note})")
                    return CellAdmixFit(self, str(run_dir), manifest)
                if verbose:
                    print(f"Parameters changed for run {run_id} ({', '.join(diff)}); refitting")
                shutil.rmtree(run_dir)

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
            annotation_hash=annotation_hash,
            verbose=verbose,
            **kwargs,
        )
        explicit_ncv_k = kwargs.get("ncv_k")
        if explicit_ncv_k and int(explicit_ncv_k) > 0:
            n_genes = len(result.get("genes", []))
            if n_genes > 100 * int(explicit_ncv_k):
                import math
                import warnings

                warnings.warn(
                    f"ncv_k={explicit_ncv_k} is small for a {n_genes}-gene panel: "
                    "neighborhoods carry almost no gene co-occurrence signal and "
                    "NMF factors become seed-dependent. Consider the automatic "
                    f"default (omit ncv_k; ~{round(20 * math.sqrt(n_genes / 400))} "
                    "for this panel, subject to cell size) or gene subsetting.",
                    stacklevel=2,
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
