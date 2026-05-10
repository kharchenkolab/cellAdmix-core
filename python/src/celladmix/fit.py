"""Fitted cellAdmix run objects."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import _core
from .io import annotation_vectors, discover_xenium_membrane_image, discover_xenium_stain_image
from .score import CellAdmixScore


class CellAdmixFit:
    """A fitted cellAdmix run."""

    def __init__(self, dataset, run_path: str, manifest: dict):
        self.dataset = dataset
        self.run_path = Path(run_path)
        self.manifest = manifest
        self._molecule_numeric_cache = None

    def __repr__(self) -> str:
        return (
            f"CellAdmixFit(run_path={str(self.run_path)!r}, "
            f"n_factors={self.manifest.get('n_factors')}, "
            f"n_cells={self.manifest.get('n_cells')})"
        )

    def cells(self) -> pd.DataFrame:
        """Return per-cell factor fractions and coordinates."""
        return pd.read_parquet(self.manifest["paths"]["cells_parquet"])

    def cell_factors(self) -> pd.DataFrame:
        """Alias for ``cells()`` for parity with the R API naming."""
        return self.cells()

    def molecules(self, columns=None) -> pd.DataFrame:
        """Return per-molecule labels and coordinates."""
        frame = pd.read_parquet(self.manifest["paths"]["molecules_parquet"], columns=columns)
        return frame

    def region(self, *, bbox=None, cache: bool = True) -> pd.DataFrame:
        """Return molecule records in a physical-coordinate bbox.

        The persisted run stores gene and cell identifiers as integer indexes;
        this helper maps them back to strings after spatial filtering.
        """
        columns = ["x", "y", "gene_idx", "cell_idx", "factor_label"]
        if cache:
            if self._molecule_numeric_cache is None:
                self._molecule_numeric_cache = pd.read_parquet(
                    self.manifest["paths"]["molecules_parquet"], columns=columns
                )
            frame = self._molecule_numeric_cache
        else:
            frame = pd.read_parquet(self.manifest["paths"]["molecules_parquet"], columns=columns)
        if bbox is not None:
            xmin, xmax, ymin, ymax = bbox
            frame = frame[
                frame["x"].between(xmin, xmax)
                & frame["y"].between(ymin, ymax)
            ].copy()
        else:
            frame = frame.copy()
        genes = self._gene_index()
        cells = self._cell_index()
        frame["gene"] = frame["gene_idx"].map(genes)
        frame["cell_id"] = frame["cell_idx"].map(cells)
        frame["factor_label"] = frame["factor_label"].astype(int) + 1
        return frame

    def _gene_index(self) -> dict[int, str]:
        path = self.dataset.store_dir / "genes.parquet"
        genes = pd.read_parquet(path)
        return dict(zip(genes["gene_idx"].astype(int), genes["gene"].astype(str)))

    def _cell_index(self) -> dict[int, str]:
        path = self.dataset.store_dir / "cells.parquet"
        cells = pd.read_parquet(path)
        return dict(zip(cells["cell_idx"].astype(int), cells["cell_id"].astype(str)))

    def factor_loadings(self) -> pd.DataFrame:
        """Return an approximate gene-by-factor loading matrix."""
        h = self.manifest.get("h")
        if h is None:
            # Older or externally written runs may only expose loadings parquet.
            factors = pd.read_parquet(self.manifest["paths"]["factors_parquet"])
            if {"factor_id", "gene", "loading"}.issubset(factors.columns):
                wide = factors.pivot(index="gene", columns="factor_id", values="loading")
                wide = wide.rename(columns={col: f"F{int(col)}" for col in wide.columns})
                return wide
            return factors
        data = h["data"]
        rows = int(h["rows"])
        cols = int(h["cols"])
        matrix = pd.DataFrame(
            [data[i * cols : (i + 1) * cols] for i in range(rows)],
            index=[f"F{i + 1}" for i in range(rows)],
            columns=self.manifest.get("genes", [f"gene_{i + 1}" for i in range(cols)]),
        )
        return matrix.T

    def counts(self):
        """Return run counts as a scipy CSC matrix with genes x cells orientation."""
        from scipy import sparse
        import numpy as np

        raw = _core.collect_counts(str(self.run_path))
        data = np.asarray(raw["data"], dtype=float)
        indices = np.asarray(raw["indices"], dtype=np.int32)
        indptr = np.asarray(raw["indptr"], dtype=np.int64)
        genes = pd.Index(list(raw["genes"]), name="gene")
        cells = pd.Index(list(raw["cells"]), name="cell_id")
        matrix = sparse.csc_matrix((data, indices, indptr), shape=(len(genes), len(cells)))
        return matrix, genes, cells

    def score_factor_sources(self, *, annotation=None, counts=None, **kwargs):
        """Score factor source cell types using marker-weighted gene content."""
        from .factor_sources import score_factor_sources

        return score_factor_sources(self, annotation=annotation, counts=counts, **kwargs)

    def score_membrane(self, *, image_path=None, pixel_size=None, num_threads=None, verbose=True, **kwargs):
        """Score admixture using the membrane-stain heuristic."""
        if image_path is None or pixel_size is None:
            # Resolve Xenium image metadata lazily so callers can override
            # either the path or pixel size without rebuilding the dataset.
            discovered = discover_xenium_membrane_image(self.dataset.source)
            image_path = image_path or discovered["image_path"]
            pixel_size = pixel_size or discovered["pixel_size"]
            kwargs.setdefault("x_offset", discovered["x_offset"])
            kwargs.setdefault("y_offset", discovered["y_offset"])
        cell_ids, labels = annotation_vectors(self.dataset.annotation)
        # Annotation vectors are passed explicitly because stored runs can be
        # scored against different cell-type labels without refitting.
        result = _core.score_membrane(
            str(self.run_path),
            str(image_path),
            float(pixel_size),
            num_threads=int(num_threads or self.dataset.num_threads),
            annotation_cell_ids=cell_ids,
            annotation_labels=labels,
            verbose=verbose,
            **kwargs,
        )
        return CellAdmixScore(self, "membrane", result)

    def score_bridge(self, *, num_threads=None, verbose=True, **kwargs):
        """Score admixture using the molecular bridge heuristic."""
        cell_ids, labels = annotation_vectors(self.dataset.annotation)
        # Keep the bridge wrapper parallel to membrane scoring: same run,
        # current annotation, and method-specific native options.
        result = _core.score_bridge(
            str(self.run_path),
            num_threads=int(num_threads or self.dataset.num_threads),
            annotation_cell_ids=cell_ids,
            annotation_labels=labels,
            verbose=verbose,
            **kwargs,
        )
        return CellAdmixScore(self, "bridge", result)

    def correct(self, rules: pd.DataFrame, *, name: str = "clean"):
        """Remove molecules matching factor/target-cell-type correction rules."""
        from .correction import CellAdmixCorrection

        if rules.empty:
            raise ValueError("rules table is empty")
        cell_ids, labels = annotation_vectors(self.dataset.annotation)
        out_dir = self.run_path / "corrected" / name
        result = _core.correct_run(
            str(self.run_path),
            str(out_dir),
            factors=rules["factor"].astype(int).tolist(),
            target_cell_types=rules["target_cell_type"].astype(str).tolist(),
            annotation_cell_ids=cell_ids,
            annotation_labels=labels,
        )
        return CellAdmixCorrection(str(out_dir), result)

    def plot_loadings(self, *, n_genes: int = 8, **kwargs):
        from .plotting import plot_loadings

        return plot_loadings(self.factor_loadings(), n_genes=n_genes, **kwargs)

    def plot_stability(self, **kwargs):
        """Plot NMF restart stability against factor importance."""
        from .plotting import plot_stability

        return plot_stability(self, **kwargs)

    def stain(self, stain: str = "membrane", **kwargs) -> dict:
        """Resolve Xenium stain metadata without loading image pixels."""
        return discover_xenium_stain_image(self.dataset.source, stain=stain, **kwargs)

    def stain_image(self, stain: str = "membrane", **kwargs) -> dict:
        """Alias for ``stain()`` for API parity with R."""
        return self.stain(stain, **kwargs)

    def stain_crop(self, image: dict, *, bbox, max_pixels: int = 512) -> dict:
        """Read a small crop from a stain image descriptor."""
        from .examples import read_stain_crop

        return read_stain_crop(image, bbox, max_pixels=max_pixels)

    def cell_boundaries(self, *, boundary_path=None, cells=None, bbox=None) -> pd.DataFrame:
        """Read Xenium cell-boundary vertices for plotting."""
        from .examples import discover_cell_boundaries, read_cell_boundaries

        path = discover_cell_boundaries(self.dataset.source, boundary_path=boundary_path)
        if path is None:
            return pd.DataFrame(columns=["cell_id", "x", "y"])
        return read_cell_boundaries(path, cells=cells, bbox=bbox)

    def prepare_cell_example(self, example, **kwargs) -> dict:
        """Prepare one score-selected example cell for plotting."""
        from .examples import prepare_cell_example

        return prepare_cell_example(self, example, **kwargs)

    def plot_cell_example(self, example, **kwargs):
        """Plot one score-selected example cell."""
        from .examples import plot_cell_example, prepare_cell_example

        if not isinstance(example, dict) or "molecules" not in example:
            example = prepare_cell_example(self, example, **kwargs)
            kwargs = {}
        return plot_cell_example(example, **kwargs)
