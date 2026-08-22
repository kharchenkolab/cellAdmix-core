"""Fitted cellAdmix run objects."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import _core
from .io import annotation_vectors, discover_xenium_membrane_image, discover_xenium_stain_image
from .score import CellAdmixScore

# Minimum per-type molecule count before the over-removal warning applies
# (kept module-level so tests can lower it for small fixtures).
OVERREMOVAL_MIN_MOLECULES = 10000


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

    @classmethod
    def load(
        cls,
        run_path,
        *,
        source=None,
        annotation=None,
        annotation_col=None,
        cell_id_col: str = "cell_id",
        num_threads=None,
    ):
        """Attach to a persisted run directory without refitting.

        Parameters
        ----------
        run_path:
            Path to a run directory or its ``run.json`` manifest.
        source:
            Optional original Xenium bundle. If omitted, the loader uses the
            source path recorded in the run manifest when available. Methods
            that need bundle-side files, such as membrane-image discovery, will
            ask for ``source`` if no usable path is attached.
        annotation:
            Optional cell annotation used by scoring/correction methods.
        """
        from .dataset import CellAdmix

        run_path = Path(run_path)
        manifest = _core.read_run_manifest(str(run_path))
        run_dir = Path(manifest["paths"]["root_dir"])
        output_dir = _output_dir_from_run_dir(run_dir)
        source = source or _source_from_run_manifest(manifest)
        dataset = CellAdmix.attach_existing(
            output_dir,
            source=source,
            annotation=annotation,
            annotation_col=annotation_col,
            cell_id_col=cell_id_col,
            num_threads=num_threads,
        )
        return cls(dataset, str(run_dir), manifest)

    def cells(self) -> pd.DataFrame:
        """Return per-cell factor fractions and coordinates."""
        return pd.read_parquet(self.manifest["paths"]["cells_parquet"])

    def cell_factors(self) -> pd.DataFrame:
        """Alias for ``cells()`` for parity with the R API naming."""
        return self.cells()

    def molecules(self, columns=None, *, raw: bool = False) -> pd.DataFrame:
        """Return per-molecule labels and coordinates.

        Public Python output uses one-based factors: ``factor`` is an integer
        in ``1..K`` and ``factor_label`` is ``F1..FK``. The persisted parquet
        file stores the implementation label as a zero-based integer; pass
        ``raw=True`` to read those columns without conversion.
        """
        if raw:
            return pd.read_parquet(self.manifest["paths"]["molecules_parquet"], columns=columns)

        requested = _normalize_columns(columns)
        parquet_columns = _molecule_parquet_columns(requested)
        frame = pd.read_parquet(self.manifest["paths"]["molecules_parquet"], columns=parquet_columns)
        return _public_molecule_factors(frame, requested)

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
        return _public_molecule_factors(frame, requested=None)

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

    def score_molecules(
        self,
        molecules: pd.DataFrame,
        *,
        gene_col: str = "gene",
        cell_col: str = "cell_id",
        x_col: str = "x",
        y_col: str = "y",
        return_scores: bool = False,
        smooth: bool = False,
        graph_k: int | None = None,
        same_label_ratio: float | None = None,
        max_iterations: int = 20,
    ) -> pd.DataFrame:
        """Score new same-fit molecules against the learned factor loadings.

        This is a lightweight same-bundle/synthetic helper. It uses the current
        gene-loading molecule potential, so no NMF refit is performed. Unknown
        genes receive a uniform factor score. When ``smooth=True``, ``cell_id``,
        ``x``, and ``y`` are required and a Python implementation of the same
        within-cell ICM smoothing used by the core pipeline is applied.
        """
        original_index = molecules.index
        frame = molecules.copy().reset_index(drop=True)
        if gene_col not in frame.columns:
            raise ValueError(f"score_molecules() requires a gene column {gene_col!r}")

        loadings = self.factor_loadings()
        factor_labels = list(loadings.columns.astype(str))
        n_factors = len(factor_labels)
        if n_factors == 0:
            raise ValueError("fit does not contain any factor loadings")
        h = loadings.to_numpy(dtype=float)
        h = np.maximum(h, 0.0)
        gene_to_row = {str(gene): idx for idx, gene in enumerate(loadings.index.astype(str))}

        scores = np.full((len(frame), n_factors), 1.0 / n_factors, dtype=float)
        genes = frame[gene_col].astype(str).to_numpy()
        for row, gene in enumerate(genes):
            idx = gene_to_row.get(gene)
            if idx is None:
                continue
            values = h[idx, :]
            total = float(np.sum(values))
            if total > np.finfo(float).eps:
                scores[row, :] = values / total

        factors = np.argmax(scores, axis=1) + 1
        if smooth:
            for required in (cell_col, x_col, y_col):
                if required not in frame.columns:
                    raise ValueError(f"smooth=True requires column {required!r}")
            factors = _smooth_molecule_scores(
                frame,
                scores,
                cell_col=cell_col,
                x_col=x_col,
                y_col=y_col,
                graph_k=graph_k if graph_k is not None else int(self.manifest.get("pipeline_options", {}).get("graph_k", 10)),
                same_label_ratio=same_label_ratio
                if same_label_ratio is not None
                else float(self.manifest.get("pipeline_options", {}).get("same_label_ratio", 5.0)),
                max_iterations=max_iterations,
            )

        sorted_scores = np.sort(scores, axis=1)
        margins = sorted_scores[:, -1] - (sorted_scores[:, -2] if n_factors > 1 else 0.0)
        out = frame.copy()
        out["factor"] = factors.astype(int)
        out["factor_label"] = ["F" + str(int(x)) for x in factors]
        out["factor_margin"] = margins
        if return_scores:
            for idx, label in enumerate(factor_labels):
                out[label] = scores[:, idx]
        out.index = original_index
        return out

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

    def score_neighbor_enrichment(self, **kwargs):
        """Score admixture by source-cell neighborhood enrichment."""
        from .neighbor import score_neighbor_enrichment

        return score_neighbor_enrichment(self, **kwargs)

    def audit_admixture(self, **kwargs):
        """Estimate per-cell-type-pair admixture from spatial exposure.

        Computed from this run's counts and cell table; the audit involves
        no factorization, and ``dataset.audit_admixture()`` gives the same
        measurement straight from the input store.
        """
        from .audit import CellAdmixAudit

        return CellAdmixAudit.from_fit(self, **kwargs)

    def score_factor_sources(self, *, annotation=None, counts=None, **kwargs):
        """Score factor source cell types using marker-weighted gene content."""
        from .factor_sources import score_factor_sources

        return score_factor_sources(self, annotation=annotation, counts=counts, **kwargs)

    def score_membrane(self, *, image_path=None, pixel_size=None, num_threads=None, verbose=True, **kwargs):
        """Score admixture using the membrane-stain heuristic."""
        if image_path is None or pixel_size is None:
            # Resolve Xenium image metadata lazily so callers can override
            # either the path or pixel size without rebuilding the dataset.
            if self.dataset.source is None:
                raise ValueError(
                    "score_membrane() needs the original Xenium bundle path. "
                    "Pass source=... to CellAdmix.attach_existing() or CellAdmixFit.load()."
                )
            discovered = discover_xenium_membrane_image(self.dataset.source)
            image_path = image_path or discovered["image_path"]
            pixel_size = pixel_size or discovered["pixel_size"]
            kwargs.setdefault("x_offset", discovered["x_offset"])
            kwargs.setdefault("y_offset", discovered["y_offset"])
        cell_ids, labels = annotation_vectors(self.dataset.annotation)
        # Annotation vectors are passed explicitly because stored runs can be
        # scored against different cell-type labels without refitting.
        ensemble_member = int(kwargs.pop("ensemble_member", -1))
        result = _core.score_membrane(
            str(self.run_path),
            str(image_path),
            float(pixel_size),
            num_threads=int(num_threads or self.dataset.num_threads),
            annotation_cell_ids=cell_ids,
            annotation_labels=labels,
            verbose=verbose,
            ensemble_member=ensemble_member,
            **kwargs,
        )
        params = dict(kwargs)
        params.update(
            image_path=str(image_path),
            pixel_size=float(pixel_size),
            num_threads=int(num_threads or self.dataset.num_threads),
            verbose=False,
        )
        return CellAdmixScore(self, "membrane", result, params=params)

    def score_bridge(self, *, num_threads=None, verbose=True, **kwargs):
        """Score admixture using the molecular bridge heuristic."""
        cell_ids, labels = annotation_vectors(self.dataset.annotation)
        # Keep the bridge wrapper parallel to membrane scoring: same run,
        # current annotation, and method-specific native options.
        ensemble_member = int(kwargs.pop("ensemble_member", -1))
        result = _core.score_bridge(
            str(self.run_path),
            num_threads=int(num_threads or self.dataset.num_threads),
            annotation_cell_ids=cell_ids,
            annotation_labels=labels,
            verbose=verbose,
            ensemble_member=ensemble_member,
            **kwargs,
        )
        params = dict(kwargs)
        params.update(
            num_threads=int(num_threads or self.dataset.num_threads),
            verbose=False,
        )
        return CellAdmixScore(self, "bridge", result, params=params)

    def score_neighbor_frequency(self, **kwargs):
        """Alias for :meth:`score_neighbor_enrichment`."""
        return self.score_neighbor_enrichment(**kwargs)

    def correct(self, rules: pd.DataFrame, *, name: str = "clean",
                rule_members=None, min_votes: int = 1):
        """Remove molecules matching factor/target-cell-type correction rules.

        With ``rule_members``, each rule row votes with the labeling of one
        ensemble member (see :meth:`CellAdmixScore.correct`); a molecule is
        removed when at least ``min_votes`` members remove it.
        """
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
            rule_members=[int(m) for m in (rule_members or [])],
            min_votes=int(min_votes),
        )
        correction = CellAdmixCorrection(str(out_dir), result, rules=rules)
        # Over-removal guard: a rule set whose factors jointly cover a cell
        # type's whole molecule content (typically because the factorization
        # has no native factor for a small type) erases the type rather than
        # cleaning it. Surface that immediately, not only at audit time.
        import warnings

        try:
            removal = correction.summary()
        except Exception:
            removal = None
        if removal is not None and len(removal):
            heavy = removal[(removal["cell_type"] != "all")
                & (removal["molecules_before"] >= OVERREMOVAL_MIN_MOLECULES)
                & (removal["fraction_molecules_removed"] > 0.6)]
            for _, row in heavy.iterrows():
                warnings.warn(
                    f"Correction removed {100 * row['fraction_molecules_removed']:.0f}% "
                    f"of all molecules from {row['cell_type']} - this is likely "
                    "erasing native expression (does the factorization have a "
                    "native factor for this type?)", stacklevel=2)
        return correction

    def plot_loadings(self, *, n_genes: int = 8, **kwargs):
        from .plotting import plot_loadings

        return plot_loadings(self.factor_loadings(), n_genes=n_genes, **kwargs)

    def plot_stability(self, **kwargs):
        """Plot NMF restart stability against factor importance."""
        from .plotting import plot_stability

        return plot_stability(self, **kwargs)

    def stain(self, stain: str = "membrane", **kwargs) -> dict:
        """Resolve Xenium stain metadata without loading image pixels."""
        if self.dataset.source is None:
            raise ValueError(
                "stain() needs the original Xenium bundle path. "
                "Pass source=... to CellAdmix.attach_existing() or CellAdmixFit.load()."
            )
        return discover_xenium_stain_image(self.dataset.source, stain=stain, **kwargs)

    def stain_image(self, stain: str = "membrane", **kwargs) -> dict:
        """Alias for ``stain()`` for API parity with R."""
        return self.stain(stain, **kwargs)

    def stains(self, names=("dapi", "membrane")) -> dict:
        """Resolve all available stain images, empty for non-Xenium sources."""
        from .io import discover_xenium_stain_images

        return discover_xenium_stain_images(self.dataset.source, names)

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


def _output_dir_from_run_dir(run_dir: Path) -> Path:
    run_dir = Path(run_dir)
    return run_dir.parent.parent if run_dir.parent.name == "runs" else run_dir.parent


def _source_from_run_manifest(manifest: dict):
    source = manifest.get("source") or {}
    path = source.get("path")
    return path or None


def _molecule_parquet_columns(requested):
    if requested is None:
        return None
    parquet_columns = []
    for column in requested:
        source = "factor_label" if column in {"factor", "factor_label"} else column
        if source not in parquet_columns:
            parquet_columns.append(source)
    return parquet_columns


def _normalize_columns(columns):
    if columns is None:
        return None
    if isinstance(columns, str):
        return [columns]
    return list(columns)


def _public_molecule_factors(frame: pd.DataFrame, requested) -> pd.DataFrame:
    if "factor_label" not in frame.columns:
        return frame if requested is None else frame[[c for c in requested if c in frame.columns]]

    raw = pd.to_numeric(frame["factor_label"], errors="coerce")
    factor = raw.astype("Int64") + 1
    invalid = raw.isna() | (raw < 0)
    factor = factor.mask(invalid)
    labels = pd.Series(pd.NA, index=frame.index, dtype="object")
    valid = factor.notna()
    labels.loc[valid] = "F" + factor.loc[valid].astype(int).astype(str)

    out = frame.drop(columns=["factor_label"]).copy()
    out["factor"] = factor
    out["factor_label"] = labels
    if requested is None:
        return out
    keep = [column for column in requested if column in out.columns]
    return out.loc[:, keep]


def _smooth_molecule_scores(
    frame: pd.DataFrame,
    scores: np.ndarray,
    *,
    cell_col: str,
    x_col: str,
    y_col: str,
    graph_k: int,
    same_label_ratio: float,
    max_iterations: int,
) -> np.ndarray:
    from scipy.spatial import cKDTree

    labels = np.argmax(scores, axis=1).astype(int)
    n_factors = scores.shape[1]
    smoothness = np.log(max(float(same_label_ratio), 1.0))
    if graph_k <= 0 or max_iterations <= 0 or smoothness <= 0.0 or len(frame) == 0:
        return labels + 1

    log_scores = np.log(np.maximum(scores, 1e-12))
    for _, index in frame.groupby(cell_col, sort=False).groups.items():
        idx = np.asarray(index, dtype=int)
        if len(idx) <= 1:
            continue
        xy = frame.loc[idx, [x_col, y_col]].to_numpy(dtype=float)
        finite = np.isfinite(xy).all(axis=1)
        if finite.sum() <= 1:
            continue
        active = idx[finite]
        coords = xy[finite]
        k = min(int(graph_k) + 1, len(active))
        _, neighbors = cKDTree(coords).query(coords, k=k)
        if k == 1:
            continue
        if neighbors.ndim == 1:
            neighbors = neighbors[:, None]
        adjacency = [set() for _ in range(len(active))]
        for local, row_neighbors in enumerate(neighbors):
            for neighbor in row_neighbors:
                neighbor = int(neighbor)
                if neighbor == local:
                    continue
                adjacency[local].add(neighbor)
                adjacency[neighbor].add(local)
        local_labels = labels[active].copy()
        local_scores = log_scores[active, :]
        counts = np.zeros(n_factors, dtype=int)
        for _ in range(int(max_iterations)):
            changed = 0
            for local, neighbor_set in enumerate(adjacency):
                counts.fill(0)
                for neighbor in neighbor_set:
                    counts[local_labels[neighbor]] += 1
                candidate = local_scores[local, :] + smoothness * counts
                best = int(np.argmax(candidate))
                if best != local_labels[local]:
                    local_labels[local] = best
                    changed += 1
            if changed == 0:
                break
        labels[active] = local_labels
    return labels + 1
