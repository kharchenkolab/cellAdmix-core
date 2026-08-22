"""Generative admixture model.

Decomposes each cell's counts into its own expression, contamination from
each detected source type, ambient background, and neighborhood-induced
expression. The fitted model exposes the per-cell decomposition and the
retained induced genes, and derives corrections (with or without
retaining the induced expression) without refitting. The model, its
validation, and the limits of the separation are described in
docs/generative.md. The fit runs in the C++ core; this module marshals
the audit's measurements in and wraps the results.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp


class GenerativeCorrection:
    """Corrected counts derived from a fitted generative model."""

    def __init__(self, matrix, genes, cells, rules, name, retained_induced,
                 cell_types=None):
        self._matrix = matrix
        self._genes = genes
        self._cells = cells
        self.rules = rules
        self.name = name
        self.retained_induced = retained_induced
        self._cell_types = cell_types

    def counts(self):
        """Corrected counts as ``(matrix, genes, cells)``."""
        return self._matrix, list(self._genes), list(self._cells)

    def cell_state_umap(self, *, annotation=None, cells_max=5000,
                        min_molecules=10, min_genes=5, n_variable_genes=1000,
                        pca_dims=30, graph_k=15, umap_neighbors=15,
                        umap_epochs=200, normalization_scale=5000.0, seed=1,
                        num_threads=1):
        """Cell-state UMAP of the corrected counts, matching the other
        corrections' ``cell_state_umap``."""
        from . import _core
        from .state import clustering_result_to_frame

        m = self._matrix.tocsc()
        result = _core.cluster_counts_matrix(
            indptr=m.indptr.tolist(), indices=m.indices.tolist(),
            values=m.data.astype(float).tolist(),
            genes=[str(g) for g in self._genes],
            cells=[str(c) for c in self._cells],
            min_molecules=min_molecules, min_genes=min_genes,
            cells_max=int(cells_max), n_variable_genes=n_variable_genes,
            pca_dims=pca_dims, graph_k=graph_k,
            umap_neighbors=umap_neighbors, umap_epochs=umap_epochs,
            num_threads=int(num_threads),
            normalization_scale=normalization_scale, seed=int(seed))
        if annotation is None:
            annotation = self._cell_types
        return clustering_result_to_frame(result, annotation=annotation)


class GenerativeModel:
    """A fitted generative admixture model.

    Attributes: ``pairs`` (per-pair dose, expected removal and induced
    molecule totals), ``induced`` (the retained induced genes with their
    disproportionality statistics), ``n_programs`` (expression programs
    fitted per type), ``whole_cell_profiles`` (True when the run carries
    no nucleus flag and whole-cell transfer profiles were used).
    """

    def __init__(self, audit, matrix, genes, cells, removed, removed_strict,
                 pairs, induced, composition, ambient_scale, n_programs,
                 whole_cell_profiles, init):
        self._cell_types = audit._ctypes
        self._audit = audit
        self._matrix = matrix
        self._genes = genes
        self._cells = cells
        self._removed = removed
        self._removed_strict = removed_strict
        self.pairs = pairs
        self.induced = induced
        self._composition = composition
        self._ambient_scale = ambient_scale
        self.n_programs = n_programs
        self.whole_cell_profiles = whole_cell_profiles
        self.init = init

    def composition(self, source=None, target=None):
        """Per-cell decomposition for one pair (or all pairs).

        Columns: cell_id, source, target, dose (the exposure-derived prior
        contamination fraction), contamination (the fitted fraction of the
        cell's molecules attributed to the source), induced_activity (the
        per-cell multiplier of the pair's induced term), ambient (the
        fitted ambient fraction of the cell).
        """
        d = self._composition
        if source is not None:
            d = d[d["source"] == source]
        if target is not None:
            d = d[d["target"] == target]
        return d.reset_index(drop=True)

    def correct(self, name="generative", retain_induced=True):
        """Derive a correction from the fitted model, without refitting.

        With ``retain_induced=False`` the induced share of each count is
        removed along with the contamination and ambient shares.
        """
        removed = self._removed if retain_induced else self._removed_strict
        corrected = sp.csc_matrix(
            (self._matrix.data - removed, self._matrix.indices.copy(),
             self._matrix.indptr.copy()), shape=self._matrix.shape)
        corrected.data = np.maximum(corrected.data, 0.0)
        rules = self.pairs[["source", "target"]].rename(columns={
            "source": "source_cell_type", "target": "target_cell_type"})
        return GenerativeCorrection(corrected, self._genes, self._cells,
                                    rules, name, retain_induced,
                                    cell_types=self._cell_types)

    def __repr__(self):
        removed = float(self._removed.sum())
        return (f"GenerativeModel(pairs={len(self.pairs)}, "
                f"induced_genes={len(self.induced)}, "
                f"removed_expected={removed:,.0f}, init={self.init!r})")


def fit_generative(audit, *, init=None, n_programs=4, num_threads=None,
                   seed=1, **options):
    """Fit the generative model on the audit's detected pairs.

    ``init`` selects the expression-program initialization: an NMF fit
    object uses its factor-labeled molecules (the default for a
    fit-derived audit — the recommended configuration); ``"clusters"``
    derives programs by clustering each type's cells, weighted toward
    lightly dosed ones (the default for a dataset-derived audit);
    ``"pseudobulk"`` uses one pooled profile per type; a dict mapping
    cell-type names to
    profile matrices (programs x genes, columns aligned with the count
    matrix genes) supplies explicit programs, e.g. from an external
    reference. Extra keyword arguments override the model's constants
    (see ``GenerativeOptions`` in the C++ core).
    """
    from . import _core

    matrix = audit._matrix.tocsc()
    genes = audit._genes
    cells = audit._cells
    types = audit._types
    type_of = {t: i for i, t in enumerate(types)}

    cell_tbl = audit._cells_table
    order = pd.Index(cell_tbl["cell_id"].astype(str))
    pos = order.get_indexer(cells)
    if (pos < 0).any():
        raise ValueError("audit cells missing from the fit cell table")
    x = cell_tbl["x"].to_numpy()[pos]
    y = cell_tbl["y"].to_numpy()[pos]
    ctypes = audit._ctypes
    type_codes = np.array([type_of.get(ctypes.get(c, None), -1)
                           for c in cells], dtype=int)

    pair_specs = []
    pair_names = []
    for key, d in audit._pairs.items():
        if not d["detected"]:
            continue
        pair_specs.append(dict(
            source_type=type_of[d["source"]],
            target_type=type_of[d["target"]],
            pool=[int(g) for g in d["pool"]],
            strict=[int(g) for g in d["strict"]],
        ))
        pair_names.append((d["source"], d["target"]))

    # Resolve the initialization.
    molecules_parquet = ""
    cells_parquet = ""
    programs_flat: list[float] = []
    program_type: list[int] = []
    if init is None:
        init = audit._default_init
    if isinstance(init, str):
        if init not in ("clusters", "pseudobulk"):
            raise ValueError("init must be an NMF fit, 'clusters', "
                             "'pseudobulk', or a dict of program matrices")
        init_mode = init
        init_label = init
    elif isinstance(init, dict):
        init_mode = "pseudobulk"  # fallback for types without programs
        init_label = "explicit"
        for t, mat in init.items():
            if t not in type_of:
                raise ValueError(f"unknown cell type in init: {t}")
            arr = np.asarray(mat, dtype=float)
            if arr.ndim == 1:
                arr = arr[None, :]
            if arr.shape[1] != len(genes):
                raise ValueError("init programs must have one column per gene")
            for row in arr:
                programs_flat.extend(float(v) for v in row)
                program_type.append(type_of[t])
    else:
        init_mode = "factors"
        init_label = "nmf_factors"
        paths = init.manifest["paths"]
        molecules_parquet = str(paths["molecules_parquet"])
        cells_parquet = str(paths["cells_parquet"])

    opts = dict(options)
    opts.setdefault("neighbor_k", int(audit.params["neighbor_k"]))
    opts["init_mode"] = init_mode
    opts["n_programs"] = int(n_programs)
    opts["seed"] = int(seed)
    if num_threads is not None:
        opts["num_threads"] = int(num_threads)

    res = _core.fit_generative(
        counts_indptr=matrix.indptr.astype(np.int64).tolist(),
        counts_indices=matrix.indices.astype(np.int64).tolist(),
        counts_values=matrix.data.astype(np.float64).tolist(),
        n_genes=len(genes),
        cell_ids=[str(c) for c in cells],
        x=np.asarray(x, dtype=float).tolist(),
        y=np.asarray(y, dtype=float).tolist(),
        type_codes=type_codes.tolist(),
        n_types=len(types),
        molecules_parquet=molecules_parquet,
        cells_parquet=cells_parquet,
        pair_specs=pair_specs,
        factor_to_type=[],
        programs_flat=programs_flat,
        program_type=program_type,
        options=opts,
    )

    removed = np.asarray(res["removed"], dtype=float)
    removed_strict = np.asarray(res["removed_without_retention"], dtype=float)
    ambient = np.asarray(res["ambient_scale"], dtype=float)
    comp_rows = []
    for j, (source, target) in enumerate(pair_names):
        cols = np.asarray(res["pair_cells"][j], dtype=int)
        if not len(cols):
            continue
        comp_rows.append(pd.DataFrame(dict(
            cell_id=np.asarray(cells)[cols],
            source=source, target=target,
            dose=np.asarray(res["pair_dose"][j], dtype=float),
            contamination=np.asarray(res["pair_alpha"][j], dtype=float),
            induced_activity=np.asarray(res["pair_rho"][j], dtype=float),
            induced=np.asarray(res["pair_induced"][j], dtype=float),
            ambient=ambient[cols])))
    composition = pd.concat(comp_rows, ignore_index=True) if comp_rows else \
        pd.DataFrame(columns=["cell_id", "source", "target", "dose",
                              "contamination", "induced_activity", "ambient"])

    induced = pd.DataFrame([
        dict(source=pair_names[r["pair"]][0], target=pair_names[r["pair"]][1],
             gene=genes[r["gene"]], excess=r["excess"],
             expected=r["expected"], z=r["z"],
             fold=r["excess"] / max(r["expected"], 1e-9))
        for r in res["induced"]])

    summaries = pd.DataFrame([
        dict(source=pair_names[r["pair"]][0], target=pair_names[r["pair"]][1],
             estimated_dose_molecules=r["prior_molecules"],
             removed_molecules_expected=r["posterior_molecules"],
             induced_molecules=r["induced_molecules"],
             mean_dose_exposed=r["mean_dose_exposed"])
        for r in res["pairs"]])

    n_programs_used = pd.Series(res["n_programs_used"], index=types)

    return GenerativeModel(audit, matrix, genes, cells, removed,
                           removed_strict, summaries, induced, composition,
                           ambient, n_programs_used,
                           bool(res["whole_cell_profiles"]), init_label)


def correct_generative(audit, *, name="generative", **kwargs):
    """Fit the generative model and derive its correction in one call."""
    return fit_generative(audit, **kwargs).correct(name=name)
