"""Generative admixture correction.

Decomposes each cell's counts into its own expression, contamination from
each detected source type, ambient background, and neighborhood-induced
expression, then removes the contamination and ambient shares while
retaining the induced ones. The model, its validation, and the limits of
the separation are described in docs/generative.md. The fit runs in the
C++ core; this module marshals the audit's measurements in and wraps the
result as a correction object that ``audit.evaluate`` accepts.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp


class GenerativeCorrection:
    """Corrected counts plus the model's per-cell decomposition."""

    def __init__(self, audit, matrix, genes, cells, pairs, induced,
                 composition, ambient_scale, whole_cell_profiles, name):
        self._matrix = matrix
        self._genes = genes
        self._cells = cells
        self.name = name
        self.pairs = pairs
        self.induced = induced
        self.whole_cell_profiles = whole_cell_profiles
        self._composition = composition
        self._ambient_scale = ambient_scale
        # rule-style pair listing for symmetry with rule-based corrections
        self.rules = pairs[["source", "target"]].rename(columns={
            "source": "source_cell_type", "target": "target_cell_type"})

    def counts(self):
        """Corrected counts as ``(matrix, genes, cells)``."""
        return self._matrix, list(self._genes), list(self._cells)

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


def correct_generative(audit, *, name="generative", num_threads=None,
                       **options):
    """Fit the generative model on the audit's detected pairs.

    Extra keyword arguments override the model's constants (see
    ``GenerativeOptions`` in the C++ core; e.g. ``use_induced=False`` or
    ``outer_rounds=1``).
    """
    from . import _core

    fit = audit.fit
    matrix = audit._matrix.tocsc()
    genes = audit._genes
    cells = audit._cells
    types = audit._types
    type_of = {t: i for i, t in enumerate(types)}

    cell_tbl = fit.cell_factors()
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

    opts = dict(options)
    opts.setdefault("neighbor_k", int(audit.params["neighbor_k"]))
    if num_threads is not None:
        opts["num_threads"] = int(num_threads)

    paths = fit.manifest["paths"]
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
        molecules_parquet=str(paths["molecules_parquet"]),
        cells_parquet=str(paths["cells_parquet"]),
        pair_specs=pair_specs,
        factor_to_type=[],
        options=opts,
    )

    removed = np.asarray(res["removed"], dtype=float)
    corrected = sp.csc_matrix(
        (matrix.data - removed, matrix.indices.copy(), matrix.indptr.copy()),
        shape=matrix.shape)
    corrected.data = np.maximum(corrected.data, 0.0)

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

    return GenerativeCorrection(
        audit, corrected, genes, cells, summaries, induced, composition,
        ambient, bool(res["whole_cell_profiles"]), name)
