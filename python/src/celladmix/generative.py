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
                 whole_cell_profiles, init, components=None):
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
        self._components = components or {}
        self._expectations = {}

    def expectation(self, target):
        """The model's expectation per cell and gene for the cells of one
        target type, in molecules, as dense arrays (cells x genes):
        ``own`` (the cell's expression programs), ``admixed`` (a dict by
        source type), ``ambient``, ``induced`` (a dict by source type,
        only sources with induced genes) and ``floor`` (the uniform
        floor), with ``cells`` the column indices of the cells and
        ``totals`` their molecule counts."""
        if target in self._expectations:
            return self._expectations[target]
        comp = self._components
        if target not in comp.get("programs", {}):
            raise ValueError(f"no fitted pairs into {target!r}")
        n_genes = len(self._genes)
        cols = comp["cells_of"][target]
        totals = np.asarray(self._matrix[:, cols].sum(axis=0)).ravel()
        cell_ids = pd.Index(np.asarray(self._cells)[cols])
        own = (comp["weights"][target] @ comp["programs"][target]) * totals[:, None]
        admixed, induced = {}, {}
        for (source, t), psi in comp["pair_profile"].items():
            if t != target:
                continue
            d = self._composition[(self._composition["source"] == source)
                                  & (self._composition["target"] == t)]
            d = d.set_index("cell_id").reindex(cell_ids)
            alpha = np.nan_to_num(d["contamination"].to_numpy(dtype=float))
            admixed[source] = np.outer(alpha * totals, psi)
            mi = comp["pair_induced_profile"][(source, t)]
            if mi.sum() > 0:
                ind = np.nan_to_num(d["induced"].to_numpy(dtype=float))
                induced[source] = np.outer(ind * totals, mi)
        ambient = np.outer(self._ambient_scale[cols] * totals,
                           comp["ambient_profile"].get(target, np.zeros(n_genes)))
        floor = np.full(own.shape, comp["floor_total"] / n_genes) * totals[:, None]
        out = dict(cells=cols, totals=totals, own=own, admixed=admixed,
                   ambient=ambient, induced=induced, floor=floor)
        self._expectations[target] = out
        return out

    def detect_responses(self, source=None, target=None, z_min=6.0,
                         min_share=0.0, n_iter=200, copy_scale="free",
                         verbose=True):
        """Detect neighborhood-induced responses with a contact term on the
        fitted expectation. For each pair (or the pairs selected by
        ``source`` and ``target``) and each gene, the counts of the target
        cells are regressed, as a Poisson model with non-negative scales,
        on the cell's own expectation, the source's admixed expectation and
        the number of source cells among the cell's nearest neighbors (the
        exposure count), with the other sources' admixed and induced
        expectations, the ambient expectation and the floor held fixed; a
        gene responds when the contact term's z reaches ``z_min`` and the
        molecules it explains make at least ``min_share`` of the gene's
        content. The induced term the fit itself carries for the pair is
        left out, so that the test reads the response from the counts
        rather than from the screen; ``screen_flagged`` says what the
        screen made of the gene. With ``copy_scale="free"`` the admixed
        expectation enters with a free per-gene scale, which absorbs a
        contact-proportional response on genes the source profile carries;
        with ``"fixed"`` it enters as fitted, so such a response loads on
        the contact term and a profile error does too. Returns one row
        per pair and gene."""
        import time
        from ._score_utils import source_exposure_counts

        audit = self._audit
        cells_table = audit._cells_table
        ann = audit._ctypes[audit._ctypes != "None"]
        exposure, types = source_exposure_counts(
            cells_table, ann, neighbor_k=self._components["neighbor_k"])
        pos = pd.Index(cells_table["cell_id"].astype(str)).get_indexer(
            [str(c) for c in self._cells])
        genes = np.asarray(self._genes)
        flagged = {(r.source, r.target, r.gene) for r in self.induced.itertuples()} \
            if len(self.induced.columns) else set()
        rows = []
        for (s, t) in self._components["pair_profile"]:
            if (source is not None and s != source) or (target is not None and t != target):
                continue
            E = self.expectation(t)
            cols = E["cells"]
            Y = np.asarray(self._matrix[:, cols].todense()).T.astype(float)
            contact = exposure[pos[cols], types.index(s)].astype(float)
            fixed = E["ambient"] + E["floor"]
            for s2, A in E["admixed"].items():
                if s2 != s:
                    fixed = fixed + A
            for s2, I in E["induced"].items():
                if s2 != s:
                    fixed = fixed + I
            t0 = time.time()
            if copy_scale == "fixed":
                scales, se = _poisson_scales(Y, fixed + E["admixed"][s], E["own"],
                                             np.zeros_like(E["own"]), contact, n_iter)
                a, b = scales[0], scales[2]
                gamma = np.ones(len(genes))
            else:
                scales, se = _poisson_scales(Y, fixed, E["own"], E["admixed"][s], contact, n_iter)
                a, gamma, b = scales
            if verbose:
                print(f"[detect_responses] {s} -> {t}: {len(cols)} cells, "
                      f"{time.time() - t0:.0f}s", flush=True)
            z = np.where(se[2] > 0, b / np.where(se[2] > 0, se[2], 1.0), np.nan)
            observed = Y.sum(axis=0)
            response = b * contact.sum()
            share = response / np.maximum(observed, 1e-9)
            for g in range(len(genes)):
                rows.append(dict(
                    source=s, target=t, gene=str(genes[g]), observed=float(observed[g]),
                    own_fit=float(a[g] * E["own"][:, g].sum()),
                    admixed_fit=float(gamma[g] * E["admixed"][s][:, g].sum()),
                    contact_rate=float(b[g]), contact_se=float(se[2][g]), z=float(z[g]),
                    response_molecules=float(response[g]), response_share=float(share[g]),
                    screen_flagged=(s, t, str(genes[g])) in flagged,
                    called=bool(np.isfinite(z[g]) and z[g] >= z_min and share[g] >= min_share)))
        return pd.DataFrame(rows)

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

    def realize(self, fit, name="generative", seed=0, retain_induced=True):
        """A molecule-level correction: the model's removal in every cell
        and gene, rounded stochastically, drawn at random from the cell's
        molecules of that gene and written as a corrected store under the
        fit's run (``corrected/<name>``, the layout of the factor-based
        cleanup). ``fit`` is the NMF fit the audit was built from. Within
        one gene in one cell nothing tells a transferred molecule from an
        own one on stacked pairs, and position adds under a point of
        precision on lateral ones, so the draw is random."""
        from .correction import realize_correction

        return realize_correction(fit, self, name=name, seed=seed,
                                  retain_induced=retain_induced)

    def __repr__(self):
        removed = float(self._removed.sum())
        return (f"GenerativeModel(pairs={len(self.pairs)}, "
                f"induced_genes={len(self.induced)}, "
                f"removed_expected={removed:,.0f}, init={self.init!r})")


def _poisson_scales(Y, fixed, own, copy, contact, n_iter=200, tol=1e-6):
    """Non-negative Poisson regression of every gene at once:
    E[y] = fixed + a * own + gamma * copy + b * contact, by multiplicative
    updates; returns the scales (3 x genes) and their standard errors from
    the Fisher information."""
    n_genes = Y.shape[1]
    a = np.full(n_genes, 0.1)
    gamma = np.full(n_genes, 0.1)
    b = np.full(n_genes, 0.1)
    col_own = own.sum(axis=0)
    col_copy = copy.sum(axis=0)
    col_contact = float(contact.sum())
    contact_col = contact[:, None]
    for _ in range(n_iter):
        mu = fixed + own * a + copy * gamma + contact_col * b
        ratio = Y / np.maximum(mu, 1e-12)
        a_new = np.where(col_own > 0, a * (own * ratio).sum(axis=0) / np.maximum(col_own, 1e-12), 0.0)
        g_new = np.where(col_copy > 0, gamma * (copy * ratio).sum(axis=0) / np.maximum(col_copy, 1e-12), 0.0)
        b_new = (b * (contact @ ratio) / col_contact) if col_contact > 0 else np.zeros(n_genes)
        change = max(np.max(np.abs(a_new - a) / np.maximum(a, 1e-8)),
                     np.max(np.abs(g_new - gamma) / np.maximum(gamma, 1e-8)),
                     np.max(np.abs(b_new - b) / np.maximum(b, 1e-8)))
        a, gamma, b = a_new, g_new, b_new
        if change < tol:
            break
    mu = fixed + own * a + copy * gamma + contact_col * b
    w = 1.0 / np.maximum(mu, 1e-12)
    info = np.zeros((n_genes, 3, 3))
    info[:, 0, 0] = (own * own * w).sum(axis=0)
    info[:, 0, 1] = info[:, 1, 0] = (own * copy * w).sum(axis=0)
    info[:, 0, 2] = info[:, 2, 0] = (own * contact_col * w).sum(axis=0)
    info[:, 1, 1] = (copy * copy * w).sum(axis=0)
    info[:, 1, 2] = info[:, 2, 1] = (copy * contact_col * w).sum(axis=0)
    info[:, 2, 2] = (contact_col * contact_col * w).sum(axis=0)
    cov = np.linalg.pinv(info)
    se = np.sqrt(np.maximum(np.stack([cov[:, 0, 0], cov[:, 1, 1], cov[:, 2, 2]]), 0.0))
    return np.stack([a, gamma, b]), se


def fit_generative(audit, *, init=None, n_programs=4, num_threads=None,
                   seed=1, verbose=True, **options):
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
    reference. ``verbose`` prints timestamped ``[INFO ...]`` progress
    lines to standard output while the fit runs: the stages (spatial
    covariates, molecule-table preparation, program initialization,
    interface-local source profiles), then per outer round and target
    type the cell, source, and count sizes, the dose prior, the
    induced-expression screen, each expectation-maximization pass with
    its iteration count, log-likelihood, and relative change, each
    refinement of the dose prior, and the molecules attributed to
    contamination. Within long loops a line is printed at most every
    ``progress_interval_sec`` seconds (default 5). Each line is flushed
    as it is written, so a redirected log follows the fit. Extra keyword
    arguments override the model's constants (see ``GenerativeOptions``
    in the C++ core).
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
    opts["verbose"] = bool(verbose)
    opts["type_names"] = [str(t) for t in types]
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
             mean_dose_exposed=r["mean_dose_exposed"],
             prior_unexposed_molecules=r["prior_unexposed_molecules"],
             cross_attributed_molecules=r["cross_attributed_molecules"])
        for r in res["pairs"]])

    n_programs_used = pd.Series(res["n_programs_used"], index=types)

    # the fitted components, for the model's expectation per cell and gene
    n_genes = len(genes)
    programs, weights, cells_of, ambient_profile = {}, {}, {}, {}
    for ti, t in enumerate(types):
        progs = np.asarray(res["type_programs"][ti], dtype=float)
        if progs.size == 0:
            continue
        programs[t] = progs.reshape(-1, n_genes)
        cols = np.flatnonzero(type_codes == ti)
        cells_of[t] = cols
        weights[t] = np.array([res["cell_weights"][c] for c in cols], dtype=float).reshape(len(cols), -1)
        ambient_profile[t] = np.asarray(res["type_ambient_profile"][ti], dtype=float)
    components = dict(
        programs=programs, weights=weights, cells_of=cells_of,
        ambient_profile=ambient_profile,
        pair_profile={pair_names[j]: np.asarray(res["pair_profile"][j], dtype=float)
                      for j in range(len(pair_names)) if len(res["pair_profile"][j])},
        pair_induced_profile={pair_names[j]: np.asarray(res["pair_induced_profile"][j], dtype=float)
                              for j in range(len(pair_names)) if len(res["pair_induced_profile"][j])},
        floor_total=float(opts.get("floor_total", 0.002)),
        neighbor_k=int(opts["neighbor_k"]))

    return GenerativeModel(audit, matrix, genes, cells, removed,
                           removed_strict, summaries, induced, composition,
                           ambient, n_programs_used,
                           bool(res["whole_cell_profiles"]), init_label,
                           components=components)


def correct_generative(audit, *, name="generative", **kwargs):
    """Fit the generative model and derive its correction in one call."""
    return fit_generative(audit, **kwargs).correct(name=name)
