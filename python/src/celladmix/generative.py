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
                         min_share=0.0, n_iter=200,
                         kernels="contact", range_um=(20.0, 60.0),
                         verbose=True):
        """Detect associations between the presence of neighboring cell
        types and the expression of a target type, on the fitted
        expectation. For each target type (or ``target``) and each gene,
        the counts of the target cells are regressed, as a Poisson model
        with non-negative scales, on the cell's own expectation, the
        admixed expectation of every source the fit carries for the type,
        each with its own free per-gene scale, and, for every other
        labeled type, the number of its cells among the cell's nearest
        neighbors (the exposure count); the ambient expectation and the
        floor are held fixed, and the fit's induced terms are left out, so
        that the test reads the response from the counts rather than from
        the screen (``screen_flagged`` says what the screen made of the
        gene).

        The free admixed scales absorb an error of a source profile and a
        remainder of admixture the fit left in the cells, which follows
        that source's neighbors; with every type's neighbors in one
        regression, such a remainder is not credited to another type that
        happens to sit near the source. The fit's admixed profiles are
        zeroed on the target's own genes; here each source's admixed
        expectation is extended to the target-owned genes the source
        expresses on its own (at least a tenth of the target's rate,
        measured in cells away from each other), so that admixture on
        genes both types express is not read as a response. The free
        scale also absorbs a response proportional to a source's dose on
        genes that source carries.

        Each row is a neighboring type (``source``), target and gene: the
        type's rate per neighbor and its z, and the gene is called for
        that type when the z reaches ``z_min`` and the molecules the
        type's terms explain make at least ``min_share`` of the gene's
        content. ``joint_z`` tests all neighboring types together (a Wald
        test over every neighbor term except those of the gene's
        ``owner``, the type with the highest rate, whose association with
        its own gene cannot be told from admixture), as a z of equal tail
        probability; ``joint_called`` marks a gene associated with the
        neighborhood, including where types that co-occur share the
        association and none reaches ``z_min`` alone.

        With ``kernels="both"`` each type also enters with its cells
        between ``range_um[0]`` and ``range_um[1]`` micrometers of the
        cell per 1000 square micrometers, so that a response to cells
        beyond the nearest neighbors loads on it; the table then carries
        ``range_rate``, ``range_se``, ``range_z`` and ``range_share`` (the
        range term's share of the type's response), and a gene is called
        for a type on either term's z."""
        import time
        from scipy import stats

        exposure, types, pos = self._exposure()
        both = kernels == "both"
        if both:
            density, _, _ = self._range(tuple(range_um))
        elif kernels != "contact":
            raise ValueError("kernels must be \"contact\" or \"both\"")
        comp = self._components
        genes = np.asarray(self._genes).astype(str)
        flagged = {(r.source, r.target, r.gene) for r in self.induced.itertuples()} \
            if len(self.induced.columns) else set()
        matrix = self._matrix.tocsc()
        col_type = self._cell_types.reindex(pd.Index(np.asarray(self._cells).astype(str))) \
            .astype(object).to_numpy()
        labeled = [x for x in types if (col_type == x).any()]
        sums = np.column_stack([np.asarray(matrix[:, col_type == x].sum(axis=1)).ravel() for x in labeled])
        owner = np.array(labeled)[np.argmax(sums / np.maximum(sums.sum(axis=0), 1.0), axis=1)]
        sources_of = {}
        for (s, t) in comp["pair_profile"]:
            sources_of.setdefault(t, []).append(s)

        def rate(cols):
            x = np.asarray(matrix[:, cols].sum(axis=1)).ravel()
            return x / max(float(x.sum()), 1.0)

        rows = []
        for t in comp["cells_of"]:
            if target is not None and t != target:
                continue
            E = self.expectation(t)
            cols = E["cells"]
            p = pos[cols]
            Y = np.asarray(matrix[:, cols].todense()).T.astype(float)
            srcs = [s for s in sources_of.get(t, []) if s in E["admixed"]]
            terms = [E["own"]]
            for s in srcs:
                psi = comp["pair_profile"][(s, t)]
                admixed = E["admixed"][s]
                s_cols = np.flatnonzero(col_type == s)
                s_far = s_cols[exposure[pos[s_cols], types.index(t)] == 0]
                t_far = cols[exposure[p, types.index(s)] == 0]
                if len(s_far) >= 50 and len(t_far) >= 50 and psi.sum() > 0:
                    r_s, r_t = rate(s_far), rate(t_far)
                    shared = (owner == t) & (psi <= 0) & (r_s > 0) & (r_s >= _SHARED_MIN * r_t)
                    if shared.any():
                        full = psi.copy()
                        full[shared] = r_s[shared]
                        admixed = np.outer(admixed.sum(axis=1) / psi.sum(), full)
                terms.append(admixed)
            names, C = [], []
            for x in types:
                if x == t:
                    continue
                for kind, values in [("contact", exposure[p, types.index(x)].astype(float))] + \
                        ([("range", density[p, types.index(x)])] if both else []):
                    if values.sum() > 0:
                        names.append((x, kind))
                        C.append(values)
            if not C:
                continue
            C = np.column_stack(C)
            t0 = time.time()
            coef, cov = _poisson_terms(Y, E["ambient"] + E["floor"], terms, C, n_iter)
            K = len(terms)
            if verbose:
                print(f"[detect_responses] {t}: {len(cols)} cells, {len(srcs)} admixed terms, "
                      f"{C.shape[1]} neighbor terms, {time.time() - t0:.0f}s", flush=True)
            se = np.sqrt(np.maximum(np.stack([cov[:, k, k] for k in range(K + C.shape[1])]), 0.0))
            b, b_se = coef[K:], se[K:]
            z = np.where(b_se > 0, b / np.where(b_se > 0, b_se, 1.0), np.nan)
            response = b * C.sum(axis=0)[:, None]
            observed = Y.sum(axis=0)
            of_type = np.array([x for x, _ in names])
            joint_stat = np.zeros(len(genes))
            joint_df = np.zeros(len(genes), dtype=int)
            joint_response = np.zeros(len(genes))
            for g in range(len(genes)):
                sel = np.flatnonzero(of_type != owner[g])
                if len(sel):
                    joint_stat[g] = float(b[sel, g] @ np.linalg.pinv(cov[g][np.ix_(K + sel, K + sel)]) @ b[sel, g])
                    joint_df[g] = len(sel)
                    joint_response[g] = response[sel, g].sum()
            p_joint = np.where(joint_df > 0, stats.chi2.sf(joint_stat, np.maximum(joint_df, 1)), 1.0)
            joint_z = stats.norm.isf(np.clip(p_joint, 1e-300, 1.0))
            joint_called = (joint_df > 0) & (joint_z >= z_min) & \
                (joint_response / np.maximum(observed, 1e-9) >= min_share)
            for x in dict.fromkeys(of_type):
                if source is not None and x != source:
                    continue
                ic = names.index((x, "contact")) if (x, "contact") in names else None
                ir = names.index((x, "range")) if (x, "range") in names else None
                k_adm = 1 + srcs.index(x) if x in srcs else None
                for g in range(len(genes)):
                    r_c = response[ic, g] if ic is not None else 0.0
                    r_r = response[ir, g] if ir is not None else 0.0
                    total = r_c + r_r
                    share = total / max(observed[g], 1e-9)
                    z_c = z[ic, g] if ic is not None else np.nan
                    row = dict(
                        source=x, target=t, gene=genes[g], observed=float(observed[g]),
                        own_fit=float(coef[0, g] * E["own"][:, g].sum()),
                        admixed_fit=float(coef[k_adm, g] * terms[k_adm][:, g].sum()) if k_adm else 0.0,
                        contact_rate=float(b[ic, g]) if ic is not None else 0.0,
                        contact_se=float(b_se[ic, g]) if ic is not None else np.nan, z=float(z_c),
                        response_molecules=float(total), response_share=float(share),
                        screen_flagged=(x, t, genes[g]) in flagged)
                    significant = bool(np.isfinite(z_c) and z_c >= z_min)
                    if both:
                        z_r = z[ir, g] if ir is not None else np.nan
                        row.update(range_rate=float(b[ir, g]) if ir is not None else 0.0,
                                   range_se=float(b_se[ir, g]) if ir is not None else np.nan,
                                   range_z=float(z_r), range_share=float(r_r / max(total, 1e-9)),
                                   range_in_um=float(range_um[0]), range_out_um=float(range_um[1]))
                        significant = significant or bool(np.isfinite(z_r) and z_r >= z_min)
                    row.update(called=significant and share >= min_share, owner=str(owner[g]),
                               joint_z=float(joint_z[g]), joint_called=bool(joint_called[g]))
                    rows.append(row)
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

    def retain_responses(self, responses):
        """Keep the called responses: for each called neighboring type,
        target and gene in ``responses`` (a table from
        :meth:`detect_responses`), the response molecules in every target
        cell, the fitted rate per neighbor times the cell's count of that
        type among its nearest neighbors (and the range term likewise),
        are taken out of the cell's removal of that gene, and nothing else
        changes.
        Returns a new model whose corrections carry the retained
        responses; the removal without retention is unchanged."""
        from scipy import sparse

        called = responses[responses["called"]] if "called" in responses.columns else responses
        removed = self._removed.copy()
        matrix = self._matrix
        position = sparse.csc_matrix(
            (np.arange(matrix.nnz) + 1, matrix.indices, matrix.indptr), shape=matrix.shape)
        gene_index = pd.Index(np.asarray(self._genes).astype(str))
        exposure, types, pos = self._exposure()
        with_range = "range_rate" in called.columns
        for (s, t), rows in called.groupby(["source", "target"]):
            if t not in self._components.get("cells_of", {}) or s not in types:
                continue
            cols = self._components["cells_of"][t]
            contact = exposure[pos[cols], types.index(s)].astype(float)
            if with_range:
                r0 = rows.iloc[0]
                density, _, _ = self._range((float(r0.range_in_um), float(r0.range_out_um)))
                at_range = density[pos[cols], types.index(s)]
            for r in rows.itertuples():
                g = gene_index.get_loc(str(r.gene))
                response = r.contact_rate * contact
                if with_range:
                    response = response + r.range_rate * at_range
                p = np.asarray(position[g, cols].todense()).ravel() - 1
                ok = p >= 0
                removed[p[ok]] = np.maximum(removed[p[ok]] - response[ok], 0.0)
        model = GenerativeModel(self._audit, self._matrix, self._genes, self._cells,
                                removed, self._removed_strict, self.pairs, self.induced,
                                self._composition, self._ambient_scale, self.n_programs,
                                self.whole_cell_profiles, self.init,
                                components=self._components)
        model.retained_responses = called.reset_index(drop=True)
        return model

    def _range(self, range_um):
        """Source cells between range_um[0] and range_um[1] micrometers of
        each cell per 1000 square micrometers, per column of the count
        matrix: the densities, the types and the row of each column in the
        cell table."""
        from ._score_utils import source_range_density

        cache = getattr(self, "_range_cache", {})
        key = (float(range_um[0]), float(range_um[1]))
        if key not in cache:
            audit = self._audit
            cells_table = audit._cells_table
            ann = audit._ctypes[audit._ctypes != "None"]
            density, types = source_range_density(cells_table, ann, r_in=key[0], r_out=key[1])
            pos = pd.Index(cells_table["cell_id"].astype(str)).get_indexer(
                [str(c) for c in self._cells])
            cache[key] = (density, types, pos)
            self._range_cache = cache
        return cache[key]

    def _exposure(self):
        """Source cells among each cell's nearest neighbors (the exposure
        count), per column of the count matrix: the counts, the types and
        the row of each column in the cell table."""
        from ._score_utils import source_exposure_counts

        if not hasattr(self, "_exposure_cache"):
            audit = self._audit
            cells_table = audit._cells_table
            ann = audit._ctypes[audit._ctypes != "None"]
            exposure, types = source_exposure_counts(
                cells_table, ann, neighbor_k=self._components["neighbor_k"])
            pos = pd.Index(cells_table["cell_id"].astype(str)).get_indexer(
                [str(c) for c in self._cells])
            self._exposure_cache = (exposure, types, pos)
        return self._exposure_cache

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


# Target-owned genes enter a source's admixed expectation in detection when
# the source expresses them at least at this fraction of the target's rate.
_SHARED_MIN = 0.1


def _poisson_terms(Y, fixed, terms, covariates, n_iter=200, tol=1e-6):
    """Non-negative Poisson regression of every gene at once:
    E[y] = fixed + sum_k c_k * terms[k] + sum_j b_j * covariate_j, by
    multiplicative updates; ``terms`` are cells x genes expectations, each
    with its own scale per gene, and ``covariates`` a cells x m matrix of
    per-cell terms shared by all genes. Returns the scales ((k + m) x
    genes) and their covariance per gene (genes x (k + m) x (k + m)) from
    the Fisher information."""
    n_genes = Y.shape[1]
    C = np.asarray(covariates, dtype=float)
    if C.ndim == 1:
        C = C[:, None]
    K, m = len(terms), C.shape[1]
    coef = np.full((K + m, n_genes), 0.1)
    col_t = [T.sum(axis=0) for T in terms]
    col_c = C.sum(axis=0)

    def predict(coef):
        mu = fixed + C @ coef[K:]
        for k, T in enumerate(terms):
            mu = mu + T * coef[k]
        return mu

    for _ in range(n_iter):
        ratio = Y / np.maximum(predict(coef), 1e-12)
        new = np.zeros_like(coef)
        for k, T in enumerate(terms):
            new[k] = np.where(col_t[k] > 0, coef[k] * (T * ratio).sum(axis=0) / np.maximum(col_t[k], 1e-12), 0.0)
        new[K:] = np.where(col_c[:, None] > 0, coef[K:] * (C.T @ ratio) / np.maximum(col_c[:, None], 1e-12), 0.0)
        change = np.max(np.abs(new - coef) / np.maximum(coef, 1e-8))
        coef = new
        if change < tol:
            break
    w = 1.0 / np.maximum(predict(coef), 1e-12)
    k = K + m
    info = np.zeros((n_genes, k, k))
    for a in range(k):
        for b in range(a, k):
            if b < K:
                v = (terms[a] * terms[b] * w).sum(axis=0)
            elif a < K:
                v = (terms[a] * (C[:, b - K:b - K + 1] * w)).sum(axis=0)
            else:
                v = (C[:, a - K] * C[:, b - K]) @ w
            info[:, a, b] = info[:, b, a] = v
    return coef, np.linalg.pinv(info)


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
             cross_attributed_molecules=r["cross_attributed_molecules"],
             prior_strength=r["prior_strength"])
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
