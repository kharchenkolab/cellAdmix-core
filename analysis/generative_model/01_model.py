"""Generative admixture model: per-target-type Poisson mixture over genes.

For each target cell c of annotated type T, the molecule counts over genes
are modeled as a Poisson mixture of four components:

1. own expression: a mixture of state-specific expression programs F_k
   (initialized from the profiles of factor-labeled molecules of the type's
   own factors, refined by expectation-maximization);
2. contamination from each source type S with a detected pair S -> T:
   proportional to S's cytoplasmic expression profile psi_S (estimated from
   source-cell molecules outside the nucleus, iteratively decontaminated),
   scaled by a per-cell contamination fraction alpha_cS. The prior mean of
   alpha_cS is an isotonic (monotone) dose-response of guide-gene content on
   the number of S cells among the target's 15 nearest neighbors, converted
   to a whole-profile fraction; each cell's expression updates its own alpha
   through a Gamma-posterior step, bounded to a 3-fold factor of the prior
   mean so that removal remains anchored to the measured exposure signal;
3. an ambient background restricted to the strict genes (near-zero native
   baseline in the target), with profile and scale estimated from target
   cells far (>250 um) from any source cell (production configuration only);
4. a sparse per-gene induced-expression term: genes whose exposure-linked
   excess significantly exceeds the source-profile-proportional expectation
   (a per-gene z test, the sparsity mechanism) get a per-gene rate m_Sg with
   a per-cell activity multiplier, and that content is RETAINED rather than
   removed.

A small uniform floor keeps posterior fractions well-defined where own
expression is zero. The observed count Y_cg is then split by posterior
component responsibility; the corrected matrix keeps the own + induced +
floor share of each entry.

Arms (--arm):
  validation      dose evidence restricted to the A halves of the marker
                  pools and to genes outside the held-out B halves; no
                  ambient component; induced support restricted to non-B
                  genes. Honest gene-split power test.
  production      full pools for dose, ambient on, own-program support
                  zeroed on strict genes, induced term on.
  production_noind  production with the induced term disabled (comparison
                  against the regression corrector, which cannot retain
                  induced genes).
  shuffled        production dose machinery with each pair's exposure values
                  permuted across target cells; ambient off (any
                  exposure-independent removal would register as power by
                  construction, matching the audit's control convention).

Outputs: a MatrixMarket matrix of removed molecule mass (gene x cell, same
dimensions as the exported count matrix) plus per-pair diagnostic tables.
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scipy.io
import scipy.sparse as sp
from scipy.spatial import cKDTree

GM = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
DATASET = os.environ.get("GM_DATASET", "pancreas")
_SUF = "" if DATASET == "pancreas" else f"_{DATASET}"
PREFIX = "" if DATASET == "pancreas" else f"{DATASET}_"
DATA = os.path.join(GM, "data" + _SUF)
RESULTS = os.path.join(GM, "results")
_EX = "/home/pkharchenko/cellAdmix/cellAdmix-core/examples"
RUN = {
    "pancreas": f"{_EX}/xenium_pancreas_membrane_377_full/out/runs/"
                "fit_manual_rank9_invsqrt_kl",
    "nsclc": f"{_EX}/cosmx_nsclc_giotto/out/runs/fit_manual_rank8_ls_nmf",
    "breast": f"{_EX}/xenium_breast_membrane_5k_full/out/runs/"
              "fit_manual_rank8_invsqrt_kl",
}[DATASET]

# ---- model hyperparameters (fixed before evaluation; see REPORT.md) ----
C0_ALPHA = 300.0      # prior strength for alpha (pseudo-molecules)
ALPHA_CAP = 5.0       # per-cell alpha may move at most 5-fold above its prior mean
TOPUP_PASSES = 3      # residual-gradient top-up passes after the main EM
C0_AMBIENT = 300.0    # prior strength for the ambient scale
AMBIENT_CAP = 5.0     # per-cell ambient scale cap (fold over the far-field prior)
EPS_TOTAL = 0.002     # uniform floor: fraction of cell content spread over genes
LAMBDA_MAX = 0.6      # cap on the prior contamination fraction
IND_Z = 8.0           # threshold on the overdispersed z of the induced screen
IND_MIN_EXCESS = 50.0 # minimum exposure-linked excess (molecules) to consider
PROF_CV = 0.15        # multiplicative uncertainty of the source profile: the
                      # screen's variance includes (PROF_CV * proportional
                      # expectation)^2, so mild disproportionality on large
                      # channels reads as profile noise, not induction
NEAR_UM = 30.0        # source cells within this distance of a target-type
                      # cell define the pair-level contamination profile
NEAR_MIN_MOL = 30000  # minimum cytoplasmic molecules for the near profile
DOSE_W = 20.0         # own-program learning downweights dosed cells by
                      # 1 / (1 + DOSE_W * prior contamination fraction)
RHO_SHAPE = 2.0       # Gamma prior shape for the per-cell induced activity
RHO_CAP = 30.0
N_EM = 40             # EM iterations per target type
N_OUTER = 2           # outer rounds of source-profile decontamination
FAR_UM = 250.0        # far-field distance for the ambient reference
MIN_FAR_MOL = 20000   # minimum far-field molecules to trust the estimate


def pava_increasing(y, w):
    """Weighted pool-adjacent-violators: monotone non-decreasing fit."""
    y = list(map(float, y)); w = list(map(float, w))
    val, wt, idx = [], [], []
    for i, (yi, wi) in enumerate(zip(y, w)):
        val.append(yi); wt.append(wi); idx.append([i])
        while len(val) > 1 and val[-2] > val[-1] + 1e-15:
            v = (val[-2] * wt[-2] + val[-1] * wt[-1]) / (wt[-2] + wt[-1])
            wt[-2] += wt[-1]; val[-2] = v
            idx[-2].extend(idx[-1])
            val.pop(); wt.pop(); idx.pop()
    out = np.empty(len(y))
    for v, ii in zip(val, idx):
        out[ii] = v
    return out


class Inputs:
    """Load everything exported by 00_export.R plus molecule-level tables."""

    def __init__(self):
        self.genes = [l.strip() for l in open(os.path.join(DATA, "genes.txt"))]
        self.cells = [l.strip() for l in open(os.path.join(DATA, "cells.txt"))]
        self.G = len(self.genes)
        self.gene_of = {g: i for i, g in enumerate(self.genes)}
        self.col_of = {c: i for i, c in enumerate(self.cells)}
        print("reading counts ...", flush=True)
        self.counts = sp.csc_matrix(scipy.io.mmread(os.path.join(DATA, "counts.mtx")))
        assert self.counts.shape == (self.G, len(self.cells))
        self.totals = np.asarray(self.counts.sum(axis=0)).ravel()

        run = json.load(open(os.path.join(RUN, "run.json")))
        assert run["genes"] == self.genes, "gene order mismatch run.json vs export"

        ct = pd.read_csv(os.path.join(DATA, "cell_types.csv"))
        self.cell_type = dict(zip(ct.cell_id, ct.cell_type))
        self.pairs = pd.read_csv(os.path.join(DATA, "pairs.csv"))
        self.expo = pd.read_csv(os.path.join(DATA, "exposure.csv.gz"))
        prof = pd.read_csv(os.path.join(DATA, "type_profiles_cpm.csv"))
        self.type_names = [c for c in prof.columns if c != "gene"]
        P = prof[self.type_names].to_numpy()
        order = [self.gene_of[g] for g in prof.gene]
        self.profiles_cpm = np.zeros((self.G, len(self.type_names)))
        self.profiles_cpm[order, :] = P
        self.top_type = np.array(self.type_names)[np.argmax(self.profiles_cpm, 1)]

        cells_tbl = pq.read_table(os.path.join(RUN, "cells.parquet"),
            columns=["cell_idx", "cell_id", "x", "y"]).to_pandas()
        self.cell_xy = np.full((len(self.cells), 2), np.nan)
        keep = cells_tbl.cell_id.isin(self.col_of)
        cc = cells_tbl[keep]
        self.cell_xy[[self.col_of[c] for c in cc.cell_id]] = cc[["x", "y"]].to_numpy()
        self.idx_of_run_cell = dict(zip(cells_tbl.cell_idx, cells_tbl.cell_id))

        print("reading molecules ...", flush=True)
        mol_path = os.path.join(RUN, "molecules.parquet")
        mol_cols = pq.ParquetFile(mol_path).schema_arrow.names
        has_nucleus = "overlaps_nucleus" in mol_cols
        mol = pq.read_table(mol_path,
            columns=["gene_idx", "cell_idx", "factor_label"]
                    + (["overlaps_nucleus"] if has_nucleus else [])).to_pandas()
        run_cell_type = np.array([
            self.cell_type.get(self.idx_of_run_cell.get(i), None) or "NA"
            for i in range(int(cells_tbl.cell_idx.max()) + 1)])
        mol_type = run_cell_type[mol.cell_idx.to_numpy()]

        # Cytoplasmic (non-nuclear) molecule counts per gene and cell, as a
        # sparse matrix: any cell subset's cytoplasmic profile is a column
        # slice away. Falls back to all molecules when the platform provides
        # no nucleus-overlap flag (e.g. this CosMx export).
        if has_nucleus:
            cyto = mol[mol.overlaps_nucleus.to_numpy() == 0]
        else:
            print("no nucleus-overlap flag: whole-cell profiles used")
            cyto = mol
        cyto_cols = np.array([self.col_of.get(self.idx_of_run_cell.get(i), -1)
                              for i in range(int(cells_tbl.cell_idx.max()) + 1)])
        cc = cyto_cols[cyto.cell_idx.to_numpy()]
        ok = cc >= 0
        self.cyto_counts = sp.coo_matrix(
            (np.ones(ok.sum()), (cyto.gene_idx.to_numpy()[ok], cc[ok])),
            shape=(self.G, len(self.cells))).tocsc()

        # Raw cytoplasmic profiles per type.
        cyto_type = run_cell_type[cyto.cell_idx.to_numpy()]
        self.psi_raw = {}
        for t in self.type_names:
            cnt = np.bincount(cyto.gene_idx.to_numpy()[cyto_type == t],
                              minlength=self.G).astype(float)
            self.psi_raw[t] = cnt / max(cnt.sum(), 1.0)
        self._near_profile_cache = {}

        # Own-program initializations: profiles of factor-labeled molecules
        # within each type, for the type's own factors (aligned factor plus
        # unaligned factors carrying >=5% of the type's content). 0-based
        # factor labels; factor_types.csv is 1-based.
        ft_path = ("/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/"
                   f"audit_guided/results/factor_types{_SUF}.csv")
        if os.path.exists(ft_path):
            ft = pd.read_csv(ft_path)
            aligned = {int(r.factor) - 1: r.cell_type for r in ft.itertuples()}
        else:
            # Derive the factor-to-type alignment from the data: each
            # factor's gene profile (from its labeled molecules) is compared
            # with the type pseudobulk profiles; a factor aligns to the best
            # match when it is clearly ahead of the second best.
            nf = int(mol.factor_label.max()) + 1
            aligned = {}
            tp = np.sqrt(self.profiles_cpm
                         / np.maximum(self.profiles_cpm.sum(0), 1.0))
            for k in range(nf):
                cnt = np.bincount(
                    mol.gene_idx.to_numpy()[mol.factor_label.to_numpy() == k],
                    minlength=self.G).astype(float)
                if cnt.sum() < 5000:
                    continue
                fk = np.sqrt(cnt / cnt.sum())
                cos = fk @ tp / max(np.linalg.norm(fk), 1e-9)                     / np.maximum(np.linalg.norm(tp, axis=0), 1e-9)
                o = np.argsort(-cos)
                if cos[o[0]] >= 0.3 and cos[o[0]] >= 1.15 * cos[o[1]]:
                    aligned[k] = self.type_names[o[0]]
            pd.DataFrame(dict(factor=[k + 1 for k in aligned],
                cell_type=list(aligned.values()))).to_csv(ft_path, index=False)
            print("derived factor alignment:",
                  {k + 1: v for k, v in aligned.items()})
        own_factors = {t: [k for k, tt in aligned.items() if tt == t]
                       for t in self.type_names}
        n_fac = int(mol.factor_label.max()) + 1
        for t in self.type_names:
            sub = mol[mol_type == t]
            shares = np.bincount(sub.factor_label.to_numpy(),
                                 minlength=n_fac).astype(float)
            shares /= max(shares.sum(), 1.0)
            for k in range(n_fac):
                if k not in aligned and shares[k] >= 0.05:
                    own_factors[t].append(k)
        self.program_init = {}
        for t in self.type_names:
            sub = mol[mol_type == t]
            progs = []
            for k in own_factors[t]:
                cnt = np.bincount(
                    sub.gene_idx.to_numpy()[sub.factor_label.to_numpy() == k],
                    minlength=self.G).astype(float)
                if cnt.sum() >= 5000:
                    progs.append(cnt / cnt.sum())
            if not progs:  # e.g. Mural/pericyte has no own factor
                cnt = np.bincount(sub.gene_idx.to_numpy(),
                                  minlength=self.G).astype(float)
                progs = [cnt / max(cnt.sum(), 1.0)]
            self.program_init[t] = np.vstack(progs)
        del mol, cyto

        # Per-pair parsed gene sets and per-cell exposure.
        self.pair_info = {}
        for r in self.pairs.itertuples():
            e = self.expo[self.expo.pair == r.pair]
            cols = np.array([self.col_of[c] for c in e.cell_id])
            self.pair_info[r.pair] = dict(
                S=r.source, T=r.target,
                pool=[g for g in r.pool.split(";") if g],
                A=[g for g in r.A.split(";") if g],
                B=[g for g in r.B.split(";") if g],
                strict=[g for g in str(r.strict).split(";") if g and g != "nan"],
                cells=cols, exposure=e.exposure.to_numpy().astype(int))

        self.type_cols = {t: np.array([self.col_of[c] for c in self.cells
                                       if self.cell_type.get(c) == t])
                          for t in self.type_names}

    def gidx(self, names):
        return np.array([self.gene_of[g] for g in names if g in self.gene_of],
                        dtype=int)

    def near_source_profile(self, S, T):
        """Cytoplasmic profile of the S cells that actually border T cells
        (within NEAR_UM of the nearest T cell), so that spatial composition
        shifts within the source type (e.g. activated fibroblasts near
        tumor) are part of the contamination profile rather than appearing
        as induced expression in the target. Falls back to a wider radius
        and then to all S cells when the near subset is too small."""
        key = (S, T)
        if key in self._near_profile_cache:
            return self._near_profile_cache[key]
        s_cols = self.type_cols[S]
        t_xy = self.cell_xy[self.type_cols[T]]
        t_xy = t_xy[np.isfinite(t_xy).all(axis=1)]
        s_xy = self.cell_xy[s_cols]
        oks = np.isfinite(s_xy).all(axis=1)
        d = np.full(len(s_cols), np.inf)
        if len(t_xy):
            d[oks] = cKDTree(t_xy).query(s_xy[oks])[0]
        prof = None
        for radius in (NEAR_UM, 100.0, np.inf):
            sel = s_cols[d <= radius]
            cnt = np.asarray(self.cyto_counts[:, sel].sum(axis=1),
                             dtype=float).ravel()
            if cnt.sum() >= NEAR_MIN_MOL or radius == np.inf:
                prof = cnt / max(cnt.sum(), 1.0)
                break
        self._near_profile_cache[key] = prof
        return prof


def dose_response(Y_t, tot_t, cell_pos, pinfo, guide_idx, share_guide,
                  exposure, verbose_rows=None, pair=None):
    """Isotonic dose-response of guide-gene content on exposure, converted to
    a whole-profile contamination fraction per exposure stratum.

    Y_t: dense counts (cells_of_type x genes); cell_pos: positions of the
    pair's target cells within Y_t; returns lambda per stratum dict."""
    g = np.asarray(Y_t[cell_pos][:, guide_idx].sum(axis=1)).ravel()
    t = tot_t[cell_pos]
    strata = np.unique(exposure)
    m_s = np.array([g[exposure == s].sum() for s in strata])
    t_s = np.array([t[exposure == s].sum() for s in strata])
    rate = pava_increasing(m_s / np.maximum(t_s, 1.0), t_s)
    if strata[0] != 0:
        return {}
    excess = np.maximum(rate - rate[0], 0.0)
    lam = np.minimum(excess / max(share_guide, 1e-6), LAMBDA_MAX)
    if verbose_rows is not None:
        for s, r, ex, l, ts in zip(strata, rate, excess, lam, t_s):
            verbose_rows.append(dict(pair=pair, stratum=int(s), totals=ts,
                rate=r, excess_rate=ex, lam=l))
    return dict(zip(strata.tolist(), lam))


def induced_screen(Y_t, tot_t, pinfo_list, inp, psi, allowed_mask):
    """Per-gene proportionality screen for one target type.

    For each pair S -> T, compares the per-gene exposure-linked excess in
    exposed target cells against the best proportional fit to the source
    profile; genes with standardized residual > IND_Z and excess >
    IND_MIN_EXCESS form the induced-term support. Returns per-pair dict of
    (support gene indices, initial m values, diagnostics rows)."""
    out = {}
    rows = []
    for pinfo in pinfo_list:
        S = pinfo["S"]
        e = pinfo["exposure"]
        pos = pinfo["pos"]
        exp_pos = pos[e > 0]
        un_pos = pos[e == 0]
        if len(exp_pos) < 50 or len(un_pos) < 50:
            out[pinfo["pair"]] = (np.array([], int), np.array([]))
            continue
        c_exp = np.asarray(Y_t[exp_pos].sum(axis=0)).ravel()
        t_exp = tot_t[exp_pos].sum()
        c_un = np.asarray(Y_t[un_pos].sum(axis=0)).ravel()
        t_un = tot_t[un_pos].sum()
        r0 = c_un / max(t_un, 1.0)
        excess = c_exp - r0 * t_exp
        var = c_exp + (t_exp / max(t_un, 1.0)) ** 2 * c_un + 1.0
        gset = np.flatnonzero((inp.top_type == S) & allowed_mask)
        p = psi[gset]
        w = 1.0 / var[gset]
        slope = max((w * excess[gset] * p).sum() /
                    max((w * p * p).sum(), 1e-12), 0.0)
        # Overdispersed residual test: counting noise plus a multiplicative
        # profile-uncertainty term. This is the sparsity mechanism - only
        # genes whose excess deviates from proportionality beyond both
        # Poisson noise AND plausible profile error load the induced term.
        resid = excess[gset] - slope * p
        z = resid / np.sqrt(var[gset] + (PROF_CV * slope * p) ** 2)
        sel = (z > IND_Z) & (excess[gset] > IND_MIN_EXCESS)
        sup = gset[sel]
        # initial per-gene induced rate per unit dose: the non-proportional
        # share of the excess divided by the pair's total dose molecules.
        dose_tot = max((tot_t[pos] * pinfo["lam_cell"]).sum(), 1.0)
        m0 = np.maximum(excess[sup] - slope * psi[sup], 0.0) / dose_tot
        out[pinfo["pair"]] = (sup, m0)
        for gi, zi in zip(gset[sel], z[sel]):
            rows.append(dict(pair=pinfo["pair"], gene=inp.genes[gi],
                resid_z=float(zi), excess=float(excess[gi]),
                proportional_expected=float(slope * psi[gi])))
    return out, rows


def fit_target(inp, T, own_frac_by_type, arm, rng, diag, return_rates=False):
    """Fit the mixture for one target type; returns (removed sparse block in
    global gene x cell coordinates, own-fraction per gene for profile
    decontamination, diagnostics). own_frac_by_type carries the previous
    outer round's per-gene own-expression fractions of each type (None on
    the first round) used to decontaminate the source profiles.

    The EM operates on the nonzero count pattern: the posterior weight W is
    zero wherever the observed count is zero, so every update is a sparse
    matrix product over the nonzero entries, and the expected-count totals
    needed for the likelihood are computed in closed form. This changes
    nothing about the model - only the cost, from (cells x genes) dense to
    the number of observed molecules."""
    cols = inp.type_cols[T]
    if len(cols) == 0:
        return None, None
    Ys = sp.csr_matrix(inp.counts[:, cols].T)
    t = np.asarray(Ys.sum(axis=1)).ravel()
    n, G = Ys.shape
    y_nz = Ys.data.astype(np.float64)
    gcols = Ys.indices
    rows_nz = np.repeat(np.arange(n), np.diff(Ys.indptr))
    t_nz = t[rows_nz]
    pos_of_col = {c: i for i, c in enumerate(cols)}

    plist = []
    for pair, pinfo in inp.pair_info.items():
        if pinfo["T"] != T:
            continue
        d = dict(pinfo)
        d["pair"] = pair
        d["pos"] = np.array([pos_of_col[c] for c in pinfo["cells"]])
        plist.append(d)
    if not plist:
        return None, None
    sources = [d["S"] for d in plist]
    S_n = len(plist)

    # Validation arm: genes whose information may not guide the dose.
    maskB = np.zeros(G, bool)
    if arm == "validation":
        for d in plist:
            maskB[inp.gidx(d["B"])] = True
    allowed = ~maskB

    # Per-pair contamination profiles: cytoplasmic profile of the source
    # cells bordering this target type, decontaminated by the source type's
    # fitted own-expression fraction (previous outer round), and zeroed on
    # the target's own genes - contamination there is unidentifiable from
    # the target's own expression and is deliberately left in place, which
    # also guarantees zero removal of the target's own markers.
    t_owned = inp.top_type == T
    Psi = np.zeros((S_n, G))
    for j, d in enumerate(plist):
        p = inp.near_source_profile(d["S"], T).copy()
        of = own_frac_by_type.get(d["S"]) if own_frac_by_type else None
        if of is not None:
            p = p * of
            p /= max(p.sum(), 1e-12)
        p[t_owned] = 0.0
        Psi[j] = p

    # Per-pair dose-response -> per-cell prior contamination fraction.
    Lam = np.zeros((n, S_n))
    dose_rows = diag.setdefault("dose", [])
    for j, d in enumerate(plist):
        guide = d["A"] if arm == "validation" else d["pool"]
        gidx = inp.gidx(guide)
        share_guide = Psi[j, gidx].sum()
        e = d["exposure"].copy()
        if arm == "shuffled":
            e = rng.permutation(e)
        lam_of = dose_response(Ys, t, d["pos"], d, gidx, share_guide, e,
                               verbose_rows=dose_rows, pair=d["pair"])
        lam_cell = np.array([lam_of.get(int(s), 0.0) for s in e])
        Lam[d["pos"], j] = lam_cell
        d["lam_cell"] = lam_cell
        d["e_used"] = e
        d["guide_idx"] = gidx
        d["share_guide"] = share_guide

    # Ambient component (production only): strict-union genes, profile and
    # scale from target cells far from any source cell.
    strictU = np.zeros(G, bool)
    for d in plist:
        strictU[inp.gidx(d["strict"])] = True
    a_prof = np.zeros(G)
    lam_a = 0.0
    use_ambient = arm in ("production", "production_noind")
    if use_ambient and strictU.any():
        xy = inp.cell_xy[cols]
        dmin = np.full(n, np.inf)
        for S in set(sources):
            sxy = inp.cell_xy[inp.type_cols[S]]
            sxy = sxy[np.isfinite(sxy).all(axis=1)]
            if len(sxy):
                dmin = np.minimum(dmin, cKDTree(sxy).query(xy)[0])
        far = dmin > FAR_UM
        far_mol = t[far].sum()
        src_est = "far_target_cells"
        if far_mol >= MIN_FAR_MOL:
            rates = np.asarray(Ys[far].sum(axis=0)).ravel() / far_mol
        else:
            # fallback: half of the pooled exposure-0 strict rate (the audit
            # measured the far-field level at roughly 0.5-0.7 of the pooled
            # zero-exposure baseline).
            src_est = "half_exposure0_baseline"
            e0 = np.zeros(n, bool)
            for d in plist:
                e0[d["pos"][d["exposure"] == 0]] = True
            rates = 0.5 * np.asarray(Ys[e0].sum(axis=0)).ravel() \
                / max(t[e0].sum(), 1.0)
        lam_a = float(rates[strictU].sum())
        if lam_a > 0:
            a_prof[strictU] = rates[strictU] / lam_a
        diag.setdefault("ambient", []).append(dict(target=T, lam_a=lam_a,
            n_far=int(far.sum()) if far_mol >= MIN_FAR_MOL else 0,
            estimate=src_est))

    # Own programs.
    F = inp.program_init[T].copy()
    K = F.shape[0]
    if arm in ("production", "production_noind") and strictU.any():
        F[:, strictU] = 0.0
        F /= np.maximum(F.sum(axis=1, keepdims=True), 1e-12)

    # Induced-expression support (per pair).
    M = np.zeros((S_n, G))
    Rho = np.ones((n, S_n))
    use_induced = arm in ("validation", "production")
    ind_rows = []
    if use_induced:
        for j, d in enumerate(plist):
            sc, rows = induced_screen(Ys, t, [d], inp, Psi[j], allowed)
            sup, m0 = sc[d["pair"]]
            M[j, sup] = m0
            ind_rows.extend(rows)
    diag.setdefault("induced", []).extend(ind_rows)

    # Initial cell-level parameters.
    Alpha = Lam.copy()
    beta = np.full(n, lam_a)
    lam_tot = Lam.sum(axis=1)
    own0 = np.maximum(1.0 - lam_tot - lam_a, 0.2)
    colsum = np.asarray(Ys.sum(axis=0)).ravel()
    prog_share = np.maximum(F @ (colsum / max(colsum.sum(), 1.0)), 1e-3)
    prog_share /= prog_share.sum()
    Theta = own0[:, None] * prog_share[None, :]
    eps_g = EPS_TOTAL / G

    share_allowed = np.maximum((Psi * allowed[None, :]).sum(axis=1), 1e-6)
    PsiA = Psi * allowed[None, :]

    # Own-expression programs are learned predominantly from lightly dosed
    # cells, so that admixed content in heavily exposed cells cannot be
    # absorbed into a "cell state".
    w_dose = 1.0 / (1.0 + DOSE_W * Lam.sum(axis=1))

    # The induced-term dose scale stays at the initial prior (the induced
    # rate is defined per unit of the measured dose); the contamination
    # prior Lam itself grows through the top-up passes below.
    U = Lam.copy()
    ll_state = [-np.inf]

    def nz_rates():
        """Component rates at the nonzero count entries."""
        r_own = np.zeros(y_nz.shape[0])
        for k in range(K):
            r_own += Theta[rows_nz, k] * F[k, gcols]
        r_bg = np.zeros(y_nz.shape[0])
        for j in range(S_n):
            r_bg += Alpha[rows_nz, j] * Psi[j, gcols]
        if lam_a > 0:
            r_bg += beta[rows_nz] * a_prof[gcols]
        r_ind = np.zeros(y_nz.shape[0])
        if use_induced and M.any():
            UR = U * Rho
            for j in range(S_n):
                if M[j].any():
                    r_ind += UR[rows_nz, j] * M[j, gcols]
        return r_own, r_bg, r_ind

    def expected_total():
        """Sum of expected counts over every (cell, gene) entry, in closed
        form: each component factorizes into a cell-side and gene-side sum."""
        tot = float((t[:, None] * Theta).sum(axis=0) @ F.sum(axis=1))
        tot += float((t[:, None] * Alpha).sum(axis=0) @ Psi.sum(axis=1))
        if lam_a > 0:
            tot += float((t * beta).sum() * a_prof.sum())
        if use_induced and M.any():
            tot += float((t[:, None] * (U * Rho)).sum(axis=0) @ M.sum(axis=1))
        return tot + float(t.sum()) * G * eps_g

    def em_pass(n_iter, update_F=True):
        nonlocal Theta, Alpha, beta, Rho, M, F
        for it in range(n_iter):
            r_own, r_bg, r_ind = nz_rates()
            rate_nz = r_own + r_bg + r_ind + eps_g
            R_nz = t_nz * rate_nz
            check_ll = (it + 1) % 10 == 0 or it == n_iter - 1
            exp_tot = expected_total() if check_ll else None
            W = sp.csr_matrix((y_nz / np.maximum(R_nz, 1e-300),
                               gcols, Ys.indptr), shape=(n, G))
            Wt = W.T

            # own programs
            Theta = Theta * (W @ F.T)
            # contamination fractions: Gamma-posterior update on allowed
            # genes, scaled by the source-profile share of allowed genes,
            # bounded to ALPHA_CAP-fold above the prior mean.
            ev = Alpha * (W @ PsiA.T) * t[:, None]       # evidence molecules
            Alpha = (ev / np.maximum(share_allowed, 1e-6)[None, :]
                     + C0_ALPHA * Lam) / (t[:, None] + C0_ALPHA)
            Alpha = np.minimum(Alpha, ALPHA_CAP * Lam)
            # ambient scale
            if use_ambient and lam_a > 0:
                evb = beta * (W @ a_prof) * t
                beta = (evb + C0_AMBIENT * lam_a) / (t + C0_AMBIENT)
                beta = np.minimum(beta, AMBIENT_CAP * lam_a)
            # induced term
            if use_induced and M.any():
                zind = (U * Rho) * (W @ M.T) * t[:, None]
                denom = t[:, None] * U * M.sum(axis=1)[None, :]
                Rho = np.minimum((RHO_SHAPE + zind) / (RHO_SHAPE + denom),
                                 RHO_CAP)
                wsum = (Wt @ (t[:, None] * U * Rho)).T       # S_n x G
                dsum = (t[:, None] * U * Rho).sum(axis=0)
                M *= wsum / np.maximum(dsum[:, None], 1e-9)
            # own program profiles (dose-downweighted cells)
            if update_F:
                Fnum = F * (Wt @ (Theta * (t * w_dose)[:, None])).T
                if arm in ("production", "production_noind") and strictU.any():
                    Fnum[:, strictU] = 0.0
                F = Fnum + 1e-8
                F /= np.maximum(F.sum(axis=1, keepdims=True), 1e-12)

            if check_ll:
                ll = float((y_nz * np.log(np.maximum(R_nz, 1e-300))).sum()
                           - exp_tot)
                gain = ll - ll_state[0]
                ll_state[0] = ll
                if it > 10 and abs(gain) < 1e-7 * abs(ll):
                    break

    def keep_at_nz():
        r_own, r_bg, r_ind = nz_rates()
        tot = r_own + r_bg + r_ind + eps_g
        keep = (r_own + r_ind + eps_g) / np.maximum(tot, 1e-300)
        return keep, r_own, tot

    em_pass(N_EM)

    # Top-up passes (mirroring the regression corrector): re-measure the
    # residual exposure gradient of the guide genes on the corrected counts
    # (excluding induced-support genes, whose gradient is retained by
    # design), add it to the prior dose, and continue the EM.
    for topup in range(TOPUP_PASSES):
        keep_nz, _, _ = keep_at_nz()
        Ycorr = sp.csr_matrix((y_nz * keep_nz, gcols, Ys.indptr),
                              shape=(n, G))
        added = 0.0
        for j, d in enumerate(plist):
            gidx2 = np.array([g for g in d["guide_idx"] if M[j, g] <= 0],
                             dtype=int)
            if len(gidx2) < 3:
                continue
            share2 = Psi[j, gidx2].sum()
            if share2 <= 1e-6:
                continue
            lam_of = dose_response(Ycorr, t, d["pos"], d, gidx2, share2,
                                   d["e_used"])
            lam_extra = np.array([lam_of.get(int(s), 0.0)
                                  for s in d["e_used"]])
            Lam[d["pos"], j] = np.minimum(Lam[d["pos"], j] + lam_extra,
                                          LAMBDA_MAX)
            added += float((lam_extra * t[d["pos"]]).sum())
        if added < 1e-4 * t.sum():
            break
        Alpha = np.maximum(Alpha, Lam)
        # own programs stay frozen during top-up: they were estimated from
        # lightly dosed cells, and the top-up only recalibrates delivery of
        # the measured residual gradient.
        em_pass(15, update_F=False)

    # Final rates and the retained (kept) share of each observed count.
    keep_nz, r_own_nz, tot_nz = keep_at_nz()
    rem = sp.csr_matrix((y_nz * (1.0 - keep_nz), gcols, Ys.indptr),
                        shape=(n, G)).T.tocsc()  # G x n block

    # own fraction per gene (for source-profile decontamination)
    own_num = np.bincount(gcols, weights=y_nz
                          * (r_own_nz / np.maximum(tot_nz, 1e-300)),
                          minlength=G)
    own_den = np.bincount(gcols, weights=y_nz, minlength=G)
    own_frac_g = own_num / np.maximum(own_den, 1.0)

    # per-pair diagnostics: prior vs posterior contamination molecules
    for j, d in enumerate(plist):
        pos = d["pos"]
        diag.setdefault("alpha", []).append(dict(pair=d["pair"],
            prior_molecules=float((Lam[pos, j] * t[pos]).sum()),
            posterior_molecules=float((Alpha[pos, j] * t[pos]).sum()),
            induced_molecules=float((U[pos, j] * Rho[pos, j] * t[pos]
                                     * M[j].sum()).sum()),
            mean_lambda_exposed=float(Lam[pos, j][d["e_used"] > 0].mean()
                                      if (d["e_used"] > 0).any() else 0)))

    if return_rates:
        # Expected molecule counts per entry, by component (for the
        # molecule-level realization): own expression, contamination,
        # ambient, induced, floor. Dense, but only requested for the one
        # spike-in target type.
        rates = dict(
            own=t[:, None] * (Theta @ F),
            cont=t[:, None] * (Alpha @ Psi),
            amb=t[:, None] * (beta[:, None] * a_prof[None, :]),
            ind=t[:, None] * ((U * Rho) @ M),
            eps=t[:, None] * np.full((1, G), eps_g),
            cols=cols, strictU=strictU,
            pair_names=[d["pair"] for d in plist],
            sources=[d["S"] for d in plist], Psi=Psi)
        return (cols, rem), own_frac_g, rates
    return (cols, rem), own_frac_g


_PAR_INP = None


def _fit_one(args):
    """Worker for the per-type parallel fit: types within an outer round are
    independent, so they run as forked processes sharing the loaded inputs."""
    T, own_frac, arm, seed = args
    rng = np.random.default_rng(seed)
    d = {}
    out = fit_target(_PAR_INP, T, own_frac, arm, rng, d)
    return T, out, d


def run_arm(arm, seed=1):
    global _PAR_INP
    os.makedirs(RESULTS, exist_ok=True)
    inp = Inputs()
    n_workers = min(int(os.environ.get("GM_WORKERS", "8")),
                    len(inp.type_names))
    own_frac = None
    n_outer = 1 if arm == "shuffled" else N_OUTER
    diag = {}
    removed_blocks = None
    for rnd in range(n_outer):
        diag = {}
        removed_blocks = {}
        own_frac_new = {}
        args = [(T, own_frac, arm, seed * 100003 + ti)
                for ti, T in enumerate(inp.type_names)]
        if n_workers > 1:
            import multiprocessing as mp
            _PAR_INP = inp
            with mp.get_context("fork").Pool(n_workers) as pool:
                results = pool.map(_fit_one, args)
        else:
            _PAR_INP = inp
            results = [_fit_one(a) for a in args]
        for T, out, dlocal in results:
            res, of = out
            for k, v in dlocal.items():
                diag.setdefault(k, []).extend(v)
            if res is None:
                continue
            removed_blocks[T] = res
            own_frac_new[T] = of
            tot_rem = float(res[1].sum())
            print(f"  [{arm} round {rnd + 1}] {T}: removed {tot_rem:,.0f} "
                  f"molecules from {len(res[0])} cells", flush=True)
        # source-profile decontamination for the next round
        own_frac = own_frac_new

    # assemble the full removed matrix (gene x all cells)
    rows, cols_ix, vals = [], [], []
    for T, (cols, rem) in removed_blocks.items():
        co = rem.tocoo()
        rows.append(co.row)
        cols_ix.append(cols[co.col])
        vals.append(co.data)
    full = sp.coo_matrix((np.concatenate(vals),
        (np.concatenate(rows), np.concatenate(cols_ix))),
        shape=(inp.G, len(inp.cells)))
    out_mtx = os.path.join(DATA, f"gm_removed_{arm}.mtx")
    scipy.io.mmwrite(out_mtx, full)
    print(f"wrote {out_mtx}: {full.sum():,.0f} molecules removed total")

    for key, fname in [("dose", f"{PREFIX}gm_{arm}_pair_dose.csv"),
                       ("induced", f"{PREFIX}gm_{arm}_induced.csv"),
                       ("alpha", f"{PREFIX}gm_{arm}_alpha.csv"),
                       ("ambient", f"{PREFIX}gm_{arm}_ambient.csv")]:
        if key in diag and diag[key]:
            pd.DataFrame(diag[key]).to_csv(os.path.join(RESULTS, fname),
                                           index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True,
        choices=["validation", "production", "production_noind", "shuffled"])
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    run_arm(args.arm, args.seed)
