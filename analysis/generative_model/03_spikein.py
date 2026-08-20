"""Spike-in retention/removal test for the generative model (gate 3).

Reuses the four-class spike-in design of analysis/audit_guided/14_t7_spikeins.py
(same pair, same seed, same placement rules): into real exposed target cells of
the pair Exocrine epithelial -> Endothelial, inject labeled synthetic
molecules of four classes:

  inplane   - multi-gene fragments copied from a nearby source cell, placed
              at the target edge facing the source (true admixture);
  outplane  - the same fragments placed over the target footprint with z
              pushed toward the slab surfaces (out-of-plane admixture);
  induction - molecules of three source-marker genes placed like the
              target's own molecules, counts scaled with exposure
              (genuine induced expression - must be RETAINED);
  ambient   - single molecules of random pool genes (must be removed).

The spiked molecules are added to the count matrix, the generative model is
refit in its production configuration (all types, one profile-refinement
round, then the target refit), and each planted molecule is scored by the
posterior removal fraction of its (cell, gene) entry. Gate targets:
induction retention >= 0.8; in-plane and out-of-plane removal >= 0.8;
ambient removed.
"""

import importlib.util
import os
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scipy.sparse as sp
from scipy.spatial import cKDTree

GM = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
spec = importlib.util.spec_from_file_location("gm_model",
    os.path.join(GM, "01_model.py"))
gm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gm)

RUN = gm.RUN
PAIR = "Exocrine epithelial -> Endothelial"
N_CELLS_PER_CLASS = 400
rng_spike = np.random.default_rng(29)

print("loading inputs ...", flush=True)
inp = gm.Inputs()

# ---- molecule-level tables for fragment sampling (as in 14_t7_spikeins.py) --
mol = pq.read_table(os.path.join(RUN, "molecules.parquet"),
    columns=["x", "y", "z", "gene_idx", "cell_idx"]).to_pandas()
cells_tbl = pq.read_table(os.path.join(RUN, "cells.parquet"),
    columns=["cell_idx", "cell_id", "x", "y"]).to_pandas()
idx_of_cell = dict(zip(cells_tbl.cell_id, cells_tbl.cell_idx))
mol = mol[mol.cell_idx >= 0]
mol = mol.iloc[np.argsort(mol.cell_idx.to_numpy(), kind="stable")]
cell_arr = mol.cell_idx.to_numpy()
starts = np.searchsorted(cell_arr, np.arange(cell_arr.max() + 2))
xyz = mol[["x", "y", "z"]].to_numpy()
garr = mol.gene_idx.to_numpy()

pinfo = inp.pair_info[PAIR]
S, T = pinfo["S"], pinfo["T"]
pool_ids = inp.gidx(pinfo["pool"])
expo_cells = np.array(inp.cells)[pinfo["cells"]]
exposed = [c for c, e in zip(expo_cells, pinfo["exposure"]) if e > 0]
e_of = dict(zip(expo_cells, pinfo["exposure"]))
s_cells = [c for c in np.array(inp.cells)[inp.type_cols[S]] if c in idx_of_cell]
s_xy = cells_tbl.set_index("cell_id").loc[s_cells, ["x", "y"]].to_numpy()
s_tree = cKDTree(s_xy)
zlo, zhi = np.percentile(xyz[:, 2], [10, 90])


def cell_mols(cid):
    i = idx_of_cell[cid]
    sl = slice(starts[i], starts[i + 1])
    return xyz[sl], garr[sl]


def sample_fragment():
    for _ in range(50):
        cid = s_cells[rng_spike.integers(len(s_cells))]
        c, g = cell_mols(cid)
        if len(g) < 8:
            continue
        seed = rng_spike.integers(len(g))
        d = np.linalg.norm(c - c[seed], axis=1)
        sel = np.flatnonzero(d <= 1.5)
        if len(sel) >= 4:
            sel = sel[:12]
            return c[sel] - c[sel].mean(0), g[sel]
    return None, None


records = []
target_cells = rng_spike.permutation(exposed)
groups = np.array_split(target_cells[:4 * N_CELLS_PER_CLASS], 4)
# Primary variant (GM_IND_RANK=0, the audit's T7 design): induction on the
# pair's top-3 pool genes - the largest contamination channels. The
# supplementary variant (GM_IND_RANK=9) plants the same molecules on
# mid-rank pool genes, where the planted excess is large relative to the
# proportional contamination expectation.
IND_RANK = int(os.environ.get("GM_IND_RANK", "0"))
SUFFIX = "" if IND_RANK == 0 else f"_rank{IND_RANK}"
induction_genes = pool_ids[IND_RANK:IND_RANK + 3]
print("planted induction genes:", [inp.genes[g] for g in induction_genes])

for cls, cell_group in zip(["inplane", "outplane", "induction", "ambient"],
                           groups):
    for cid in cell_group:
        c, g = cell_mols(cid)
        if len(g) < 10:
            continue
        cx = c.mean(0)
        if cls in ("inplane", "outplane"):
            fc, fg = sample_fragment()
            if fc is None:
                continue
            if cls == "inplane":
                _, si = s_tree.query(cx[:2])
                to_s = s_xy[si] - cx[:2]
                to_s = to_s / max(np.linalg.norm(to_s), 1e-9)
                r90 = np.percentile(np.linalg.norm(c[:, :2] - cx[:2], axis=1), 90)
                anchor = np.array([cx[0] + to_s[0] * r90,
                                   cx[1] + to_s[1] * r90, cx[2]])
            else:
                ang = rng_spike.uniform(0, 2 * np.pi)
                rad = rng_spike.uniform(0, np.percentile(
                    np.linalg.norm(c[:, :2] - cx[:2], axis=1), 70))
                zpos = zlo if rng_spike.random() < 0.5 else zhi
                anchor = np.array([cx[0] + rad * np.cos(ang),
                                   cx[1] + rad * np.sin(ang), zpos])
            pos = fc + anchor
            for p, gg in zip(pos, fg):
                records.append((cid, cls, *p, int(gg)))
        elif cls == "induction":
            n_add = 2 + 2 * min(int(e_of.get(cid, 1)), 4)
            sel = rng_spike.choice(len(c), n_add)
            jit = rng_spike.normal(0, 0.5, (n_add, 3))
            for p, gg in zip(c[sel] + jit,
                             rng_spike.choice(induction_genes, n_add)):
                records.append((cid, cls, *p, int(gg)))
        else:
            n_add = 2
            sel = rng_spike.choice(len(c), n_add)
            jit = rng_spike.normal(0, 3.0, (n_add, 3))
            for p, gg in zip(c[sel] + jit, rng_spike.choice(pool_ids, n_add)):
                records.append((cid, cls, *p, int(gg)))

spk = pd.DataFrame(records,
    columns=["cell_id", "cls", "x", "y", "z", "gene_idx"])
print("injected:", spk.groupby("cls").size().to_dict())
del mol  # keep xyz/garr/starts: the molecule-level scorer needs positions

# ---- add the spikes to the count matrix ------------------------------------
col_ids = np.array([inp.col_of[c] for c in spk.cell_id])
add = sp.coo_matrix((np.ones(len(spk)),
    (spk.gene_idx.to_numpy(), col_ids)),
    shape=inp.counts.shape).tocsc()
inp.counts = (inp.counts + add).tocsc()
inp.totals = np.asarray(inp.counts.sum(axis=0)).ravel()

# ---- refit, production configuration, with one profile-refinement round ----
arm = "production"
CACHE = os.path.join(gm.DATA, f"gm_spikein_rates_cache{SUFFIX}.npz")
if os.path.exists(CACHE):
    print("using cached spiked-refit rates:", CACHE, flush=True)
    z = np.load(CACHE)
    rates_T = dict(own=z["own"], cont=z["cont"], amb=z["amb"], ind=z["ind"],
                   eps=z["eps"], cols=z["cols"])
    cols = z["cols"]
    Yt = np.asarray(inp.counts[:, cols].todense(), dtype=float).T
    tot_rate = rates_T["own"] + rates_T["cont"] + rates_T["amb"] \
        + rates_T["ind"] + rates_T["eps"]
    keep = (rates_T["own"] + rates_T["ind"] + rates_T["eps"]) \
        / np.maximum(tot_rate, 1e-300)
    removed_T = (cols, sp.csc_matrix((Yt * (1 - keep)).T))
else:
    rng = np.random.default_rng(1)
    own_frac = None
    removed_T = None
    rates_T = None
    diag_last = None
    for rnd in range(gm.N_OUTER):
        diag = {}
        own_frac_new = {}
        for Ttype in inp.type_names:
            if Ttype == T:
                res, of, rates_T = gm.fit_target(inp, Ttype, own_frac, arm,
                                                 rng, diag, return_rates=True)
            else:
                res, of = gm.fit_target(inp, Ttype, own_frac, arm, rng, diag)
            if res is None:
                continue
            own_frac_new[Ttype] = of
            if Ttype == T:
                removed_T = res
            print(f"  [spike round {rnd + 1}] {Ttype}: removed "
                  f"{float(res[1].sum()):,.0f}", flush=True)
        own_frac = own_frac_new
        diag_last = diag
    pd.DataFrame(diag_last.get("induced", [])).to_csv(
        os.path.join(gm.RESULTS, f"gm_spikein_induced_screen{SUFFIX}.csv"), index=False)
    np.savez_compressed(CACHE, own=rates_T["own"], cont=rates_T["cont"],
                        amb=rates_T["amb"], ind=rates_T["ind"],
                        eps=rates_T["eps"], cols=rates_T["cols"])

cols_T, rem_T = removed_T           # rem_T: genes x cells-of-T block
pos_of_col = {c: i for i, c in enumerate(cols_T)}
Yspk = inp.counts[:, cols_T]        # spiked counts for target cells
rem_T = rem_T.tocsr()
Yspk = Yspk.tocsr()

rows = []
for (cid, cls, g), grp in spk.groupby(
        [spk.cell_id, spk.cls, spk.gene_idx]):
    j = pos_of_col[inp.col_of[cid]]
    y = Yspk[g, j]
    r = rem_T[g, j]
    frac_removed = float(r / y) if y > 0 else 0.0
    rows.append(dict(cls=cls, gene=inp.genes[g], n=len(grp),
                     removal_prob=frac_removed))
per_mol = pd.DataFrame(rows)
cls_stats = per_mol.groupby("cls").apply(
    lambda d: pd.Series(dict(
        molecules=int(d.n.sum()),
        mean_removal=float((d.removal_prob * d.n).sum() / d.n.sum()),
        frac_removed_over_half=float(
            ((d.removal_prob > 0.5) * d.n).sum() / d.n.sum()))),
    include_groups=False).reset_index()
cls_stats["mean_retention"] = 1 - cls_stats["mean_removal"]
print("\nspike-in gate results (count-level posterior):")
print(cls_stats.to_string(index=False))
cls_stats.to_csv(os.path.join(gm.RESULTS, f"gm_spikein_classes{SUFFIX}.csv"),
                 index=False)
per_mol.to_csv(os.path.join(gm.RESULTS, f"gm_spikein_per_gene{SUFFIX}.csv"),
               index=False)

ind = cls_stats.set_index("cls")
print("\ncount-level (rate-share) summary:")
for k, v in [("induction_retention", float(ind.loc["induction", "mean_retention"])),
             ("inplane_removal", float(ind.loc["inplane", "mean_removal"])),
             ("outplane_removal", float(ind.loc["outplane", "mean_removal"])),
             ("ambient_removal", float(ind.loc["ambient", "mean_removal"]))]:
    print(f"  {k}: {v:.3f}")
spk.to_csv(os.path.join(gm.RESULTS, f"gm_spikein_molecules{SUFFIX}.csv.gz"),
           index=False)

# ===========================================================================
# Molecule-level realization: per-molecule posterior origin.
#
# For each count-matrix entry (cell, gene) of the target, the observed y
# molecules are split between a removable origin (contamination + ambient,
# a Poisson count process with mean mu_c from the fitted cell-level model)
# and a retained origin (own + induced + floor, mean mu_r), extended by a
# spike-and-slab term: with small prior probability PI0 an entry carries an
# unmodeled induced burst whose extra-molecule count follows a geometric
# distribution with mean NU0. Each molecule i carries a spatial-coherence
# feature x_i - the number of DISTINCT other source-owned genes among the
# cell's molecules within PATCH_R um - modeled as Poisson with mean mu_pc
# for admixed fragments and mu_po for native-placed molecules (both means
# estimated from reference molecule sets in the same cells: real strict-gene
# molecules in exposed cells for mu_pc, real native-marker molecules for
# mu_po). The exact posterior over which molecules are removable combines
# the count prior with the per-molecule feature likelihood ratios through
# elementary symmetric polynomials, yielding a removal probability per
# molecule. Gene identity thus sets the prior; spatial coherence moves it.
# ===========================================================================
from scipy.stats import poisson

PATCH_R = 1.0  # the audit found 1 um informative; larger radii saturate
PI0 = 0.02   # prior probability of an unmodeled induced burst at an entry
NU0 = 4.0    # mean size of such a burst (geometric)

cols_T = removed_T[0]
pos_of_col = {c: i for i, c in enumerate(cols_T)}
mu_rem_M = rates_T["cont"] + rates_T["amb"]          # n x G removable mean
mu_ret_M = rates_T["own"] + rates_T["ind"] + rates_T["eps"]
src_owned = np.flatnonzero(inp.top_type == S)
src_owned_set = set(src_owned.tolist())
native = pd.read_csv(os.path.join(gm.DATA, "native_markers.csv"))
native_T = set(inp.gidx(native[native.cell_type == T].gene.tolist()).tolist())
strict_T = set(inp.gidx(pinfo["strict"]).tolist())
if not strict_T:
    strict_T = set(pool_ids.tolist())

def geom_pmf(j, nu):
    p = 1.0 / (1.0 + nu)
    return p * (1 - p) ** j

def entry_posterior(y_feats, mu_c, mu_r):
    """y_feats: list of per-molecule feature pairs (x_src, x_nat): the
    number of distinct other source-owned genes and the number of
    target-native-marker molecules within PATCH_R um. Both are modeled as
    class-conditional Poisson (admixed vs native origin), with means
    estimated from reference molecule sets. Returns per-molecule removal
    probabilities (posterior probability of a removable origin), exact
    over subsets via symmetric polynomials."""
    y = len(y_feats)
    lr = np.array([
        (poisson.pmf(xs, mu_src_c) / max(poisson.pmf(xs, mu_src_o), 1e-12))
        * (poisson.pmf(xn, mu_nat_c) / max(poisson.pmf(xn, mu_nat_o), 1e-12))
        for xs, xn in y_feats])
    lr = np.clip(lr, 1e-3, 1e3)
    # elementary symmetric polynomials e_0..e_y of lr
    e = np.zeros(y + 1); e[0] = 1.0
    for w in lr:
        e[1:] = e[1:] + w * e[:-1]
    # count prior over k removable molecules
    base = np.zeros(y + 1)
    for k in range(y + 1):
        ret = (1 - PI0) * poisson.pmf(y - k, mu_r) \
            + PI0 * sum(poisson.pmf(j, mu_r) * geom_pmf(y - k - j, NU0)
                        for j in range(y - k + 1))
        base[k] = poisson.pmf(k, mu_c) * ret
    from math import comb
    pk = np.array([base[k] * e[k] / comb(y, k) for k in range(y + 1)])
    if pk.sum() <= 0:
        return np.full(y, mu_c / max(mu_c + mu_r, 1e-12))
    pk /= pk.sum()
    # per-molecule marginal: P(i removable) = sum_k pk[k]/e_k * w_i * e_{k-1}(without i)
    probs = np.zeros(y)
    for i in range(y):
        # e without molecule i (deflation)
        ei = np.zeros(y); ei[0] = 1.0
        for j, w in enumerate(lr):
            if j == i:
                continue
            ei[1:] = ei[1:] + w * ei[:-1]
        acc = 0.0
        for k in range(1, y + 1):
            if e[k] > 0 and pk[k] > 0:
                # fraction of size-k subsets (weighted) containing i
                acc += pk[k] * lr[i] * ei[k - 1] / e[k]
        probs[i] = min(acc, 1.0)
    return probs

# --- estimate patch-feature means from reference molecule sets -------------
def mol_features(allc, allg, anchors):
    """(x_src, x_nat) for each anchor index: distinct other source-owned
    genes and native-marker molecule count within PATCH_R um."""
    tree = cKDTree(allc)
    out = []
    for i in anchors:
        nb = tree.query_ball_point(allc[i], PATCH_R)
        gs = {int(allg[j]) for j in nb if j != i
              and int(allg[j]) in src_owned_set and allg[j] != allg[i]}
        xn = sum(1 for j in nb if j != i and int(allg[j]) in native_T)
        out.append((len(gs), xn))
    return out

def patch_counts_for(cell_ids, gene_filter, max_mols=4000):
    out = []
    for cid in cell_ids:
        c0, g0 = cell_mols(cid)
        s = spk[spk.cell_id == cid]
        allc = np.vstack([c0, s[["x", "y", "z"]].to_numpy()]) if len(s) else c0
        allg = np.concatenate([g0, s.gene_idx.to_numpy()]) if len(s) else g0
        anchors = [i for i in range(len(g0)) if g0[i] in gene_filter]
        out.extend(mol_features(allc, allg, anchors))
        if len(out) >= max_mols:
            break
    return np.array(out[:max_mols]) if out else np.zeros((0, 2))

exposed_set = [c for c in exposed if c in idx_of_cell]
rng2 = np.random.default_rng(7)
samp = list(rng2.permutation(exposed_set)[:800])
xc = patch_counts_for(samp, strict_T)
xo = patch_counts_for(samp, native_T)
mu_src_c = max(xc[:, 0].mean(), 0.05) if len(xc) else 1.0
mu_src_o = max(xo[:, 0].mean(), 0.02) if len(xo) else 0.1
mu_nat_c = max(xc[:, 1].mean(), 0.02) if len(xc) else 0.1
mu_nat_o = max(xo[:, 1].mean(), 0.05) if len(xo) else 1.0
print(f"\npatch-feature means (radius {PATCH_R} um):")
print(f"  distinct other source genes: admixed ref {mu_src_c:.2f} vs "
      f"native ref {mu_src_o:.2f}  ({len(xc)}/{len(xo)} molecules)")
print(f"  native-marker neighbors:     admixed ref {mu_nat_c:.2f} vs "
      f"native ref {mu_nat_o:.2f}")

# --- score every planted molecule ------------------------------------------
Ycsr = inp.counts.tocsr()
rows_m = []
for cid, cell_spk in spk.groupby("cell_id"):
    j = pos_of_col.get(inp.col_of[cid])
    if j is None:
        continue
    c0, g0 = cell_mols(cid)
    allc = np.vstack([c0, cell_spk[["x", "y", "z"]].to_numpy()])
    allg = np.concatenate([g0, cell_spk.gene_idx.to_numpy()])
    inj_off = len(g0)
    feats = {}
    for g in np.unique(cell_spk.gene_idx):
        idxs = np.flatnonzero(allg == g)
        xs = mol_features(allc, allg, idxs)
        mu_c = float(mu_rem_M[j, g])
        mu_r = float(mu_ret_M[j, g])
        probs = entry_posterior(xs, mu_c, mu_r)
        for i, pr in zip(idxs, probs):
            if i >= inj_off:
                feats[i] = pr
    for local, (_, srow) in enumerate(cell_spk.iterrows()):
        i = inj_off + local
        if i in feats:
            rows_m.append(dict(cls=srow.cls, gene=inp.genes[int(srow.gene_idx)],
                               removal_prob=float(feats[i])))

mol_df = pd.DataFrame(rows_m)
mstats = mol_df.groupby("cls").removal_prob.agg(["mean", "count"]).reset_index()
mstats["retention"] = 1 - mstats["mean"]
print("\nmolecule-level spike-in results (primary, gate 3):")
print(mstats.rename(columns={"mean": "mean_removal",
                             "count": "molecules"}).to_string(index=False))
mstats.to_csv(os.path.join(gm.RESULTS, f"gm_spikein_classes_molecule{SUFFIX}.csv"),
              index=False)
mol_df.to_csv(os.path.join(gm.RESULTS, f"gm_spikein_per_molecule{SUFFIX}.csv.gz"),
              index=False)
mi = mstats.set_index("cls")
print("\ngate 3 (molecule level): induction retention "
      f"{1 - mi.loc['induction', 'mean']:.3f} (target >= 0.8); "
      f"inplane removal {mi.loc['inplane', 'mean']:.3f}, "
      f"outplane removal {mi.loc['outplane', 'mean']:.3f} (targets >= 0.8); "
      f"ambient removal {mi.loc['ambient', 'mean']:.3f}")
print("SPIKE DONE")
