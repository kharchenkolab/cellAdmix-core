"""T1/T2 molecule-level diagnostics: multi-gene patchiness and baseline
decomposition.

T1 (exposed target cells): for each pair, the fraction of source-pool
molecules with a molecule of a DIFFERENT pool gene within RADIUS um (3D),
against a within-cell permutation null (pool marks reassigned to random
molecules of the same cell, pool-gene multiset preserved). Admixed
fragments are gene-diverse patches; induction/native bursts are not.

T2 (zero-exposure target cells):
  a. same cross-gene statistic on pool molecules (out-of-plane admixture
     should be patchy there too; ambient should not),
  b. strict-gene rate versus lateral distance to the nearest source-type
     cell (out-of-plane leakage decays; ambient is flat),
  c. z-position of strict molecules within their cell's z-distribution
     (out-of-plane material biases toward slab surfaces).

Usage: run from the pancreas example directory.
"""
import gzip
import json
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.spatial import cKDTree

AG = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/audit_guided"
RUN = "out/runs/fit_manual_rank9_invsqrt_kl"
RADIUS = 2.5
N_PERM = 10
rng = np.random.default_rng(11)

genes = json.load(open(f"{RUN}/run.json"))["genes"]
gene_idx = {g: i for i, g in enumerate(genes)}

mol = pq.read_table(f"{RUN}/molecules.parquet",
    columns=["x", "y", "z", "gene_idx", "cell_idx"]).to_pandas()
cells = pq.read_table(f"{RUN}/cells.parquet",
    columns=["cell_idx", "cell_id", "x", "y", "transcript_count"]).to_pandas()
cell_id_of = dict(zip(cells.cell_idx, cells.cell_id))
idx_of_cell = dict(zip(cells.cell_id, cells.cell_idx))
ctypes = pd.read_csv(f"{AG}/results/cell_types.csv")
type_of = dict(zip(ctypes.cell_id, ctypes.cell_type))

mol = mol[mol.cell_idx >= 0].reset_index(drop=True)
order = np.argsort(mol.cell_idx.to_numpy(), kind="stable")
mol = mol.iloc[order].reset_index(drop=True)
cell_arr = mol.cell_idx.to_numpy()
starts = np.searchsorted(cell_arr, np.arange(cell_arr.max() + 2))
xyz = mol[["x", "y", "z"]].to_numpy()
garr = mol.gene_idx.to_numpy()

meta = pd.read_csv(f"{AG}/results/pair_meta.csv")
expo = pd.read_csv(f"{AG}/results/pair_exposure.csv.gz")

def cell_slice(ci):
    return slice(starts[ci], starts[ci + 1])

def crossgene_stat(coords, gids):
    """Fraction of molecules with a different-gene molecule within RADIUS."""
    if len(gids) < 2 or len(set(gids)) < 2:
        return 0.0, len(gids)
    tree = cKDTree(coords)
    pairs = tree.query_pairs(RADIUS, output_type="ndarray")
    hit = np.zeros(len(gids), dtype=bool)
    if len(pairs):
        diff = gids[pairs[:, 0]] != gids[pairs[:, 1]]
        hit[pairs[diff, 0]] = True
        hit[pairs[diff, 1]] = True
    return float(hit.sum()), len(gids)

def patch_test(cell_ids, pool_ids):
    """Observed vs permuted cross-gene co-location over a set of cells."""
    obs_hit = obs_n = 0.0
    null_hit = np.zeros(N_PERM)
    null_n = 0.0
    for cid in cell_ids:
        ci = idx_of_cell.get(cid)
        if ci is None:
            continue
        sl = cell_slice(ci)
        g = garr[sl]
        mask = np.isin(g, pool_ids)
        n_pool = int(mask.sum())
        if n_pool < 2 or len(set(g[mask])) < 2:
            continue
        c = xyz[sl]
        h, n = crossgene_stat(c[mask], g[mask])
        obs_hit += h; obs_n += n
        pool_genes_here = g[mask]
        for p in range(N_PERM):
            sel = rng.choice(len(g), n_pool, replace=False)
            h0, _ = crossgene_stat(c[sel], rng.permutation(pool_genes_here))
            null_hit[p] += h0
        null_n += n_pool
    if obs_n == 0 or null_n == 0:
        return None
    obs = obs_hit / obs_n
    null = null_hit / null_n
    return dict(observed=obs, null_mean=float(null.mean()),
        null_sd=float(null.std()), n_molecules=int(obs_n),
        enrichment=obs / max(null.mean(), 1e-9))

rows_t1, rows_t2a, rows_t2b, rows_t2c = [], [], [], []
src_centroids = {}
for t in set(type_of.values()):
    ids = [c for c, ty in type_of.items() if ty == t]
    sub = cells[cells.cell_id.isin(ids)]
    src_centroids[t] = cKDTree(sub[["x", "y"]].to_numpy())

for _, m in meta.iterrows():
    pair = m["pair"]; S = m["source"]; T = m["target"]
    pool = [g for g in m["pool"].split(";") if g]
    strict = [g for g in str(m["strict"]).split(";") if g and g != "nan"]
    pool_ids = np.array([gene_idx[g] for g in pool if g in gene_idx])
    strict_ids = np.array([gene_idx[g] for g in strict if g in gene_idx])
    ex = expo[expo.pair == pair]
    exposed = ex.cell_id[ex.exposure > 0].tolist()
    unexposed = ex.cell_id[ex.exposure == 0].tolist()

    r = patch_test(exposed, pool_ids)
    if r: rows_t1.append(dict(pair=pair, **r))

    r = patch_test(unexposed, pool_ids)
    if r: rows_t2a.append(dict(pair=pair, **r))

    # T2b: strict rate vs lateral distance to nearest S cell (e=0 cells).
    if len(strict_ids) >= 2 and len(unexposed) > 200:
        sub = cells[cells.cell_id.isin(unexposed)]
        d, _ = src_centroids[S].query(sub[["x", "y"]].to_numpy())
        sr = []
        for cid in sub.cell_id:
            ci = idx_of_cell[cid]
            g = garr[cell_slice(ci)]
            sr.append(int(np.isin(g, strict_ids).sum()))
        sub = sub.assign(strict=sr, dist=d)
        bins = [0, 30, 60, 120, 250, np.inf]
        lab = ["<30", "30-60", "60-120", "120-250", ">250"]
        grp = sub.groupby(pd.cut(sub.dist, bins, labels=lab), observed=True)
        prof = (grp.strict.sum() / grp.transcript_count.sum() * 1e3)
        rows_t2b.append(dict(pair=pair, **{f"rate_{k}": float(v)
            for k, v in prof.items()}, n_cells=len(sub)))

    # T2c: z-percentile of strict molecules within their cell (e=0 cells).
    if len(strict_ids):
        dev_obs, dev_null = [], []
        for cid in unexposed:
            ci = idx_of_cell.get(cid)
            if ci is None: continue
            sl = cell_slice(ci)
            g = garr[sl]; z = xyz[sl, 2]
            mask = np.isin(g, strict_ids)
            if not mask.any() or len(z) < 5:
                continue
            pct = (np.argsort(np.argsort(z)) + 0.5) / len(z)
            dev_obs.extend(np.abs(pct[mask] - 0.5))
            sel = rng.choice(len(z), int(mask.sum()), replace=False)
            dev_null.extend(np.abs(pct[sel] - 0.5))
        if len(dev_obs) >= 50:
            rows_t2c.append(dict(pair=pair, mean_absdev_strict=float(np.mean(dev_obs)),
                mean_absdev_null=float(np.mean(dev_null)), n=len(dev_obs)))

pd.DataFrame(rows_t1).to_csv(f"{AG}/results/t1_patch_exposed.csv", index=False)
pd.DataFrame(rows_t2a).to_csv(f"{AG}/results/t2a_patch_unexposed.csv", index=False)
pd.DataFrame(rows_t2b).to_csv(f"{AG}/results/t2b_distance_decay.csv", index=False)
pd.DataFrame(rows_t2c).to_csv(f"{AG}/results/t2c_z_position.csv", index=False)

for name, rows in [("T1 exposed", rows_t1), ("T2a unexposed", rows_t2a)]:
    df = pd.DataFrame(rows)
    if len(df):
        w = df.n_molecules
        print(f"{name}: pairs={len(df)} obs={np.average(df.observed, weights=w):.3f} "
              f"null={np.average(df.null_mean, weights=w):.3f} "
              f"median enrichment={df.enrichment.median():.2f}")
print("T2b pairs:", len(rows_t2b), " T2c pairs:", len(rows_t2c))
print("PATCH TESTS DONE")
