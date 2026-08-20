"""T1 v2 (non-saturating patchiness) and far-field reference quantification.

T1v2: distance from each pool molecule to the nearest DIFFERENT-pool-gene
molecule in the same cell, observed vs within-cell permutation null,
restricted to cells with 2-10 pool molecules (the saturation regime is
excluded). Reported as median distance ratio (null/observed; >1 = patchy)
and 1.0um co-location enrichment.

Far-field: per pair, the pool-gene rate in zero-exposure cells stratified
by lateral distance to the nearest source cell; the >250um stratum is the
candidate true ambient reference. Reports the implied increase in the
pair's excess estimate if the far-field reference replaced the pooled
zero-exposure reference.

Run from the pancreas example directory.
"""
import json
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.spatial import cKDTree

AG = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/audit_guided"
RUN = "out/runs/fit_manual_rank9_invsqrt_kl"
rng = np.random.default_rng(13)
N_PERM = 10

genes = json.load(open(f"{RUN}/run.json"))["genes"]
gene_idx = {g: i for i, g in enumerate(genes)}
mol = pq.read_table(f"{RUN}/molecules.parquet",
    columns=["x", "y", "z", "gene_idx", "cell_idx"]).to_pandas()
cells = pq.read_table(f"{RUN}/cells.parquet",
    columns=["cell_idx", "cell_id", "x", "y", "transcript_count"]).to_pandas()
idx_of_cell = dict(zip(cells.cell_id, cells.cell_idx))
ctypes = pd.read_csv(f"{AG}/results/cell_types.csv")
type_of = dict(zip(ctypes.cell_id, ctypes.cell_type))

mol = mol[mol.cell_idx >= 0]
mol = mol.iloc[np.argsort(mol.cell_idx.to_numpy(), kind="stable")].reset_index(drop=True)
cell_arr = mol.cell_idx.to_numpy()
starts = np.searchsorted(cell_arr, np.arange(cell_arr.max() + 2))
xyz = mol[["x", "y", "z"]].to_numpy()
garr = mol.gene_idx.to_numpy()

meta = pd.read_csv(f"{AG}/results/pair_meta.csv")
expo = pd.read_csv(f"{AG}/results/pair_exposure.csv.gz")

src_trees = {}
for t in set(type_of.values()):
    ids = [c for c, ty in type_of.items() if ty == t]
    sub = cells[cells.cell_id.isin(ids)]
    src_trees[t] = cKDTree(sub[["x", "y"]].to_numpy())

def nn_diff_stats(coords, gids):
    """Per-molecule distance to nearest different-gene molecule."""
    tree = cKDTree(coords)
    d, idx = tree.query(coords, k=min(len(coords), 6))
    out = np.full(len(coords), np.inf)
    for i in range(len(coords)):
        for j in range(1, d.shape[1]):
            if gids[idx[i, j]] != gids[i]:
                out[i] = d[i, j]
                break
    return out

rows = []
for _, m in meta.iterrows():
    pair = m["pair"]; S = m["source"]
    pool = [g for g in m["pool"].split(";") if g]
    pool_ids = np.array([gene_idx[g] for g in pool if g in gene_idx])
    ex = expo[expo.pair == m["pair"]]
    for state, cell_list in [("exposed", ex.cell_id[ex.exposure > 0]),
                             ("unexposed", ex.cell_id[ex.exposure == 0])]:
        obs_d, null_d, obs_hit, null_hit, n_tot = [], [], 0.0, 0.0, 0
        for cid in cell_list:
            ci = idx_of_cell.get(cid)
            if ci is None: continue
            sl = slice(starts[ci], starts[ci + 1])
            g = garr[sl]
            mask = np.isin(g, pool_ids)
            n_pool = int(mask.sum())
            if n_pool < 2 or n_pool > 10 or len(set(g[mask])) < 2:
                continue
            c = xyz[sl]
            d_obs = nn_diff_stats(c[mask], g[mask])
            obs_d.extend(d_obs[np.isfinite(d_obs)])
            obs_hit += float((d_obs <= 1.0).sum()); n_tot += n_pool
            for _ in range(N_PERM):
                sel = rng.choice(len(g), n_pool, replace=False)
                d0 = nn_diff_stats(c[sel], rng.permutation(g[mask]))
                null_d.extend(d0[np.isfinite(d0)])
                null_hit += float((d0 <= 1.0).sum())
        if len(obs_d) < 100:
            continue
        rows.append(dict(pair=pair, state=state,
            median_nn_obs=float(np.median(obs_d)),
            median_nn_null=float(np.median(null_d)),
            dist_ratio=float(np.median(null_d) / max(np.median(obs_d), 1e-9)),
            coloc1um_obs=obs_hit / n_tot,
            coloc1um_null=null_hit / (n_tot * N_PERM),
            coloc_enrichment=(obs_hit / n_tot) /
                max(null_hit / (n_tot * N_PERM), 1e-9),
            n_molecules=n_tot))
t1v2 = pd.DataFrame(rows)
t1v2.to_csv(f"{AG}/results/t1v2_patch.csv", index=False)

# --- far-field reference impact ---
ff_rows = []
for _, m in meta.iterrows():
    pair = m["pair"]; S = m["source"]
    pool_ids = np.array([gene_idx[g] for g in m["pool"].split(";") if g in gene_idx])
    ex = expo[expo.pair == pair]
    e_map = dict(zip(ex.cell_id, ex.exposure))
    sub = cells[cells.cell_id.isin(ex.cell_id)].copy()
    d, _ = src_trees[S].query(sub[["x", "y"]].to_numpy())
    sub["dist"] = d
    sub["e"] = sub.cell_id.map(e_map)
    pool_counts = []
    for cid in sub.cell_id:
        ci = idx_of_cell[cid]
        g = garr[slice(starts[ci], starts[ci + 1])]
        pool_counts.append(int(np.isin(g, pool_ids).sum()))
    sub["pool"] = pool_counts
    e0 = sub[sub.e == 0]
    far = e0[e0.dist > 250]
    r0_pooled = e0.pool.sum() / max(e0.transcript_count.sum(), 1)
    r0_far = (far.pool.sum() / max(far.transcript_count.sum(), 1)
              if len(far) >= 100 else np.nan)
    exposed = sub[sub.e > 0]
    ex_pooled = max(exposed.pool.sum() - r0_pooled * exposed.transcript_count.sum(), 0)
    ex_far = (max(exposed.pool.sum() - r0_far * exposed.transcript_count.sum(), 0)
              + max((e0.pool.sum() - r0_far * e0.transcript_count.sum()), 0)
              if np.isfinite(r0_far) else np.nan)
    ff_rows.append(dict(pair=pair, r0_pooled=r0_pooled * 1e3,
        r0_far=(r0_far * 1e3 if np.isfinite(r0_far) else np.nan),
        excess_current=ex_pooled, excess_farfield=ex_far,
        inflation=(ex_far / ex_pooled if ex_pooled > 0 and np.isfinite(ex_far)
                   else np.nan),
        n_far_cells=len(far)))
ff = pd.DataFrame(ff_rows)
ff.to_csv(f"{AG}/results/farfield_reference.csv", index=False)

print("T1v2 summary (exposed):")
sub = t1v2[t1v2.state == "exposed"]
print(f"  pairs={len(sub)} median dist_ratio={sub.dist_ratio.median():.2f} "
      f"median coloc enrichment={sub.coloc_enrichment.median():.2f}")
sub = t1v2[t1v2.state == "unexposed"]
print("T1v2 summary (unexposed):")
print(f"  pairs={len(sub)} median dist_ratio={sub.dist_ratio.median():.2f} "
      f"median coloc enrichment={sub.coloc_enrichment.median():.2f}")
ok = ff[np.isfinite(ff.inflation)]
print(f"far-field: pairs with >=100 far cells: {len(ok)}; "
      f"median excess inflation {ok.inflation.median():.2f}x; "
      f"total excess {ok.excess_current.sum():.0f} -> {ok.excess_farfield.sum():.0f}")
print("V2 DONE")
