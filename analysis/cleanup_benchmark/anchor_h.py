"""Anchor-based (separable) factor recovery from a persisted run's training NCV.

Rebuilds the run's exact training neighborhood-composition matrix, forms the
gene-gene co-occurrence table, selects k anchor genes by successive projection
(SPA), recovers topic-gene loadings by simplex regression of each gene's
conditional co-occurrence row onto the anchor rows, and writes H as CSV
(k rows x genes columns, header = gene names) for injection via
fit(nmf_fixed_h = ...).

Usage: anchor_h.py <run_dir> <k> <out_csv> [min_count]
"""
import sys
import json
import os
import numpy as np
import pyarrow.parquet as pq
from scipy.sparse import csr_matrix, save_npz, load_npz
from scipy.spatial import cKDTree
from scipy.optimize import nnls

run_dir, k, out_csv = sys.argv[1], int(sys.argv[2]), sys.argv[3]
min_count = float(sys.argv[4]) if len(sys.argv) > 4 else None

manifest = json.load(open(f"{run_dir}/run.json"))
genes = manifest["genes"]
G = len(genes)
ncv_k = int(manifest.get("pipeline_options", {}).get("ncv_k", 20)) or 20

cache_dir = os.path.join(os.path.dirname(os.path.abspath(out_csv)), "cache")
os.makedirs(cache_dir, exist_ok=True)
cache = os.path.join(cache_dir, os.path.basename(os.path.normpath(run_dir)) + "_ncv.npz")

if os.path.exists(cache):
    X = load_npz(cache)
    print(f"NCV loaded from cache: {X.shape}", flush=True)
else:
    tr = pq.read_table(f"{run_dir}/training_rows.parquet").to_pandas().sort_values("training_rank")
    train_ids = tr["obs_id"].to_numpy()
    mol = pq.read_table(f"{run_dir}/molecules.parquet",
        columns=["obs_id", "x", "y", "z", "gene_idx", "cell_idx"]).to_pandas()
    obs = mol["obs_id"].to_numpy()
    order = np.argsort(obs, kind="stable")
    pos = order[np.searchsorted(obs[order], train_ids)]
    assert np.array_equal(obs[pos], train_ids)
    cell = mol["cell_idx"].to_numpy()
    train_cells = np.unique(cell[pos]); train_cells = train_cells[train_cells >= 0]
    sub = mol.loc[np.isin(cell, train_cells)].reset_index(drop=True)
    sub_ids = sub["obs_id"].to_numpy()
    sorder = np.argsort(sub_ids, kind="stable")
    spos = sorder[np.searchsorted(sub_ids[sorder], train_ids)]
    gene_arr = sub["gene_idx"].to_numpy()
    xyz = sub[["x", "y", "z"]].to_numpy()
    cell_arr = sub["cell_idx"].to_numpy()
    rows, cols, vals = [], [], []
    for c in np.unique(cell_arr[spos]):
        cmask = cell_arr == c
        cidx = np.flatnonzero(cmask)
        tlocal = np.flatnonzero(np.isin(cidx, spos))
        if not len(tlocal):
            continue
        tree = cKDTree(xyz[cmask])
        k_eff = min(ncv_k, int(cmask.sum()))
        _, nn = tree.query(xyz[cmask][tlocal], k=k_eff)
        nn = np.atleast_2d(nn)
        cg = gene_arr[cmask]
        for irow, tl in enumerate(tlocal):
            u, cnt = np.unique(cg[nn[irow]], return_counts=True)
            rows.extend([int(cidx[tl])] * len(u)); cols.extend(u.tolist()); vals.extend(cnt.tolist())
    row_map = {v: i for i, v in enumerate(sorted(set(rows)))}
    X = csr_matrix((vals, ([row_map[r] for r in rows], cols)), shape=(len(row_map), G))
    save_npz(cache, X)
    print(f"NCV rebuilt: {X.shape}", flush=True)

Xd = np.asarray(X.todense(), dtype=float)
tot = Xd.sum(axis=0)
if min_count is None:
    min_count = max(50.0, 0.25 * tot[tot > 0].mean())
ok = tot >= min_count
print(f"anchor candidates: {int(ok.sum())}/{G} genes (min_count {min_count:.0f})", flush=True)

Q = Xd[:, ok].T @ Xd[:, ok]
np.fill_diagonal(Q, 0)
Qn = Q / np.maximum(Q.sum(axis=1, keepdims=True), 1e-12)
idx = np.flatnonzero(ok)
R = Qn.copy()
anchors = []
for _ in range(k):
    j = int(np.argmax((R ** 2).sum(axis=1)))
    anchors.append(j)
    u = R[j] / max(np.linalg.norm(R[j]), 1e-12)
    R = R - np.outer(R @ u, u)
print("anchors:", ", ".join(genes[idx[a]] for a in anchors), flush=True)

A = Qn[anchors]
C = np.zeros((Qn.shape[0], k))
for i in range(Qn.shape[0]):
    c, _ = nnls(A.T, Qn[i])
    s = c.sum()
    C[i] = c / s if s > 0 else 0
pg = Q.sum(axis=1) / max(Q.sum(), 1e-12)
H = np.zeros((k, G))
H[:, idx] = (C * pg[:, None]).T

with open(out_csv, "w") as f:
    f.write(",".join(genes) + "\n")
    for frow in H:
        f.write(",".join(f"{v:.8g}" for v in frow) + "\n")
print(f"wrote {out_csv}", flush=True)
