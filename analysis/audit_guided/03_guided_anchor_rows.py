"""Audit-guided soft anchor rows from a run's training NCV co-occurrence.

Instead of SPA-selecting anchor genes (anchor_h.py), the anchor gene sets are
supplied - the guide halves of the audit's source-marker pools. Each labeled
gene set becomes one H row: the occurrence-weighted average of its genes'
shrunk conditional co-occurrence rows, prevalence-scaled like anchor_h.py.
Loadings therefore spread to co-occurring genes the guide never named - the
soft-anchor generalization the A/B validation tests.

Usage: 03_guided_anchor_rows.py <run_dir> <guides.json> <out_csv>
  guides.json: {"label": ["GENE1", ...], ...}
"""
import sys
import json
import os
import hashlib
import numpy as np
import pyarrow.parquet as pq
from scipy.sparse import csr_matrix, save_npz, load_npz
from scipy.spatial import cKDTree

run_dir, guides_json, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]
guides = json.load(open(guides_json))

manifest = json.load(open(f"{run_dir}/run.json"))
genes = manifest["genes"]
G = len(genes)
gene_idx = {g: i for i, g in enumerate(genes)}
ncv_k = int(manifest.get("pipeline_options", {}).get("ncv_k", 20)) or 20

cache_dir = os.path.join(os.path.dirname(os.path.abspath(out_csv)), "cache")
os.makedirs(cache_dir, exist_ok=True)
run_key = hashlib.sha1(os.path.abspath(run_dir).encode()).hexdigest()[:10]
cache = os.path.join(cache_dir,
    os.path.basename(os.path.normpath(run_dir)) + "_" + run_key + "_ncv.npz")

if os.path.exists(cache):
    X = load_npz(cache)
    assert X.shape[1] == G
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
Q = Xd.T @ Xd
np.fill_diagonal(Q, 0)
Qn = Q / np.maximum(Q.sum(axis=1, keepdims=True), 1e-12)
p_global = Q.sum(axis=0) / max(Q.sum(), 1e-12)
tau = 300.0
lam = (tot / (tot + tau))[:, None]
Qn = lam * Qn + (1.0 - lam) * p_global[None, :]

rows_out = {}
for label, gset in guides.items():
    idx = [gene_idx[g] for g in gset if g in gene_idx and tot[gene_idx[g]] > 0]
    if not idx:
        print(f"WARNING: no usable guide genes for {label}", flush=True)
        continue
    w = tot[idx] / tot[idx].sum()
    profile = (w[:, None] * Qn[idx]).sum(axis=0)
    # Guide genes barely co-occur with themselves in Qn (diagonal removed);
    # restore their own presence at the level of their occurrence share so
    # the anchor factor can claim guide-gene molecules directly too.
    own = tot[idx] / max(tot.sum(), 1e-12)
    profile[idx] = np.maximum(profile[idx], own * profile.sum())
    rows_out[label] = profile / max(profile.sum(), 1e-12)
    top = np.argsort(-rows_out[label])[:12]
    print(f"{label}: {len(idx)} guide genes; top loadings: "
          + ", ".join(genes[t] for t in top), flush=True)

with open(out_csv, "w") as f:
    f.write("label," + ",".join(genes) + "\n")
    for label, row in rows_out.items():
        f.write(label + "," + ",".join(f"{v:.8g}" for v in row) + "\n")
print(f"wrote {out_csv}", flush=True)
