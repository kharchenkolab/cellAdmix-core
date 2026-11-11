# For each dataset: the panel size, the automatic neighborhood size, the
# median molecules per cell, and how far a k-nearest-neighbor neighborhood
# reaches into a cell: the median k-th within-cell neighbor distance as a
# share of the cell's equivalent radius (2D, from the molecule spread), and
# k as a share of the cell's molecules.
import json, sys
import numpy as np
import pyarrow.parquet as pq
from scipy.spatial import cKDTree

runs = {
    "breast rep1 (321 genes)": "/home/pkharchenko/cellAdmix/multisample_runs/breast/rep1/runs/fit_manual_rank23_ls_nmf",
    "pancreas left (377 genes)": "/home/pkharchenko/cellAdmix/multisample_runs/pancreas/left/runs/fit_manual_rank9_ls_nmf",
    "breast 5k crop (5101 genes)": "/home/pkharchenko/cellAdmix/cellAdmix-core/examples/xenium_breast_membrane_5k_full/out_medium_crop/runs/bench_seed1_ls_nmf",
}
ks = (10, 20, 50, 100)
rng = np.random.default_rng(1)
print("dataset | auto k | median molecules per cell | " + " | ".join(f"k={k}: reach/R, k/N" for k in ks))
for name, run in runs.items():
    m = json.load(open(run + "/run.json"))
    k_auto = m["pipeline_options"]["ncv_k"]
    cells = pq.read_table(run + "/cells.parquet", columns=["cell_idx", "transcript_count"]).to_pandas()
    med = int(cells.transcript_count.median())
    pool = cells[cells.transcript_count >= 30]
    samp = pool.sample(min(1500, len(pool)), random_state=1)
    mols = pq.read_table(run + "/molecules.parquet", columns=["x", "y", "cell_idx"],
                         filters=[("cell_idx", "in", samp.cell_idx.tolist())]).to_pandas()
    reach = {k: [] for k in ks}; frac = {k: [] for k in ks}
    for cid, g in mols.groupby("cell_idx"):
        P = g[["x", "y"]].to_numpy(); n = len(P)
        c = P.mean(0); R = np.sqrt(2 * np.mean(((P - c) ** 2).sum(1)))
        if R <= 0: continue
        tree = cKDTree(P)
        kmax = min(max(ks), n - 1)
        d, _ = tree.query(P, k=kmax + 1)
        for k in ks:
            if n > k:
                reach[k].append(np.median(d[:, k]) / R); frac[k].append(k / n)
    cols = [f"{np.median(reach[k]):.2f}, {np.median(frac[k]):.2f} (n={len(reach[k])})" for k in ks]
    print(f"{name} | {k_auto} | {med} | " + " | ".join(cols))
