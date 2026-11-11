# Label-free spatial scale of admixture on breast rep1: for a target cell
# type, the molecules of genes its isolated cells essentially never express
# ("foreign markers", chosen from the pseudobulk ratio against the source
# type) are admixture with near certainty. Per target cell: how many such
# molecules it holds, their radial position, their clustering within the
# cell (share of foreign markers among a foreign molecule's k nearest
# within-cell neighbors, against the cell's share), and whether they sit on
# the side of the cell facing the nearest source cell.
import json, sys
import numpy as np, pandas as pd
import pyarrow.parquet as pq
from scipy.spatial import cKDTree

run = "/home/pkharchenko/cellAdmix/multisample_runs/breast/rep1/runs/fit_manual_rank23_ls_nmf"
ann = pd.read_csv("/home/pkharchenko/cellAdmix/wt-dev_multisample/examples/xenium_breast_replicates/annotations/rep1_annotation.csv.gz")
genes = np.array(json.load(open(run + "/run.json"))["genes"])
cells = pq.read_table(run + "/cells.parquet", columns=["cell_idx", "cell_id", "x", "y", "transcript_count"]).to_pandas()
cells["cell_id"] = cells.cell_id.astype(str); ann["cell_id"] = ann.cell_id.astype(str)
cells = cells.merge(ann, on="cell_id", how="left")
type_of = dict(zip(cells.cell_idx, cells.cell_type))
mols = pq.read_table(run + "/molecules.parquet", columns=["x", "y", "gene_idx", "cell_idx"]).to_pandas()
mols = mols[mols.cell_idx >= 0]
mols["cell_type"] = mols.cell_idx.map(type_of)
pb = mols.groupby(["cell_type", "gene_idx"]).size().unstack(fill_value=0)
pb = pb.div(pb.sum(1), axis=0)  # per-type molecule share per gene
pairs = [("Stromal", "Invasive_Tumor"), ("Invasive_Tumor", "Stromal"), ("Endothelial", "Invasive_Tumor"), ("Stromal", "DCIS_1"), ("Myoepi_ACTA2+", "Invasive_Tumor")]
ks = (5, 10, 20, 50)
rng = np.random.default_rng(1)
cell_xy = cells.set_index("cell_idx")[["x", "y"]]
for source, target in pairs:
    ratio = pb.loc[source] / pb.loc[target].clip(lower=1e-7)
    markers = [g for g in pb.columns if ratio[g] >= 30 and pb.loc[source, g] >= 5e-4]
    if len(markers) < 3:
        print(f"{source} -> {target}: fewer than 3 markers"); continue
    tcells = cells[cells.cell_type == target]
    scells = cells[cells.cell_type == source]
    stree = cKDTree(scells[["x", "y"]].to_numpy())
    tm = mols[mols.cell_idx.isin(tcells.cell_idx)]
    tm = tm.assign(foreign=tm.gene_idx.isin(markers))
    per_cell = tm.groupby("cell_idx").foreign.agg(["sum", "size"])
    per_cell = per_cell[per_cell["size"] >= 30]
    dist_src, _ = stree.query(cell_xy.loc[per_cell.index].to_numpy(), k=1)
    near = per_cell[dist_src <= 20]  # target cells within 20 um of a source centroid
    far = per_cell[dist_src > 60]
    print(f"\n{source} -> {target}: {len(markers)} foreign markers ({', '.join(genes[markers][:8])}); target cells with >= 30 molecules: {len(per_cell):,}")
    print(f"  foreign molecules per cell: near source (<= 20 um, n={len(near):,}) median {near['sum'].median():.0f}, share {np.median(near['sum']/near['size']):.3f}; "
          f"far (> 60 um, n={len(far):,}) median {far['sum'].median():.0f}, share {np.median(far['sum']/far['size']):.3f}")
    use = near[near["sum"] >= 3].index.to_numpy()
    use = rng.choice(use, size=min(2000, len(use)), replace=False)
    sub = tm[tm.cell_idx.isin(use)]
    nb = {k: [] for k in ks}; base = []; rad_f = []; rad_a = []; facing = []
    for cid, g in sub.groupby("cell_idx"):
        P = g[["x", "y"]].to_numpy(); f = g.foreign.to_numpy(); n = len(P)
        c = P.mean(0); d = np.sqrt(((P - c) ** 2).sum(1)); R = np.sqrt(2 * np.mean(d ** 2))
        if R <= 0: continue
        rad_f.append(np.median(d[f]) / R); rad_a.append(np.median(d) / R); base.append(f.mean())
        _, j = stree.query(c, k=1); s = scells[["x", "y"]].to_numpy()[j]
        u = (s - c) / max(np.linalg.norm(s - c), 1e-9)
        proj = (P - c) @ u  # position along the axis toward the source cell
        facing.append((proj[f] > 0).mean() - (proj > 0).mean())
        tree = cKDTree(P); kmax = min(max(ks), n - 1)
        _, idx = tree.query(P[f], k=kmax + 1)
        for k in ks:
            nb[k].append(f[idx[:, 1:min(k, kmax) + 1]].mean())
    base = np.array(base)
    print(f"  cells used {len(base)}: median foreign share {np.median(base):.3f}; radial position / R: foreign {np.median(rad_f):.2f}, all {np.median(rad_a):.2f}; "
          f"excess share of foreign molecules on the half facing the source: {np.median(facing):+.2f}")
    print("  share of a foreign molecule's k nearest within-cell neighbors that are foreign (median), enrichment over the cell's share:")
    for k in ks:
        v = np.array(nb[k]); print(f"    k={k:>2}: {np.median(v):.2f}  enrichment {np.median(v / base):.1f}x")
