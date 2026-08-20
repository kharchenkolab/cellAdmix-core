"""T7: four-class spike-in ground truth for the discriminators.

Into real exposed target cells of one substrate pair, inject labeled
synthetic molecules:
  inplane   - multi-gene fragments sampled from a nearby source cell,
              placed at the target edge facing it (z preserved);
  outplane  - the same fragments placed anywhere over the target footprint
              with z pushed toward the slab surfaces;
  induction - three source-marker genes only, positions sampled from the
              target's own molecules (native placement), counts scaled
              with exposure;
  ambient   - single molecules of random pool genes, uniform placement.

Then score each class with the T1v2 patch statistic and run the T3
gene-level proportionality test on the spiked data to check that the
induction genes - and only they - are flagged as outliers.

Run from the pancreas example directory.
"""
import json
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.spatial import cKDTree

AG = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/audit_guided"
RUN = "out/runs/fit_manual_rank9_invsqrt_kl"
PAIR = "Exocrine epithelial -> Endothelial"
N_CELLS_PER_CLASS = 400
rng = np.random.default_rng(29)

genes = json.load(open(f"{RUN}/run.json"))["genes"]
G = len(genes)
gene_idx = {g: i for i, g in enumerate(genes)}
mol = pq.read_table(f"{RUN}/molecules.parquet",
    columns=["x", "y", "z", "gene_idx", "cell_idx"]).to_pandas()
cells = pq.read_table(f"{RUN}/cells.parquet",
    columns=["cell_idx", "cell_id", "x", "y", "transcript_count"]).to_pandas()
idx_of_cell = dict(zip(cells.cell_id, cells.cell_idx))
ctypes = pd.read_csv(f"{AG}/results/cell_types.csv")
type_of = dict(zip(ctypes.cell_id, ctypes.cell_type))
meta = pd.read_csv(f"{AG}/results/pair_meta.csv")
expo = pd.read_csv(f"{AG}/results/pair_exposure.csv.gz")

mol = mol[mol.cell_idx >= 0]
mol = mol.iloc[np.argsort(mol.cell_idx.to_numpy(), kind="stable")].reset_index(drop=True)
cell_arr = mol.cell_idx.to_numpy()
starts = np.searchsorted(cell_arr, np.arange(cell_arr.max() + 2))
xyz = mol[["x", "y", "z"]].to_numpy()
garr = mol.gene_idx.to_numpy()

m = meta[meta.pair == PAIR].iloc[0]
S, T = m.source, m.target
pool = [g for g in m.pool.split(";") if g in gene_idx]
pool_ids = np.array([gene_idx[g] for g in pool])
ex = expo[expo.pair == PAIR]
exposed = [c for c in ex.cell_id[ex.exposure > 0] if c in idx_of_cell]
e_of = dict(zip(ex.cell_id, ex.exposure))
s_cells = [c for c, t in type_of.items() if t == S and c in idx_of_cell]
s_xy = cells.set_index("cell_id").loc[s_cells, ["x", "y"]].to_numpy()
s_tree = cKDTree(s_xy)
zlo, zhi = np.percentile(xyz[:, 2], [10, 90])

def cell_mols(cid):
    sl = slice(starts[idx_of_cell[cid]], starts[idx_of_cell[cid]] +
        (starts[idx_of_cell[cid] + 1] - starts[idx_of_cell[cid]]))
    return xyz[sl], garr[sl]

def sample_fragment():
    """A 1.5um ball of molecules from a random source cell."""
    for _ in range(50):
        cid = s_cells[rng.integers(len(s_cells))]
        c, g = cell_mols(cid)
        if len(g) < 8:
            continue
        seed = rng.integers(len(g))
        d = np.linalg.norm(c - c[seed], axis=1)
        sel = np.flatnonzero(d <= 1.5)
        if len(sel) >= 4:
            sel = sel[:12]
            return c[sel] - c[sel].mean(0), g[sel]
    return None, None

records = []
target_cells = rng.permutation(exposed)
groups = np.array_split(target_cells[:4 * N_CELLS_PER_CLASS], 4)
induction_genes = pool_ids[:3]

for cls, cell_group in zip(["inplane", "outplane", "induction", "ambient"], groups):
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
                anchor = np.array([cx[0] + to_s[0] * r90, cx[1] + to_s[1] * r90,
                    cx[2]])
            else:
                ang = rng.uniform(0, 2 * np.pi)
                rad = rng.uniform(0, np.percentile(
                    np.linalg.norm(c[:, :2] - cx[:2], axis=1), 70))
                zpos = zlo if rng.random() < 0.5 else zhi
                anchor = np.array([cx[0] + rad * np.cos(ang),
                    cx[1] + rad * np.sin(ang), zpos])
            pos = fc + anchor
            for p, gg in zip(pos, fg):
                records.append((cid, cls, *p, gg))
        elif cls == "induction":
            n_add = 2 + 2 * min(int(e_of.get(cid, 1)), 4)
            sel = rng.choice(len(c), n_add)
            jit = rng.normal(0, 0.5, (n_add, 3))
            for p, gg in zip(c[sel] + jit,
                    rng.choice(induction_genes, n_add)):
                records.append((cid, cls, *p, gg))
        else:
            n_add = 2
            sel = rng.choice(len(c), n_add)
            jit = rng.normal(0, 3.0, (n_add, 3))
            for p, gg in zip(c[sel] + jit, rng.choice(pool_ids, n_add)):
                records.append((cid, cls, *p, gg))

spk = pd.DataFrame(records, columns=["cell_id", "cls", "x", "y", "z", "gene_idx"])
print("injected:", spk.groupby("cls").size().to_dict())

# --- discriminator 1: patch statistic per class (1um co-location with a
# different pool gene, among real+injected molecules of the cell) ---
res = []
for cls, sub in spk.groupby("cls"):
    hits = n = 0
    for cid, cm in sub.groupby("cell_id"):
        c0, g0 = cell_mols(cid)
        allc = np.vstack([c0, cm[["x", "y", "z"]].to_numpy()])
        allg = np.concatenate([g0, cm.gene_idx.to_numpy()])
        inj = np.arange(len(g0), len(allg))
        pm = np.flatnonzero(np.isin(allg, pool_ids))
        tree = cKDTree(allc[pm])
        gsel = allg[pm]
        pos_of = {int(p): k for k, p in enumerate(pm)}
        for i in inj:
            k = pos_of.get(int(i))
            if k is None:
                continue
            nb = tree.query_ball_point(allc[i], 1.0)
            n += 1
            if any(gsel[j] != allg[i] for j in nb if j != k):
                hits += 1
    res.append(dict(cls=cls, coloc1um=hits / max(n, 1), n=n))
d1 = pd.DataFrame(res)
print("\npatch statistic by class:")
print(d1.to_string(index=False))

# --- discriminator 2: T3 proportionality on the spiked pair ---
own_counts = np.zeros(G)
for cid in s_cells[:2000]:
    _, g = cell_mols(cid)
    own_counts += np.bincount(g, minlength=G)
psi = own_counts / own_counts.sum()

unexposed = [c for c in ex.cell_id[ex.exposure == 0] if c in idx_of_cell]
def gcounts(cell_list, add_spike=False):
    cnt = np.zeros(G); tot = 0
    for cid in cell_list:
        _, g = cell_mols(cid)
        cnt += np.bincount(g, minlength=G); tot += len(g)
    if add_spike:
        s = spk[spk.cell_id.isin(cell_list)]
        cnt += np.bincount(s.gene_idx.to_numpy(), minlength=G)
        tot += len(s)
    return cnt, tot

c_exp, t_exp = gcounts(exposed, add_spike=True)
c_un, t_un = gcounts(unexposed)
r0 = c_un / t_un
excess = c_exp - r0 * t_exp
var = c_exp + (t_exp / t_un) ** 2 * c_un + 1
gsel = np.flatnonzero((own_counts > 200) & ((c_exp + c_un) >= 20))
w = 1 / var[gsel]
p = psi[gsel]
slope = max((w * excess[gsel] * p).sum() / max((w * p * p).sum(), 1e-12), 0)
resid_z = (excess[gsel] - slope * p) / np.sqrt(var[gsel])
flag = [(genes[gsel[i]], round(float(resid_z[i]), 1))
    for i in np.argsort(-resid_z)[:6] if resid_z[i] > 3]
print("\nT3 outliers on spiked pair (expect the 3 induction genes):")
print("  injected induction genes:", [genes[g] for g in induction_genes])
print("  flagged:", flag)
d1.to_csv(f"{AG}/results/t7_spikein_patch.csv", index=False)
pd.DataFrame(flag, columns=["gene", "resid_z"]).to_csv(
    f"{AG}/results/t7_spikein_outliers.csv", index=False)
print("T7 DONE")
