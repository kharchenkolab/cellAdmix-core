"""T3: the primary admixture/induction discriminator - profile proportionality.

Admixture samples the source's cytoplasmic transcriptome, so the per-gene
exposure-linked excess in target cells should be proportional to the
source's cytoplasmic profile psi_cyto (one shared constant per pair).
Induction is gene-selective: excess far above the proportional share on
specific genes.

Level A (gene): per pair, per-gene excess in exposed target cells versus
psi variants (cytoplasmic / whole-cell / nuclear) over source-owned genes;
reports which profile fits best, the proportional fraction, and outliers.

Level B (patch): gene composition of source-factor-labeled molecules
inside exposed target cells versus psi_cyto.

Run from the pancreas example directory. Requires results/factor_types.csv
(factor -> best type map).
"""
import json
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

AG = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/audit_guided"
RUN = "out/runs/fit_manual_rank9_invsqrt_kl"

genes = json.load(open(f"{RUN}/run.json"))["genes"]
G = len(genes)
mol = pq.read_table(f"{RUN}/molecules.parquet",
    columns=["gene_idx", "cell_idx", "overlaps_nucleus", "factor_label"]).to_pandas()
cells = pq.read_table(f"{RUN}/cells.parquet",
    columns=["cell_idx", "cell_id", "transcript_count"]).to_pandas()
idx_of_cell = dict(zip(cells.cell_id, cells.cell_idx))
ctypes = pd.read_csv(f"{AG}/results/cell_types.csv")
type_of = dict(zip(ctypes.cell_id, ctypes.cell_type))
cells["cell_type"] = cells.cell_id.map(type_of)
mol = mol[mol.cell_idx >= 0]
mol = mol.merge(cells[["cell_idx", "cell_type"]], on="cell_idx", how="left")

meta = pd.read_csv(f"{AG}/results/pair_meta.csv")
expo = pd.read_csv(f"{AG}/results/pair_exposure.csv.gz")
fmap = pd.read_csv(f"{AG}/results/factor_types.csv")
factor_of_type = dict(zip(fmap.cell_type, fmap.factor))

def profile(sub):
    cnt = np.bincount(sub.gene_idx.to_numpy(), minlength=G).astype(float)
    return cnt / max(cnt.sum(), 1)

profiles = {}
for t in set(type_of.values()):
    s = mol[mol.cell_type == t]
    profiles[t] = dict(
        whole=profile(s),
        cyto=profile(s[s.overlaps_nucleus == 0]),
        nuc=profile(s[s.overlaps_nucleus == 1]))

own_type = np.array([max(profiles, key=lambda t: profiles[t]["whole"][g])
    for g in range(G)])

gene_counts_by_cellset = {}
def counts_in(cell_ids):
    key = id(cell_ids)
    sub = mol[mol.cell_idx.isin([idx_of_cell[c] for c in cell_ids
        if c in idx_of_cell])]
    return np.bincount(sub.gene_idx.to_numpy(), minlength=G).astype(float), sub

rows_a, rows_b, outlier_rows = [], [], []
for _, m in meta.iterrows():
    pair = m["pair"]; S = m["source"]; T = m["target"]
    ex = expo[expo.pair == pair]
    exposed = ex.cell_id[ex.exposure > 0].tolist()
    unexposed = ex.cell_id[ex.exposure == 0].tolist()
    c_exp, sub_exp = counts_in(exposed)
    c_un, _ = counts_in(unexposed)
    tot_exp = c_exp.sum(); tot_un = c_un.sum()
    if tot_exp < 1e4 or tot_un < 1e4:
        continue
    r0 = c_un / tot_un
    excess = c_exp - r0 * tot_exp
    var = c_exp + (tot_exp / tot_un) ** 2 * c_un + 1
    zscore = excess / np.sqrt(var)
    gsel = np.flatnonzero((own_type == S) & ((c_exp + c_un) >= 20))
    if len(gsel) < 8:
        continue
    res = dict(pair=pair, n_genes=len(gsel),
        excess_total=float(excess[gsel].sum()))
    ex_v = excess[gsel]
    for pv in ("cyto", "whole", "nuc"):
        psi = profiles[S][pv][gsel]
        w = 1 / var[gsel]
        slope = max((w * ex_v * psi).sum() / max((w * psi * psi).sum(), 1e-12), 0)
        pred = slope * psi
        ss_res = ((ex_v - pred) ** 2 / var[gsel]).sum()
        ss_tot = ((ex_v - ex_v.mean()) ** 2 / var[gsel]).sum()
        res[f"r2_{pv}"] = 1 - ss_res / max(ss_tot, 1e-12)
        res[f"corr_{pv}"] = float(np.corrcoef(ex_v, psi)[0, 1])
        if pv == "cyto":
            resid_z = (ex_v - pred) / np.sqrt(var[gsel])
            prop_frac = float(np.minimum(ex_v, pred).clip(0).sum() /
                max(ex_v.clip(0).sum(), 1e-9))
            res["proportional_fraction"] = prop_frac
            for i in np.argsort(-resid_z)[:3]:
                if resid_z[i] > 3 and ex_v[i] > 50:
                    outlier_rows.append(dict(pair=pair, gene=genes[gsel[i]],
                        excess=float(ex_v[i]), expected=float(pred[i]),
                        resid_z=float(resid_z[i])))
    rows_a.append(res)

    f_src = factor_of_type.get(S)
    if f_src is not None:
        patch = sub_exp[sub_exp.factor_label == f_src - 1] \
            if (sub_exp.factor_label.min() == 0) else \
            sub_exp[sub_exp.factor_label == f_src]
        if len(patch) >= 200:
            comp = profile(patch)
            gs = np.flatnonzero(comp + profiles[S]["cyto"] > 0)
            rows_b.append(dict(pair=pair, n_patch_molecules=len(patch),
                cos_cyto=float(np.dot(comp[gs], profiles[S]["cyto"][gs]) /
                    max(np.linalg.norm(comp[gs]) *
                        np.linalg.norm(profiles[S]["cyto"][gs]), 1e-12)),
                cos_target=float(np.dot(comp[gs], profiles[T]["whole"][gs]) /
                    max(np.linalg.norm(comp[gs]) *
                        np.linalg.norm(profiles[T]["whole"][gs]), 1e-12))))

a = pd.DataFrame(rows_a); a.to_csv(f"{AG}/results/t3_gene_proportionality.csv", index=False)
b = pd.DataFrame(rows_b); b.to_csv(f"{AG}/results/t3_patch_composition.csv", index=False)
o = pd.DataFrame(outlier_rows); o.to_csv(f"{AG}/results/t3_outlier_genes.csv", index=False)

print(f"T3A pairs={len(a)}  corr_cyto med={a.corr_cyto.median():.3f} "
      f"corr_whole med={a.corr_whole.median():.3f} corr_nuc med={a.corr_nuc.median():.3f}")
print(f"    proportional_fraction med={a.proportional_fraction.median():.3f}")
print(f"T3B pairs={len(b)}  cos(source cyto) med={b.cos_cyto.median():.3f} "
      f"cos(target) med={b.cos_target.median():.3f}")
print(f"outlier genes flagged: {len(o)}")
print(o.to_string(index=False) if len(o) else "  (none)")
print("T3 DONE")
