"""Compare the generative model's expression-program initializations.

For each dataset, the packaged model is fitted with each initialization -
the NMF fit's factor-labeled molecules, dose-weighted within-type
clustering, and a single pseudobulk profile per type - and scored on the
audit's own yardstick: excess-weighted sensitivity over the detected
pairs, own-marker false removal, total removal, entry-level agreement
with the factor-initialized fit, and retention of the induced genes the
factor-initialized model flags.

Usage: 09_init_comparison.py <pancreas|nsclc|breast>
Appends one row per initialization to results/init_comparison.csv.
"""
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, "/home/pkharchenko/cellAdmix/cellAdmix-core/python/src")
from celladmix import CellAdmix, read_annotation

DS = sys.argv[1]
EX = "/home/pkharchenko/cellAdmix/cellAdmix-core/examples"
OUT = ("/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model/"
       "results/init_comparison.csv")

if DS == "pancreas":
    os.chdir(f"{EX}/xenium_pancreas_membrane_377_full")
    tbl = pd.read_csv("annotations/annotation.csv.gz")
    ann = read_annotation(tbl, annotation_col="merged_annotation")
    ds = CellAdmix("data", output_dir="out", annotation=ann, num_threads=10)
    fit = ds.fit(nmf_variant="invsqrt_kl", verbose=False)
elif DS == "nsclc":
    os.chdir(f"{EX}/cosmx_nsclc_giotto")
    meta = pd.read_csv("prepared/cell_metadata_all.csv.gz")
    ann = pd.Series(meta.cell_type_coarse.values,
                    index=meta.cell.astype(str))
    # the tabular store was built by the R workflow; attach and load the
    # cached run rather than rebuilding
    ds = CellAdmix.attach_existing("out", annotation=ann, num_threads=10)
    fit = ds.load_fit("fit_manual_rank8_ls_nmf")
elif DS == "breast":
    os.chdir(f"{EX}/xenium_breast_membrane_5k_full")
    tbl = pd.read_csv("annotations/annotation.csv.gz")
    tbl = tbl[tbl["merged_annotation"] != "Ambiguous / low-quality"]
    ann = read_annotation(tbl, annotation_col="merged_annotation")
    ds = CellAdmix("data", output_dir="out_python", annotation=ann,
                   num_threads=10)
    fit = ds.fit(nmf_variant="invsqrt_kl", verbose=False)
else:
    raise SystemExit("unknown dataset")

audit = fit.audit_admixture()
total = float(audit._matrix.sum())

results = {}
rows = []
for label, init in [("nmf_factors", fit), ("clusters", "clusters"),
                    ("pseudobulk", "pseudobulk")]:
    t0 = time.time()
    model = audit.fit_generative(init=init, num_threads=12)
    t_fit = time.time() - t0
    correction = model.correct()
    matrix, _, _ = correction.counts()
    removed = np.asarray(audit._matrix.tocsc().data
                         - matrix.tocsc().data)
    results[label] = dict(model=model, removed=removed)
    report = audit.evaluate(correction, warn_uncovered=False)
    p = report.pairs()
    ok = p["sensitivity"].notna()
    sens = float((p.loc[ok, "sensitivity"] * p.loc[ok, "excess"]).sum()
                 / p.loc[ok, "excess"].sum())
    fr = report.false_removal()
    rows.append(dict(dataset=DS, init=label,
        removed=float(removed.sum()),
        removed_frac=float(removed.sum()) / total,
        sensitivity_weighted=sens,
        own_marker_false_removal=float(
            (fr["false_removal"] * fr["own_marker_molecules"]).sum()
            / fr["own_marker_molecules"].sum()),
        n_induced=len(model.induced),
        seconds=t_fit))
    print(f"{DS} {label}: removed {removed.sum():,.0f}, sensitivity {sens:.3f}, "
          f"{t_fit:.0f}s", flush=True)

# entry-level agreement and induced-gene retention relative to the
# factor-initialized fit
ref = results["nmf_factors"]["removed"]
ref_model = results["nmf_factors"]["model"]
genes = list(audit._genes)
gi = {g: i for i, g in enumerate(genes)}
csc = audit._matrix.tocsc()
ctypes = audit._ctypes

def flagged_retention(removed):
    kept = 0.0
    tot_x = 0.0
    from celladmix._score_utils import source_exposure_counts
    cells_xy = audit.fit.cell_factors()
    expo_all, types_e = source_exposure_counts(
        cells_xy, audit.fit.dataset.annotation, neighbor_k=15)
    expo_all = pd.DataFrame(expo_all, columns=types_e,
                            index=cells_xy["cell_id"].astype(str))
    after = csc.copy()
    after.data = csc.data - removed
    for r in ref_model.induced.itertuples():
        t_cells = [c for c in audit._cells if ctypes.get(c) == r.target]
        e = expo_all.loc[t_cells, r.source].to_numpy()
        cols = pd.Index(audit._cells).get_indexer(t_cells)
        tt = np.asarray(csc[:, cols].sum(axis=0)).ravel()
        g = gi[r.gene]
        def ex(m):
            v = np.asarray(m[g][:, cols].todense()).ravel()
            return (v[e > 0].sum()
                    - v[e == 0].sum() / max(tt[e == 0].sum(), 1.0)
                    * tt[e > 0].sum())
        b = ex(csc)
        a = ex(after)
        if b > 0:
            kept += a
            tot_x += b
    return kept / max(tot_x, 1.0)

for row in rows:
    removed = results[row["init"]]["removed"]
    keep = (removed > 0) | (ref > 0)
    row["entry_corr_vs_factors"] = float(
        np.corrcoef(removed[keep], ref[keep])[0, 1]) if keep.any() else 1.0
    row["flagged_retention"] = flagged_retention(removed)

out = pd.DataFrame(rows)
out.to_csv(OUT, mode="a", header=not os.path.exists(OUT), index=False)
print(out.to_string(index=False))
print("INIT COMPARISON DONE")
