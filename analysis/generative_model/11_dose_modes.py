"""Real-data comparison of the dose prior's estimation modes.

Dataset from GM_DATASET (pancreas, nsclc, breast); arms from the command
line (default: all arms on pancreas, the production arms plus two
validation arms elsewhere). The audit-based scoring needs the dataset's
cell x gene matrix in dense form and is skipped on breast 5K.

Fits the packaged C++ model on the harness problem specification (the same
counts, pools, strict tiers and factor alignment as 08_cpp_parity.py) under
the combinations of dose_mode (marginal / joint), dose_ladder (off / 30-60-
120-240) and the trend guard, in the production configuration (dose guided
by the full pools) and the validation configuration (guided by the A halves
only, scored on the held-out B halves by 02_eval_gates.R), plus the built-in
covariate-permutation control.

For every arm it writes the removed-molecule matrix for the gate harness
(data/gm_removed_dose_<arm>.mtx), and scores the correction with the audit:
per-pair sensitivity, own-marker false removal, the exposure-linked excess
removed in the unexposed bin (content above the ambient reference in cells
with no source among their 15 nearest), and the retention of the flagged
induced genes' excess.

Outputs: results/gm_dose_arms_summary.csv, results/gm_dose_arms_pairs.csv,
figures/dose_modes_*.png.
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scipy.io
import scipy.sparse as sp

from celladmix import _core

import json

GM = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
DATASET = os.environ.get("GM_DATASET", "pancreas")
SUF = "" if DATASET == "pancreas" else f"_{DATASET}"
PREFIX = "" if DATASET == "pancreas" else f"{DATASET}_"
DATA = f"{GM}/data{SUF}"
_EX = "/home/pkharchenko/cellAdmix/cellAdmix-core/examples"
EX = {"pancreas": f"{_EX}/xenium_pancreas_membrane_377_full",
      "nsclc": f"{_EX}/cosmx_nsclc_giotto",
      "breast": f"{_EX}/xenium_breast_membrane_5k_full"}[DATASET]
RUN = {"pancreas": f"{EX}/out/runs/fit_manual_rank9_invsqrt_kl",
       "nsclc": f"{EX}/out/runs/fit_manual_rank8_ls_nmf",
       "breast": f"{EX}/out/runs/fit_manual_rank8_invsqrt_kl"}[DATASET]
DO_AUDIT = DATASET != "breast"
LADDER = [30, 60, 120, 240]
THREADS = 12

# ---- harness inputs -----------------------------------------------------------
genes = [l.strip() for l in open(f"{DATA}/genes.txt")]
cells = [l.strip() for l in open(f"{DATA}/cells.txt")]
gene_of = {g: i for i, g in enumerate(genes)}
counts = sp.csc_matrix(scipy.io.mmread(f"{DATA}/counts.mtx"))
ct = pd.read_csv(f"{DATA}/cell_types.csv")
type_of_cell = dict(zip(ct.cell_id, ct.cell_type))
types = sorted(set(ct.cell_type.dropna()))
tcode = {t: i for i, t in enumerate(types)}
type_codes = [tcode.get(type_of_cell.get(c), -1) for c in cells]
cells_tbl = pq.read_table(f"{RUN}/cells.parquet",
                          columns=["cell_id", "x", "y"]).to_pandas()
xy = cells_tbl.set_index("cell_id").reindex(cells)
pairs_df = pd.read_csv(f"{DATA}/pairs.csv")
ft = pd.read_csv("/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/"
                 f"audit_guided/results/factor_types{SUF}.csv")
alignment = [-1] * int(json.load(open(f"{RUN}/run.json"))["n_factors"])
for r in ft.itertuples():
    alignment[int(r.factor) - 1] = tcode[r.cell_type]


def pair_specs(validation):
    out = []
    for r in pairs_df.itertuples():
        d = dict(source_type=tcode[r.source], target_type=tcode[r.target],
                 pool=[gene_of[g] for g in r.pool.split(";") if g],
                 strict=[gene_of[g] for g in str(r.strict).split(";")
                         if g and g != "nan"])
        if validation:
            d["guide"] = [gene_of[g] for g in r.A.split(";") if g]
        out.append(d)
    return out


# every arm sets the three dose options explicitly, so the table does not
# depend on the package defaults
ARMS = {
    "marginal_z0": dict(dose_mode="marginal", dose_ladder=[], dose_min_z=0.0),
    "marginal": dict(dose_mode="marginal", dose_ladder=[]),
    "joint": dict(dose_mode="joint", dose_ladder=[]),
    "ladder": dict(dose_mode="marginal", dose_ladder=LADDER),
    "joint_ladder": dict(dose_mode="joint", dose_ladder=LADDER),
    "joint_ladder_control": dict(dose_mode="joint", dose_ladder=LADDER,
                                 use_ambient=False, permute_covariates_seed=11),
    "v_marginal_z0": dict(dose_mode="marginal", dose_ladder=[], dose_min_z=0.0),
    "v_marginal": dict(dose_mode="marginal", dose_ladder=[]),
    "v_joint": dict(dose_mode="joint", dose_ladder=[]),
    "v_ladder": dict(dose_mode="marginal", dose_ladder=LADDER),
    "v_joint_ladder": dict(dose_mode="joint", dose_ladder=LADDER),
}
only = sys.argv[1:]
if not only and DATASET != "pancreas":
    only = ["marginal", "joint", "ladder", "joint_ladder", "joint_ladder_control",
            "v_marginal", "v_joint_ladder"]

base_args = dict(
    counts_indptr=[int(v) for v in counts.indptr],
    counts_indices=[int(v) for v in counts.indices],
    counts_values=[float(v) for v in counts.data],
    n_genes=len(genes), cell_ids=cells,
    x=[float(v) for v in xy["x"].to_numpy()],
    y=[float(v) for v in xy["y"].to_numpy()],
    type_codes=[int(v) for v in type_codes], n_types=len(types),
    molecules_parquet=f"{RUN}/molecules.parquet",
    cells_parquet=f"{RUN}/cells.parquet",
    factor_to_type=alignment, programs_flat=[], program_type=[])

fits = {}
pair_rows = []
for arm, opts in ARMS.items():
    if only and arm not in only:
        continue
    validation = arm.startswith("v_")
    t0 = time.time()
    res = _core.fit_generative(
        **base_args, pair_specs=pair_specs(validation),
        options=dict(num_threads=THREADS, neighbor_k=15, **opts))
    seconds = time.time() - t0
    removed = np.asarray(res["removed"], dtype=float)
    rem = sp.csc_matrix((removed, counts.indices.copy(), counts.indptr.copy()),
                        shape=counts.shape)
    scipy.io.mmwrite(f"{DATA}/gm_removed_dose_{arm}.mtx", rem)
    fits[arm] = dict(removed=removed, seconds=seconds,
                     n_induced=len(res["induced"]), induced=res["induced"])
    for s in res["pairs"]:
        r = pairs_df.iloc[s["pair"]]
        pair_rows.append(dict(arm=arm, pair=r.pair, source=r.source, target=r.target,
            prior=s["prior_molecules"], posterior=s["posterior_molecules"],
            induced=s["induced_molecules"],
            prior_unexposed=s["prior_unexposed_molecules"],
            cross_attributed=s["cross_attributed_molecules"]))
    print(f"{arm}: {seconds:.0f}s, removed {removed.sum():,.0f}, "
          f"induced genes {len(res['induced'])}", flush=True)

# ---- audit-based scoring ------------------------------------------------------
sys.path.insert(0, "/home/pkharchenko/cellAdmix/cellAdmix-core/python/src")
from celladmix import CellAdmix, read_annotation          # noqa: E402
from celladmix.audit import _bin_rates, _excess_vs_ref  # noqa: E402
from celladmix._score_utils import source_exposure_counts  # noqa: E402

if not DO_AUDIT:
    summary = pd.DataFrame([dict(arm=a, seconds=f["seconds"],
        removed=float(f["removed"].sum()), n_induced=f["n_induced"])
        for a, f in fits.items()])
    summary.to_csv(f"{GM}/results/{PREFIX}gm_dose_arms_summary.csv", index=False)
    pd.DataFrame(pair_rows).to_csv(f"{GM}/results/{PREFIX}gm_dose_arms_pairs.csv", index=False)
    print(summary.to_string(index=False))
    print("DOSE MODES DONE (no audit scoring on this dataset)")
    sys.exit(0)

os.chdir(EX)
if DATASET == "pancreas":
    tbl = pd.read_csv("annotations/annotation.csv.gz")
    ann = read_annotation(tbl, annotation_col="merged_annotation")
    ds = CellAdmix("data", output_dir="out", annotation=ann, num_threads=10)
    fit = ds.fit(nmf_variant="invsqrt_kl", verbose=False)
else:
    meta = pd.read_csv("prepared/cell_metadata_all.csv.gz")
    ann = pd.Series(meta.cell_type_coarse.values, index=meta.cell.astype(str))
    ds = CellAdmix.attach_existing("out", annotation=ann, num_threads=10)
    fit = ds.load_fit("fit_manual_rank8_ls_nmf")
audit = fit.audit_admixture()
agenes = pd.Index(audit._genes)
acells = pd.Index(audit._cells)
gi = agenes.get_indexer(genes)
ci = acells.get_indexer(cells)
assert (gi >= 0).all() and (ci >= 0).all()


class Shim:
    def __init__(self, matrix):
        self._m = matrix

    def counts(self):
        return self._m, list(genes), list(cells)


# exposure of every audit cell to every source (for induced retention)
cf = fit.cell_factors()
expo_all, types_e = source_exposure_counts(cf, ds.annotation, neighbor_k=15)
expo_all = pd.DataFrame(expo_all, columns=types_e, index=cf["cell_id"].astype(str))
csc_before = audit._matrix.tocsc()
tt_all = np.asarray(csc_before.sum(axis=0)).ravel()
ctypes = audit._ctypes

ref_fit = fits.get("marginal_z0") or fits.get("marginal") or (next(iter(fits.values())) if fits else None)
ref_induced = ref_fit["induced"] if ref_fit else []
flagged = [(pairs_df.iloc[r["pair"]].source, pairs_df.iloc[r["pair"]].target,
            genes[r["gene"]]) for r in ref_induced]


def flagged_retention(after):
    kept = tot_x = 0.0
    for src, tgt, gene in flagged:
        t_cells = [c for c in audit._cells if ctypes.get(c) == tgt]
        e = expo_all.loc[t_cells, src].to_numpy()
        cols = acells.get_indexer(t_cells)
        tt = tt_all[cols]
        g = agenes.get_loc(gene)
        def ex(m):
            v = np.asarray(m[g][:, cols].todense()).ravel()
            return (v[e > 0].sum() - v[e == 0].sum() / max(tt[e == 0].sum(), 1.0)
                    * tt[e > 0].sum())
        b, a = ex(csc_before), ex(after)
        if b > 0:
            kept += a
            tot_x += b
    return kept / max(tot_x, 1.0)


summary_rows = []
audit_pair_rows = []
for arm, f in fits.items():
    rem_h = sp.csc_matrix((f["removed"], counts.indices.copy(), counts.indptr.copy()),
                          shape=counts.shape)
    # removed matrix in the audit's gene/cell order
    rem_a = sp.csc_matrix((len(agenes), len(acells)))
    rem_a = sp.lil_matrix(rem_a.shape)
    rem_a[np.ix_(gi, ci)] = rem_h.todense()
    rem_a = rem_a.tocsc()
    after = csc_before - rem_a
    after.data = np.maximum(after.data, 0.0)
    report = audit.evaluate(Shim(after[gi][:, ci]), warn_uncovered=False)
    s = report.summary()
    # unexposed-bin excess above the ambient reference, before and after
    un_before = un_after = ex_before = ex_after = 0.0
    for d in audit._pairs.values():
        if not d["detected"]:
            continue
        cols = d["cols"]
        tot = audit._totals[cols]
        mk_a = np.asarray(after[d["pool"]][:, cols].sum(axis=0)).ravel()
        rates_a = _bin_rates(mk_a, tot, d["bins"])
        rb = d["rates"]
        ref = d["reference_rate"]
        b0 = rb[rb["bin"] == "0"].iloc[0]
        a0 = rates_a[rates_a["bin"] == "0"].iloc[0]
        ub = max(b0["rate"] - ref, 0.0) * b0["totals"]
        ua = max(a0["rate"] - ref, 0.0) * a0["totals"]
        eb = _excess_vs_ref(rb[rb["bin"] != "0"], ref)
        ea = _excess_vs_ref(rates_a[rates_a["bin"] != "0"], ref)
        un_before += ub
        un_after += ua
        ex_before += eb
        ex_after += ea
        audit_pair_rows.append(dict(arm=arm, source=d["source"], target=d["target"],
            unexposed_excess_before=ub, unexposed_excess_after=ua,
            exposed_excess_before=eb, exposed_excess_after=ea))
    summary_rows.append(dict(arm=arm, seconds=f["seconds"],
        removed=float(f["removed"].sum()), n_induced=f["n_induced"],
        leakage_removed_overall=s["leakage_removed_overall"],
        median_pair_sensitivity=s["median_pair_sensitivity"],
        own_marker_false_removal=s["own_marker_false_removal"],
        exposed_excess_removed=1 - ex_after / max(ex_before, 1.0),
        unexposed_excess_removed=1 - un_after / max(un_before, 1.0),
        unexposed_excess_before=un_before,
        flagged_retention=flagged_retention(after)))
    print(f"{arm}: leakage removed {s['leakage_removed_overall']:.3f}, "
          f"unexposed-bin excess removed {summary_rows[-1]['unexposed_excess_removed']:.3f}, "
          f"flagged retention {summary_rows[-1]['flagged_retention']:.3f}, "
          f"own-marker {s['own_marker_false_removal']:.4f}", flush=True)

summary = pd.DataFrame(summary_rows)
pairs_out = pd.DataFrame(pair_rows)
audit_pairs = pd.DataFrame(audit_pair_rows)
summary.to_csv(f"{GM}/results/{PREFIX}gm_dose_arms_summary.csv", index=False)
pairs_out.to_csv(f"{GM}/results/{PREFIX}gm_dose_arms_pairs.csv", index=False)
audit_pairs.to_csv(f"{GM}/results/{PREFIX}gm_dose_arms_audit_pairs.csv", index=False)
print(summary.to_string(index=False))
print("DOSE MODES DONE")
