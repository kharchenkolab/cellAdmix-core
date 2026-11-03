"""Figure 2 of docs/generative.md: the model fitted on the CosMx NSCLC
section (fine labels), shown on admixture from malignant cells into
fibroblasts. (a) The admixture fraction of fibroblast molecules by
exposure stratum: the malignant marker content observed, the model's
prior (the dose-response) and the fitted fraction. (b-d) Three fibroblast
genes, observed and retained after correction: a malignant marker that is
admixture only (KRT19), a malignant gene induced in fibroblasts above its
admixture expectation (MMP7), and a fibroblast gene repressed near
malignant cells (DCN).

Strata: fibroblasts with malignant cells among their 15 nearest neighbors
are binned by that number (1, 2, 3, 4+); fibroblasts with none are placed
by how far the neighborhood must be widened before a malignant cell
appears (none among 15, 30, 60, 120, 240 nearest).

The fit (package defaults, initialized from the example's NMF fit) is
saved under data_nsclc_fine/ on the first run and reused afterwards.
Writes docs/figures/generative_fig2.png.
"""
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

HERE = Path(__file__).resolve().parent
CORE = HERE.parents[1]
sys.path.insert(0, str(CORE / "python" / "src"))
from celladmix import CellAdmix  # noqa: E402
from celladmix._score_utils import source_exposure_counts, source_nearest_distance  # noqa: E402
from celladmix.audit import _interface_profile  # noqa: E402

DATA = HERE / "data_nsclc_fine"
OUT = CORE / "docs" / "figures" / "generative_fig2.png"
SOURCE, TARGET = "malignant", "fibroblast"
GENES = [("KRT19", "admixture only"), ("MMP7", "induced above admixture"), ("DCN", "repressed")]
LADDER = [15, 30, 60, 120, 240]
RUNGS = ["0/240", "0/120", "0/60", "0/30", "0/15", "1", "2", "3", "4+"]
RED, BLUE, DARK = "#c0392b", "#2e86c1", "#2c3e50"
THREADS = int(os.environ.get("CELLADMIX_THREADS", "8"))


def fit_and_save():
    """Fit the model on the NSCLC example with its fine labels (T and B cell
    subsets merged) and save the counts before and after correction, the
    per-cell composition and the cell table."""
    os.chdir(CORE / "examples" / "cosmx_nsclc_giotto")
    meta = pd.read_csv("prepared/cell_metadata_all.csv.gz")
    lab = meta.cell_type.replace({"CD4+ T cells": "T cell", "CD8+ T cells": "T cell",
                                  "Treg": "T cell", "B cells": "B cell"})
    ann = pd.Series(lab.values, index=meta.cell.astype(str))
    ds = CellAdmix.attach_existing("out", annotation=ann, num_threads=THREADS)
    fit = ds.load_fit("fit_manual_rank8_ls_nmf")
    model = fit.audit_admixture().fit_generative(init=fit, num_threads=THREADS)
    before, genes, cells = fit.counts()
    after, _, _ = model.correct().counts()
    DATA.mkdir(exist_ok=True)
    sp.save_npz(DATA / "before.npz", sp.csc_matrix(before))
    sp.save_npz(DATA / "after.npz", sp.csc_matrix(after))
    (DATA / "genes.txt").write_text("\n".join(map(str, genes)) + "\n")
    model.composition(SOURCE, TARGET).to_csv(DATA / "composition.csv.gz", index=False)
    model.induced.to_csv(DATA / "induced.csv", index=False)
    ct = pd.DataFrame(dict(cell_id=[str(c) for c in cells]))
    xy = fit.cells()[["cell_id", "x", "y"]].astype({"cell_id": str})
    ct = ct.merge(xy, on="cell_id", how="left")
    ct["annotation"] = ann.reindex(ct.cell_id).to_numpy()
    ct.to_csv(DATA / "cells.csv.gz", index=False)


if not (DATA / "after.npz").exists():
    fit_and_save()
before = sp.load_npz(DATA / "before.npz").tocsc()
after = sp.load_npz(DATA / "after.npz").tocsc()
genes = pd.Index((DATA / "genes.txt").read_text().split())
ct = pd.read_csv(DATA / "cells.csv.gz", dtype={"cell_id": str})
ann = ct.annotation.to_numpy()
annotation = pd.Series(ann, index=ct.cell_id).dropna()
cells = ct[["cell_id", "x", "y"]]

# exposure strata of the target cells
t_cols = np.flatnonzero(ann == TARGET)
expo = {}
for K in LADDER:
    counts, types = source_exposure_counts(cells, annotation, neighbor_k=K)
    expo[K] = np.asarray(counts)[t_cols, [str(t) for t in types].index(SOURCE)]
rung = np.full(len(t_cols), "0/15", dtype=object)
for n, lab in ((1, "1"), (2, "2"), (3, "3")):
    rung[expo[15] == n] = lab
rung[expo[15] >= 4] = "4+"
for K in LADDER[1:]:
    rung[(expo[15] == 0) & (expo[K] == 0)] = f"0/{K}"
totals = np.asarray(before.sum(axis=0)).ravel()[t_cols]
shown = [r for r in RUNGS if (rung == r).sum() >= 100]
xb = np.arange(len(shown))
ticks = [r.replace("0/", "none\nin ") for r in shown]
short_ticks = [r.replace("0/", "0 in ") for r in shown]
n_cells = np.array([(rung == r).sum() for r in shown])
T = np.array([totals[rung == r].sum() for r in shown])


def by_stratum(per_cell):
    """Molecule-weighted mean of a per-cell fraction in each stratum."""
    return np.array([(per_cell * totals)[rung == r].sum() / T[i] for i, r in enumerate(shown)])


# the model's prior (dose) and fitted admixture fraction for the pair
comp = pd.read_csv(DATA / "composition.csv.gz", dtype={"cell_id": str}).set_index("cell_id")
comp = comp.loc[ct.cell_id.to_numpy()[t_cols]]
prior, fitted = by_stratum(comp.dose.to_numpy()), by_stratum(comp.contamination.to_numpy())
overall = float((comp.contamination.to_numpy() * totals).sum() / totals.sum())

# malignant marker content of each stratum, read as a fraction of the whole malignant profile:
# the stratum's share of molecules on malignant-specific genes over those genes' share of the
# profile of the malignant cells bordering fibroblasts (what the prior is built from)
type_names = sorted(set(annotation))
prof = np.column_stack([np.asarray(before[:, np.flatnonzero(ann == t)].sum(axis=1)).ravel() for t in type_names])
prof = prof / prof.sum(axis=0) * 1e4
s_i = type_names.index(SOURCE)
spec = prof[:, s_i] / np.maximum(np.delete(prof, s_i, axis=1).max(axis=1), 0.1)
markers = np.flatnonzero((spec >= 5.0) & (prof[:, s_i] >= 0.5))
markers = markers[np.argsort(-prof[markers, s_i])][:20]
dist, dtypes = source_nearest_distance(cells, annotation)
s_cols = np.flatnonzero(ann == SOURCE)
near_t = np.asarray(dist)[:, [str(t) for t in dtypes].index(TARGET)]
psi = np.asarray(_interface_profile(before, s_cols, near_t[s_cols]))
M = np.asarray(before[markers][:, t_cols].sum(axis=0)).ravel()
observed = np.array([M[rung == r].sum() / T[i] for i, r in enumerate(shown)]) / float(psi[markers].sum())


def rates(matrix, gene):
    v = np.asarray(matrix[genes.get_loc(gene)][:, t_cols].todense()).ravel()
    s = np.array([v[rung == r].sum() for r in shown])
    return pd.DataFrame(dict(rate=s / T * 1e4, lo=stats.gamma.ppf(0.025, np.maximum(s, 1e-9)) / T * 1e4,
                             hi=stats.gamma.ppf(0.975, s + 1) / T * 1e4))


plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5, "axes.spines.top": False,
                     "axes.spines.right": False})
fig = plt.figure(figsize=(7.4, 5.0))

ax = fig.add_axes([0.08, 0.60, 0.60, 0.33])
ax.plot(xb, observed * 100, "s-", color=DARK, ms=4, lw=1, label="observed malignant marker content\n(scaled to the whole malignant profile)")
ax.plot(xb, prior * 100, "o--", color=RED, ms=4, lw=1, mfc="white", label="prior: observed minus the\ndeepest stratum")
ax.bar(xb, fitted * 100, width=0.6, color=RED, alpha=0.35, label="fitted admixture fraction")
ax.axhline(overall * 100, color=RED, lw=0.8, ls=":")
ax.text(2.5, overall * 100, f"all fibroblasts: {overall * 100:.1f}%", color=RED, fontsize=6.5, ha="center", va="bottom")
for i, n in enumerate(n_cells):
    ax.text(i, -0.2, f"{n:,}", fontsize=5.5, ha="center", va="top", color="#555555", transform=ax.get_xaxis_transform())
ax.set_xticks(xb, ticks)
ax.set_ylim(0, observed.max() * 100 * 1.15)
ax.set_xlabel("malignant cells among the nearest neighbors (grey: fibroblasts in each stratum)", labelpad=14)
ax.set_ylabel("share of fibroblast molecules (%)")
ax.set_title("admixture fraction of fibroblast molecules by exposure stratum")
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))

for k, (gene, kind) in enumerate(GENES):
    ax = fig.add_axes([0.08 + 0.335 * k, 0.13, 0.23, 0.26])
    ob, re = rates(before, gene), rates(after, gene)
    ax.fill_between(xb, re.rate, ob.rate, color=RED, alpha=0.12)
    ax.axhline(ob.rate[0], color=RED, lw=0.8, ls=":")   # the observed rate at the deepest stratum
    ax.errorbar(xb, ob.rate, yerr=[ob.rate - ob.lo, ob.hi - ob.rate], fmt="s-", color=RED, ms=3.5, lw=1, capsize=2, label="observed")
    ax.errorbar(xb, re.rate, yerr=[re.rate - re.lo, re.hi - re.rate], fmt="o-", color=BLUE, ms=3.5, lw=1, capsize=2, label="retained")
    ax.set_xticks(xb, short_ticks, fontsize=6, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_ylim(bottom=0)
    ax.set_title(f"{gene} in fibroblasts\n{kind}")
    if k == 1:
        ax.set_xlabel("malignant cells among the nearest neighbors (0 in K: none among the K nearest)", fontsize=7, labelpad=4)
    if k == 0:
        ax.set_ylabel("molecules per 10,000")
        ax.legend(frameon=False, loc="upper left")
for letter, (x, y) in zip("abcd", [(0.01, 0.96), (0.01, 0.45), (0.345, 0.45), (0.68, 0.45)]):
    fig.text(x, y, letter, fontsize=11, fontweight="bold")
fig.savefig(OUT, dpi=200)

print(pd.DataFrame(dict(stratum=shown, cells=n_cells, observed=observed, prior=prior, fitted=fitted)).round(4).to_string(index=False))
print(f"markers {len(markers)}; fitted over all fibroblasts {overall:.4f}")
for gene, _ in GENES:
    print(gene, "observed", rates(before, gene).rate.round(1).tolist(), "retained", rates(after, gene).rate.round(1).tolist())
