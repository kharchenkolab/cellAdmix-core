"""Figure for docs/generative.md: (a) retention of planted induced
expression as a function of its disproportionality to the transfer
expectation (spike-in sweep, pancreas), with the real flagged genes'
retained fractions overlaid at their measured disproportionality;
(b) a worked example - CFTR in ductal cells by exocrine exposure, with
the observed rate split into the model's retained and removed shares."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

GM = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
OUT = "/home/pkharchenko/cellAdmix/cellAdmix-core/docs/figures"

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.2))

# ---- (a) retention versus disproportionality --------------------------------
sw = pd.read_csv(os.path.join(GM, "results", "gm_spikein_sweep.csv"))
# The planted molecules add to the channel's excess: the disproportionality
# the screen faces is (planted + pre-existing) / proportional expectation.
sw["ratio"] = (sw.planted + sw.preexisting_excess) / sw.preexisting_excess
ax.scatter(sw.ratio, sw.retention_molecule, s=28, color="#2980b9",
           zorder=3, label="planted induction (spike-in sweep)")
led = pd.read_csv(os.path.join(GM, "results",
                               "gm_production_induced_retention.csv"))
scr = pd.read_csv(os.path.join(GM, "results", "gm_production_induced.csv"))
led = led.merge(scr[["pair", "gene", "excess", "proportional_expected"]],
                on=["pair", "gene"])
led["ratio"] = led.excess / led.proportional_expected
ax.scatter(led.ratio, led.retained_frac.clip(upper=1), s=22, marker="s",
           facecolors="none", edgecolors="#c0392b", zorder=3,
           label="real flagged genes (pancreas)")
ax.axhline(0.8, color="grey", lw=0.7, ls=":")
ax.set_xscale("log")
ax.set_xlabel("excess relative to the transfer expectation (fold)")
ax.set_ylabel("fraction of induced excess retained")
ax.set_title("(a) retention versus disproportionality", fontsize=10)
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.annotate("indistinguishable\nfrom transfer", (1.1, 0.32), fontsize=8,
            color="grey", ha="left")

# ---- (b) CFTR in ductal cells -----------------------------------------------
DATA = os.path.join(GM, "data")
genes = [l.strip() for l in open(os.path.join(DATA, "genes.txt"))]
cells = [l.strip() for l in open(os.path.join(DATA, "cells.txt"))]
col_of = {c: i for i, c in enumerate(cells)}
counts = sp.csc_matrix(scipy.io.mmread(os.path.join(DATA, "counts.mtx")))
removed = sp.csc_matrix(scipy.io.mmread(
    os.path.join(DATA, "gm_removed_production.mtx")))
totals = np.asarray(counts.sum(axis=0)).ravel()

expo = pd.read_csv(os.path.join(DATA, "exposure.csv.gz"))
PAIR = "Exocrine epithelial -> Ductal/tumor epithelial"
e = expo[expo.pair == PAIR]
cols = np.array([col_of[c] for c in e.cell_id])
ev = e.exposure.to_numpy()
g = genes.index("CFTR")
obs = np.asarray(counts[g, cols].todense()).ravel()
rem = np.asarray(removed[g, cols].todense()).ravel()
tot = totals[cols]

labels = ["0", "1", "2", "3+"]
bins = np.clip(ev, 0, 3)
x = np.arange(4)
r_obs = [obs[bins == b].sum() / tot[bins == b].sum() * 1e3 for b in range(4)]
r_keep = [(obs - rem)[bins == b].sum() / tot[bins == b].sum() * 1e3
          for b in range(4)]
ax2.plot(x, r_obs, "s-", color="#c0392b", label="observed")
ax2.plot(x, r_keep, "o-", color="#2980b9",
         label="retained (native + induced)")
ax2.fill_between(x, r_keep, r_obs, color="#c0392b", alpha=0.12,
                 label="removed (transferred)")
ax2.axhline(r_obs[0], color="grey", lw=0.7, ls=":")
ax2.set_xticks(x, labels)
ax2.set_xlabel("exocrine cells among 15 nearest neighbors")
ax2.set_ylabel("CFTR rate (per 1,000 molecules)")
ax2.set_title("(b) CFTR in ductal cells: the split", fontsize=10)
ax2.legend(frameon=False, fontsize=8, loc="upper left")

fig.tight_layout()
fig.savefig(os.path.join(OUT, "generative_fig1.png"), dpi=150,
            bbox_inches="tight")
print("FIGURE DONE")
