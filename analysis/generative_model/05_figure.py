"""Figures for docs/generative.md.

Figure 1 - the model's components and a worked example (pancreas,
exocrine -> ductal): (a) schematic of the components of a target cell's
observed counts; (b) the fitted composition of ductal-cell content as a
function of exocrine exposure; (c, d) two genes contrasted - AMY2A, whose
exposure gradient is transferred material and is removed, and CFTR, whose
gradient is an induced duct-cell program and is retained.

Figure 2 - retention of planted induced expression as a function of its
disproportionality to the transfer expectation (spike-in sweep), with the
real flagged genes overlaid.
"""
import importlib.util
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
import numpy as np
import pandas as pd
import scipy.sparse as sp

GM = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
OUT = "/home/pkharchenko/cellAdmix/cellAdmix-core/docs/figures"
TARGET = "Ductal/tumor epithelial"
PAIR = "Exocrine epithelial -> Ductal/tumor epithelial"

spec = importlib.util.spec_from_file_location("gm_model",
    os.path.join(GM, "01_model.py"))
gm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gm)

# ---- fit the production model and keep the target type's component rates ---
CACHE = os.path.join(gm.DATA, "gm_fig_rates_cache.npz")
inp = gm.Inputs()
if os.path.exists(CACHE):
    z = np.load(CACHE)
    rates = dict(own=z["own"], cont=z["cont"], amb=z["amb"], ind=z["ind"],
                 eps=z["eps"], cols=z["cols"])
else:
    rng = np.random.default_rng(1)
    own_frac = None
    rates = None
    for rnd in range(gm.N_OUTER):
        diag = {}
        own_frac_new = {}
        for T in inp.type_names:
            if T == TARGET and rnd == gm.N_OUTER - 1:
                res, of, rates = gm.fit_target(inp, T, own_frac, "production",
                                               rng, diag, return_rates=True)
            else:
                out = gm.fit_target(inp, T, own_frac, "production", rng, diag)
                res, of = out
            if of is not None:
                own_frac_new[T] = of
        own_frac = own_frac_new
    np.savez_compressed(CACHE, own=rates["own"], cont=rates["cont"],
                        amb=rates["amb"], ind=rates["ind"], eps=rates["eps"],
                        cols=rates["cols"])

cols = rates["cols"]
Yt = np.asarray(inp.counts[:, cols].todense(), dtype=float).T   # n x G
tot = Yt.sum(axis=1)
total_rate = rates["own"] + rates["cont"] + rates["amb"] + rates["ind"] \
    + rates["eps"]
keep = (rates["own"] + rates["ind"] + rates["eps"]) \
    / np.maximum(total_rate, 1e-300)

pinfo = inp.pair_info[PAIR]
pos_of_col = {c: i for i, c in enumerate(cols)}
pos = np.array([pos_of_col[c] for c in pinfo["cells"]])
expo = pinfo["exposure"]
bins = np.clip(expo, 0, 3)
labels = ["0", "1", "2", "3+"]

fig = plt.figure(figsize=(14.5, 4.0))
gs = fig.add_gridspec(1, 4, width_ratios=[1.25, 1, 0.9, 0.9], wspace=0.32)

# ---- (a) schematic ---------------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
ax.set_xlim(0, 1.05)
ax.set_ylim(0.03, 1.02)
ax.set_aspect("equal")
ax.axis("off")
rs = np.random.default_rng(4)

src = Circle((0.17, 0.50), 0.155, fc="#fdecea", ec="#c0392b", lw=1.4)
tgt = Circle((0.63, 0.48), 0.27, fc="#f4f9fd", ec="#2c3e50", lw=1.6)
ax.add_patch(src)
ax.add_patch(tgt)
ax.add_patch(Circle((0.70, 0.50), 0.085, fc="#d6e4f0", ec="#7f9db9", lw=1.0))

def dots(cx, cy, r, n, color, size=7, cluster=1.0):
    ang = rs.uniform(0, 2 * np.pi, n)
    rad = r * np.sqrt(rs.uniform(0, 1, n)) * cluster
    ax.scatter(cx + rad * np.cos(ang), cy + rad * np.sin(ang),
               s=size, color=color, zorder=5, lw=0)

dots(0.17, 0.50, 0.13, 60, "#c0392b")                 # source content
dots(0.63, 0.48, 0.24, 55, "#2980b9")                 # own expression
dots(0.415, 0.46, 0.045, 16, "#c0392b")               # transferred patch
dots(0.63, 0.48, 0.24, 8, "#8a8a8a", size=6)          # ambient inside
for cx, cy in [(0.36, 0.85), (0.95, 0.20), (0.06, 0.14), (0.93, 0.83)]:
    dots(cx, cy, 0.03, 3, "#8a8a8a", size=6)          # ambient outside
ang = np.linspace(0.6, 1.9, 9)                         # induced, near nucleus
ax.scatter(0.70 + 0.105 * np.cos(ang), 0.50 + 0.105 * np.sin(ang),
           s=9, color="#1e8449", zorder=6, lw=0)
ax.add_patch(FancyArrowPatch((0.31, 0.44), (0.385, 0.45), lw=1.2,
             arrowstyle="-|>", mutation_scale=11, color="#c0392b"))

ax.text(0.17, 0.71, "source cell", ha="center", fontsize=8.5)
ax.text(0.63, 0.80, "target cell", ha="center", fontsize=8.5)
ax.text(0.365, 0.30, "transferred\n$\\alpha_{S}\\,\\psi_{S}$", ha="center",
        fontsize=8, color="#c0392b")
ax.text(0.80, 0.13, "own programs\n$\\sum_k \\theta_{k} F_{k}$", ha="center",
        fontsize=8, color="#2980b9")
ax.text(0.885, 0.68, "induced\n$u\\,\\rho\\,m$", ha="center", fontsize=8,
        color="#1e8449")
ax.text(0.10, 0.90, "ambient $\\beta\\,a$", fontsize=8, color="#6a6a6a")
ax.set_title("(a) components of a cell's counts", fontsize=10)

# ---- (b) fitted composition by exposure ------------------------------------
ax = fig.add_subplot(gs[0, 1])
comp_names = [("own", "#2980b9", "own expression"),
              ("cont", "#c0392b", "contamination"),
              ("amb", "#8a8a8a", "ambient"),
              ("ind", "#1e8449", "induced")]
x = np.arange(4)
bottom = np.zeros(4)
for key, color, label in comp_names:
    vals = np.array([rates[key][pos[bins == b]].sum()
                     / max(tot[pos[bins == b]].sum(), 1.0) for b in range(4)])
    ax.bar(x, vals, bottom=bottom, color=color, width=0.7, label=label,
           alpha=0.9 if key != "own" else 0.75)
    bottom += vals
ax.set_ylim(0, 1.02)
ax.set_xticks(x, labels)
ax.set_xlabel("exocrine cells among 15 nearest")
ax.set_ylabel("share of ductal-cell content")
ax.legend(frameon=False, fontsize=7.5, loc="lower left")
ax.set_title("(b) fitted composition of ductal cells", fontsize=10)

# ---- (c, d) two genes: transferred vs induced ------------------------------
for panel, gene, subtitle in [(2, "AMY2A", "(c) AMY2A: transferred, removed"),
                              (3, "CFTR", "(d) CFTR: induced, retained")]:
    ax = fig.add_subplot(gs[0, panel])
    g = inp.gene_of[gene]
    obs = Yt[pos, g]
    ret = obs * keep[pos, g]
    tt = tot[pos]
    r_obs = [obs[bins == b].sum() / tt[bins == b].sum() * 1e3
             for b in range(4)]
    r_ret = [ret[bins == b].sum() / tt[bins == b].sum() * 1e3
             for b in range(4)]
    ax.plot(x, r_obs, "s-", color="#c0392b", label="observed")
    ax.plot(x, r_ret, "o-", color="#2980b9", label="retained")
    ax.fill_between(x, r_ret, r_obs, color="#c0392b", alpha=0.12)
    ax.axhline(r_obs[0], color="grey", lw=0.7, ls=":")
    ax.set_xticks(x, labels)
    ax.set_xlabel("exocrine cells among 15 nearest")
    ax.set_ylabel(f"{gene} rate (per 1,000 molecules)")
    if panel == 2:
        ax.legend(frameon=False, fontsize=8)
    ax.set_title(subtitle, fontsize=10)

fig.savefig(os.path.join(OUT, "generative_fig1.png"), dpi=150,
            bbox_inches="tight")
print("fig1 done")

# ---- Figure 2: the retention boundary --------------------------------------
fig2, ax = plt.subplots(figsize=(5.8, 4.4))
sw = pd.read_csv(os.path.join(GM, "results", "gm_spikein_sweep.csv"))
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
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.annotate("indistinguishable\nfrom transfer", (1.1, 0.32), fontsize=8,
            color="grey", ha="left")
fig2.tight_layout()
fig2.savefig(os.path.join(OUT, "generative_fig2.png"), dpi=150,
             bbox_inches="tight")
print("fig2 done")
