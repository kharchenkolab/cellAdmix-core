"""Figures for the real-data comparison of the dose prior's modes
(pancreas), from the tables written by 11_dose_modes.py and, when present,
the gate-harness per-pair results of 02_eval_gates.R.

Output: figures/dose_modes_pancreas.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GM = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
R = f"{GM}/results"
pairs = pd.read_csv(f"{R}/gm_dose_arms_pairs.csv")
audit = pd.read_csv(f"{R}/gm_dose_arms_audit_pairs.csv")
short = {"Ductal/tumor epithelial": "ductal", "Exocrine epithelial": "exocrine",
         "Fibroblast / CAF": "fibroblast", "Mural / pericyte": "mural",
         "Endothelial": "endothelial", "Endocrine": "endocrine", "Immune": "immune"}
pairs["label"] = pairs.source.map(short) + " → " + pairs.target.map(short)

def wide(df, col):
    return df.pivot_table(index="pair", columns="arm", values=col)

fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.4))

# (a) attributed molecules per pair: marginal versus joint
ax = axes[0]
w = wide(pairs, "posterior")
lab = pairs.drop_duplicates("pair").set_index("pair")["label"]
x = w["marginal"].clip(lower=1)
y = w["joint"].clip(lower=1)
ax.scatter(x, y, s=18, color="#2980b9", zorder=3)
lim = [max(x.min(), y.min()) * 0.7, max(x.max(), y.max()) * 1.4]
ax.plot(lim, lim, color="grey", lw=0.7)
ratio = (y / x)
for p in ratio.sort_values().index[:4]:
    ax.annotate(lab[p], (x[p], y[p]), fontsize=7, xytext=(4, -8),
                textcoords="offset points")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("marginal fit: contamination molecules attributed to the pair")
ax.set_ylabel("joint fit")
ax.set_title("(a) attribution per pair, marginal vs joint", fontsize=10)

# (b) unexposed-cell prior (ladder) vs the audit's unexposed-bin excess
ax = axes[1]
pu = wide(pairs, "prior_unexposed")["joint_ladder"]
ub = audit[audit.arm == "marginal"].assign(
    pair=lambda d: d.source + " -> " + d.target).set_index("pair")["unexposed_excess_before"]
common = pu.index.intersection(ub.index)
xx = ub[common].clip(lower=1); yy = pu[common].clip(lower=1)
ax.scatter(xx, yy, s=18, color="#1e8449", zorder=3)
lim = [min(xx.min(), yy.min()) * 0.7, max(xx.max(), yy.max()) * 1.4]
ax.plot(lim, lim, color="grey", lw=0.7)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("audit: unexposed-bin excess above the ambient reference (molecules)")
ax.set_ylabel("ladder prior in unexposed cells (molecules)")
ax.set_title("(b) the ladder prior against the audit's measurement", fontsize=10)

# (c) per-pair sensitivity on exposed cells, by arm
ax = axes[2]
a = audit.assign(pair=lambda d: d.source + " -> " + d.target,
                 sens=lambda d: 1 - d.exposed_excess_after / d.exposed_excess_before.clip(lower=1))
arms = ["marginal", "joint", "ladder", "joint_ladder"]
colors = {"marginal": "#7f8c8d", "joint": "#2980b9", "ladder": "#e67e22",
          "joint_ladder": "#1e8449"}
order = a[a.arm == "marginal"].sort_values("sens")["pair"].tolist()
for arm in arms:
    d = a[a.arm == arm].set_index("pair").reindex(order)
    ax.plot(d["sens"].to_numpy(), np.arange(len(order)), "o", ms=3.5,
            color=colors[arm], label=arm, alpha=0.9)
ax.set_yticks(np.arange(len(order)))
ax.set_yticklabels([f"{short[p.split(' -> ')[0]]} → {short[p.split(' -> ')[1]]}"
                    for p in order], fontsize=6)
ax.set_xlabel("exposure-linked excess removed (audit sensitivity)")
ax.set_title("(c) per-pair sensitivity by estimation mode", fontsize=10)
ax.legend(frameon=False, fontsize=7.5, loc="lower left")
ax.axvline(0.8, color="grey", lw=0.6, ls=":")
fig.tight_layout()
fig.savefig(f"{GM}/figures/dose_modes_pancreas.png", dpi=150, bbox_inches="tight")
print("figure written")
