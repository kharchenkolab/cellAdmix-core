"""Figures for the real-data comparison of the dose prior's modes, from the
tables written by 11_dose_modes.py (per-pair attribution and, where the
audit scoring ran, the audit's per-pair measurements).

For each dataset: (a) contamination molecules attributed to each pair by
the marginal and the joint fit; (b) the ladder's prior contamination in
unexposed cells against the audit's own measurement of the unexposed
bin's excess above the ambient reference; (c) per-pair sensitivity on the
exposed cells under the audit's reference, by estimation mode. Without
audit scoring (breast 5K) panel (b) shows the attribution with and without
the ladder instead, and (c) is omitted.

Output: figures/dose_modes_<dataset>.png for pancreas, nsclc and breast
(those whose tables are present).
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GM = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
R = f"{GM}/results"

SHORT = {
    "pancreas": {"Ductal/tumor epithelial": "ductal", "Exocrine epithelial": "exocrine",
                 "Fibroblast / CAF": "fibroblast", "Mural / pericyte": "mural",
                 "Endothelial": "endothelial", "Endocrine": "endocrine", "Immune": "immune"},
    "nsclc": {},
    "breast": {"Fibroblast / CAF": "fibroblast", "B / plasma": "B/plasma", "T / NK": "T/NK",
               "Epithelial": "epithelial", "Endothelial": "endothelial", "Myeloid": "myeloid"},
}


def wide(df, col):
    return df.pivot_table(index="pair", columns="arm", values=col)


def scatter_pairs(ax, x, y, lab, color, xlabel, ylabel, title, annotate="low", n=4):
    x = x.clip(lower=1)
    y = y.clip(lower=1)
    ax.scatter(x, y, s=18, color=color, zorder=3)
    lim = [min(x.min(), y.min()) * 0.7, max(x.max(), y.max()) * 1.4]
    ax.plot(lim, lim, color="grey", lw=0.7)
    ratio = (y / x).sort_values()
    picks = ratio.index[:n] if annotate == "low" else ratio.index[-n:]
    offsets = [(4, -8), (4, 4), (4, -14), (4, 10)]
    for k, p in enumerate(picks):
        ax.annotate(lab[p], (x[p], y[p]), fontsize=7, xytext=offsets[k % len(offsets)],
                    textcoords="offset points")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)


for dataset in ["pancreas", "nsclc", "breast"]:
    prefix = "" if dataset == "pancreas" else f"{dataset}_"
    pairs_path = f"{R}/{prefix}gm_dose_arms_pairs.csv"
    if not os.path.exists(pairs_path):
        continue
    pairs = pd.read_csv(pairs_path)
    audit_path = f"{R}/{prefix}gm_dose_arms_audit_pairs.csv"
    audit = pd.read_csv(audit_path) if os.path.exists(audit_path) else None
    short = SHORT[dataset]
    sh = lambda s: short.get(s, s)  # noqa: E731
    pairs["label"] = pairs.source.map(sh) + " → " + pairs.target.map(sh)
    lab = pairs.drop_duplicates("pair").set_index("pair")["label"]
    w = wide(pairs, "posterior")

    n_panels = 3 if audit is not None else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(4.9 * n_panels, 4.4))

    scatter_pairs(axes[0], w["marginal"], w["joint"], lab, "#2980b9",
                  "marginal fit: contamination molecules attributed to the pair",
                  "joint fit", "(a) attribution per pair, marginal vs joint")

    if audit is not None:
        pu = wide(pairs, "prior_unexposed")["joint_ladder"]
        ub = audit[audit.arm == "marginal"].assign(
            pair=lambda d: d.source + " -> " + d.target).set_index("pair")["unexposed_excess_before"]
        common = pu.index.intersection(ub.index)
        ax = axes[1]
        xx = ub[common].clip(lower=1)
        yy = pu[common].clip(lower=1)
        ax.scatter(xx, yy, s=18, color="#1e8449", zorder=3)
        lim = [min(xx.min(), yy.min()) * 0.7, max(xx.max(), yy.max()) * 1.4]
        ax.plot(lim, lim, color="grey", lw=0.7)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("audit: unexposed-bin excess above the ambient reference (molecules)")
        ax.set_ylabel("ladder prior in unexposed cells (molecules)")
        ax.set_title("(b) the ladder prior against the audit's measurement", fontsize=10)

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
        ax.set_yticklabels([f"{sh(p.split(' -> ')[0])} → {sh(p.split(' -> ')[1])}"
                            for p in order], fontsize=6)
        ax.set_xlabel("exposure-linked excess removed (audit sensitivity)")
        ax.set_title("(c) per-pair sensitivity by estimation mode", fontsize=10)
        ax.legend(frameon=False, fontsize=7.5, loc="lower left")
        ax.axvline(0.8, color="grey", lw=0.6, ls=":")
    else:
        scatter_pairs(axes[1], w["marginal"], w["ladder"], lab, "#e67e22",
                      "marginal fit, single unexposed stratum (molecules)",
                      "marginal fit with the ladder",
                      "(b) attribution per pair, with and without the ladder",
                      annotate="high")

    fig.tight_layout()
    fig.savefig(f"{GM}/figures/dose_modes_{dataset}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"figure written: dose_modes_{dataset}.png")
