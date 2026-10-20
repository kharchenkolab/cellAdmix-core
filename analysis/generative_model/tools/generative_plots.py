"""Reusable analysis and plotting helpers for the generative admixture model.

The Python counterparts of generative_plots.R: the figure designs of
docs/generative.md factored out of the example notebooks, kept outside the
package until the model's place in it settles. Import from a notebook with
the repository root on the path, for example::

    import sys; sys.path.insert(0, "../../analysis/generative_model/tools")
    from generative_plots import *

Conventions: ``before``/``after`` are gene x cell sparse matrices with the
gene and cell lists from ``ds.counts()``/``correction.counts()``;
``annotation`` is a pandas Series of cell types indexed by cell id;
``cells_xy`` is ``ds.cells()``; ``exposure`` is a Series of source-type
neighbor counts indexed by cell id (``generative_exposure``).
"""
import numpy as np
import pandas as pd
from scipy import stats

PALETTE = {"own expression": "#2980b9", "contamination": "#c0392b",
           "ambient": "#8a8a8a", "induced": "#1e8449"}
BIN_LABELS = ["0", "1", "2", "3+"]


def generative_exposure(cells_xy, annotation, source, neighbor_k=15):
    """Number of source-type cells among each cell's nearest neighbors."""
    from celladmix._score_utils import source_exposure_counts

    counts, types = source_exposure_counts(cells_xy, annotation, neighbor_k=neighbor_k)
    return pd.Series(counts[:, list(types).index(source)],
                     index=cells_xy["cell_id"].astype(str))


def generative_bins(exposure):
    return np.minimum(np.asarray(exposure, dtype=float), 3).astype(int)


def generative_bin_rates(matrix, gene_index, cols, bins, totals):
    """Per-bin rate of one gene (per 1,000 molecules) with exact Poisson intervals."""
    v = np.asarray(matrix[gene_index][:, cols].todense()).ravel()
    rows = []
    for b in range(4):
        m = v[bins == b].sum()
        t = max(totals[bins == b].sum(), 1.0)
        rows.append(dict(bin=BIN_LABELS[b], rate=m / t * 1e3,
                         lo=stats.gamma.ppf(0.025, max(m, 1e-9)) / t * 1e3,
                         hi=stats.gamma.ppf(0.975, m + 1) / t * 1e3))
    return pd.DataFrame(rows)


def _short(name):
    return str(name).split("/")[0].strip().lower()


# ---- admixture pattern ----------------------------------------------------------

def generative_plot_burden(before, after, cells, annotation, ax=None):
    """Share of each type's molecules removed by a correction."""
    import matplotlib.pyplot as plt

    types = annotation.reindex([str(c) for c in cells]).to_numpy()
    tb = np.asarray(before.sum(axis=0)).ravel()
    ta = np.asarray(after.sum(axis=0)).ravel()
    d = pd.DataFrame(dict(type=types, removed=tb - ta, total=tb)).dropna()
    g = d.groupby("type").sum()
    share = (g.removed / g.total).sort_values()
    if ax is None:
        _, ax = plt.subplots(figsize=(6.5, 3.2))
    ax.barh(share.index, share.to_numpy(), color=PALETTE["contamination"], alpha=0.85)
    ax.set_xlabel("share of the type's molecules removed as admixture")
    ax.xaxis.set_major_formatter(lambda v, _: f"{100 * v:.0f}%")
    return ax


def generative_plot_composition(model, target, exposure, x_label, ax=None):
    """Mean composition of the target type's content by exposure bin
    (docs/generative.md Figure 1b)."""
    import matplotlib.pyplot as plt

    comp = model.composition(target=target)
    per = comp.groupby("cell_id").agg(contamination=("contamination", "sum"),
                                      induced=("induced", "sum"),
                                      ambient=("ambient", "max"))
    per["own expression"] = (1 - per.contamination - per.induced - per.ambient).clip(lower=0)
    per["bin"] = generative_bins(exposure.reindex(per.index).to_numpy())
    agg = per.groupby("bin")[["own expression", "contamination", "ambient", "induced"]].mean()
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 3.8))
    bottom = np.zeros(len(agg))
    for part in ["own expression", "contamination", "ambient", "induced"]:
        ax.bar(agg.index, agg[part], bottom=bottom, color=PALETTE[part], width=0.72, label=part)
        bottom += agg[part].to_numpy()
    ax.set_xticks(range(4), BIN_LABELS)
    ax.set_xlabel(x_label)
    ax.set_ylabel(f"mean share of {_short(target)} content")
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    return ax


def generative_plot_composition_map(model, source, target, cells_xy, cap=0.6, ax=None):
    """Fitted contamination fraction of one pair drawn over the tissue."""
    import matplotlib.pyplot as plt

    comp = model.composition(source, target)
    xy = cells_xy.set_index(cells_xy["cell_id"].astype(str)).reindex(comp["cell_id"])
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 4))
    sc = ax.scatter(xy["x"], xy["y"], c=np.minimum(comp["contamination"], cap), s=0.5,
                    cmap="inferno", rasterized=True)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"{str(source).split('/')[0].strip()} contamination of {_short(target)} cells",
                 fontsize=10, loc="left")
    ax.figure.colorbar(sc, ax=ax, fraction=0.04,
                       label=f"{_short(source)} share of\n{_short(target)} content")
    return ax


def generative_plot_dose_vs_fitted(model, source, target, exposure, legend, ax=None):
    """Fitted contamination fraction against the exposure-derived prior."""
    import matplotlib.pyplot as plt

    comp = model.composition(source, target)
    e = np.minimum(exposure.reindex(comp["cell_id"]).to_numpy(), 3)
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    sc = ax.scatter(comp["dose"], comp["contamination"], c=e, s=3, alpha=0.5,
                    cmap="viridis", vmin=0, vmax=3.3)
    lim = [0, max(comp["dose"].max(), comp["contamination"].max()) * 1.05]
    ax.plot(lim, lim, color="grey", lw=0.6)
    ax.set_xlabel("dose: prior contamination fraction from exposure")
    ax.set_ylabel("fitted contamination fraction")
    ax.figure.colorbar(sc, ax=ax, fraction=0.04, label=legend, ticks=[0, 1, 2, 3])
    return ax


# ---- induced expression ----------------------------------------------------------

def generative_induced_changes(before, after, genes, cells, model, cells_xy, annotation,
                               target, neighbor_k=15, exposure_cap=4, min_counts=100):
    """Neighborhood-associated expression changes on the corrected counts:
    every gene of the target type fitted with a Poisson regression on all
    sources' exposures jointly (quasi errors), with the two guards of
    docs/generative.md Figure 4. One row per gene and source."""
    from celladmix.audit import _interface_profile
    from celladmix._score_utils import source_nearest_distance

    pairs = model.pairs[model.pairs["target"] == target]
    sources = list(pairs["source"])
    tcells = model.composition(sources[0], target)["cell_id"].astype(str)
    cindex = pd.Index([str(c) for c in cells])
    cols = cindex.get_indexer(tcells)
    E = np.column_stack([np.minimum(generative_exposure(cells_xy, annotation, s, neighbor_k)
                                    .reindex(tcells).to_numpy(), exposure_cap) for s in sources])
    before_t = before[:, cols].tocsc()
    after_t = after[:, cols].tocsc()
    removed_frac = 1 - np.asarray(after_t.sum(axis=1)).ravel() / np.maximum(
        np.asarray(before_t.sum(axis=1)).ravel(), 1.0)
    t_kept = np.asarray(after_t.sum(axis=0)).ravel()
    ok = t_kept > 10
    K = after_t[:, ok].tocsr()
    E = E[ok]
    t_kept = t_kept[ok]
    X = np.column_stack([np.ones(len(t_kept)), E])

    def fit_gene(y):
        beta = np.zeros(X.shape[1]); beta[0] = np.log(max(y.sum(), 0.5) / t_kept.sum())
        for _ in range(25):
            mu = t_kept * np.exp(np.clip(X @ beta, -30, 5))
            try:
                step = np.linalg.solve(X.T @ (X * mu[:, None]) + 1e-8 * np.eye(X.shape[1]),
                                       X.T @ (y - mu))
            except np.linalg.LinAlgError:
                return None
            beta += step
            if np.abs(step).max() < 1e-6:
                break
        mu = t_kept * np.exp(np.clip(X @ beta, -30, 5))
        disp = max(((y - mu) ** 2 / np.maximum(mu, 1e-9)).sum() / max(len(y) - X.shape[1], 1), 1.0)
        cov = np.linalg.inv(X.T @ (X * mu[:, None]) + 1e-8 * np.eye(X.shape[1]))
        se = np.sqrt(np.maximum(np.diag(cov), 1e-300) * disp)
        return beta[1:], beta[1:] / se[1:]

    # gene ownership and interface-local profiles with target-owned genes zeroed
    gindex = pd.Index([str(g) for g in genes])
    before_csc = before.tocsc()
    types = sorted(set(annotation.dropna()))
    prof = np.zeros((len(genes), len(types)))
    ann_cells = annotation.reindex(cindex).to_numpy()
    for j, t in enumerate(types):
        tc = np.flatnonzero(ann_cells == t)
        prof[:, j] = np.asarray(before_csc[:, tc].sum(axis=1)).ravel()
    prof = prof / np.maximum(prof.sum(axis=0), 1.0) * 1e6
    owner = np.array(types)[prof.argmax(axis=1)]
    near_dist, near_types = source_nearest_distance(cells_xy, annotation)
    near_df = pd.DataFrame(near_dist, columns=near_types,
                           index=cells_xy["cell_id"].astype(str))
    psi = []
    for s in sources:
        s_cells = [str(c) for c in annotation.index[annotation == s]]
        s_cells = [c for c in s_cells if c in cindex]
        s_cols = cindex.get_indexer(s_cells)
        p = np.asarray(_interface_profile(before_csc, s_cols,
                                          near_df.reindex(s_cells)[target].to_numpy()),
                       dtype=float).copy()
        p[owner == target] = 0.0
        psi.append(p)
    rows = []
    for gi in range(K.shape[0]):
        y = np.asarray(K[gi].todense()).ravel()
        if y.sum() < min_counts:
            continue
        f = fit_gene(y)
        if f is None:
            continue
        beta, z = f
        for j, s in enumerate(sources):
            rows.append(dict(gene=gindex[gi], source=s, log_fold_per_neighbor=float(beta[j]),
                             z=float(z[j]), removed_frac=float(removed_frac[gi]),
                             psi_share_e4=float(psi[j][gi] * 1e4)))
    out = pd.DataFrame(rows)
    out["clean"] = (out.removed_frac < 0.10) & (out.psi_share_e4 < 20)
    return out


def generative_plot_volcano(changes, source, z_cap=60, label_n=7, ax=None):
    """Volcano of the changes for one source (docs/generative.md Figure 4)."""
    import matplotlib.pyplot as plt

    d = changes[(changes.source == source) & np.isfinite(changes.z)].copy()
    src = _short(source)
    d["lfc"] = d.log_fold_per_neighbor / np.log(2)
    d["zc"] = np.minimum(d.z.abs(), z_cap)
    if ax is None:
        _, ax = plt.subplots(figsize=(6.8, 5))
    ent = d[~d.clean]
    ax.scatter(ent.lfc, ent.zc, s=12, facecolors="none", edgecolors="#b6b6b6", lw=0.8,
               label="transfer-entangled")
    up = d[d.clean & (d.z > 0)]; dn = d[d.clean & (d.z < 0)]
    ax.scatter(up.lfc, up.zc, s=14, color="#1e8449", lw=0, label=f"increased near {src}")
    ax.scatter(dn.lfc, dn.zc, s=14, color="#7d3c98", lw=0, label=f"decreased near {src}")
    ax.axhline(6, color="grey", lw=0.7, ls=":")
    for r in d[d.clean].sort_values("zc", ascending=False).head(label_n).itertuples():
        ax.annotate(r.gene, (r.lfc, r.zc), fontsize=7.5, xytext=(4, 2), textcoords="offset points")
    ax.set_xlabel(f"fold change per {src} neighbor (log2)")
    ax.set_ylabel(f"significance |z| (capped at {z_cap})")
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    return ax


def generative_plot_gene_contrast(before, after, genes_all, cells_all, cells, exposure,
                                  genes, x_label, y_label, after_alt=None,
                                  alt_label="induced share also removed", axes=None):
    """Per-bin rates of chosen genes in the target cells, observed and
    retained, with exact intervals, the removed band and an optional third
    series (docs/generative.md Figure 1c/d)."""
    import matplotlib.pyplot as plt

    gindex = pd.Index([str(g) for g in genes_all])
    cindex = pd.Index([str(c) for c in cells_all])
    cols = cindex.get_indexer([str(c) for c in cells])
    totals = np.asarray(before[:, cols].sum(axis=0)).ravel()
    bins = generative_bins(exposure.reindex([str(c) for c in cells]).to_numpy())
    series = [("observed", before, "#c0392b", "s", "-"), ("retained", after, "#2980b9", "o", "-")]
    if after_alt is not None:
        series.append((alt_label, after_alt, "#7d3c98", "o", "--"))
    if axes is None:
        _, axes = plt.subplots(1, len(genes), figsize=(4.3 * len(genes), 3.8))
    for ax, gene in zip(np.atleast_1d(axes), genes):
        gi = gindex.get_loc(str(gene))
        tabs = {lab: generative_bin_rates(m, gi, cols, bins, totals) for lab, m, *_ in series}
        x = np.arange(4)
        ax.fill_between(x, tabs["retained"].rate, tabs["observed"].rate, color="#c0392b", alpha=0.12)
        ax.axhline(tabs["observed"].rate[0], color="grey", lw=0.6, ls=":")
        for lab, _, color, marker, ls in series:
            t = tabs[lab]
            ax.errorbar(x, t.rate, yerr=[t.rate - t.lo, t.hi - t.rate], fmt=marker + ls,
                        color=color, ms=4, lw=1, capsize=2, label=lab)
        ax.set_title(gene, fontsize=10)
        ax.set_xticks(x, BIN_LABELS)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
    np.atleast_1d(axes)[0].legend(frameon=False, fontsize=8)
    return axes


def generative_plot_state_umaps(original_state, corrected_state, annotation,
                                titles=("Original counts", "Generative correction")):
    """Side-by-side cell-state UMAPs with one shared legend."""
    import matplotlib.pyplot as plt

    levels = sorted(set(annotation.dropna()))
    cmap = plt.get_cmap("tab10" if len(levels) <= 10 else "tab20")
    colors = {t: cmap(i) for i, t in enumerate(levels)}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, df, title in zip(axes, (original_state, corrected_state), titles):
        for t in levels:
            d = df[df.cell_type == t]
            ax.scatter(d.umap_1, d.umap_2, s=2, alpha=0.7, color=colors[t], label=t, lw=0)
        ax.set_title(title, fontsize=10, loc="left")
        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    axes[1].legend(frameon=False, fontsize=8, markerscale=4, bbox_to_anchor=(1.02, 1),
                   loc="upper left")
    fig.tight_layout()
    return axes
