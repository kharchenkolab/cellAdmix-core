"""Admixture audit: exposure-based estimates of per-cell-type-pair
admixture rates and admixed-molecule counts, and verification of
corrections against them.

The measure exploits the spatial structure of segmentation-driven admixture:
source-marker content in target cells rises with the number of source-type
neighbor cells, while unexposed target cells form an internal negative
control. See docs/benchmarks.md for the methodology and validation.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

BIN_EDGES = [-0.5, 0.5, 1.5, 2.5, np.inf]
BIN_LABELS = ["0", "1", "2", "3+"]

# Minimum own-marker molecule count before the severe-over-removal warning
# applies (module-level so tests can lower it for small fixtures).
OWN_MARKER_WARN_MIN = 1000


def _pseudobulk(matrix, cell_index, cells, scale=1e6):
    cols = cell_index.get_indexer(cells)
    cols = cols[cols >= 0]
    if not len(cols):
        return np.full(matrix.shape[0], np.nan)
    sums = np.asarray(matrix[:, cols].sum(axis=1)).ravel()
    return sums / max(sums.sum(), 1) * scale


def _marker_pool(profiles, types, source, baseline, n_pool=20, min_source_cpm=50.0,
                 baseline_floor=20.0):
    top = np.array(types)[np.argmax(profiles, axis=1)]
    s_idx = types.index(source)
    ok = (top == source) & (profiles[:, s_idx] >= min_source_cpm)
    cand = np.flatnonzero(ok)
    contrast = profiles[cand, s_idx] / np.maximum(baseline[cand], baseline_floor)
    order = np.lexsort((-profiles[cand, s_idx], -contrast))
    return cand[order][:n_pool]


def _bin_rates(marker_counts, totals, bins):
    from scipy import stats as sps

    rows = []
    for b in BIN_LABELS:
        mask = bins == b
        m = float(marker_counts[mask].sum())
        M = max(float(totals[mask].sum()), 1.0)
        rows.append(dict(bin=b, markers=m, totals=M, n_cells=int(mask.sum()),
            rate=m / M,
            rate_lo=sps.gamma.ppf(0.025, max(m, 1e-9)) / M,
            rate_hi=sps.gamma.ppf(0.975, m + 1) / M))
    return pd.DataFrame(rows)


def _excess(rates):
    from scipy import stats as sps

    r0 = float(rates.loc[rates["bin"] == "0", "rate"].iloc[0])
    exposed = rates[rates["bin"] != "0"]
    excess = float((np.maximum(exposed["rate"] - r0, 0) * exposed["totals"]).sum())
    expected = r0 * float(exposed["totals"].sum())
    m = float(exposed["markers"].sum())
    p = 1.0 if (expected <= 0 or m <= expected) else float(
        sps.poisson.sf(m - 1, expected))
    return excess, p


def _power(rates_before, rates_after):
    r0b = float(rates_before.loc[rates_before["bin"] == "0", "rate"].iloc[0])
    r0a = float(rates_after.loc[rates_after["bin"] == "0", "rate"].iloc[0])
    eb = np.maximum(rates_before.loc[rates_before["bin"] != "0", "rate"] - r0b, 0) \
        * rates_before.loc[rates_before["bin"] != "0", "totals"]
    ea = np.maximum(rates_after.loc[rates_after["bin"] != "0", "rate"] - r0a, 0) \
        * rates_before.loc[rates_before["bin"] != "0", "totals"]
    return np.nan if eb.sum() <= 0 else 1 - float(ea.sum()) / float(eb.sum())


def _pava_dec(y, w):
    """Weighted monotone-decreasing fit (pool-adjacent violators):
    stabilizes each rung's estimate with the whole ladder, without
    extrapolating beyond it."""
    val = list(map(float, y))
    wt = list(map(float, w))
    idx = [[i] for i in range(len(val))]
    i = 0
    while i < len(val) - 1:
        if val[i] < val[i + 1] - 1e-15:
            mw = wt[i] + wt[i + 1]
            mv = (val[i] * wt[i] + val[i + 1] * wt[i + 1]) / mw
            val[i], wt[i] = mv, mw
            idx[i].extend(idx[i + 1])
            del val[i + 1], wt[i + 1], idx[i + 1]
            i = max(0, i - 1)
        else:
            i += 1
    out = np.empty(sum(len(g) for g in idx))
    for j, group in enumerate(idx):
        for g in group:
            out[g] = val[j]
    return out


def _reference_rate(marker_counts, totals, zero_by_K, min_tail_totals=2e4):
    """Ambient reference: the pooled panel rate among target cells with zero
    source-type cells among their K nearest neighbor cells, for the largest
    K that retains enough reference molecules.

    A cell with zero source neighbors among its 15 nearest can still sit in
    source-rich surroundings - including source cells above or below the
    section plane - and carry their material; requiring zero source cells
    among progressively more neighbors selects cells in progressively
    source-free surroundings while spanning only a ~100 um lateral radius,
    so the comparison never leaves the local tissue neighborhood.
    ``zero_by_K``: ordered mapping of label -> boolean mask over the target
    cells, ascending in K; the first entry is the base exposure definition.
    """
    labels = list(zero_by_K)
    m_k = np.array([float(marker_counts[z].sum()) for z in zero_by_K.values()])
    M_k = np.array([float(totals[z].sum()) for z in zero_by_K.values()])
    rate_k = m_k / np.maximum(M_k, 1.0)
    # Monotone fit across all rungs; the reference is the fitted value at
    # the deepest well-populated rung - every rung contributes, nothing
    # beyond the measured ladder is assumed.
    fitted = _pava_dec(rate_k, np.maximum(M_k, 1.0))
    pooled = float(rate_k[0])
    ok = np.flatnonzero(M_k >= min_tail_totals)
    pick = int(ok.max()) if len(ok) else 0
    return (min(float(fitted[pick]), pooled), labels[pick], pooled,
            fitted)


def _reference_trend(d):
    lad = d["ref_ladder"]
    K = int(d["reference_kind"][1:])
    idx = lad.index[lad["K"] == K]
    if not len(idx) or idx[0] == 0:
        return np.nan
    prev = lad["rate_fitted"].iloc[idx[0] - 1]
    return float((prev - lad["rate_fitted"].iloc[idx[0]]) / max(prev, 1e-12))


def _excess_vs_ref(rates, ref_rate):
    """Leakage above an external reference rate, summed over all exposure
    bins - including the unexposed bin, whose content above the ambient
    reference is contamination from unobserved (out-of-section) neighbors."""
    return float((np.maximum(rates["rate"] - ref_rate, 0) * rates["totals"]).sum())


def _interface_profile(matrix, s_cols, dist_T, radii=(30.0, 100.0, np.inf),
                       min_molecules=3e4):
    """Expression profile of the source cells that actually border the
    target type. Cell states are not uniform across a tissue: source cells
    at an interface can express activation genes several-fold above the
    source average, and material transferred from them carries that
    elevated share. Measuring transfer proportionality against the
    interface-local profile keeps such genes attributed to transfer;
    against the global profile they would look induced in the target. The
    radius widens when the near subset carries too few molecules."""
    s_cols = np.asarray(s_cols)
    d = np.nan_to_num(np.asarray(dist_T, dtype=float), nan=np.inf)
    for r in radii:
        sel = s_cols[d <= r]
        if not len(sel):
            continue
        cnt = np.asarray(matrix[:, sel].sum(axis=1), dtype=float).ravel()
        if cnt.sum() >= min_molecules or not np.isfinite(r):
            return cnt / max(cnt.sum(), 1.0)
    cnt = np.asarray(matrix[:, s_cols].sum(axis=1), dtype=float).ravel()
    return cnt / max(cnt.sum(), 1.0)


def _screen_panel(matrix, candidates, cols_exposed, cols_unexposed,
                  source_profile, n_pool=20, z_min=8.0,
                  excess_min=50.0, profile_cv=0.15):
    """Exclude likely induced genes from the marker panel.

    A gene whose exposure-linked excess in target cells is
    disproportionate to the interface-local source profile is more likely
    induced by proximity than transferred; it is excluded and the panel is
    filled from the remaining candidates in rank order. The proportional
    expectation is a weighted regression of per-gene excess on the
    profile, and the residual is standardized by counting noise plus a
    multiplicative profile-uncertainty term - without the latter, a small
    relative deviation on a large transfer channel reaches an arbitrary
    significance simply because the channel is large.
    """
    empty_stats = pd.DataFrame(columns=["gene", "excess", "expected",
                                        "fold", "z"])
    gset = np.asarray(list(candidates), dtype=int)
    t_exp = float(matrix[:, cols_exposed].sum())
    t_un = float(matrix[:, cols_unexposed].sum())
    if t_exp <= 0 or t_un <= 0 or not len(gset):
        return (np.asarray(candidates[:n_pool]), np.array([], dtype=int),
                empty_stats)
    m_exp = np.asarray(matrix[gset][:, cols_exposed].sum(axis=1)).ravel()
    m_un = np.asarray(matrix[gset][:, cols_unexposed].sum(axis=1)).ravel()
    excess = m_exp - m_un / t_un * t_exp
    var_count = m_exp + (t_exp / t_un) ** 2 * m_un + 1
    psi = np.maximum(source_profile[gset], 0.0)
    w = 1.0 / var_count
    slope = max(float((w * excess * psi).sum())
                / max(float((w * psi * psi).sum()), 1e-12), 0.0)
    expected = slope * psi
    z = (excess - expected) / np.sqrt(var_count + (profile_cv * expected) ** 2)
    sel = (z > z_min) & (excess > excess_min)
    flagged = gset[sel]
    stats = pd.DataFrame(dict(gene=flagged, excess=excess[sel],
        expected=expected[sel], fold=excess[sel] / np.maximum(expected[sel], 1e-9),
        z=z[sel]))
    pool = [g for g in gset.tolist() if g not in set(flagged.tolist())]
    return np.asarray(pool[:n_pool]), flagged, stats


class CellAdmixAudit:
    """Per-cell-type-pair admixture estimates and cleanup verification."""

    def __init__(self, fit, *, neighbor_k=15, n_pool=20, q_thresh=0.01,
                 min_excess=200, min_target_cells=200, min_reference_cells=100):
        from ._score_utils import source_exposure_counts, source_nearest_distance

        annotation = getattr(getattr(fit, "dataset", None), "annotation", None)
        if annotation is None:
            raise ValueError("audit_admixture requires a cell-type annotation")
        self.fit = fit
        self.params = dict(neighbor_k=neighbor_k, n_pool=n_pool,
            q_thresh=q_thresh, min_excess=min_excess)
        cells = fit.cell_factors()
        matrix, genes, cell_ids = fit.counts()
        self._matrix = matrix.tocsc()
        self._genes = pd.Index(genes)
        self._cells = pd.Index([str(c) for c in cell_ids])
        self._totals = np.asarray(matrix.sum(axis=0)).ravel()
        exp_counts, types = source_exposure_counts(cells, annotation,
            neighbor_k=neighbor_k)
        exp_df = pd.DataFrame(exp_counts, columns=types,
            index=cells["cell_id"].astype(str))
        # Exposure at growing neighborhood sizes, used only to pick each
        # pair's reference cells (zero source neighbors among the K nearest).
        ladder_K = sorted({int(neighbor_k), 30, 60, 120, 240})
        ladder_K = [K for K in ladder_K if K >= int(neighbor_k)]
        expo_ladder = {}
        for K in ladder_K:
            cts, _ = source_exposure_counts(cells, annotation, neighbor_k=K)
            expo_ladder[f"k{K}"] = pd.DataFrame(cts, columns=types,
                index=cells["cell_id"].astype(str))
        ctypes = pd.Series(
            [str(annotation.get(str(c), None)) for c in cells["cell_id"].astype(str)],
            index=cells["cell_id"].astype(str))
        self._types = list(types)
        # Distance from every cell to the nearest cell of each type, for
        # selecting the interface-local source cells of each pair's screen.
        near_dist, near_types = source_nearest_distance(cells, annotation)
        near_df = pd.DataFrame(near_dist, columns=near_types,
            index=cells["cell_id"].astype(str))
        profiles = np.column_stack([
            _pseudobulk(self._matrix, self._cells,
                ctypes.index[ctypes == t]) for t in types])
        self._profiles = profiles
        top = np.array(types)[np.argmax(profiles, axis=1)]
        self._native_markers = {}
        for ti, t in enumerate(types):
            cand = np.flatnonzero((top == t) & (profiles[:, ti] >= 200))
            others = [i for i in range(len(types)) if i != ti]
            max_other = profiles[cand][:, others].max(axis=1) if len(cand) else np.array([])
            spec = profiles[cand, ti] / np.maximum(profiles[cand, ti] + max_other, 1e-9)
            self._native_markers[t] = cand[np.argsort(-spec)][:10]
        self._ctypes = ctypes

        self._pairs = {}
        for source in types:
            s_idx = types.index(source)
            for target in types:
                if target == source:
                    continue
                T_cells = ctypes.index[ctypes == target]
                T_cells = T_cells[self._cells.get_indexer(T_cells) >= 0]
                if len(T_cells) < min_target_cells:
                    continue
                expo = exp_df.loc[T_cells, source].to_numpy()
                if (expo == 0).sum() < min_reference_cells:
                    continue
                baseline = _pseudobulk(self._matrix, self._cells, T_cells[expo == 0])
                bins = pd.cut(expo, BIN_EDGES, labels=BIN_LABELS).astype(str)
                cols = self._cells.get_indexer(T_cells)
                tot = self._totals[cols]
                candidates = _marker_pool(profiles, self._types, source,
                    np.nan_to_num(baseline), n_pool=10 * n_pool)
                S_cells = ctypes.index[ctypes == source]
                S_cells = S_cells[self._cells.get_indexer(S_cells) >= 0]
                psi_interface = _interface_profile(self._matrix,
                    self._cells.get_indexer(S_cells),
                    near_df.loc[S_cells, target].to_numpy())
                pool, induced, induced_stats = _screen_panel(self._matrix,
                    list(candidates), cols[expo > 0], cols[expo == 0],
                    psi_interface, n_pool=n_pool)
                if len(pool) < 3:
                    continue
                strict = pool[np.nan_to_num(baseline)[pool] <
                    0.05 * profiles[pool, s_idx]]
                mk_by_cell = np.asarray(
                    self._matrix[pool][:, cols].sum(axis=0)).ravel()
                rates = _bin_rates(mk_by_cell, tot, bins)
                excess_gradient, p = _excess(rates)
                # The reported leakage is measured against the ambient
                # reference (zero source neighbors among the largest
                # sufficiently populated neighborhood); detection stays
                # gradient-based.
                zero_by_K = {label: df.loc[T_cells, source].to_numpy() == 0
                             for label, df in expo_ladder.items()}
                ref_rate, ref_kind, ref_unexposed, ref_fitted = _reference_rate(
                    mk_by_cell, tot, zero_by_K)
                ref_mask = zero_by_K[ref_kind]
                ref_ladder = pd.DataFrame(dict(
                    K=[int(k[1:]) for k in zero_by_K],
                    rate=[float(mk_by_cell[z].sum()) / max(float(tot[z].sum()), 1.0)
                          for z in zero_by_K.values()],
                    rate_fitted=ref_fitted,
                    totals=[float(tot[z].sum()) for z in zero_by_K.values()]))
                excess = _excess_vs_ref(rates, ref_rate)
                if len(strict) >= 2:
                    mks = np.asarray(self._matrix[strict][:, cols].sum(axis=0)).ravel()
                    rates_strict = _bin_rates(mks, tot, bins)
                    ref_s, _, _, _ = _reference_rate(mks, tot, zero_by_K)
                    excess_strict = _excess_vs_ref(rates_strict, ref_s)
                else:
                    rates_strict, excess_strict = None, np.nan
                # marker-pool share of the source transcriptome: the
                # extrapolation factor from pool leakage to total admixture
                coverage = float(profiles[pool, s_idx].sum()) / max(
                    float(profiles[:, s_idx].sum()), 1e-9)
                self._pairs[(source, target)] = dict(
                    source=source, target=target, cols=cols, bins=bins,
                    ref_mask=ref_mask, ref_ladder=ref_ladder,
                    pool=pool, strict=strict, induced=induced,
                    induced_stats=induced_stats,
                    rates=rates, rates_strict=rates_strict,
                    excess=excess, excess_gradient=excess_gradient, p=p,
                    reference_rate=ref_rate, reference_kind=ref_kind,
                    reference_rate_unexposed=ref_unexposed,
                    excess_strict=excess_strict, coverage=coverage,
                    target_molecules=float(tot.sum()))
        pvals = np.array([d["p"] for d in self._pairs.values()])
        if len(pvals):
            order = np.argsort(pvals)
            q = np.empty_like(pvals)
            ranked = pvals[order] * len(pvals) / (np.arange(len(pvals)) + 1)
            q[order] = np.minimum.accumulate(ranked[::-1])[::-1]
            for d, qv in zip(self._pairs.values(), q):
                d["q"] = float(qv)
                d["detected"] = (qv < q_thresh
                    and d["excess_gradient"] >= min_excess)
        induced_all = sorted({int(g) for d in self._pairs.values()
                              for g in d["induced"]})
        if induced_all:
            names = [str(self._genes[g]) for g in induced_all[:8]]
            print(f"Excluded {len(induced_all)} likely induced gene"
                  f"{'s' if len(induced_all) > 1 else ''} from marker panels "
                  "(exposure-linked excess far above the source-profile "
                  f"expectation): {', '.join(names)}")

    # ---- tables ----
    def pairs(self, detected_only=False):
        rows = []
        for d in self._pairs.values():
            exposed = d["cols"][np.asarray(d["bins"]) != "0"]
            rows.append(dict(source=d["source"], target=d["target"],
                rate=d["excess"] / (d["coverage"] * max(d["target_molecules"], 1)),
                admixed_molecules=round(d["excess"] / d["coverage"]),
                excess=round(d["excess"]),
                excess_strict=None if np.isnan(d["excess_strict"]) else round(d["excess_strict"]),
                coverage=d["coverage"],
                q_value=d["q"], detected=d["detected"],
                reference_kind=d["reference_kind"],
                # ratio of the pooled unexposed rate to the ambient
                # reference: values above 1 quantify contamination from
                # out-of-section neighbors present in the unexposed cells
                reference_inflation=d["reference_rate_unexposed"]
                    / max(d["reference_rate"], 1e-12),
                # fractional decline over the ladder's final step: values
                # well above 0 mean the reference had not flattened at the
                # chosen neighborhood, so the estimate stays conservative
                reference_trend=_reference_trend(d),
                n_exposed=len(exposed), n_reference=int((np.asarray(d["bins"]) == "0").sum()),
                n_markers=len(d["pool"]), n_strict=len(d["strict"]),
                n_induced=len(d["induced"])))
        df = pd.DataFrame(rows)
        return df[df["detected"]] if detected_only else df

    def correct_generative(self, **kwargs):
        """Fit the generative admixture model on the detected pairs.

        Decomposes each target cell's counts into own expression,
        per-source contamination, ambient background, and induced
        expression, removing the contamination and ambient shares while
        retaining the induced ones; see docs/generative.md. Returns a
        correction object accepted by ``evaluate``.
        """
        from .generative import correct_generative

        return correct_generative(self, **kwargs)

    def markers(self, source, target):
        d = self._pair(source, target)
        stats = d["induced_stats"].copy()
        if len(stats):
            stats["gene"] = self._genes[stats["gene"].to_numpy(dtype=int)]
        return dict(pool=list(self._genes[d["pool"]]),
            strict=list(self._genes[d["strict"]]),
            induced=list(self._genes[d["induced"]]),
            induced_stats=stats)

    def _pair(self, source, target):
        d = self._pairs.get((source, target))
        if d is None:
            raise KeyError(f"No audited pair {source} -> {target}")
        return d

    # ---- plots ----
    def plot_map(self, value="rate", detected_only=True, ax=None):
        import matplotlib.pyplot as plt

        df = self.pairs(detected_only=detected_only)
        if df.empty:
            raise ValueError("No detected admixture pairs")
        sources = sorted(df["source"].unique())
        targets = sorted(df["target"].unique())
        grid = np.full((len(sources), len(targets)), np.nan)
        text = {}
        for _, r in df.iterrows():
            i, j = sources.index(r["source"]), targets.index(r["target"])
            if value == "rate":
                grid[i, j] = 100 * r["rate"]
                text[(i, j)] = f"{100 * r['rate']:.1f}"
            else:
                grid[i, j] = r["admixed_molecules"]
                text[(i, j)] = (f"{r['admixed_molecules']/1000:.0f}k"
                    if r["admixed_molecules"] >= 1000
                    else f"{int(r['admixed_molecules'])}")
        if ax is None:
            _, ax = plt.subplots(figsize=(6.5, 5))
        im = ax.imshow(grid, cmap="Oranges", aspect="auto")
        for (i, j), t in text.items():
            ax.text(j, i, t, ha="center", va="center", fontsize=8)
        ax.set_xticks(range(len(targets)), targets, rotation=40, ha="right")
        ax.set_yticks(range(len(sources)), sources)
        ax.set_xlabel("target cell type")
        ax.set_ylabel("source cell type")
        ax.set_title("Estimated admixture by cell-type pair\n"
            + ("estimated admixture rate: % of the target type's molecules"
               if value == "rate" else "cell text: estimated admixed molecules"),
            fontsize=10)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("grey")
        plt.colorbar(im, ax=ax, fraction=0.04,
            label="% of target-type molecules" if value == "rate"
            else "admixed molecules")
        return ax

    def _state_rates(self, matrix, source=None, target=None, strict=True):
        if source is None:
            # Pool each pair's excess over its own ambient reference.
            # Pooling raw rates would mix pairs with very different
            # baselines, and the shifting pair composition across bins can
            # then produce spurious non-monotone curves (Simpson's paradox).
            excess = np.zeros(len(BIN_LABELS))
            var_tot = np.zeros(len(BIN_LABELS))
            totals = np.zeros(len(BIN_LABELS))
            n_pairs = 0
            for d in self._pairs.values():
                if not d["detected"]:
                    continue
                n_pairs += 1
                genes = d["strict"] if strict and len(d["strict"]) >= 2 else d["pool"]
                mk = np.asarray(matrix[genes][:, d["cols"]].sum(axis=0)).ravel()
                rates = _bin_rates(mk, self._totals[d["cols"]], np.asarray(d["bins"]))
                # each pair's excess over its own ambient reference, so the
                # zero-neighbor bin shows its structured contamination too
                m_ref = float(mk[d["ref_mask"]].sum())
                M_ref = max(float(self._totals[d["cols"]][d["ref_mask"]].sum()), 1.0)
                r_ref = m_ref / M_ref
                for bi, b in enumerate(BIN_LABELS):
                    row = rates[rates["bin"] == b].iloc[0]
                    excess[bi] += max(row["rate"] - r_ref, 0.0) * row["totals"]
                    var_tot[bi] += row["markers"] + (row["totals"] / M_ref) ** 2 * m_ref
                    totals[bi] += row["totals"]
            if n_pairs == 0:
                raise ValueError("No detected pairs to pool")
            sd = np.sqrt(var_tot)
            denom = np.maximum(totals, 1.0)
            return pd.DataFrame(dict(bin=BIN_LABELS,
                rate=excess / denom,
                rate_lo=np.maximum(excess - 1.96 * sd, 0.0) / denom,
                rate_hi=(excess + 1.96 * sd) / denom))
        d = self._pair(source, target)
        genes = d["strict"] if strict and len(d["strict"]) >= 2 else d["pool"]
        mk = np.asarray(matrix[genes][:, d["cols"]].sum(axis=0)).ravel()
        return _bin_rates(mk, self._totals[d["cols"]], np.asarray(d["bins"]))

    def plot_exposure(self, source=None, target=None, correction=None,
                      strict=True, ax=None):
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(4.6, 3.6))
        states = [("before", self._matrix, "#c0392b")]
        if correction is not None:
            after, genes, cells = correction.counts()
            after = _align_counts(after, pd.Index(genes),
                pd.Index([str(c) for c in cells]), self._genes, self._cells)
            states.append(("after cleanup", after, "#2980b9"))
        pooled = source is None
        for name, mat, color in states:
            rates = self._state_rates(mat, source, target, strict)
            x = np.arange(len(BIN_LABELS))
            ax.errorbar(x, rates["rate"] * 1e3,
                yerr=[(rates["rate"] - rates["rate_lo"]) * 1e3,
                      (rates["rate_hi"] - rates["rate"]) * 1e3],
                fmt="o-", ms=4, capsize=2, color=color, label=name)
            if name == "before" and not pooled:
                ax.axhline(rates["rate"].iloc[0] * 1e3, ls=":", c="grey", lw=0.8)
                d = self._pair(source, target)
                genes = d["strict"] if strict and len(d["strict"]) >= 2 else d["pool"]
                mk0 = np.asarray(self._matrix[genes][:, d["cols"]].sum(axis=0)).ravel()
                r_ref = float(mk0[d["ref_mask"]].sum()) / max(
                    float(self._totals[d["cols"]][d["ref_mask"]].sum()), 1.0)
                ax.axhline(r_ref * 1e3, ls="--", c="0.25", lw=0.9)
        ax.set_xticks(range(len(BIN_LABELS)), BIN_LABELS)
        ax.set_xlabel("source-type cells among nearest neighbors")
        ax.set_ylabel("excess admixture-marker rate over pair ambient reference\n(per 1,000 molecules)"
            if pooled else "admixture-marker rate (per 1,000 molecules)\n"
            "(dotted: zero-neighbor rate; dashed: ambient reference)")
        title = ("Admixture exposure profile (pooled over detected pairs)"
            if source is None else f"{source} → {target}")
        ax.set_title(title, fontsize=10)
        if correction is not None:
            ax.legend(frameon=False, fontsize=8)
        return ax

    def plot_reference(self, source, target, ax=None):
        import matplotlib.pyplot as plt

        d = self._pair(source, target)
        lad = d["ref_ladder"]
        if ax is None:
            _, ax = plt.subplots(figsize=(4.2, 3.2))
        ax.plot(lad["K"], lad["rate_fitted"] * 1e3, "-", color="0.6", zorder=1)
        chosen = [f"k{int(k)}" == d["reference_kind"] for k in lad["K"]]
        usable = lad["totals"] >= 2e4
        for i in range(len(lad)):
            ax.plot(lad["K"].iloc[i], lad["rate"].iloc[i] * 1e3,
                "o" if usable.iloc[i] else "o",
                mfc=("#c0392b" if chosen[i] else "#34495e") if usable.iloc[i]
                    else "none",
                mec="#c0392b" if chosen[i] else "#34495e", ms=6, zorder=2)
        ax.axhline(d["reference_rate"] * 1e3, ls="--", c="0.25", lw=0.9)
        ax.set_xscale("log", base=2)
        ax.set_xticks(lad["K"].tolist(), [str(int(k)) for k in lad["K"]])
        ax.set_xlabel("neighborhood size K (zero source cells among K nearest)")
        ax.set_ylabel("admixture-marker rate in reference cells\n(per 1,000 molecules)")
        ax.set_title(f"{source} → {target}: ambient reference", fontsize=10)
        return ax

    def plot_remaining(self, corrections=None, ax=None):
        import matplotlib.pyplot as plt

        corrections = corrections or {}
        states = {"uncorrected": self._matrix}
        for name, corr in corrections.items():
            after, genes, cells = corr.counts()
            states[name] = _align_counts(after, pd.Index(genes),
                pd.Index([str(c) for c in cells]), self._genes, self._cells)
        denom = max(self._totals.sum(), 1)
        ys, los, his = [], [], []
        for mat in states.values():
            tot, var = 0.0, 0.0
            for d in self._pairs.values():
                if not d["detected"]:
                    continue
                mk = np.asarray(mat[d["pool"]][:, d["cols"]].sum(axis=0)).ravel()
                rates = _bin_rates(mk, self._totals[d["cols"]], np.asarray(d["bins"]))
                r0 = float(rates.loc[rates["bin"] == "0", "rate"].iloc[0])
                exposed = rates[rates["bin"] != "0"]
                # extrapolate each pair's pool excess by its marker coverage
                tot += float((np.maximum(exposed["rate"] - r0, 0) *
                    exposed["totals"]).sum()) / d["coverage"]
                M_exp = float(exposed["totals"].sum())
                M_ref = max(float(rates.loc[rates["bin"] == "0", "totals"].iloc[0]), 1)
                var += (float(exposed["markers"].sum()) + (M_exp / M_ref) ** 2 *
                    float(rates.loc[rates["bin"] == "0", "markers"].iloc[0])) / \
                    d["coverage"] ** 2
            ys.append(100 * tot / denom)
            los.append(100 * max(tot - 1.96 * np.sqrt(var), 0) / denom)
            his.append(100 * (tot + 1.96 * np.sqrt(var)) / denom)
        if ax is None:
            _, ax = plt.subplots(figsize=(4.4, 3.4))
        x = np.arange(len(states))
        ax.bar(x, ys, width=0.6, color="#34495e", alpha=0.9)
        ax.errorbar(x, ys, yerr=[np.array(ys) - los, np.array(his) - ys],
            fmt="none", ecolor="black", capsize=3, lw=0.8)
        ax.set_xticks(x, list(states))
        ax.set_ylabel("estimated admixture (% of all molecules)")
        ax.set_title("Remaining admixture by correction\n"
            "estimate over detected pairs; 95% intervals", fontsize=9.5)
        return ax

    # ---- verification ----
    def evaluate(self, correction, warn_uncovered=True):
        after, genes, cells = correction.counts()
        after = _align_counts(after, pd.Index(genes),
            pd.Index([str(c) for c in cells]), self._genes, self._cells)
        rows = []
        for d in self._pairs.values():
            if not d["detected"]:
                continue
            cols = d["cols"]
            tot = self._totals[cols]
            mask = d["ref_mask"]

            def ref_of(counts_g):
                return float(counts_g[mask].sum()) / max(
                    float(tot[mask].sum()), 1.0)

            def power_vs_ref(genes, rates_before, ref_before):
                mk = np.asarray(after[genes][:, cols].sum(axis=0)).ravel()
                rates_a = _bin_rates(mk, tot, d["bins"])
                eb = _excess_vs_ref(rates_before, ref_before)
                ea = _excess_vs_ref(rates_a, min(ref_of(mk), ref_before))
                return np.nan if eb <= 0 else 1 - ea / eb

            sens = power_vs_ref(d["pool"], d["rates"], d["reference_rate"])
            sens_strict = np.nan
            if d["rates_strict"] is not None:
                mks_b = np.asarray(
                    self._matrix[d["strict"]][:, cols].sum(axis=0)).ravel()
                sens_strict = power_vs_ref(d["strict"], d["rates_strict"],
                    ref_of(mks_b))
            rows.append(dict(source=d["source"], target=d["target"],
                admixed_molecules=round(d["excess"] / d["coverage"]),
                excess=round(d["excess"]), sensitivity=sens,
                sensitivity_strict=sens_strict))
        pairs = pd.DataFrame(rows)
        fr_rows = []
        for t, gene_idx in self._native_markers.items():
            tc = self._ctypes.index[self._ctypes == t]
            cols = self._cells.get_indexer(tc)
            cols = cols[cols >= 0]
            before = float(self._matrix[gene_idx][:, cols].sum())
            post = float(after[gene_idx][:, cols].sum())
            fr_rows.append(dict(cell_type=t, own_marker_molecules=before,
                false_removal=1 - post / max(before, 1)))
        false_removal = pd.DataFrame(fr_rows)
        severe = false_removal[(false_removal["false_removal"] > 0.25)
            & (false_removal["own_marker_molecules"] >= OWN_MARKER_WARN_MIN)]
        for _, row in severe.iterrows():
            warnings.warn(
                f"Correction removed {100 * row['false_removal']:.0f}% of "
                f"{row['cell_type']}'s own-marker molecules - severe "
                "over-removal of near-surely-genuine expression", stacklevel=2)
        if warn_uncovered and len(pairs):
            # Warn from the measured removal, not from the rule list:
            # ensemble members can remove molecules for pairs absent from
            # the primary rules, and rules can exist yet remove nothing.
            for _, r in pairs.iterrows():
                if (np.isfinite(r["sensitivity"]) and r["sensitivity"] < 0.2
                        and r["excess"] >= self.params["min_excess"] * 5):
                    warnings.warn(
                        f"Correction removed only {100 * max(r['sensitivity'], 0):.0f}% "
                        f"of the estimated ~{r['admixed_molecules']:,} admixed "
                        f"molecules from {r['source']} into {r['target']}",
                        stacklevel=2)
        return CellAdmixCleanupReport(self, correction, pairs, false_removal)


def _align_counts(matrix, genes, cells, ref_genes, ref_cells):
    """Reindex a counts matrix onto the audit's gene/cell order."""
    from scipy import sparse

    matrix = matrix.tocsc()
    gi = ref_genes.get_indexer(genes)
    ci = ref_cells.get_indexer([str(c) for c in cells])
    out = sparse.lil_matrix((len(ref_genes), len(ref_cells)))
    gmask = gi >= 0
    cmask = ci >= 0
    sub = matrix[np.flatnonzero(gmask)][:, np.flatnonzero(cmask)]
    out[np.ix_(gi[gmask], ci[cmask])] = sub.todense()
    return out.tocsc()


class CellAdmixCleanupReport:
    """Per-pair cleanup verification produced by ``audit.evaluate()``."""

    def __init__(self, audit, correction, pairs, false_removal):
        self.audit = audit
        self.correction = correction
        self._pairs = pairs
        self._false_removal = false_removal

    def pairs(self):
        return self._pairs

    def false_removal(self):
        return self._false_removal

    def summary(self):
        p, fr = self._pairs, self._false_removal
        total = p["admixed_molecules"].sum()
        removed = (p["admixed_molecules"] * p["sensitivity"].clip(lower=0)).sum()
        return dict(detected_pairs=len(p),
            estimated_admixed_molecules=int(total),
            leakage_removed_overall=removed / max(total, 1),
            median_pair_sensitivity=float(p["sensitivity"].median()),
            own_marker_false_removal=float(
                (fr["false_removal"] * fr["own_marker_molecules"]).sum() /
                max(fr["own_marker_molecules"].sum(), 1)))

    def plot_cleanup(self, ax=None):
        import matplotlib.pyplot as plt

        p = self._pairs.sort_values("admixed_molecules")
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 0.22 * len(p) + 1.5))
        y = np.arange(len(p))
        ax.barh(y, p["sensitivity"].clip(lower=0), color="#34495e", alpha=0.9)
        ax.scatter(np.full(len(p), 1.06), y, s=np.sqrt(p["admixed_molecules"]) / 3,
            color="#e67e22", alpha=0.85)
        ax.set_yticks(y, [f"{s} → {t}" for s, t in zip(p["source"], p["target"])],
            fontsize=7)
        ax.set_xlim(0, 1.12)
        s = self.summary()
        ax.set_xlabel("estimated cleanup sensitivity")
        ax.set_title("Cleanup verification by cell-type pair\n"
            f"leakage removed overall: {100*s['leakage_removed_overall']:.0f}%; "
            f"own-marker false removal: {100*s['own_marker_false_removal']:.2f}%",
            fontsize=9.5)
        return ax
