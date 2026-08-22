"""Robustness and circularity checks for the generative model (pancreas).

A. Ownership under correction: gene ownership (the type with the highest
   pseudobulk rate) is computed from observed, contaminated counts, and
   ownership decides which genes are protected and screened. Recompute it
   from the model's retained counts and count the flips.
C. Marginal versus joint dose: each pair's contamination dose is measured
   marginally (guide-gene content against that source's exposure). Compare
   with a joint Poisson regression on all sources' exposures, which
   controls for correlated neighborhoods, and report per-pair divergence.
E. Detection against the raw counts: rerun the general induced-change
   test on the observed counts (no removal) and compare hit sets - a
   change reported only on the corrected counts depends on the transfer
   split; one reported on both does not.

Requires the per-target rate caches written by 06_induced_general.py.
"""
import importlib.util
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp

GM = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
spec = importlib.util.spec_from_file_location("gm_model",
    os.path.join(GM, "01_model.py"))
gm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gm)
inp = gm.Inputs()

def rates_of(T):
    cache = os.path.join(gm.DATA,
        f"gm_rates_cache_{T.split('/')[0].split()[0].lower()}.npz")
    z = np.load(cache)
    return {k: z[k] for k in ["own", "cont", "amb", "ind", "eps", "cols"]}

# ---- A: ownership stability ------------------------------------------------
print("== A. gene ownership recomputed on corrected counts ==")
G = inp.G
prof_corr = np.zeros((G, len(inp.type_names)))
kept_by_type = {}
for ti, T in enumerate(inp.type_names):
    r = rates_of(T)
    cols = r["cols"]
    Yt = np.asarray(inp.counts[:, cols].todense(), dtype=float).T
    tot = r["own"] + r["cont"] + r["amb"] + r["ind"] + r["eps"]
    keep = (r["own"] + r["ind"] + r["eps"]) / np.maximum(tot, 1e-300)
    K = Yt * keep
    kept_by_type[T] = (cols, K)
    csum = K.sum(axis=0)
    prof_corr[:, ti] = csum / max(csum.sum(), 1.0) * 1e6
top_corr = np.array(inp.type_names)[np.argmax(prof_corr, axis=1)]
flips = np.flatnonzero(top_corr != inp.top_type)
print(f"flips: {len(flips)} of {G} genes")
for g in flips:
    o, c = inp.top_type[g], top_corr[g]
    cpm_o = inp.profiles_cpm[g, inp.type_names.index(o)]
    cpm_c = inp.profiles_cpm[g, inp.type_names.index(c)]
    print(f"  {inp.genes[g]:10s} {o.split('/')[0][:10]:10s} -> "
          f"{c.split('/')[0][:10]:10s}  observed cpm {cpm_o:8.0f} vs {cpm_c:8.0f}")

# ---- C: marginal versus joint dose -----------------------------------------
print("\n== C. marginal versus joint dose per pair ==")
rows_c = []
for T in inp.type_names:
    plist = [(p, pi) for p, pi in inp.pair_info.items() if pi["T"] == T]
    if not plist:
        continue
    cols, K = kept_by_type[T]
    pos_of_col = {c: i for i, c in enumerate(cols)}
    n = len(cols)
    Yt = np.asarray(inp.counts[:, cols].todense(), dtype=float).T
    t = Yt.sum(axis=1)
    E = np.zeros((n, len(plist)))
    for j, (p, pi) in enumerate(plist):
        E[[pos_of_col[c] for c in pi["cells"]], j] = np.minimum(pi["exposure"], 4)
    for j, (p, pi) in enumerate(plist):
        gidx = inp.gidx(pi["pool"])
        y = Yt[:, gidx].sum(axis=1)
        # marginal: isotonic excess over the zero-exposure stratum
        e = E[:, j]
        strata = np.unique(e)
        m_s = np.array([y[e == s].sum() for s in strata])
        t_s = np.array([t[e == s].sum() for s in strata])
        rate = gm.pava_increasing(m_s / np.maximum(t_s, 1.0), t_s)
        marg = float(((rate - rate[0]) * t_s)[1:].sum())
        # joint: Poisson regression on all exposures, dose = predicted
        # minus the all-zero-exposure prediction
        X = np.column_stack([np.ones(n), E])
        beta = np.zeros(X.shape[1]); beta[0] = np.log(max(y.sum(), .5) / t.sum())
        for _ in range(25):
            mu = t * np.exp(np.clip(X @ beta, -30, 5))
            try:
                step = np.linalg.solve(X.T @ (X * mu[:, None])
                                       + 1e-8 * np.eye(X.shape[1]),
                                       X.T @ (y - mu))
            except np.linalg.LinAlgError:
                break
            beta += step
            if np.abs(step).max() < 1e-6:
                break
        X0 = X.copy(); X0[:, 1 + j] = 0.0
        mu_full = t * np.exp(np.clip(X @ beta, -30, 5))
        mu_wo = t * np.exp(np.clip(X0 @ beta, -30, 5))
        joint = float((mu_full - mu_wo).sum())
        rows_c.append(dict(pair=p, marginal=marg, joint=joint,
                           ratio=joint / max(marg, 1.0)))
dc = pd.DataFrame(rows_c)
dc.to_csv(os.path.join(gm.RESULTS, "gm_dose_marginal_vs_joint.csv"),
          index=False)
print(f"pairs: {len(dc)}; joint/marginal ratio median "
      f"{dc.ratio.median():.2f}, quartiles "
      f"{dc.ratio.quantile(.25):.2f}-{dc.ratio.quantile(.75):.2f}")
big = dc[(dc.ratio < 0.67) | (dc.ratio > 1.5)]
print(f"pairs diverging beyond 1.5-fold: {len(big)}")
print(big.sort_values("ratio").to_string(index=False))

# ---- E: detection on raw versus corrected counts ---------------------------
print("\n== E. general screen on raw versus corrected counts ==")
d_corr = pd.read_csv(os.path.join(gm.RESULTS,
    "gm_induced_general_pancreas.csv"))
hits_corr = d_corr[d_corr.hit]
def glm_all(K, E, t_kept):
    n, S_n = E.shape
    X = np.column_stack([np.ones(n), E])
    res = {}
    for gi in range(K.shape[1]):
        y = K[:, gi]
        if y.sum() < 100:
            continue
        beta = np.zeros(S_n + 1)
        beta[0] = np.log(max(y.sum(), .5) / t_kept.sum())
        bad = False
        for _ in range(25):
            mu = t_kept * np.exp(np.clip(X @ beta, -30, 5))
            try:
                step = np.linalg.solve(X.T @ (X * mu[:, None])
                                       + 1e-8 * np.eye(S_n + 1),
                                       X.T @ (y - mu))
            except np.linalg.LinAlgError:
                bad = True; break
            beta += step
            if np.abs(step).max() < 1e-6:
                break
        if bad:
            continue
        mu = t_kept * np.exp(np.clip(X @ beta, -30, 5))
        disp = max(((y - mu) ** 2 / np.maximum(mu, 1e-9)).sum()
                   / max(n - S_n - 1, 1), 1.0)
        cov = np.linalg.inv(X.T @ (X * mu[:, None]) + 1e-8 * np.eye(S_n + 1))
        se = np.sqrt(np.maximum(np.diag(cov), 1e-300) * disp)
        res[gi] = (beta, se)
    return res

rows_e = []
for T in inp.type_names:
    cols, K = kept_by_type[T]
    Yt = np.asarray(inp.counts[:, cols].todense(), dtype=float).T
    pos_of_col = {c: i for i, c in enumerate(cols)}
    plist = [(p, pi) for p, pi in inp.pair_info.items() if pi["T"] == T]
    if not plist:
        continue
    E = np.zeros((len(cols), len(plist)))
    for j, (p, pi) in enumerate(plist):
        E[[pos_of_col[c] for c in pi["cells"]], j] = np.minimum(pi["exposure"], 4)
    t_raw = Yt.sum(axis=1)
    okc = t_raw > 10
    raw = glm_all(Yt[okc], E[okc], t_raw[okc])
    for gi, (beta, se) in raw.items():
        for j, (p, pi) in enumerate(plist):
            z = beta[j + 1] / se[j + 1]
            hit = abs(z) > 6 and abs(beta[j + 1]) > 0.05
            rows_e.append(dict(target=T, source=pi["S"],
                gene=inp.genes[gi], z_raw=float(z),
                lfc_raw=float(beta[j + 1]), hit_raw=bool(hit)))
de = pd.DataFrame(rows_e)
m = hits_corr.merge(de, on=["target", "source", "gene"], how="left")
m["clean"] = (m.removed_frac < 0.10) & (m.psi_share_e4 < 20)
m["sign_consistent"] = np.sign(m.lfc_raw.fillna(0)) == np.sign(
    m.log_fold_per_neighbor)
both = m.hit_raw.fillna(False)
print(f"corrected-count hits: {len(m)}; also present on raw counts: "
      f"{both.sum()} ({both.mean():.0%}); sign-consistent on raw: "
      f"{m.sign_consistent.mean():.0%}")
print(f"clean tier: {m[m.clean].hit_raw.fillna(False).mean():.0%} on raw, "
      f"sign-consistent {m[m.clean].sign_consistent.mean():.0%}")
only_corr = m[m.clean & ~m.hit_raw.fillna(False) & (m.z.abs() > 10)]
print(f"clean hits absent on raw at |z|>10: {len(only_corr)}")
if len(only_corr):
    print(only_corr.nlargest(6, "z")[["target", "source", "gene", "z",
        "z_raw"]].to_string(index=False))
m.to_csv(os.path.join(gm.RESULTS, "gm_induced_raw_vs_corrected.csv"),
         index=False)
print("ROBUSTNESS DONE")
