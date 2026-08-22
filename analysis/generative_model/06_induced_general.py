"""General induced-change detection on the model's retained counts.

For every target type and every gene - any ownership, both directions -
tests whether the gene's retained expression (transfer and ambient
removed by the production fit) changes with each source's proximity,
jointly across all sources' exposures so that spatially correlated
neighborhoods do not cross-attribute:

    kept_cg ~ quasi-Poisson( t_c * exp(b0_g + sum_S bS_g * e_cS) )

with e_cS the number of S cells among the cell's 15 nearest (capped at 4)
and t_c the cell's retained total. Errors use the Fisher information
scaled by the Pearson dispersion, so overdispersion beyond Poisson widens
them rather than inflating significance.

Guards reported with every hit:
  removed_frac  - fraction of the gene's observed content in the target
                  cells that the fit removed as transfer or ambient; hits
                  on genes with substantial removal can be residue of
                  imperfect removal rather than expression change.
  psi_share     - the gene's share of the source's interface-local
                  transfer profile (per 10,000); large values mark genes
                  whose transfer from this source is a plausible confound.
Calibration: the identical fit is rerun with the exposure rows permuted
across cells (one joint permutation, preserving the correlation between
sources), and the number of hits at the reporting threshold is recorded.

Outputs: results/gm_induced_general_pancreas.csv (all gene-source rows)
and results/gm_induced_general_calibration.csv (hit counts per target,
observed versus permuted).
"""
import importlib.util
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp

GM = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
Z_HIT = 6.0
LFC_HIT = 0.05

spec = importlib.util.spec_from_file_location("gm_model",
    os.path.join(GM, "01_model.py"))
gm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gm)

inp = gm.Inputs()

# ---- production fit, keeping each target's component rates (cached) --------
def target_rates(T):
    cache = os.path.join(gm.DATA,
        f"gm_rates_cache_{T.split('/')[0].split()[0].lower()}.npz")
    if os.path.exists(cache):
        z = np.load(cache)
        return {k: z[k] for k in ["own", "cont", "amb", "ind", "eps", "cols"]}
    return None

missing = [T for T in inp.type_names if target_rates(T) is None]
if missing:
    print("fitting production model for:", missing)
    rng = np.random.default_rng(1)
    own_frac = None
    for rnd in range(gm.N_OUTER):
        diag = {}
        own_frac_new = {}
        for T in inp.type_names:
            if rnd == gm.N_OUTER - 1:
                out = gm.fit_target(inp, T, own_frac, "production", rng, diag,
                                    return_rates=True)
                if out[0] is None:
                    continue
                res, of, rates = out
                cache = os.path.join(gm.DATA,
                    f"gm_rates_cache_{T.split('/')[0].split()[0].lower()}.npz")
                np.savez_compressed(cache, own=rates["own"],
                    cont=rates["cont"], amb=rates["amb"], ind=rates["ind"],
                    eps=rates["eps"], cols=rates["cols"])
            else:
                out = gm.fit_target(inp, T, own_frac, "production", rng, diag)
                if out[0] is None:
                    continue
                res, of = out
            own_frac_new[T] = of
        own_frac = own_frac_new

def glm_genes(K, E, t_kept):
    """Per-gene Poisson IRLS with offset log(t_kept), quasi errors."""
    n, S_n = E.shape
    X = np.column_stack([np.ones(n), E])
    out = []
    for gi in range(K.shape[1]):
        y = K[:, gi]
        if y.sum() < 100:
            out.append(None)
            continue
        beta = np.zeros(S_n + 1)
        beta[0] = np.log(max(y.sum(), 0.5) / t_kept.sum())
        okfit = True
        for _ in range(25):
            mu = t_kept * np.exp(np.clip(X @ beta, -30, 5))
            H = X.T @ (X * mu[:, None])
            g = X.T @ (y - mu)
            try:
                step = np.linalg.solve(H + 1e-8 * np.eye(S_n + 1), g)
            except np.linalg.LinAlgError:
                okfit = False
                break
            beta = beta + step
            if np.abs(step).max() < 1e-6:
                break
        if not okfit:
            out.append(None)
            continue
        mu = t_kept * np.exp(np.clip(X @ beta, -30, 5))
        disp = max(((y - mu) ** 2 / np.maximum(mu, 1e-9)).sum()
                   / max(n - S_n - 1, 1), 1.0)
        try:
            cov = np.linalg.inv(X.T @ (X * mu[:, None])
                                + 1e-8 * np.eye(S_n + 1))
        except np.linalg.LinAlgError:
            out.append(None)
            continue
        se = np.sqrt(np.maximum(np.diag(cov), 1e-300) * disp)
        out.append((beta, se, disp))
    return out

rows, calib = [], []
rng_perm = np.random.default_rng(11)
for T in inp.type_names:
    rates = target_rates(T)
    if rates is None:
        continue
    cols = rates["cols"]
    Yt = np.asarray(inp.counts[:, cols].todense(), dtype=float).T
    total_rate = rates["own"] + rates["cont"] + rates["amb"] + rates["ind"] \
        + rates["eps"]
    keep = (rates["own"] + rates["ind"] + rates["eps"]) \
        / np.maximum(total_rate, 1e-300)
    K = Yt * keep
    removed_frac_g = 1.0 - K.sum(axis=0) / np.maximum(Yt.sum(axis=0), 1.0)

    pos_of_col = {c: i for i, c in enumerate(cols)}
    sources, E, psi_list = [], [], []
    t_owned = inp.top_type == T
    for pair, pi in inp.pair_info.items():
        if pi["T"] != T:
            continue
        e = np.zeros(len(cols))
        e[[pos_of_col[c] for c in pi["cells"]]] = np.minimum(pi["exposure"], 4)
        sources.append(pi["S"])
        E.append(e)
        p = inp.near_source_profile(pi["S"], T).copy()
        p[t_owned] = 0.0
        psi_list.append(p)
    if not sources:
        continue
    E = np.column_stack(E)
    t_kept = K.sum(axis=1)
    okc = t_kept > 10
    K, E, t_kept = K[okc], E[okc], t_kept[okc]
    print(f"{T}: {okc.sum()} cells, sources {len(sources)}", flush=True)

    fits = glm_genes(K, E, t_kept)
    E_perm = E[rng_perm.permutation(E.shape[0])]
    fits_perm = glm_genes(K, E_perm, t_kept)

    n_obs = n_perm = 0
    for gi in range(K.shape[1]):
        for tag, f in [("obs", fits[gi]), ("perm", fits_perm[gi])]:
            if f is None:
                continue
            beta, se, disp = f
            for j, S in enumerate(sources):
                zj = beta[j + 1] / se[j + 1]
                hit = abs(zj) > Z_HIT and abs(beta[j + 1]) > LFC_HIT
                if tag == "perm":
                    n_perm += hit
                    continue
                n_obs += hit
                rows.append(dict(target=T, source=S, gene=inp.genes[gi],
                    owner=str(inp.top_type[gi]),
                    retained_molecules=float(K[:, gi].sum()),
                    log_fold_per_neighbor=float(beta[j + 1]),
                    z=float(zj), dispersion=float(disp),
                    removed_frac=float(removed_frac_g[gi]),
                    psi_share_e4=float(psi_list[j][gi] * 1e4),
                    hit=bool(hit)))
    calib.append(dict(target=T, tests=int(sum(f is not None for f in fits))
                      * len(sources), hits_observed=n_obs,
                      hits_permuted=n_perm))
    print(f"  hits: observed {n_obs}, permuted {n_perm}", flush=True)

d = pd.DataFrame(rows)
d.to_csv(os.path.join(gm.RESULTS, "gm_induced_general_pancreas.csv"),
         index=False)
pd.DataFrame(calib).to_csv(
    os.path.join(gm.RESULTS, "gm_induced_general_calibration.csv"),
    index=False)
print("\ncalibration (hits at |z|>%.0f, |log fold|>%.2f):" % (Z_HIT, LFC_HIT))
print(pd.DataFrame(calib).to_string(index=False))
print("GENERAL SCREEN DONE")
