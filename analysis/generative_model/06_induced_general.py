"""General induced-gene detection: for one target type, test every gene -
any ownership, both directions - for expression change associated with each
source's proximity, on the model's retained counts (transfer and ambient
removed), jointly across all sources' exposures.

For each gene g of target cells, the retained count is modeled as
    kept_cg ~ quasi-Poisson( t_c * exp(b0_g + sum_S bS_g * e_cS) )
where e_cS is the number of S cells among the cell's 15 nearest (capped at
4) and t_c the cell's retained total. The joint fit attributes a change to
the source whose own exposure carries it, so spatially correlated
neighborhoods do not cross-attribute. Significance uses the Fisher
information scaled by the Pearson dispersion (quasi-likelihood), so
overdispersion beyond Poisson widens the errors instead of inflating z.

Run after 05_figure.py (uses its cached production-fit rates for the
target). Output: results/gm_induced_general_<target>.csv.
"""
import importlib.util
import os

import numpy as np
import pandas as pd

GM = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
TARGET = os.environ.get("GM_IND_TARGET", "Ductal/tumor epithelial")

spec = importlib.util.spec_from_file_location("gm_model",
    os.path.join(GM, "01_model.py"))
gm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gm)

inp = gm.Inputs()
z = np.load(os.path.join(gm.DATA, "gm_fig_rates_cache.npz"))
rates = {k: z[k] for k in ["own", "cont", "amb", "ind", "eps", "cols"]}
cols = rates["cols"]
import scipy.sparse as sp
Yt = np.asarray(inp.counts[:, cols].todense(), dtype=float).T   # n x G
total_rate = rates["own"] + rates["cont"] + rates["amb"] + rates["ind"] \
    + rates["eps"]
keep = (rates["own"] + rates["ind"] + rates["eps"]) \
    / np.maximum(total_rate, 1e-300)
K = Yt * keep                                                   # retained
n, G = K.shape

# Exposure covariates: one per source with a detected pair into the target.
pos_of_col = {c: i for i, c in enumerate(cols)}
sources, E = [], []
for pair, pi in inp.pair_info.items():
    if pi["T"] != TARGET:
        continue
    e = np.zeros(n)
    e[[pos_of_col[c] for c in pi["cells"]]] = np.minimum(pi["exposure"], 4)
    sources.append(pi["S"])
    E.append(e)
E = np.column_stack(E)                                          # n x S
S_n = len(sources)
print(f"target {TARGET}: {n} cells, {G} genes, sources: {sources}")

t_kept = K.sum(axis=1)
okc = t_kept > 10
K, E, t_kept = K[okc], E[okc], t_kept[okc]
n = okc.sum()
X = np.column_stack([np.ones(n), E])                            # intercept

def glm_gene(y):
    """Poisson IRLS with offset log(t_kept); quasi-likelihood errors."""
    beta = np.zeros(S_n + 1)
    beta[0] = np.log(max(y.sum(), 0.5) / t_kept.sum())
    for _ in range(25):
        mu = t_kept * np.exp(np.clip(X @ beta, -30, 5))
        Wx = X * mu[:, None]
        H = X.T @ Wx
        g = X.T @ (y - mu)
        try:
            step = np.linalg.solve(H + 1e-8 * np.eye(S_n + 1), g)
        except np.linalg.LinAlgError:
            return None
        beta = beta + step
        if np.abs(step).max() < 1e-6:
            break
    mu = t_kept * np.exp(np.clip(X @ beta, -30, 5))
    disp = max(((y - mu) ** 2 / np.maximum(mu, 1e-9)).sum() / max(n - S_n - 1, 1), 1.0)
    try:
        cov = np.linalg.inv(X.T @ (X * mu[:, None]) + 1e-8 * np.eye(S_n + 1))
    except np.linalg.LinAlgError:
        return None
    se = np.sqrt(np.maximum(np.diag(cov), 1e-300) * disp)
    return beta, se, disp

rows = []
for gi in range(G):
    y = K[:, gi]
    if y.sum() < 100:
        continue
    res = glm_gene(y)
    if res is None:
        continue
    beta, se, disp = res
    for j, S in enumerate(sources):
        rows.append(dict(target=TARGET, source=S, gene=inp.genes[gi],
            retained_molecules=float(y.sum()),
            log_fold_per_neighbor=float(beta[j + 1]),
            z=float(beta[j + 1] / se[j + 1]), dispersion=float(disp),
            owner=str(inp.top_type[gi])))

d = pd.DataFrame(rows)
out = os.path.join(gm.RESULTS,
    f"gm_induced_general_{TARGET.split('/')[0].split()[0].lower()}.csv")
d.to_csv(out, index=False)
print(f"wrote {out}: {len(d)} gene-source rows")

for S in sources:
    ds = d[d.source == S]
    up = ds[(ds.z > 6) & (ds.log_fold_per_neighbor > 0.05)]
    dn = ds[(ds.z < -6) & (ds.log_fold_per_neighbor < -0.05)]
    if len(up) or len(dn):
        print(f"\n== {S} -> {TARGET}")
        for tag, t in [("up", up.nlargest(8, "z")),
                       ("down", dn.nsmallest(8, "z"))]:
            for r in t.itertuples():
                print(f"  {tag:4s} {r.gene:10s} owner={r.owner.split('/')[0][:10]:10s} "
                      f"fold/neighbor={np.exp(r.log_fold_per_neighbor):5.2f} "
                      f"z={r.z:7.1f} kept={r.retained_molecules:9.0f}")
print("GENERAL SCREEN DONE")
