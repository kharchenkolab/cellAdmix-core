"""Parity and performance: the C++ generative fit (the packaged
implementation, celladmix._core.fit_generative) versus the Python
implementation's pancreas production run, on the identical problem
specification (same pools, strict tiers, factor alignment, counts)."""
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scipy.io
import scipy.sparse as sp

from celladmix import _core

GM = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
DATA = f"{GM}/data"
RUN = ("/home/pkharchenko/cellAdmix/cellAdmix-core/examples/"
       "xenium_pancreas_membrane_377_full/out/runs/fit_manual_rank9_invsqrt_kl")

genes = [l.strip() for l in open(f"{DATA}/genes.txt")]
cells = [l.strip() for l in open(f"{DATA}/cells.txt")]
gene_of = {g: i for i, g in enumerate(genes)}
counts = sp.csc_matrix(scipy.io.mmread(f"{DATA}/counts.mtx"))

ct = pd.read_csv(f"{DATA}/cell_types.csv")
type_of_cell = dict(zip(ct.cell_id, ct.cell_type))
types = sorted(set(ct.cell_type.dropna()))
tcode = {t: i for i, t in enumerate(types)}
type_codes = [tcode.get(type_of_cell.get(c), -1) for c in cells]

cells_tbl = pq.read_table(f"{RUN}/cells.parquet",
                          columns=["cell_id", "x", "y"]).to_pandas()
xy = cells_tbl.set_index("cell_id").reindex(cells)
x = xy["x"].to_numpy()
y = xy["y"].to_numpy()

pairs_df = pd.read_csv(f"{DATA}/pairs.csv")
pair_specs = []
for r in pairs_df.itertuples():
    pair_specs.append(dict(
        source_type=tcode[r.source], target_type=tcode[r.target],
        pool=[gene_of[g] for g in r.pool.split(";") if g],
        strict=[gene_of[g] for g in str(r.strict).split(";")
                if g and g != "nan"]))

ft = pd.read_csv("/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/"
                 "audit_guided/results/factor_types.csv")
n_factors = 9
alignment = [-1] * n_factors
for r in ft.itertuples():
    alignment[int(r.factor) - 1] = tcode[r.cell_type]

t0 = time.time()
res = _core.fit_generative(
    counts_indptr=[int(v) for v in counts.indptr],
    counts_indices=[int(v) for v in counts.indices],
    counts_values=[float(v) for v in counts.data],
    n_genes=len(genes), cell_ids=cells,
    x=[float(v) for v in x], y=[float(v) for v in y],
    type_codes=[int(v) for v in type_codes], n_types=len(types),
    molecules_parquet=f"{RUN}/molecules.parquet",
    cells_parquet=f"{RUN}/cells.parquet",
    pair_specs=pair_specs, factor_to_type=alignment,
    options=dict(num_threads=12, neighbor_k=15),
)
t_cpp = time.time() - t0
removed_cpp = np.asarray(res["removed"])
print(f"C++ fit: {t_cpp:.1f}s, removed {removed_cpp.sum():,.0f}")

ref = sp.csc_matrix(scipy.io.mmread(f"{DATA}/gm_removed_production.mtx"))
removed_ref = sp.csc_matrix(
    (removed_cpp, counts.indices.copy(), counts.indptr.copy()),
    shape=counts.shape)
print(f"python reference removed {ref.sum():,.0f}")
print(f"total ratio C++/python: {removed_cpp.sum() / ref.sum():.4f}")
diff = removed_ref - ref
denom = max(float(ref.sum()), 1.0)
print(f"sum |difference| / total: {abs(diff).sum() / denom:.4f}")
a = np.asarray(ref.todense()).ravel()
b = np.asarray(removed_ref.todense()).ravel()
keep = (a > 0) | (b > 0)
print(f"entry correlation: {np.corrcoef(a[keep], b[keep])[0, 1]:.6f}")

# per-pair comparison against the python run's dose diagnostics
alpha_py = pd.read_csv(f"{GM}/results/gm_production_alpha.csv")
rows = []
for j, r in enumerate(pairs_df.itertuples()):
    s = [p for p in res["pairs"] if p["pair"] == j]
    if not s:
        continue
    ref_row = alpha_py[alpha_py.pair == r.pair]
    if not len(ref_row):
        continue
    rows.append(dict(pair=r.pair,
                     post_cpp=s[0]["posterior_molecules"],
                     post_py=float(ref_row.posterior_molecules.iloc[0])))
d = pd.DataFrame(rows)
d["ratio"] = d.post_cpp / d.post_py.clip(lower=1)
print(f"per-pair posterior-molecule ratio: median {d.ratio.median():.3f}, "
      f"quartiles {d.ratio.quantile(.25):.3f}-{d.ratio.quantile(.75):.3f}")
print(d.reindex(d.ratio.sub(1).abs().sort_values(ascending=False).index)
      .head(5).to_string(index=False))
print("PARITY DONE")
