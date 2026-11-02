"""Batch-separating factors in a joint fit over several samples.

For a joint run whose cell ids carry a sample prefix, report per factor the
share of its molecules coming from each sample next to the sample's share
of all molecules (concentration = factor share / sample share), and list
pairs of factors whose gene profiles are near-identical (cosine >= 0.9)
while their molecules come predominantly from different samples: the
signature of one program split by sample.

Usage: joint_split_factors.py <run_dir> <label>
"""
import sys

import numpy as np
import pandas as pd

run_dir, label = sys.argv[1], sys.argv[2]
profiles = pd.read_csv(f"{run_dir}/../../profiles.csv").set_index("factor")
cells = pd.read_parquet(f"{run_dir}/cells.parquet")
cells["sample"] = cells["cell_id"].astype(str).str.split("_", n=1).str[0]
fraction_cols = sorted([c for c in cells.columns if c.startswith("factor_") and c.endswith("_fraction")],
                       key=lambda c: int(c.split("_")[1]))
molecules = cells[fraction_cols].to_numpy(float) * cells["transcript_count"].to_numpy(float)[:, None]
by_sample = pd.DataFrame(molecules, columns=[f"F{i + 1}" for i in range(len(fraction_cols))]).groupby(cells["sample"].to_numpy()).sum()
sample_share = by_sample.sum(axis=1) / by_sample.values.sum()
factor_share = by_sample / by_sample.sum(axis=0)
concentration = factor_share.div(sample_share, axis=0)
H = profiles.to_numpy(float)
norms = np.linalg.norm(H, axis=1, keepdims=True)
cosine = (H / norms) @ (H / norms).T
np.fill_diagonal(cosine, 0)
pairs = []
for i in range(len(H)):
    for j in range(i + 1, len(H)):
        if cosine[i, j] >= 0.9:
            dom_i, dom_j = factor_share.iloc[:, i].idxmax(), factor_share.iloc[:, j].idxmax()
            pairs.append((f"F{i+1}", f"F{j+1}", round(float(cosine[i, j]), 3), dom_i, dom_j, dom_i != dom_j))
print(f"{label}: factors {len(H)}, max off-diagonal cosine {cosine.max():.3f}, "
      f"max sample concentration {concentration.values.max():.2f} (1 = proportional)")
print("  per-factor concentration of the smaller sample:",
      [round(v, 2) for v in concentration.loc[sample_share.idxmin()].tolist()])
print("  near-identical factor pairs:", pairs if pairs else "none")
