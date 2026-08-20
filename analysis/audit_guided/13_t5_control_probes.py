"""T5: negative-control probes as the measured ambient null.

Control probes and blank codewords have no biology: their spatial behavior
inside cells is the pure noise floor. Read them straight from the bundle's
transcript table and compute the same statistics used on strict genes:
per-cell dispersion (nearest same-class neighbor distance vs random),
z-position within the cell, and nucleus association - the empirical null
for T1/T2.

Run from the pancreas example directory.
"""
import numpy as np
import pandas as pd
import pyarrow.csv as pcsv
from scipy.spatial import cKDTree

AG = "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/audit_guided"
rng = np.random.default_rng(17)

tab = pcsv.read_csv("data/transcripts.csv.gz").to_pandas()
tab = tab[(tab.qv >= 20) & (tab.cell_id != "UNASSIGNED")]
is_ctrl = tab.feature_name.str.startswith(("NegControlProbe", "NegControlCodeword",
    "BLANK", "antisense", "UnassignedCodeword", "DeprecatedCodeword"))
print("control molecules in cells:", int(is_ctrl.sum()),
      "of", len(tab), "assigned molecules")

ctrl = tab[is_ctrl]
per_cell = ctrl.groupby("cell_id").size()
multi = per_cell[per_cell >= 2].index
print("cells with >=2 control molecules:", len(multi))

obs_d, null_d, zdev_obs, zdev_null = [], [], [], []
nuc_ctrl = float(ctrl.overlaps_nucleus.mean())
groups = tab[tab.cell_id.isin(multi)].groupby("cell_id")
for cid, cellm in groups:
    cm = cellm[["x_location", "y_location", "z_location"]].to_numpy()
    isc = cellm.feature_name.str.startswith(("NegControlProbe",
        "NegControlCodeword", "BLANK", "antisense", "UnassignedCodeword",
        "DeprecatedCodeword")).to_numpy()
    n_c = int(isc.sum())
    if n_c < 2 or len(cm) < 5:
        continue
    pts = cm[isc]
    tree = cKDTree(pts)
    d, _ = tree.query(pts, k=2)
    obs_d.extend(d[:, 1])
    sel = rng.choice(len(cm), n_c, replace=False)
    tree0 = cKDTree(cm[sel])
    d0, _ = tree0.query(cm[sel], k=2)
    null_d.extend(d0[:, 1])
    z = cm[:, 2]
    pct = (np.argsort(np.argsort(z)) + 0.5) / len(z)
    zdev_obs.extend(np.abs(pct[isc] - 0.5))
    zdev_null.extend(np.abs(pct[sel] - 0.5))

nuc_all = float(tab.overlaps_nucleus.mean())
out = dict(
    n_ctrl_in_cells=int(is_ctrl.sum()),
    cells_multi=len(multi),
    nn_dist_ctrl_median=float(np.median(obs_d)),
    nn_dist_null_median=float(np.median(null_d)),
    dispersion_ratio=float(np.median(obs_d) / max(np.median(null_d), 1e-9)),
    z_absdev_ctrl=float(np.mean(zdev_obs)),
    z_absdev_null=float(np.mean(zdev_null)),
    nucleus_overlap_ctrl=nuc_ctrl,
    nucleus_overlap_all=nuc_all)
pd.DataFrame([out]).to_csv(f"{AG}/results/t5_control_null.csv", index=False)
for k, v in out.items():
    print(f"{k}: {v}")
print("T5 DONE")
