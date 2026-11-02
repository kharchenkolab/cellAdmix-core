"""Cut one Xenium bundle into two pseudo-replicate bundles.

The section is split at the median x coordinate of the cells (or into a
checkerboard of square tiles with --tile). Each half gets its own bundle
directory with filtered transcripts.parquet and cells.parquet, a copy of
experiment.xenium, and symlinks to the morphology images and boundary
files of the source bundle. Coordinates are left in the original frame,
so the shared images stay valid for both halves.

An optional batch effect thins the molecules of the second half gene by
gene: for gene g a log-efficiency deviation eps_g ~ N(0, sigma^2) is
drawn and each molecule of g is kept with probability min(1, exp(eps_g)).
--uniform p additionally keeps every molecule with probability p.

Usage:
  make_pseudo_replicates.py <bundle_dir> <out_dir> [--sigma S] [--uniform P]
                            [--seed N] [--tile MICRONS] [--annotation CSV]
Writes <out_dir>/<name>/ for name in (left, right) or (even, odd) with
--tile, plus <out_dir>/<name>_annotation.csv.gz when --annotation is given.
"""
import argparse
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bundle_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--sigma", type=float, default=0.0)
    ap.add_argument("--uniform", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--tile", type=float, default=0.0)
    ap.add_argument("--annotation", default=None)
    args = ap.parse_args()

    src = Path(args.bundle_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cells = pq.read_table(src / "cells.parquet")
    transcripts = pq.read_table(src / "transcripts.parquet")
    cell_ids = cells.column("cell_id")
    cx = cells.column("x_centroid").to_numpy()
    cy = cells.column("y_centroid").to_numpy()

    if args.tile > 0:
        side = ((np.floor(cx / args.tile) + np.floor(cy / args.tile)) % 2).astype(int)
        names = ("even", "odd")
    else:
        side = (cx >= np.median(cx)).astype(int)
        names = ("left", "right")
    side_by_cell = dict(zip(cell_ids.to_pylist(), side.tolist()))

    tx_cells = transcripts.column("cell_id").to_pylist()
    tx_side = np.array([side_by_cell.get(c, -1) for c in tx_cells])
    rng = np.random.default_rng(args.seed)
    genes = transcripts.column("feature_name").to_pylist()
    gene_names = sorted(set(genes))
    eps = rng.normal(0.0, args.sigma, size=len(gene_names)) if args.sigma > 0 else np.zeros(len(gene_names))
    keep_prob = {g: min(1.0, float(np.exp(e))) * args.uniform for g, e in zip(gene_names, eps)}
    draw = rng.uniform(size=len(genes))
    thin_keep = np.array([draw[i] < keep_prob[g] for i, g in enumerate(genes)])

    report = {"source": str(src), "sigma": args.sigma, "uniform": args.uniform,
              "seed": args.seed, "tile": args.tile, "halves": {}}
    for s, name in enumerate(names):
        dest = out / name
        dest.mkdir(exist_ok=True)
        cell_mask = pa.array(side == s)
        tx_mask = tx_side == s
        if s == 1:
            tx_mask = tx_mask & thin_keep
        half_cells = cells.filter(cell_mask)
        half_tx = transcripts.filter(pa.array(tx_mask))
        pq.write_table(half_cells, dest / "cells.parquet")
        pq.write_table(half_tx, dest / "transcripts.parquet")
        shutil.copy(src / "experiment.xenium", dest / "experiment.xenium")
        for entry in os.listdir(src):
            if entry in ("cells.parquet", "transcripts.parquet", "experiment.xenium",
                         "cells.csv.gz", "transcripts.csv.gz"):
                continue
            target = dest / entry
            if not target.exists():
                os.symlink((src / entry).resolve(), target)
        report["halves"][name] = {"cells": half_cells.num_rows, "transcripts": half_tx.num_rows}
        if args.annotation:
            import pandas as pd
            ann = pd.read_csv(args.annotation)
            id_col = "cell_id"
            label_col = "merged_annotation" if "merged_annotation" in ann.columns else "cell_type"
            keep_ids = set(half_cells.column("cell_id").to_pylist())
            sub = ann[ann[id_col].astype(str).isin({str(c) for c in keep_ids})]
            sub[[id_col, label_col]].rename(columns={label_col: "cell_type"}).to_csv(
                out / f"{name}_annotation.csv.gz", index=False)
            report["halves"][name]["annotated"] = int(len(sub))
    if args.sigma > 0:
        report["gene_log_efficiency_sd"] = float(np.std(eps))
    (out / "pseudo_replicates.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    sys.exit(main())
