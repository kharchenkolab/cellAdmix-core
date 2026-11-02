"""Cut one Xenium bundle into a grid of pseudo-sample bundles.

Cells are binned into an nx-by-ny grid by quantiles of their centroids (x
first, then y within each x bin), so every tile holds about the same number
of cells. Each tile becomes a bundle directory with filtered
transcripts.parquet and cells.parquet, a copy of experiment.xenium, and
symlinks to the morphology images and boundary files of the source bundle;
coordinates stay in the original frame, so the shared images remain valid.
Molecules not assigned to a cell are dropped.

An optional batch effect thins the molecules of chosen tiles gene by gene:
for gene g a log-efficiency deviation eps_g ~ N(0, sigma^2) is drawn per
thinned tile and each molecule of g is kept with probability min(1,
exp(eps_g)).

Usage:
  make_tiles.py <bundle_dir> <out_dir> --prefix NAME [--grid 2x2]
                [--thin TILE:SIGMA ...] [--annotation CSV] [--seed N]

Writes <out_dir>/<prefix>_t<k>/ for k = 1..nx*ny, plus
<out_dir>/<prefix>_t<k>_annotation.csv.gz when --annotation is given, and
<out_dir>/<prefix>_tiles.json.
"""
import argparse
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bundle_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--grid", default="2x2")
    ap.add_argument("--thin", action="append", default=[], help="TILE:SIGMA, repeatable")
    ap.add_argument("--annotation", default=None)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    nx, ny = (int(v) for v in args.grid.lower().split("x"))
    thin = {int(item.split(":")[0]): float(item.split(":")[1]) for item in args.thin}

    src = Path(args.bundle_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cells = pq.read_table(src / "cells.parquet")
    cx = cells.column("x_centroid").to_numpy()
    cy = cells.column("y_centroid").to_numpy()
    x_bin = np.minimum((np.searchsorted(np.quantile(cx, np.linspace(0, 1, nx + 1)[1:-1]), cx, side="right")), nx - 1)
    tile = np.zeros(len(cx), dtype=int)
    for b in range(nx):
        rows = np.flatnonzero(x_bin == b)
        edges = np.quantile(cy[rows], np.linspace(0, 1, ny + 1)[1:-1])
        y_bin = np.minimum(np.searchsorted(edges, cy[rows], side="right"), ny - 1)
        tile[rows] = b * ny + y_bin + 1
    cell_ids = cells.column("cell_id").to_pandas()
    tile_of_cell = pd.Series(tile, index=cell_ids.astype(str).to_numpy())

    transcripts = pq.read_table(src / "transcripts.parquet")
    tx_cell = transcripts.column("cell_id").to_pandas().astype(str)
    tx_tile = tile_of_cell.reindex(tx_cell.to_numpy()).fillna(0).to_numpy().astype(int)
    genes = pc.dictionary_encode(transcripts.column("feature_name")).combine_chunks()
    gene_index = genes.indices.to_numpy(zero_copy_only=False).astype(int)
    gene_names = [str(g) for g in genes.dictionary.to_pylist()]

    report = {"source": str(src), "grid": args.grid, "seed": args.seed, "thin": thin, "tiles": {}}
    for k in range(1, nx * ny + 1):
        name = f"{args.prefix}_t{k}"
        dest = out / name
        dest.mkdir(exist_ok=True)
        tx_mask = tx_tile == k
        sigma = thin.get(k, 0.0)
        if sigma > 0:
            rng = np.random.default_rng(args.seed * 1000 + k)
            eps = rng.normal(0.0, sigma, size=len(gene_names))
            keep_prob = np.minimum(1.0, np.exp(eps))[gene_index]
            tx_mask = tx_mask & (rng.uniform(size=len(gene_index)) < keep_prob)
            report["tiles"][name] = {"gene_log_efficiency_sd": float(np.std(eps))}
        else:
            report["tiles"][name] = {}
        tile_cells = cells.filter(pa.array(tile == k))
        tile_tx = transcripts.filter(pa.array(tx_mask))
        pq.write_table(tile_cells, dest / "cells.parquet")
        pq.write_table(tile_tx, dest / "transcripts.parquet")
        shutil.copy(src / "experiment.xenium", dest / "experiment.xenium")
        for entry in os.listdir(src):
            if entry in ("cells.parquet", "transcripts.parquet", "experiment.xenium",
                         "cells.csv.gz", "transcripts.csv.gz"):
                continue
            target = dest / entry
            if not target.exists():
                os.symlink((src / entry).resolve(), target)
        report["tiles"][name].update({"cells": tile_cells.num_rows, "transcripts": tile_tx.num_rows,
                                      "x_range": [float(cx[tile == k].min()), float(cx[tile == k].max())],
                                      "y_range": [float(cy[tile == k].min()), float(cy[tile == k].max())]})
        if args.annotation:
            ann = pd.read_csv(args.annotation)
            label_col = "merged_annotation" if "merged_annotation" in ann.columns else "cell_type"
            keep_ids = set(tile_cells.column("cell_id").to_pandas().astype(str))
            sub = ann[ann["cell_id"].astype(str).isin(keep_ids)]
            sub[["cell_id", label_col]].rename(columns={label_col: "cell_type"}).to_csv(
                out / f"{name}_annotation.csv.gz", index=False)
            report["tiles"][name]["annotated"] = int(len(sub))
    (out / f"{args.prefix}_tiles.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    sys.exit(main())
