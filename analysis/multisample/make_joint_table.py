"""Write several Xenium bundles as one molecule table for the joint stand-in.

Assigned molecules of gene features are concatenated with the sample name
prefixed to every cell id and a sample_id column; the samples are placed
side by side with a large coordinate offset so that no neighborhood query
can cross from one sample to the next. The filters mirror the Xenium store
builder's defaults (unassigned molecules dropped, control features dropped,
no quality threshold). A combined annotation table is written alongside.

Usage: make_joint_table.py <out_prefix> <gap_um> <sample>=<bundle>=<annotation.csv[.gz]> ...
Writes <out_prefix>_molecules.parquet and <out_prefix>_annotation.csv.gz.
"""
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

CONTROL_PREFIXES = ("deprecatedcodeword", "unassignedcodeword", "negcontrol", "negativecontrol",
                    "negative_control", "genomiccontrol", "genomic_control", "intergenic_region", "blank")


def is_gene(name: str) -> bool:
    lower = name.lower()
    return bool(name) and not any(lower.startswith(p) for p in CONTROL_PREFIXES)


def main():
    out_prefix, gap = sys.argv[1], float(sys.argv[2])
    specs = [spec.split("=") for spec in sys.argv[3:]]
    frames, annotations = [], []
    offset = 0.0
    for sample, bundle, annotation_path in specs:
        table = pq.read_table(f"{bundle}/transcripts.parquet",
                              columns=["cell_id", "feature_name", "x_location", "y_location", "z_location", "qv"])
        df = table.to_pandas()
        df["feature_name"] = df["feature_name"].map(lambda v: v.decode() if isinstance(v, bytes) else str(v))
        df["cell_id"] = df["cell_id"].map(lambda v: v.decode() if isinstance(v, bytes) else str(v))
        unassigned = df["cell_id"].isin(["", "0", "-1", "UNASSIGNED", "cell_0"])
        keep = (~unassigned) & df["feature_name"].map(is_gene)
        df = df.loc[keep]
        shift = offset - float(df["x_location"].min())
        out = pd.DataFrame({
            "x": df["x_location"].to_numpy(dtype=float) + shift,
            "y": df["y_location"].to_numpy(dtype=float),
            "z": df["z_location"].to_numpy(dtype=float),
            "gene": df["feature_name"].to_numpy(),
            "cell": sample + "_" + df["cell_id"].to_numpy(),
            "sample_id": sample,
            "qv": df["qv"].to_numpy(dtype=float),
        })
        offset = float(out["x"].max()) + gap
        frames.append(out)
        ann = pd.read_csv(annotation_path)
        ann = ann[ann["cell_type"].notna() & (ann["cell_type"] != "Unlabeled")]
        annotations.append(pd.DataFrame({"cell_id": sample + "_" + ann["cell_id"].astype(str),
                                         "cell_type": ann["cell_type"].astype(str), "sample": sample}))
        print(sample, "molecules", len(out), "x range", round(out["x"].min()), round(out["x"].max()))
    molecules = pd.concat(frames, ignore_index=True)
    pq.write_table(pa.Table.from_pandas(molecules, preserve_index=False), f"{out_prefix}_molecules.parquet")
    pd.concat(annotations, ignore_index=True).to_csv(f"{out_prefix}_annotation.csv.gz", index=False)
    print("total molecules", len(molecules), "cells", molecules["cell"].nunique())


if __name__ == "__main__":
    main()
