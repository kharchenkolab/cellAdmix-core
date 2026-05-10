#!/usr/bin/env python3
"""CLI batch runner for the Python cellAdmix bindings.

The script mirrors the R batch entrypoint for the core Xenium workflow:
construct a dataset, fit factors, score admixture, optionally correct molecules,
write summary tables, and optionally write a compact self-contained HTML report.
"""

from __future__ import annotations

import argparse
import base64
import html
import io
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch cellAdmix processing through the Python bindings.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python scripts/celladmix_batch.py --input data --output out \\
    --annotation annotations/annotation.csv.gz --annotation-col merged_annotation \\
    --score auto --report --threads 10

  python scripts/celladmix_batch.py --input data --output out \\
    --annotation annotations/annotation.csv.gz --annotation-col merged_annotation \\
    --score bridge --no-null --threads 10

Notes:
  Python v1 supports Xenium bundle inputs. Tabular input options are accepted
  for CLI parity with the R script, but currently fail with explicit messages
  until the Python tabular store binding is added.
""",
    )

    core = parser.add_argument_group("Core options")
    core.add_argument("--input", required=True, help="Xenium bundle or input source path")
    core.add_argument("--output", required=True, help="Output/cache directory")
    core.add_argument("--format", choices=("auto", "xenium", "tabular"), default="auto")
    core.add_argument("--threads", type=int, default=max(1, os.cpu_count() or 1))
    core.add_argument("--verbose", action="store_true")
    core.add_argument("--overwrite", action="store_true", help="Recompute fit even when cached")

    ann = parser.add_argument_group("Annotation")
    ann.add_argument("--annotation", help="Annotation CSV/CSV.GZ/Parquet")
    ann.add_argument("--annotation-col", help="Annotation label column")
    ann.add_argument("--cell-id-col", default="cell_id", help="Annotation cell id column")
    ann.add_argument("--annotation-name", default="manual", help="Recorded annotation name")
    ann.add_argument("--auto-annotate", action="store_true", help="Run clustering if no annotation is supplied")
    ann.add_argument("--annotation-report", action="store_true", help="Write annotation report")
    ann.add_argument("--cluster-name", default="cluster")
    ann.add_argument("--cluster-cells-max", type=int, default=-1, help="Maximum cells for auto-clustering; -1 uses all cells")
    ann.add_argument("--cluster-min-molecules", type=int, default=10)
    ann.add_argument("--cluster-min-genes", type=int, default=5)
    ann.add_argument("--cluster-resolution", type=float, default=1.0)
    ann.add_argument("--cluster-umap", type=parse_bool, default=True)
    ann.add_argument("--cluster-markers", type=parse_bool, default=True)
    ann.add_argument("--cluster-top-markers", type=int, default=5)
    ann.add_argument("--cluster-marker-min-fraction", type=float, default=0.05)
    ann.add_argument("--cluster-marker-min-logfc", type=float, default=0.0)
    ann.add_argument("--cluster-marker-max-panels", type=int, default=25)

    tab = parser.add_argument_group("Tabular input")
    tab.add_argument("--molecules", help="Molecule table path; defaults to --input for tabular")
    tab.add_argument("--cell-metadata", help="Cell metadata CSV/CSV.GZ/Parquet")
    tab.add_argument("--x-col", default="x")
    tab.add_argument("--y-col", default="y")
    tab.add_argument("--z-col", default="z")
    tab.add_argument("--gene-col", default="gene")
    tab.add_argument("--cell-col", default="cell")
    tab.add_argument("--qv-col")
    tab.add_argument("--cell-type-col")
    tab.add_argument("--cell-metadata-cell", default="cell")
    tab.add_argument("--cell-metadata-cell-type")
    tab.add_argument("--segmentation-mask")
    tab.add_argument("--sample-id-col")
    tab.add_argument("--fov-id-col")

    fit = parser.add_argument_group("Fit")
    fit.add_argument("--rank", default="auto", help="auto or integer factor rank")
    fit.add_argument("--rank-multiplier", type=float, default=1.2)
    fit.add_argument("--rank-cap", type=int, default=30)
    fit.add_argument("--nmf-variant", choices=("invsqrt_kl", "kl", "sqrt_kl", "ls_nmf"), default="invsqrt_kl")
    fit.add_argument("--nmf-init", choices=("auto", "random", "cluster"), default="auto")
    fit.add_argument("--nmf-runs", default="auto", help="auto or integer NMF restarts")
    fit.add_argument("--nmf-iterations", type=int)
    fit.add_argument("--molecule-scoring", default="gene_loadings")
    fit.add_argument("--seed", type=int, default=1)
    fit.add_argument("--run-id")
    fit.add_argument("--min-qv", type=float, default=-1.0)
    fit.add_argument("--keep-unassigned", action="store_true")
    fit.add_argument("--keep-non-gene", action="store_true", help="Keep Xenium control/codeword/non-gene features")
    fit.add_argument("--cell-filter", help="Comma-separated cell ids to keep")
    fit.add_argument("--gene-filter", help="Comma-separated genes to keep")

    score = parser.add_argument_group("Scoring and correction")
    score.add_argument("--score", choices=("auto", "membrane", "bridge", "coherence"), default="auto")
    score.add_argument("--score-name")
    score.add_argument("--p-thresh", type=float, default=0.1)
    score.add_argument("--adjust-p", action="store_true")
    score.add_argument("--targets", help="Comma-separated target cell types")
    score.add_argument("--no-correct", action="store_true")
    score.add_argument("--correction-name")
    score.add_argument("--max-cells-per-type-pair", type=int, default=400)
    score.add_argument("--candidate-pairs-per-type-pair", type=int, default=400)
    score.add_argument("--min-factor-molecules", type=int, default=5)
    score.add_argument("--min-pairs", type=int, default=5)
    score.add_argument("--no-null", action="store_true")
    score.add_argument("--null-iterations", type=int, default=3)

    report = parser.add_argument_group("Reports")
    report.add_argument("--report", action="store_true", help="Write compact HTML report")
    report.add_argument("--report-file", help="Report output path")
    report.add_argument("--top-genes", type=int, default=10)

    return parser


def parse_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).lower() in {"1", "true", "yes", "y", "on"}


def log(message: str) -> None:
    print(f"[batch {time.strftime('%H:%M:%S')}] {message}", flush=True)


def csv_list(value: str | None) -> list[str] | None:
    if value is None or value == "":
        return None
    return [x.strip() for x in value.split(",") if x.strip()]


def infer_format(source: Path, requested: str) -> str:
    if requested != "auto":
        return requested
    if source.is_dir() and (source / "experiment.xenium").exists():
        return "xenium"
    return "tabular"


def read_optional_annotation(args) -> pd.Series | None:
    import celladmix as ca

    if args.annotation:
        return ca.read_annotation(args.annotation, annotation_col=args.annotation_col, cell_id_col=args.cell_id_col)
    if args.cell_metadata and args.cell_metadata_cell_type:
        return ca.read_annotation(
            args.cell_metadata,
            annotation_col=args.cell_metadata_cell_type,
            cell_id_col=args.cell_metadata_cell,
        )
    return None


def _cluster_label(value, prefix: str) -> str:
    try:
        number = float(value)
        if np.isfinite(number) and number.is_integer():
            return f"{prefix}_{int(number)}"
    except (TypeError, ValueError):
        pass
    return f"{prefix}_{value}"


def auto_cluster_annotation(
    ds,
    args,
    output_dir: Path,
    *,
    compute_markers: bool,
) -> tuple[pd.Series, pd.DataFrame, pd.DataFrame]:
    """Cluster cell expression and return a cell-id-indexed annotation."""
    log("Auto-annotating cells by expression clustering")
    # The native clustering path returns a cell table; the CLI records cluster
    # labels as a reusable annotation file for downstream reruns.
    frame = ds.cell_state_umap(
        cells_max=None if args.cluster_cells_max < 0 else args.cluster_cells_max,
        min_molecules=args.cluster_min_molecules,
        min_genes=args.cluster_min_genes,
        cluster_resolution=args.cluster_resolution,
        compute_umap=args.cluster_umap,
        num_threads=args.threads,
        seed=args.seed,
    )
    labels = frame["cluster"].map(lambda x: _cluster_label(x, args.cluster_name))
    labels = pd.Series(labels.to_numpy(), index=frame["cell_id"].astype(str).to_numpy(), name=args.cluster_name)
    frame = frame.copy()
    frame["annotation"] = labels.reindex(frame["cell_id"].astype(str)).to_numpy()

    annotation_out = pd.DataFrame({"cell_id": labels.index, args.cluster_name: labels.to_numpy()})
    annotation_out.to_csv(output_dir / "annotation_clusters.csv", index=False)
    markers = compute_annotation_markers(ds.store_dir, labels, args) if compute_markers else pd.DataFrame()
    if not markers.empty:
        markers.to_csv(output_dir / "annotation_cluster_markers.csv", index=False)
    return labels, frame, markers


def compute_annotation_markers(store_dir: Path, annotation: pd.Series, args) -> pd.DataFrame:
    """Compute quick one-vs-rest marker summaries from sparse input-store counts."""
    # This intentionally uses input-store parquet files directly so marker
    # reporting works before a factor fit exists.
    cells = pd.read_parquet(store_dir / "cells.parquet", columns=["cell_idx", "cell_id"])
    cells["cell_id"] = cells["cell_id"].astype(str)
    cells["label"] = annotation.reindex(cells["cell_id"]).to_numpy()
    cells = cells[cells["label"].notna() & (cells["label"].astype(str) != "")].copy()
    if cells.empty:
        return pd.DataFrame()

    counts = pd.read_parquet(store_dir / "cell_gene_counts.parquet", columns=["cell_idx", "gene_idx", "count"])
    counts = counts.merge(cells[["cell_idx", "label"]], on="cell_idx", how="inner")
    if counts.empty:
        return pd.DataFrame()

    genes = pd.read_parquet(store_dir / "genes.parquet", columns=["gene_idx", "gene"])
    gene_names = dict(zip(genes["gene_idx"].astype(int), genes["gene"].astype(str)))
    gene_ids = genes["gene_idx"].astype(int).to_list()
    labels = sorted(cells["label"].astype(str).unique())

    counts["label"] = counts["label"].astype(str)
    counts["count"] = counts["count"].astype(float)
    sum_matrix = counts.groupby(["label", "gene_idx"], observed=True)["count"].sum().unstack(fill_value=0.0)
    detected_matrix = counts.groupby(["label", "gene_idx"], observed=True)["cell_idx"].nunique().unstack(fill_value=0)
    sum_matrix = sum_matrix.reindex(index=labels, columns=gene_ids, fill_value=0.0)
    detected_matrix = detected_matrix.reindex(index=labels, columns=gene_ids, fill_value=0)

    n_cells = cells.groupby("label", observed=True)["cell_idx"].nunique().reindex(labels).astype(float)
    label_totals = sum_matrix.sum(axis=1)
    all_gene_counts = sum_matrix.sum(axis=0)
    all_gene_detected = detected_matrix.sum(axis=0)
    total_counts = float(label_totals.sum())
    total_cells = float(n_cells.sum())
    n_genes = max(1, len(gene_ids))
    pseudocount = 0.5

    rows: list[dict] = []
    for label in labels:
        cluster_counts = sum_matrix.loc[label]
        rest_counts = all_gene_counts - cluster_counts
        cluster_total = float(label_totals.loc[label])
        rest_total = max(0.0, total_counts - cluster_total)
        cluster_n = max(1.0, float(n_cells.loc[label]))
        rest_n = max(1.0, total_cells - cluster_n)
        cluster_rate = (cluster_counts + pseudocount) / (cluster_total + pseudocount * n_genes)
        rest_rate = (rest_counts + pseudocount) / (rest_total + pseudocount * n_genes)
        log2fc = np.log2(cluster_rate / rest_rate)
        fraction_in = detected_matrix.loc[label].astype(float) / cluster_n
        fraction_rest = (all_gene_detected - detected_matrix.loc[label]).astype(float) / rest_n
        marker_frame = pd.DataFrame(
            {
                "cluster": label,
                "gene_idx": gene_ids,
                "gene": [gene_names[int(gene_idx)] for gene_idx in gene_ids],
                "log2_fold_change": log2fc.to_numpy(dtype=float),
                "mean_count": (cluster_counts / cluster_n).to_numpy(dtype=float),
                "fraction_in_cluster": fraction_in.to_numpy(dtype=float),
                "fraction_rest": fraction_rest.to_numpy(dtype=float),
                "cluster_count": cluster_counts.to_numpy(dtype=float),
                "rest_count": rest_counts.to_numpy(dtype=float),
            }
        )
        marker_frame = marker_frame[
            (marker_frame["fraction_in_cluster"] >= args.cluster_marker_min_fraction)
            & (marker_frame["log2_fold_change"] >= args.cluster_marker_min_logfc)
        ]
        marker_frame = marker_frame.sort_values(
            ["log2_fold_change", "fraction_in_cluster", "mean_count"],
            ascending=[False, False, False],
        ).head(max(1, int(args.cluster_top_markers)))
        marker_frame.insert(1, "rank", range(1, len(marker_frame) + 1))
        rows.extend(marker_frame.to_dict("records"))

    return pd.DataFrame(rows)


def membrane_available(fit) -> bool:
    try:
        image = fit.stain("membrane")
    except Exception:
        return False
    return Path(image["image_path"]).exists()


def resolve_score_method(requested: str, fit) -> str:
    if requested == "auto":
        return "membrane" if membrane_available(fit) else "bridge"
    if requested == "membrane" and not membrane_available(fit):
        raise RuntimeError("--score membrane was requested, but no membrane stain could be discovered")
    if requested == "coherence":
        raise NotImplementedError("Python batch v1 does not expose coherence scoring yet")
    return requested


def factor_top_genes(fit, n: int) -> pd.DataFrame:
    loadings = fit.factor_loadings()
    rows = []
    for factor in loadings.columns:
        values = loadings[factor].astype(float).sort_values(ascending=False).head(max(1, int(n)))
        denom = float(loadings[factor].astype(float).sum())
        for rank, (gene, loading) in enumerate(values.items(), start=1):
            rows.append(
                {
                    "factor": factor,
                    "rank": rank,
                    "gene": gene,
                    "loading": float(loading),
                    "loading_fraction": float(loading / denom) if denom > 0 else float("nan"),
                }
            )
    return pd.DataFrame(rows)


def fit_summary(
    fit,
    *,
    args,
    score_method: str | None = None,
    score_name: str | None = None,
    n_rules: int | None = None,
    corrected: bool = False,
) -> dict:
    options = fit.manifest.get("pipeline_options", {}) or {}
    diagnostics = fit.manifest.get("nmf_diagnostics", {}) or {}
    return {
        "input": str(Path(args.input).resolve()),
        "output_dir": str(Path(args.output).resolve()),
        "format": args.format,
        "annotation_name": args.annotation_name,
        "rank": int(fit.manifest.get("n_factors", options.get("rank", 0)) or 0),
        "n_cells": int(fit.manifest.get("n_cells", 0) or 0),
        "n_transcripts": int(fit.manifest.get("n_transcripts", 0) or 0),
        "n_training_rows": int(fit.manifest.get("n_training_rows", 0) or 0),
        "nmf_variant": options.get("nmf_variant", args.nmf_variant),
        "nmf_n_runs": int(options.get("nmf_n_runs", 0) or 0),
        "nmf_final_objective": diagnostics.get("final_objective"),
        "score_method": score_method,
        "score_name": score_name,
        "n_rules": n_rules,
        "corrected": bool(corrected),
        "threads": int(args.threads),
    }


def write_outputs(output_dir: Path, fit, score, rules: pd.DataFrame, correction, top_genes: pd.DataFrame, summary: dict) -> None:
    """Write stable CSV/JSON outputs consumed by batch workflows."""
    pd.DataFrame([summary]).to_csv(output_dir / "batch_summary.csv", index=False)
    (output_dir / "batch_summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n", encoding="utf-8")
    pd.DataFrame([fit_summary_like_row(summary)]).to_csv(output_dir / "fit_summary.csv", index=False)
    top_genes.to_csv(output_dir / "factor_top_genes.csv", index=False)
    score.summary.to_csv(output_dir / "score_summary.csv", index=False)
    rules.to_csv(output_dir / "score_rules.csv", index=False)
    if correction is not None:
        correction.summary().to_csv(output_dir / "correction_summary.csv", index=False)


def fit_summary_like_row(summary: dict) -> dict:
    keys = [
        "rank",
        "n_cells",
        "n_transcripts",
        "n_training_rows",
        "nmf_variant",
        "nmf_n_runs",
        "nmf_final_objective",
    ]
    return {key: summary.get(key) for key in keys}


def render_report(path: Path, fit, score, rules: pd.DataFrame, correction, top_genes: pd.DataFrame, *, p_thresh: float, top_genes_n: int) -> None:
    """Render a compact self-contained HTML report with embedded figures."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sections = [
        "<h1>cellAdmix Batch Report</h1>",
        "<h2>Summary</h2>",
        table_html(pd.DataFrame([fit_summary_like_row(fit_summary_from_manifest(fit))])),
        "<h2>Score Rules</h2>",
        table_html(rules[["factor", "source_cell_type", "target_cell_type", "p_value"]].head(20) if not rules.empty else rules),
        "<h2>Correction Summary</h2>",
        table_html(correction.summary() if correction is not None else pd.DataFrame({"note": ["correction skipped"]})),
        "<h2>Top Factor Genes</h2>",
        table_html(top_genes.head(max(20, int(fit.manifest.get("n_factors", 1)) * top_genes_n))),
    ]

    figures = [
        ("Factor Loadings", fit.plot_loadings(n_genes=top_genes_n)),
        ("NMF Stability", fit.plot_stability()),
        ("Score Heatmap", score.plot_heatmap(p_thresh=p_thresh)),
        ("Score Pair Plots", score.plot_pairs(p_thresh=p_thresh)),
    ]
    if correction is not None:
        figures.append(("Molecules Removed", correction.plot_removed_molecules()))

    for title, obj in figures:
        fig = obj.figure if hasattr(obj, "figure") else obj
        sections.append(f"<h2>{html.escape(title)}</h2>")
        sections.append(f'<img src="data:image/png;base64,{fig_to_base64(fig)}" />')
        plt.close(fig)

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "<!doctype html><html><head><meta charset='utf-8'>",
                "<style>body{font-family:sans-serif;max-width:1200px;margin:24px auto;}"
                "img{max-width:100%;height:auto;} table{border-collapse:collapse;font-size:13px;}"
                "th,td{border:1px solid #ddd;padding:4px 6px;} th{background:#f5f5f5;}</style>",
                "</head><body>",
                *sections,
                "</body></html>",
            ]
        ),
        encoding="utf-8",
    )


def render_annotation_report(path: Path, frame: pd.DataFrame, markers: pd.DataFrame, *, label_col: str = "annotation") -> None:
    """Render the optional auto-annotation QC report."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    label = frame[label_col].astype(str)
    sizes = label.value_counts().rename_axis("annotation").reset_index(name="n_cells")
    sections = [
        "<h1>cellAdmix Annotation Report</h1>",
        "<h2>Cell Counts</h2>",
        table_html(sizes),
    ]
    if markers is not None and not markers.empty:
        sections.extend(
            [
                "<h2>Top Marker Genes</h2>",
                table_html(markers[["cluster", "rank", "gene", "log2_fold_change", "fraction_in_cluster"]].head(100)),
            ]
        )

    figures = [
        ("Cell-State UMAP", plot_annotation_scatter(frame, "umap_1", "umap_2", label_col, "Cell-state UMAP")),
        ("Spatial Labels", plot_annotation_scatter(frame, "x", "y", label_col, "Spatial labels")),
    ]
    for title, fig in figures:
        sections.append(f"<h2>{html.escape(title)}</h2>")
        sections.append(f'<img src="data:image/png;base64,{fig_to_base64(fig)}" />')
        plt.close(fig)

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "<!doctype html><html><head><meta charset='utf-8'>",
                "<style>body{font-family:sans-serif;max-width:1200px;margin:24px auto;}"
                "img{max-width:100%;height:auto;} table{border-collapse:collapse;font-size:13px;}"
                "th,td{border:1px solid #ddd;padding:4px 6px;} th{background:#f5f5f5;}</style>",
                "</head><body>",
                *sections,
                "</body></html>",
            ]
        ),
        encoding="utf-8",
    )


def plot_annotation_scatter(frame: pd.DataFrame, x_col: str, y_col: str, label_col: str, title: str):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    if x_col not in frame.columns or y_col not in frame.columns:
        ax.text(0.5, 0.5, f"{x_col}/{y_col} unavailable", ha="center", va="center")
        ax.set_axis_off()
        return fig
    plot_frame = frame[[x_col, y_col, label_col]].copy()
    plot_frame = plot_frame[np.isfinite(plot_frame[x_col]) & np.isfinite(plot_frame[y_col])]
    if plot_frame.empty:
        ax.text(0.5, 0.5, "no cells to plot", ha="center", va="center")
        ax.set_axis_off()
        return fig
    labels = sorted(plot_frame[label_col].astype(str).unique())
    cmap = plt.get_cmap("tab20", max(1, len(labels)))
    for i, label in enumerate(labels):
        sub = plot_frame[plot_frame[label_col].astype(str) == label]
        ax.scatter(sub[x_col], sub[y_col], s=2, alpha=0.75, color=cmap(i), linewidths=0, label=label)
    ax.set_title(title)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_aspect("equal", adjustable="box")
    if len(labels) <= 30:
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), markerscale=4, frameon=False, fontsize=7)
    return fig


def fit_summary_from_manifest(fit) -> dict:
    options = fit.manifest.get("pipeline_options", {}) or {}
    diagnostics = fit.manifest.get("nmf_diagnostics", {}) or {}
    return {
        "rank": fit.manifest.get("n_factors"),
        "n_cells": fit.manifest.get("n_cells"),
        "n_transcripts": fit.manifest.get("n_transcripts"),
        "n_training_rows": fit.manifest.get("n_training_rows"),
        "nmf_variant": options.get("nmf_variant"),
        "nmf_n_runs": options.get("nmf_n_runs"),
        "nmf_final_objective": diagnostics.get("final_objective"),
    }


def table_html(frame: pd.DataFrame) -> str:
    if frame is None or frame.empty:
        frame = pd.DataFrame({"note": ["empty"]})
    return frame.to_html(index=False, escape=True, border=0)


def fig_to_base64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    return base64.b64encode(buf.getvalue()).decode("ascii")


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    source = Path(args.input)
    fmt = infer_format(source, args.format)
    args.format = fmt
    if fmt != "xenium":
        raise NotImplementedError(
            "Python batch v1 currently supports Xenium bundle inputs. "
            "Use scripts/celladmix_batch.R for tabular inputs until the Python tabular store binding is added."
        )

    import celladmix as ca

    log("Constructing dataset")
    annotation = read_optional_annotation(args)
    ds = ca.CellAdmix(
        source,
        output_dir=output_dir,
        format=fmt,
        annotation=annotation,
        annotation_col=args.annotation_col,
        cell_id_col=args.cell_id_col,
        num_threads=args.threads,
        cell_filter=csv_list(args.cell_filter),
        gene_filter=csv_list(args.gene_filter),
        min_qv=args.min_qv,
        keep_unassigned=args.keep_unassigned,
        keep_non_gene=args.keep_non_gene,
    )
    print(ds)

    annotation_frame = None
    annotation_markers = pd.DataFrame()
    if annotation is None:
        if not args.auto_annotate:
            raise RuntimeError(
                "No annotation was supplied or inferred. Provide --annotation/--annotation-col "
                "or pass --auto-annotate to create cluster labels."
            )
        render_cluster_report = args.annotation_report or args.report
        annotation, annotation_frame, annotation_markers = auto_cluster_annotation(
            ds,
            args,
            output_dir,
            compute_markers=bool(render_cluster_report and args.cluster_markers),
        )
        ds.set_annotation(annotation)
        args.annotation_name = args.cluster_name
    elif args.annotation_report:
        log("Computing cell-state embedding for annotation report")
        annotation_frame = ds.cell_state_umap(
            cells_max=None if args.cluster_cells_max < 0 else args.cluster_cells_max,
            min_molecules=args.cluster_min_molecules,
            min_genes=args.cluster_min_genes,
            cluster_resolution=args.cluster_resolution,
            compute_umap=args.cluster_umap,
            num_threads=args.threads,
            seed=args.seed,
        )
        annotation_frame = annotation_frame.copy()
        annotation_frame["annotation"] = annotation.reindex(annotation_frame["cell_id"].astype(str)).to_numpy()
        annotation_markers = compute_annotation_markers(ds.store_dir, annotation, args) if args.cluster_markers else pd.DataFrame()
        if not annotation_markers.empty:
            annotation_markers.to_csv(output_dir / "annotation_markers.csv", index=False)

    if annotation_frame is not None and (args.annotation_report or (args.report and args.auto_annotate)):
        annotation_report_file = output_dir / "annotation_report.html"
        log(f"Writing annotation report: {annotation_report_file}")
        render_annotation_report(annotation_report_file, annotation_frame, annotation_markers)

    rank = None if args.rank == "auto" else int(args.rank)
    nmf_n_runs = None if args.nmf_runs == "auto" else int(args.nmf_runs)
    fit_kwargs = {}
    if args.nmf_iterations is not None:
        fit_kwargs["nmf_iterations"] = args.nmf_iterations

    log("Fitting factors")
    fit = ds.fit(
        rank=rank,
        rank_multiplier=args.rank_multiplier,
        rank_cap=args.rank_cap,
        nmf_variant=args.nmf_variant,
        nmf_init=args.nmf_init,
        nmf_n_runs=nmf_n_runs,
        molecule_scoring=args.molecule_scoring,
        run_id=args.run_id,
        overwrite=args.overwrite,
        num_threads=args.threads,
        verbose=args.verbose,
        seed=args.seed,
        **fit_kwargs,
    )
    print(fit)

    score_method = resolve_score_method(args.score, fit)
    score_name = args.score_name or score_method
    log(f"Scoring with method: {score_method}")
    score_kwargs = dict(
        num_threads=args.threads,
        min_factor_molecules=args.min_factor_molecules,
        min_pairs=args.min_pairs,
        max_cells_per_type_pair=args.max_cells_per_type_pair,
        verbose=args.verbose,
    )
    # Candidate-pair caps are shared by membrane and bridge scoring; null
    # options only apply to the bridge scorer.
    if score_method in {"membrane", "bridge"}:
        score_kwargs["candidate_pairs_per_type_pair"] = args.candidate_pairs_per_type_pair
    if score_method == "bridge":
        score_kwargs["compute_null"] = not args.no_null
        score_kwargs["null_iterations"] = args.null_iterations
    score = fit.score_membrane(**score_kwargs) if score_method == "membrane" else fit.score_bridge(**score_kwargs)

    targets = csv_list(args.targets)
    rules = score.rules(p_thresh=args.p_thresh, adjust_p=args.adjust_p, targets=targets)
    correction = None
    if not args.no_correct:
        if rules.empty:
            log("No correction rules were called; skipping correction")
        else:
            correction_name = args.correction_name or f"{score_name}_clean"
            log("Applying correction rules")
            correction = score.correct(rules=rules, name=correction_name)
            print(correction.summary())

    top_genes = factor_top_genes(fit, args.top_genes)
    summary = fit_summary(
        fit,
        args=args,
        score_method=score_method,
        score_name=score_name,
        n_rules=len(rules),
        corrected=correction is not None,
    )
    write_outputs(output_dir, fit, score, rules, correction, top_genes, summary)

    if args.report:
        report_file = Path(args.report_file) if args.report_file else output_dir / "report_minimal.html"
        if not report_file.is_absolute():
            report_file = output_dir / report_file
        log(f"Writing report: {report_file}")
        render_report(report_file, fit, score, rules, correction, top_genes, p_thresh=args.p_thresh, top_genes_n=args.top_genes)

    log("Batch run complete")
    print(f"Output directory: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise SystemExit(1)
