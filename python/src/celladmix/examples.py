"""Example-cell selection and rendering helpers for Python notebooks."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def factor_id(value) -> int | None:
    """Return a one-based factor id from integers or labels such as ``F3``."""
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return None
    if isinstance(value, str):
        value = value.strip().lstrip("Ff").lstrip("_")
    try:
        out = int(value)
    except (TypeError, ValueError):
        return None
    return out if out > 0 else None


def factor_label(value) -> str:
    """Return a display label for a one-based factor id."""
    out = factor_id(value)
    return f"F{out}" if out is not None else "factor"


def discover_cell_boundaries(source, boundary_path=None) -> Path | None:
    """Find the standard Xenium cell-boundary table for a source bundle."""
    if boundary_path is not None:
        path = Path(boundary_path)
        if not path.exists():
            raise FileNotFoundError(f"Cell boundary path does not exist: {path}")
        return path
    source = Path(source)
    source_dir = source if source.is_dir() else source.parent
    for name in ("cell_boundaries.parquet", "cell_boundaries.csv.gz", "cell_boundaries.csv"):
        path = source_dir / name
        if path.exists():
            return path
    return None


def read_cell_boundaries(path, *, cells: Iterable[str] | None = None, bbox=None) -> pd.DataFrame:
    """Read Xenium-style cell boundary vertices."""
    path = Path(path)
    if path.suffix == ".parquet":
        frame = pd.read_parquet(path)
    else:
        frame = pd.read_csv(path)
    cols = {
        "cell": _first_present(frame, ("cell_id", "cell", "cellID")),
        "x": _first_present(frame, ("vertex_x", "x", "X")),
        "y": _first_present(frame, ("vertex_y", "y", "Y")),
    }
    out = pd.DataFrame(
        {
            "cell_id": frame[cols["cell"]].astype(str),
            "x": pd.to_numeric(frame[cols["x"]], errors="coerce"),
            "y": pd.to_numeric(frame[cols["y"]], errors="coerce"),
        }
    )
    out = out[out["cell_id"].ne("") & np.isfinite(out["x"]) & np.isfinite(out["y"])]
    if cells is not None:
        keep = set(map(str, cells))
        out = out[out["cell_id"].isin(keep)]
    if bbox is not None and not out.empty:
        xmin, xmax, ymin, ymax = bbox
        hit = out["x"].between(xmin, xmax) & out["y"].between(ymin, ymax)
        out = out[out["cell_id"].isin(out.loc[hit, "cell_id"].unique())]
    return out.reset_index(drop=True)


def read_stain_crop(image: dict, bbox, *, max_pixels: int = 512) -> dict:
    """Read a downsampled physical-coordinate crop from a tiled Xenium OME-TIFF.

    The file is read with OME interpretation disabled: multimodal Xenium
    bundles split the focus channels across files whose shared OME-XML makes
    tifffile aggregate them into one multi-channel series, in which case every
    file would otherwise yield channel 0 (DAPI) regardless of the requested
    stain. Reading file-local pages keeps each focus file's own channel.
    """
    import tifffile
    import zarr

    xmin, xmax, ymin, ymax = map(float, bbox)
    pixel_size = float(image.get("pixel_size", 1.0))
    x_offset = float(image.get("x_offset", 0.0))
    y_offset = float(image.get("y_offset", 0.0))
    channel = image.get("channel")
    image_path = str(image["image_path"])
    with tifffile.TiffFile(image_path, is_ome=False) as tif:
        series = tif.series[0]
        levels = sorted(series.levels, key=lambda l: l.shape[-1], reverse=True)
        full_width = levels[0].shape[-1]
        width_px = max(1.0, (xmax - xmin) / pixel_size)
        height_px = max(1.0, (ymax - ymin) / pixel_size)
        level = levels[0]
        # Use the coarsest pyramid level that still satisfies the requested
        # pixel budget; this keeps notebook example rendering lightweight.
        for candidate in levels:
            scale = full_width / candidate.shape[-1]
            if max(width_px / scale, height_px / scale) <= max_pixels:
                level = candidate
                break
        scale = full_width / level.shape[-1]
        store = level.aszarr()
        try:
            data = zarr.open(store, mode="r")
            if isinstance(data, zarr.Group):
                data = data[next(iter(sorted(data.array_keys())))]
            px0 = int(np.floor((xmin - x_offset) / pixel_size / scale))
            px1 = int(np.ceil((xmax - x_offset) / pixel_size / scale))
            py0 = int(np.floor((ymin - y_offset) / pixel_size / scale))
            py1 = int(np.ceil((ymax - y_offset) / pixel_size / scale))
            px0 = max(0, min(px0, level.shape[-1] - 1))
            py0 = max(0, min(py0, level.shape[-2] - 1))
            px1 = max(px0 + 1, min(px1, level.shape[-1]))
            py1 = max(py0 + 1, min(py1, level.shape[-2]))
            if len(level.shape) == 3:
                values = np.asarray(
                    data[int(channel or 0), py0:py1, px0:px1], dtype=float)
            else:
                values = np.asarray(data[py0:py1, px0:px1], dtype=float)
        finally:
            store.close()
    return {
        "values": values,
        "bbox": (xmin, xmax, ymin, ymax),
        "pixel_size": pixel_size * scale,
        "stain": image.get("stain", "stain"),
    }


def compose_stain_background(crops: dict[str, dict] | None) -> np.ndarray | None:
    """Compose stain crops as a light RGB background."""
    if not crops:
        return None
    crops = {k: v for k, v in crops.items() if v is not None}
    if not crops:
        return None
    first = next(iter(crops.values()))["values"]
    rgb = np.ones((first.shape[0], first.shape[1], 3), dtype=float)
    colors = {
        "dapi": np.array([0.84, 0.73, 0.96]),
        "nucleus": np.array([0.84, 0.73, 0.96]),
        "membrane": np.array([0.75, 0.93, 0.74]),
        "polya": np.array([0.95, 0.82, 0.56]),
        "poly_a": np.array([0.95, 0.82, 0.56]),
    }
    for name, crop in crops.items():
        channel = _normalize_channel(crop["values"])
        color = colors.get(name.lower(), np.array([0.35, 0.35, 0.35]))
        rgb *= 1.0 - channel[:, :, None] * (1.0 - color[None, None, :])
    return np.clip(rgb, 0, 1)


_EXAMPLE_COLUMNS = [
    "target_cell", "target_cell_type", "top_admix_factor", "top_admix_source",
    "top_admix_score", "top_factor_molecules", "top_factor_fraction",
    "top_evidence", "n_candidate_factors", "total_factor_molecules",
    "transcript_count", "dominant_factor",
]


def _normalize_cells(cells) -> list[str]:
    if isinstance(cells, str):
        return [cells]
    return [str(c) for c in cells]


def _requested_cell_examples(score, cells, examples, cell_data) -> pd.DataFrame:
    """Assemble the example table for explicitly requested cells.

    Score-derived rows are kept; requested cells without score evidence come
    back as NA-evidence rows, and the requested order is preserved.
    """
    have = set() if examples is None or examples.empty else set(examples["target_cell"].astype(str))
    missing = [c for c in cells if c not in have]
    if missing:
        try:
            annotation = score.fit.dataset.annotation
        except AttributeError:
            annotation = None
        rows = pd.DataFrame({"target_cell": missing})
        rows["target_cell_type"] = (
            [annotation.get(c, np.nan) for c in missing] if annotation is not None else np.nan
        )
        for column in _EXAMPLE_COLUMNS[2:]:
            rows[column] = np.nan
        if cell_data is not None and not cell_data.empty and "cell_id" in cell_data.columns:
            lookup = cell_data.set_index(cell_data["cell_id"].astype(str))
            for column in ("transcript_count", "dominant_factor"):
                if column in lookup.columns:
                    rows[column] = [
                        lookup[column].get(c, np.nan) for c in missing
                    ]
        if examples is not None and not examples.empty:
            examples = pd.concat([examples, rows[examples.columns]], ignore_index=True)
        else:
            examples = rows
    examples = examples.set_index(examples["target_cell"].astype(str)).reindex(cells)
    return examples.dropna(subset=["target_cell"]).reset_index(drop=True)


def select_example_cells(
    score,
    *,
    rules: pd.DataFrame | None = None,
    score_annotation: dict | None = None,
    cell_data: pd.DataFrame | None = None,
    targets=None,
    n_per_target: int = 4,
    p_thresh: float = 0.1,
    adjust_p: bool = False,
    use_rules: bool = True,
    min_molecules: int = 50,
    min_factor_molecules: int = 3,
    cells=None,
) -> pd.DataFrame:
    """Select target cells with strong non-native factor evidence.

    When ``cells`` is supplied, selection is restricted to those cell IDs, the
    evidence filters are relaxed so every requested cell is returned (with NA
    factor columns when the score has no evidence for it), and the cells come
    back in the requested order.
    """
    if cells is not None:
        cells = _normalize_cells(cells)
        use_rules = False
        min_molecules = 0
        min_factor_molecules = 0
        targets = None
    pairs = score.pairs.copy()
    if cells is not None and not pairs.empty:
        pairs = pairs[pairs["target_cell"].astype(str).isin(cells)]
    if cell_data is None:
        cell_data = score.fit.cell_factors()
    if pairs.empty:
        if cells is None:
            return pd.DataFrame()
        return _requested_cell_examples(score, cells, None, cell_data)
    score_annotation = score_annotation or score.annotation(p_thresh=p_thresh, adjust_p=adjust_p)
    if use_rules and rules is None:
        rules = score.rules(p_thresh=p_thresh, adjust_p=adjust_p, targets=targets)
    if targets is not None:
        targets = list(map(str, targets))
        pairs = pairs[pairs["target_cell_type"].astype(str).isin(targets)]
    if "used_in_summary" in pairs.columns and cells is None:
        pairs = pairs[pairs["used_in_summary"].astype(bool)]
    pairs = pairs[np.isfinite(pd.to_numeric(pairs["mean_score"], errors="coerce"))]
    if use_rules and rules is not None and not rules.empty:
        keys = set(zip(rules["target_cell_type"].astype(str), rules["factor"].astype(int)))
        pairs = pairs[
            [
                (str(row.target_cell_type), int(row.factor)) in keys
                for row in pairs.itertuples(index=False)
            ]
        ]
    if "factor_count" in pairs.columns:
        pairs = pairs[pairs["factor_count"].astype(float) >= min_factor_molecules]
    if pairs.empty:
        if cells is None:
            return pd.DataFrame()
        return _requested_cell_examples(score, cells, None, cell_data)

    source_calls = score_annotation.get("source_calls", {})
    per_factor = (
        pairs.groupby(["target_cell", "target_cell_type", "factor"], as_index=False)
        .agg(
            max_score=("mean_score", "max"),
            factor_count=("factor_count", "max"),
            scored_molecules=("scored_molecules", "max")
            if "scored_molecules" in pairs.columns
            else ("factor_count", "max"),
        )
    )
    per_factor["source_cell_type"] = per_factor["factor"].astype(int).map(source_calls)
    if cell_data is not None and not cell_data.empty:
        optional = [
            col
            for col in ("cell_id", "transcript_count", "dominant_factor")
            if col in cell_data.columns
        ]
        per_factor = per_factor.merge(
            cell_data[optional],
            left_on="target_cell",
            right_on="cell_id",
            how="left",
            sort=False,
        )
    per_factor = per_factor[
        per_factor["source_cell_type"].isna()
        | (per_factor["source_cell_type"].astype(str) != per_factor["target_cell_type"].astype(str))
    ]
    if "dominant_factor" in per_factor.columns:
        per_factor = per_factor[
            per_factor["dominant_factor"].isna()
            | (per_factor["factor"].astype(int) != per_factor["dominant_factor"].astype(int))
        ]
    if "transcript_count" in per_factor.columns:
        per_factor = per_factor[
            per_factor["transcript_count"].isna()
            | (per_factor["transcript_count"].astype(float) >= min_molecules)
        ]
    if cells is None:
        per_factor = per_factor[np.isfinite(per_factor["max_score"]) & (per_factor["max_score"] > 0)]
    if per_factor.empty:
        if cells is None:
            return pd.DataFrame()
        return _requested_cell_examples(score, cells, None, cell_data)

    rows = []
    for _, group in per_factor.groupby("target_cell", sort=False):
        group = group.sort_values(["factor_count", "max_score"], ascending=False)
        top = group.iloc[0]
        transcript_count = float(top.get("transcript_count", np.nan))
        factor_count = float(top.get("factor_count", np.nan))
        rows.append(
            {
                "target_cell": top["target_cell"],
                "target_cell_type": top["target_cell_type"],
                "top_admix_factor": int(top["factor"]),
                "top_admix_source": top.get("source_cell_type", np.nan),
                "top_admix_score": float(top["max_score"]),
                "top_factor_molecules": factor_count,
                "top_factor_fraction": factor_count / transcript_count
                if np.isfinite(transcript_count) and transcript_count > 0
                else np.nan,
                "top_evidence": factor_count * max(0.0, float(top["max_score"]))
                if np.isfinite(factor_count)
                else np.nan,
                "n_candidate_factors": int(len(group)),
                "total_factor_molecules": float(group["factor_count"].sum()),
                "transcript_count": transcript_count,
                "dominant_factor": top.get("dominant_factor", np.nan),
            }
        )
    out = pd.DataFrame(rows).sort_values(
        ["top_evidence", "top_factor_molecules", "total_factor_molecules", "top_admix_score"],
        ascending=False,
    )
    if cells is not None:
        return _requested_cell_examples(score, cells, out, cell_data)
    order = targets or out["target_cell_type"].drop_duplicates().tolist()
    pieces = [out[out["target_cell_type"].astype(str) == str(t)].head(n_per_target) for t in order]
    return pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()


def resolve_example_stains(fit, stains) -> dict:
    """Resolve the stains argument of example plotting into stain descriptors.

    ``"auto"`` (the default everywhere) discovers the available stain images
    from the fit's Xenium bundle and is empty for non-Xenium sources; ``None``
    disables backgrounds; a single descriptor or a name-keyed mapping of
    descriptors is passed through.
    """
    if isinstance(stains, str):
        if stains != "auto":
            raise ValueError(f"stains must be 'auto', None, or descriptors; got {stains!r}")
        from .io import discover_xenium_stain_images

        source = getattr(getattr(fit, "dataset", None), "source", None)
        return discover_xenium_stain_images(source)
    if not stains:
        return {}
    if isinstance(stains, dict) and "image_path" in stains:
        return {str(stains.get("stain", "stain")): stains}
    return dict(stains)


def resolve_example_cell_types(fit, cell_types):
    """Resolve the cell_types argument into a cell_id -> label mapping."""
    if isinstance(cell_types, str):
        if cell_types != "auto":
            raise ValueError(f"cell_types must be 'auto', None, or a mapping; got {cell_types!r}")
        annotation = getattr(getattr(fit, "dataset", None), "annotation", None)
        return None if annotation is None else annotation
    if cell_types is None:
        return None
    if isinstance(cell_types, pd.DataFrame):
        type_col = next(
            (c for c in ("cell_type", "merged_annotation", "cluster_label") if c in cell_types.columns),
            None,
        )
        if type_col is None or "cell_id" not in cell_types.columns:
            raise ValueError("cell_types data frame must contain cell_id and cell_type columns")
        return pd.Series(
            cell_types[type_col].astype(str).to_numpy(),
            index=cell_types["cell_id"].astype(str),
        )
    return cell_types


def prepare_cell_example(
    fit,
    example,
    *,
    score_annotation: dict | None = None,
    cell_data: pd.DataFrame | None = None,
    boundaries=None,
    stains="auto",
    cell_types="auto",
    markers=None,
    padding: float = 5,
    min_side: float = 24,
    max_pixels: int = 512,
) -> dict:
    """Collect molecules, contours, and optional stain crops for one example.

    Stain backgrounds and cell boundaries are discovered automatically for
    Xenium-backed fits (``stains=None`` disables backgrounds); ``cell_types``
    labels neighbor contours for cell-type coloring and defaults to the
    dataset's annotation.
    """
    stains = resolve_example_stains(fit, stains)
    if isinstance(example, pd.DataFrame):
        row = example.iloc[0]
    elif isinstance(example, pd.Series):
        row = example
    else:
        row = pd.Series(dict(example))
    cell_data = fit.cell_factors() if cell_data is None else cell_data
    target_cell = str(row["target_cell"])
    target_type = str(row["target_cell_type"])
    center = cell_data[cell_data["cell_id"].astype(str) == target_cell]
    if center.empty:
        raise ValueError(f"Target cell is not present in cell_data: {target_cell}")

    boundary_path = None
    if boundaries is None:
        boundary_path = discover_cell_boundaries(fit.dataset.source)
        boundaries = (
            read_cell_boundaries(boundary_path, cells=[target_cell])
            if boundary_path is not None
            else pd.DataFrame()
        )
    elif isinstance(boundaries, (str, Path)):
        boundaries = read_cell_boundaries(boundaries)
    boundary = boundaries[boundaries["cell_id"].astype(str) == target_cell] if not boundaries.empty else pd.DataFrame()
    if not boundary.empty:
        bbox = square_bbox((boundary["x"].min(), boundary["x"].max()), (boundary["y"].min(), boundary["y"].max()), padding=padding, min_side=min_side)
    else:
        bbox = square_bbox((float(center["x"].iloc[0]), float(center["x"].iloc[0])), (float(center["y"].iloc[0]), float(center["y"].iloc[0])), padding=padding, min_side=min_side)
    if boundary_path is not None:
        # The discovery path read only the target cell above; add the
        # neighboring cell polygons inside the plotting window.
        neighbors = read_cell_boundaries(boundary_path, bbox=bbox)
        # Deduplicate whole cells, not rows: each ring's closing vertex repeats
        # its first vertex, and row-level deduplication would unclose the ring.
        neighbors = neighbors[~neighbors["cell_id"].isin(set(boundaries["cell_id"]))]
        boundaries = pd.concat([boundaries, neighbors], ignore_index=True)
    if not boundaries.empty:
        boundaries = boundaries[
            boundaries["cell_id"].isin(
                boundaries.loc[
                    boundaries["x"].between(bbox[0], bbox[1])
                    & boundaries["y"].between(bbox[2], bbox[3]),
                    "cell_id",
                ].unique()
            )
            | (boundaries["cell_id"].astype(str) == target_cell)
        ]

    stain_crops = {}
    for name, stain in (stains or {}).items():
        stain_crops[name] = read_stain_crop(stain, bbox, max_pixels=max_pixels)
    background = compose_stain_background(stain_crops)
    molecules = fit.region(bbox=bbox)
    molecules["inside_target"] = molecules["cell_id"].astype(str) == target_cell

    source_calls = (score_annotation or {}).get("source_calls", {})
    native = [int(f) for f, source in source_calls.items() if str(source) == target_type]
    if "dominant_factor" in center.columns and pd.notna(center["dominant_factor"].iloc[0]):
        native.append(int(center["dominant_factor"].iloc[0]))
    top_factor = factor_id(row.get("top_admix_factor"))
    # Native factors are shown separately from the top non-native factor that
    # caused this cell to be selected.
    native = sorted(set(native) - ({top_factor} if top_factor else set()))
    native_name = f"native ({','.join(factor_label(x) for x in native)})" if native else "native factors"
    top_name = f"top admixture ({factor_label(top_factor)})" if top_factor else "top admixture"
    other_name = "other admixture factors"
    molecules["role"] = other_name
    if native:
        molecules.loc[molecules["factor"].isin(native), "role"] = native_name
    if top_factor:
        molecules.loc[molecules["factor"].astype("Int64") == top_factor, "role"] = top_name
    marker_set = set(markers or [])
    molecules["is_marker"] = molecules["inside_target"] & molecules["gene"].isin(marker_set)
    contours = boundaries.copy() if boundaries is not None else pd.DataFrame()
    type_map = resolve_example_cell_types(fit, cell_types)
    if not contours.empty:
        if type_map is not None:
            contours["cell_type"] = contours["cell_id"].astype(str).map(
                lambda c: type_map.get(c, np.nan))
        else:
            contours["cell_type"] = np.nan
    return {
        "example": row.to_dict(),
        "bbox": bbox,
        "background": background,
        "molecules": molecules,
        "contours": contours,
        "target_cell": target_cell,
        "target_cell_type": target_type,
        "role_levels": [native_name, top_name, other_name],
    }


def plot_cell_example(
    example: dict,
    *,
    ax=None,
    outside_size: float = 8.0,
    inside_size: float = 18.0,
    marker_size_multiplier: float = 1.1,
    non_marker_size_multiplier: float = 0.9,
    contour_color: str = "#e85d04",
    color_cell_types: bool = True,
    cell_type_linewidth: float = 0.9,
):
    """Render one prepared example-cell overlay."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(4, 4))
    bbox = example["bbox"]
    if example.get("background") is not None:
        ax.imshow(example["background"], extent=(bbox[0], bbox[1], bbox[3], bbox[2]), origin="upper")
    contours = example.get("contours", pd.DataFrame())
    type_handles = []
    if not contours.empty:
        nearby = contours[contours["cell_id"].astype(str) != example["target_cell"]]
        target = contours[contours["cell_id"].astype(str) == example["target_cell"]]
        typed = nearby.iloc[0:0]
        if color_cell_types and "cell_type" in contours.columns and contours["cell_type"].notna().any():
            from matplotlib import colormaps
            from matplotlib.lines import Line2D

            # Only contours are colored by cell type; filled polygons would
            # occlude the stain background.
            typed = nearby[nearby["cell_type"].notna()]
            types = sorted(typed["cell_type"].astype(str).unique())
            cmap = colormaps["tab20"]
            type_colors = {t: cmap(i % 20) for i, t in enumerate(types)}
            for (cell_id, cell_type), group in typed.groupby(["cell_id", "cell_type"], sort=False):
                ax.plot(group["x"], group["y"], color=type_colors[str(cell_type)],
                        linewidth=cell_type_linewidth, alpha=0.95)
            type_handles = [
                Line2D([0], [0], color=type_colors[t], linewidth=2, label=t) for t in types
            ]
        plain = nearby if typed.empty else nearby[~nearby["cell_id"].isin(typed["cell_id"])]
        for _, group in plain.groupby("cell_id", sort=False):
            ax.plot(group["x"], group["y"], color="white", linewidth=0.35, alpha=0.35)
        for _, group in target.groupby("cell_id", sort=False):
            ax.plot(group["x"], group["y"], color="white", linewidth=1.2, alpha=0.95)
            ax.plot(group["x"], group["y"], color=contour_color, linewidth=0.75, alpha=0.98)

    molecules = example["molecules"].copy()
    colors = dict(zip(example["role_levels"], ["#2b6cb0", "#c92a2a", "#f08c00"]))
    # Draw outside-target molecules first so target-cell evidence remains
    # visible on top of the local tissue context.
    for role in example["role_levels"]:
        subset = molecules[(~molecules["inside_target"]) & (molecules["role"] == role)]
        if not subset.empty:
            ax.scatter(subset["x"], subset["y"], s=outside_size * non_marker_size_multiplier, c=colors[role], alpha=0.38, linewidths=0)
    for role in example["role_levels"]:
        subset = molecules[molecules["inside_target"] & (~molecules["is_marker"]) & (molecules["role"] == role)]
        if not subset.empty:
            ax.scatter(subset["x"], subset["y"], s=inside_size * non_marker_size_multiplier, c=colors[role], alpha=0.92, linewidths=0)
        markers = molecules[molecules["inside_target"] & molecules["is_marker"] & (molecules["role"] == role)]
        if not markers.empty:
            ax.scatter(markers["x"], markers["y"], s=inside_size * marker_size_multiplier, c=colors[role], alpha=0.92, linewidths=0)

    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", color=colors[role], label=role, markersize=5)
        for role in example["role_levels"]
    ]
    legend = ax.legend(handles=handles, loc="upper left", frameon=True, fontsize=6.5)
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_alpha(0.55)
    legend.get_frame().set_linewidth(0)
    if type_handles:
        type_legend = ax.legend(handles=type_handles, loc="lower right", frameon=True,
                                fontsize=5.8, title="Cell type", title_fontsize=6.2)
        type_legend.get_frame().set_facecolor("white")
        type_legend.get_frame().set_alpha(0.55)
        type_legend.get_frame().set_linewidth(0)
        ax.add_artist(legend)
    ax.set_xlim(bbox[0], bbox[1])
    ax.set_ylim(bbox[3], bbox[2])
    ax.set_aspect("equal")
    ax.set_title(f"{example['target_cell_type']}: {example['target_cell']}", fontsize=8, fontweight="bold")
    ax.axis("off")
    return ax


def plot_examples(examples, *, fit, score_annotation=None, cell_data=None, boundaries=None, stains="auto", cell_types="auto", markers=None, ncol: int = 2, **kwargs):
    """Render a grid of score-selected example cells.

    Stain backgrounds are discovered automatically for Xenium-backed fits;
    pass ``stains=None`` to disable or a mapping of descriptors to override.
    """
    import matplotlib.pyplot as plt

    if examples is None or len(examples) == 0:
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, "No example cells selected", ha="center", va="center")
        ax.axis("off")
        return fig
    stains = resolve_example_stains(fit, stains)
    if isinstance(boundaries, (str, Path)):
        boundaries = read_cell_boundaries(boundaries)
    elif boundaries is None:
        path = discover_cell_boundaries(fit.dataset.source)
        boundaries = read_cell_boundaries(path) if path is not None else pd.DataFrame()
    ncol = max(1, int(ncol))
    nrow = math.ceil(len(examples) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 4.35 * nrow), squeeze=False)
    for ax, (_, row) in zip(axes.ravel(), examples.iterrows()):
        prepared = prepare_cell_example(
            fit,
            row,
            score_annotation=score_annotation,
            cell_data=cell_data,
            boundaries=boundaries,
            stains=stains,
            cell_types=cell_types,
            markers=markers,
        )
        plot_cell_example(prepared, ax=ax, **kwargs)
    for ax in axes.ravel()[len(examples) :]:
        ax.axis("off")
    fig.tight_layout()
    return fig


def square_bbox(xrange, yrange, *, padding: float = 5, min_side: float = 24):
    """Return a square ``(xmin, xmax, ymin, ymax)`` bbox."""
    xmin, xmax = map(float, xrange)
    ymin, ymax = map(float, yrange)
    xmid = (xmin + xmax) / 2
    ymid = (ymin + ymax) / 2
    side = max(xmax - xmin, ymax - ymin, min_side) + 2 * padding
    half = side / 2
    return (xmid - half, xmid + half, ymid - half, ymid + half)


def _normalize_channel(values, probs=(0.01, 0.998), gamma=0.85):
    vals = np.asarray(values, dtype=float)
    finite = vals[np.isfinite(vals)]
    if finite.size == 0:
        return np.zeros_like(vals)
    lo, hi = np.quantile(finite, probs)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = np.nanmin(finite), np.nanmax(finite)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.zeros_like(vals)
    return np.clip((vals - lo) / (hi - lo), 0, 1) ** gamma


def _first_present(frame: pd.DataFrame, names) -> str:
    for name in names:
        if name in frame.columns:
            return name
    raise ValueError(f"Table is missing required columns among: {', '.join(names)}")
