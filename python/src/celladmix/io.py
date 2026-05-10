"""Small I/O helpers shared by the Python API and notebooks."""

from __future__ import annotations

import json
import math
import os
import re
from pathlib import Path
from typing import Iterable, Optional, Tuple

import pandas as pd


DEFAULT_ANNOTATION_COLUMNS = (
    "cell_type",
    "merged_annotation",
    "annotation",
    "cluster",
    "cluster_id",
)


def read_annotation(
    annotation,
    *,
    cell_id_col: str = "cell_id",
    annotation_col: Optional[str] = None,
) -> pd.Series:
    """Return a cell-id-indexed annotation series from a path, Series, dict, or DataFrame."""
    if annotation is None:
        raise ValueError("annotation is required")
    if isinstance(annotation, pd.Series):
        out = annotation.copy()
        out = out.dropna()
        out.index = out.index.astype(str)
        out = out.astype(str)
        out = out[out != ""]
        out.name = annotation_col or annotation.name or "annotation"
        return out
    if isinstance(annotation, dict):
        return pd.Series(annotation, dtype="object", name=annotation_col or "annotation").astype(str)
    if isinstance(annotation, (str, os.PathLike)):
        path = Path(annotation)
        if path.suffix == ".parquet":
            frame = pd.read_parquet(path)
        else:
            frame = pd.read_csv(path)
        return read_annotation(frame, cell_id_col=cell_id_col, annotation_col=annotation_col)
    if not isinstance(annotation, pd.DataFrame):
        raise TypeError("annotation must be a pandas Series, dict, DataFrame, or file path")

    frame = annotation.copy()
    if cell_id_col not in frame.columns:
        if frame.index.name:
            frame = frame.reset_index()
            cell_id_col = frame.columns[0]
        else:
            raise ValueError(f"annotation table does not contain cell id column {cell_id_col!r}")
    label_col = annotation_col
    if label_col is None:
        # Prefer common annotation column names so examples can pass raw
        # clustering/metadata tables without extra boilerplate.
        for candidate in DEFAULT_ANNOTATION_COLUMNS:
            if candidate in frame.columns and candidate != cell_id_col:
                label_col = candidate
                break
    if label_col is None or label_col not in frame.columns:
        raise ValueError(
            "annotation_col was not provided and no standard annotation column was found"
        )
    valid = frame[cell_id_col].notna() & frame[label_col].notna()
    frame = frame.loc[valid]
    out = pd.Series(
        frame[label_col].astype(str).to_numpy(),
        index=frame[cell_id_col].astype(str).to_numpy(),
        name=label_col,
    )
    out = out[out != ""]
    return out[~out.index.duplicated(keep="first")]


def annotation_vectors(annotation: Optional[pd.Series]) -> Tuple[list[str], list[str]]:
    """Return C++-friendly annotation vectors."""
    if annotation is None:
        return [], []
    clean = annotation.dropna().astype(str)
    clean = clean[clean != ""]
    return clean.index.astype(str).tolist(), clean.tolist()


def recommended_rank(annotation: Optional[pd.Series], *, multiplier: float = 1.2, cap: int = 30) -> int:
    """Default NMF rank policy: about 1.2 factors per annotated cell type."""
    if annotation is None:
        return 4
    n_labels = int(annotation.dropna().astype(str).nunique())
    if n_labels <= 0:
        return 4
    return max(2, min(cap, int(math.ceil(multiplier * n_labels))))


def discover_xenium_stain_image(
    bundle_dir: str | os.PathLike,
    stain: str = "membrane",
    *,
    focus_index: Optional[int] = None,
) -> dict:
    """Find a Xenium morphology focus image and pixel size.

    Xenium bundles store image paths in ``experiment.xenium``. The first focus
    image is usually DAPI/nuclear and focus index 1 is the membrane channel for
    the 10x multimodal membrane examples.
    """
    bundle = Path(bundle_dir)
    manifest_path = bundle / "experiment.xenium"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Cannot find Xenium manifest: {manifest_path}")
    with manifest_path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    images = manifest.get("images", {})
    focus = images.get("morphology_focus_filepath")
    if not focus:
        raise ValueError("Xenium manifest does not list morphology_focus_filepath")
    if focus_index is None:
        key = stain.lower().replace("-", "_")
        # The Xenium focus stack convention used by the membrane examples puts
        # membrane signal in focus plane 1 and DAPI/nuclear signal in plane 0.
        focus_index = 1 if key in {"membrane", "boundary", "cell_boundary"} else 0
    focus_path = Path(focus)
    image_name = re.sub(r"_(\d{4})(?=\.)", f"_{focus_index:04d}", focus_path.name)
    image_path = bundle / focus_path.parent / image_name
    if not image_path.exists():
        image_path = bundle / focus
    return {
        "stain": stain,
        "image_path": str(image_path),
        "pixel_size": float(manifest.get("pixel_size", 1.0)),
        "x_offset": 0.0,
        "y_offset": 0.0,
    }


def discover_xenium_membrane_image(bundle_dir: str | os.PathLike) -> dict:
    """Find the default Xenium membrane-focused image and pixel size."""
    return discover_xenium_stain_image(bundle_dir, "membrane", focus_index=1)


def ensure_list(value: Optional[Iterable[str]]) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return list(value)
