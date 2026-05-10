"""Cell-state embedding helpers backed by the C++ clustering/UMAP code."""

from __future__ import annotations

import numpy as np
import pandas as pd


def clustering_result_to_frame(result: dict, *, annotation: pd.Series | None = None) -> pd.DataFrame:
    """Convert a C++ clustering result dictionary into a tidy cell table."""
    frame = pd.DataFrame(
        {
            "cell_id": result["cell_id"],
            "cell_type": result.get("cell_type", [""] * len(result["cell_id"])),
            "cluster": result["cluster"],
            "transcript_count": result["transcript_count"],
            "detected_genes": result.get("detected_genes", [np.nan] * len(result["cell_id"])),
            "x": result["x"],
            "y": result["y"],
            "z": result["z"],
            "analysis_crop": result.get("analysis_crop", [""] * len(result["cell_id"])),
            "umap_1": result["umap_1"],
            "umap_2": result["umap_2"],
        }
    )
    frame["cell_id"] = frame["cell_id"].astype(str)
    if annotation is not None:
        # Annotation is keyed by cell id; reindex after native filtering so the
        # returned table stays aligned to the C++ result order.
        frame["cell_type"] = annotation.reindex(frame["cell_id"]).to_numpy()

    pcs = np.asarray(result.get("pcs", []), dtype=float)
    if pcs.ndim == 2 and pcs.shape[0] == len(frame):
        # Expose PCs as ordinary columns so plotting and kNN purity helpers can
        # work with plain pandas frames.
        for i in range(pcs.shape[1]):
            frame[f"pc_{i + 1}"] = pcs[:, i]
    return frame


def annotation_knn_purity(
    frame: pd.DataFrame,
    *,
    annotation_col: str = "cell_type",
    pc_prefix: str = "pc_",
    k: int = 15,
) -> dict:
    """Score local annotation separation using same-label kNN purity in PCA space."""
    from sklearn.neighbors import NearestNeighbors

    pc_cols = [col for col in frame.columns if col.startswith(pc_prefix)]
    if not pc_cols:
        raise ValueError("frame does not contain PCA columns")
    valid = frame[annotation_col].notna() & (frame[annotation_col].astype(str) != "")
    data = frame.loc[valid, pc_cols].to_numpy(dtype=float)
    labels = frame.loc[valid, annotation_col].astype(str).to_numpy()
    if data.shape[0] <= k:
        raise ValueError("not enough annotated cells for kNN purity")

    neighbors = NearestNeighbors(n_neighbors=k + 1, metric="cosine").fit(data)
    neighbor_idx = neighbors.kneighbors(data, return_distance=False)[:, 1:]
    purity = float(np.mean(labels[neighbor_idx] == labels[:, None]))
    # Normalize against the same-label rate expected from label frequencies.
    frequencies = pd.Series(labels).value_counts(normalize=True).to_numpy()
    baseline = float(np.sum(frequencies**2))
    normalized = float((purity - baseline) / max(1e-12, 1 - baseline))
    return {
        "knn_same_type_fraction": purity,
        "baseline_fraction": baseline,
        "normalized_purity": normalized,
    }
