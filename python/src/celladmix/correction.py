"""Corrected run objects."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import _core
from .state import clustering_result_to_frame


class CellAdmixCorrection:
    """A corrected cellAdmix run."""

    def __init__(self, run_path: str, manifest: dict, rules=None):
        self.run_path = Path(run_path)
        self.manifest = manifest
        self.rules = rules

    def __repr__(self) -> str:
        removed = self.manifest.get("n_removed", "unknown")
        return f"CellAdmixCorrection(run_path={str(self.run_path)!r}, n_removed={removed})"

    def summary(self) -> pd.DataFrame:
        """Return a compact per-cell-type correction summary."""
        return _compact_correction_summary(self.cell_summary())

    def ensemble(self) -> dict:
        """Return molecule-vote ensemble diagnostics for this correction."""
        histogram = self.manifest.get("vote_histogram")
        return {
            "members": self.manifest.get("n_members", 1),
            "min_votes": self.manifest.get("min_votes", 1),
            "n_removed": self.manifest.get("n_removed"),
            "vote_histogram": histogram,
        }

    def cell_summary(self) -> pd.DataFrame:
        """Return per-cell molecule removal statistics."""
        path = self.manifest.get("correction_summary_parquet")
        if path is None:
            path = self.run_path / "correction_summary.parquet"
        return pd.read_parquet(path)

    def counts(self):
        """Return corrected counts as a scipy CSC matrix with genes x cells orientation."""
        from scipy import sparse

        # The native layer returns CSC arrays directly; keep orientation explicit
        # here because AnnData expects the transpose.
        raw = _core.collect_counts(str(self.run_path))
        data = np.asarray(raw["data"], dtype=float)
        indices = np.asarray(raw["indices"], dtype=np.int32)
        indptr = np.asarray(raw["indptr"], dtype=np.int64)
        genes = list(raw["genes"])
        cells = list(raw["cells"])
        matrix = sparse.csc_matrix((data, indices, indptr), shape=(len(genes), len(cells)))
        return matrix, pd.Index(genes, name="gene"), pd.Index(cells, name="cell_id")

    def to_anndata(self):
        """Return corrected counts as an AnnData object with cells x genes orientation."""
        import anndata as ad

        matrix, genes, cells = self.counts()
        return ad.AnnData(X=matrix.T.tocsr(), obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))

    def cell_state_umap(
        self,
        *,
        annotation: pd.Series | None = None,
        cells_max: int = 5000,
        min_molecules: int = 10,
        min_genes: int = 5,
        n_variable_genes: int = 1000,
        pca_dims: int = 30,
        graph_k: int = 15,
        umap_neighbors: int = 15,
        umap_epochs: int = 200,
        normalization_scale: float = 5000.0,
        seed: int = 1,
        num_threads: int = 1,
    ) -> pd.DataFrame:
        """Compute a cell-state UMAP from corrected counts."""
        # Corrected runs are already materialized on disk, so clustering can
        # reuse the same native count-collection path as uncorrected runs.
        result = _core.cluster_run_counts(
            str(self.run_path),
            min_molecules=min_molecules,
            min_genes=min_genes,
            cells_max=int(cells_max),
            n_variable_genes=n_variable_genes,
            pca_dims=pca_dims,
            graph_k=graph_k,
            compute_umap=True,
            umap_neighbors=umap_neighbors,
            umap_epochs=umap_epochs,
            num_threads=int(num_threads),
            normalization_scale=normalization_scale,
            seed=int(seed),
        )
        return clustering_result_to_frame(result, annotation=annotation)

    def plot_removed_molecules(self):
        """Plot per-cell molecule removal counts and fractions."""
        from .plotting import plot_molecule_removal

        return plot_molecule_removal(self.cell_summary())


def _compact_correction_summary(summary: pd.DataFrame) -> pd.DataFrame:
    required = {
        "cell_type",
        "n_molecules_before",
        "n_molecules_after",
        "n_removed",
        "fraction_removed",
    }
    if summary.empty:
        return pd.DataFrame()
    missing = required.difference(summary.columns)
    if missing:
        raise ValueError(f"correction summary is missing columns: {', '.join(sorted(missing))}")
    data = summary.copy()
    data["cell_type"] = data["cell_type"].fillna("unknown").replace("", "unknown").astype(str)
    data["modified"] = data["n_removed"].astype(float) > 0

    def summarize(group: pd.DataFrame, label: str) -> dict:
        modified = group[group["modified"]]
        molecules_before = float(group["n_molecules_before"].sum())
        molecules_removed = float(group["n_removed"].sum())
        return {
            "cell_type": label,
            "n_cells": int(len(group)),
            "n_modified_cells": int(group["modified"].sum()),
            "fraction_cells_modified": float(group["modified"].sum() / max(1, len(group))),
            "molecules_before": molecules_before,
            "molecules_after": float(group["n_molecules_after"].sum()),
            "molecules_removed": molecules_removed,
            "fraction_molecules_removed": float(molecules_removed / max(1.0, molecules_before)),
            "median_removed_per_modified_cell": float(modified["n_removed"].median()) if len(modified) else 0.0,
            "median_fraction_removed_per_modified_cell": (
                float(modified["fraction_removed"].median()) if len(modified) else 0.0
            ),
        }

    rows = [summarize(data, "all")]
    for cell_type, group in data.groupby("cell_type", sort=True):
        rows.append(summarize(group, str(cell_type)))
    return pd.DataFrame(rows)


# ---- capping a molecule-level correction at a count-level amount ------------

def _store_molecule_votes(fit, member_rules, rows):
    """Votes per molecule row of the parent store: the number of ensemble
    members whose rules remove it. ``member_rules`` is the dict
    ``CellAdmixScore._ensemble_member_rules`` returns; ``rows`` restricts
    the computation to those molecule rows (the removed ones)."""
    import pyarrow.parquet as pq

    run = Path(fit.run_path)
    mol = pq.read_table(run / "molecules.parquet", columns=["cell_idx"]).column(0).to_numpy()[rows]
    cells = pq.read_table(run / "cells.parquet", columns=["cell_idx", "cell_id"]).to_pandas()
    ann = fit.dataset.annotation
    cell_type = pd.Series(
        [str(ann.get(str(c), None)) for c in cells["cell_id"].astype(str)],
        index=cells["cell_idx"].to_numpy())
    types = pd.Index(sorted(set(cell_type)))
    type_code = types.get_indexer(cell_type.reindex(mol).to_numpy())
    rules = member_rules["rules"].reset_index(drop=True)
    members = np.asarray(member_rules["member"], dtype=int)
    n_factors = int(fit.manifest["n_factors"])
    votes = np.zeros(len(rows), dtype=np.int32)
    for member in np.unique(members):
        sub = rules[members == member]
        table = np.zeros((len(types), n_factors), dtype=bool)
        for factor, target in zip(sub["factor"].astype(int), sub["target_cell_type"].astype(str)):
            t = types.get_loc(target) if target in types else -1
            if t >= 0 and 1 <= factor <= n_factors:
                table[t, factor - 1] = True
        path = run / f"ensemble_labels_m{int(member)}.parquet"
        if int(member) >= 0 and path.exists():
            labels = pq.read_table(path, columns=["label"]).column(0).to_numpy()[rows]
        else:
            labels = pq.read_table(run / "molecules.parquet",
                                   columns=["factor_label"]).column(0).to_numpy()[rows]
        ok = (type_code >= 0) & (labels >= 0) & (labels < n_factors)
        hit = np.zeros(len(rows), dtype=bool)
        hit[ok] = table[type_code[ok], labels[ok]]
        votes += hit
    return votes


def _cap_amounts(cap, fit):
    """The count-level removal amounts (genes x cells, fractional) mapped to
    the parent store's gene and cell indices, as a sparse matrix over
    (store cell_idx, store gene_idx). ``cap`` is a fitted generative model
    or a ``(removed_matrix, genes, cell_ids)`` triple with the matrix in
    genes x cells orientation."""
    import pyarrow.parquet as pq
    from scipy import sparse

    if hasattr(cap, "_removed") and hasattr(cap, "_matrix"):
        matrix = sparse.csc_matrix(
            (np.asarray(cap._removed, dtype=float), cap._matrix.indices, cap._matrix.indptr),
            shape=cap._matrix.shape)
        genes, cell_ids = cap._genes, cap._cells
    else:
        matrix, genes, cell_ids = cap
        matrix = sparse.csc_matrix(matrix, dtype=float)
    store_genes = pd.Index([str(g) for g in fit.manifest["genes"]])
    cells = pq.read_table(Path(fit.run_path) / "cells.parquet",
                          columns=["cell_idx", "cell_id"]).to_pandas()
    store_cells = pd.Index(cells["cell_id"].astype(str))
    gi = store_genes.get_indexer([str(g) for g in genes])
    ci = store_cells.get_indexer([str(c) for c in cell_ids])
    coo = matrix.tocoo()
    ok = (gi[coo.row] >= 0) & (ci[coo.col] >= 0) & (coo.data > 0)
    cell_idx = cells["cell_idx"].to_numpy()[ci[coo.col[ok]]]
    n_cells = int(cells["cell_idx"].max()) + 1
    return sparse.csr_matrix((coo.data[ok], (cell_idx, gi[coo.row[ok]])),
                             shape=(n_cells, len(store_genes)))


def cap_correction(fit, correction, cap, *, votes=None, seed=0, scale=1.0):
    """Limit a molecule-level correction to a count-level removal amount.

    For every cell and gene the number of molecules the correction removed
    is reduced to the amount ``cap`` attributes to admixture there (a fitted
    generative model, or a ``(removed_matrix, genes, cell_ids)`` triple),
    times ``scale``, rounded stochastically; the surplus molecules are restored, those with
    the fewest ensemble votes first (``votes`` per removed molecule, in the
    order of the parent store's removed rows; ties and missing votes at
    random). The corrected store is rewritten in place (molecules, cells,
    the per-cell summary and the manifest) and the correction object
    returned with its manifest updated.
    """
    import json

    import pyarrow as pa
    import pyarrow.parquet as pq

    run = Path(fit.run_path)
    out = Path(correction.run_path)
    parent = pq.read_table(run / "molecules.parquet")
    n = parent.num_rows
    obs = parent.column("obs_id").to_numpy()
    if not np.array_equal(obs, np.arange(n)):
        raise ValueError("the parent store's molecule identifiers are not its row indices")
    kept = pq.read_table(out / "molecules.parquet", columns=["obs_id"]).column(0).to_numpy()
    removed = np.ones(n, dtype=bool)
    removed[kept] = False
    rows = np.flatnonzero(removed)
    gene_idx = parent.column("gene_idx").to_numpy()[rows]
    cell_idx = parent.column("cell_idx").to_numpy()[rows]
    amounts = _cap_amounts(cap, fit)
    n_genes = amounts.shape[1]
    key = cell_idx.astype(np.int64) * n_genes + gene_idx
    uniq, inverse = np.unique(key, return_inverse=True)
    u_cell, u_gene = uniq // n_genes, uniq % n_genes
    amount = np.asarray(amounts[u_cell, u_gene]).ravel() * float(scale)
    rng = np.random.default_rng(seed)
    cap_count = np.floor(amount).astype(np.int64) + (rng.random(len(uniq)) < amount - np.floor(amount))
    v = np.zeros(len(rows)) if votes is None else -np.asarray(votes, dtype=float)
    order = np.lexsort((rng.random(len(rows)), v, inverse))       # by key, most votes first, random
    rank = np.empty(len(rows), dtype=np.int64)
    starts = np.r_[0, np.flatnonzero(np.diff(inverse[order])) + 1]
    rank[order] = np.arange(len(rows)) - np.repeat(starts, np.diff(np.r_[starts, len(rows)]))
    restore = rank >= cap_count[inverse]
    keep_new = ~removed
    keep_new[rows[restore]] = True
    # the corrected store: molecules in the writer's tile order, the cells table
    # and the per-cell summary recomputed from the kept molecules
    kept_table = parent.filter(pa.array(keep_new)).sort_by(
        [("crop_idx", "ascending"), ("tile_x", "ascending"), ("tile_y", "ascending"), ("obs_id", "ascending")])
    row_group = int(fit.manifest.get("storage_options", {}).get("parquet_row_group_size", 65536))
    pq.write_table(kept_table, out / "molecules.parquet", row_group_size=row_group, compression="none")
    n_factors = int(fit.manifest["n_factors"])
    labels = parent.column("factor_label").to_numpy()[keep_new]
    cells_kept = parent.column("cell_idx").to_numpy()[keep_new]
    cells_tbl = pq.read_table(out / "cells.parquet")
    cells_df = cells_tbl.to_pandas()
    n_cells = len(cells_df)
    pos = pd.Index(cells_df["cell_idx"].to_numpy()).get_indexer(cells_kept)
    counts = np.bincount(pos[pos >= 0], minlength=n_cells)
    frac = np.zeros((n_cells, n_factors))
    lab_ok = (pos >= 0) & (labels >= 0) & (labels < n_factors)
    np.add.at(frac, (pos[lab_ok], labels[lab_ok]), 1.0)
    frac = frac / np.maximum(counts, 1)[:, None]
    cells_df["transcript_count"] = counts.astype(np.int32)
    cells_df["dominant_factor"] = np.where(counts > 0, frac.argmax(axis=1) + 1, -1).astype(np.int32)
    cells_df["dominant_fraction"] = np.where(counts > 0, frac.max(axis=1), 0.0)
    # one fraction column per factor of the fit (the native writer keeps only
    # as many as the kept labels reach; restored molecules may need more)
    cells_df = cells_df[[c for c in cells_df.columns if not c.startswith("factor_")]]
    fields = [f for f in cells_tbl.schema if not f.name.startswith("factor_")]
    for k in range(n_factors):
        cells_df[f"factor_{k + 1}_fraction"] = frac[:, k]
        fields.append(pa.field(f"factor_{k + 1}_fraction", pa.float64()))
    pq.write_table(pa.Table.from_pandas(cells_df, schema=pa.schema(fields), preserve_index=False),
                   out / "cells.parquet", row_group_size=row_group, compression="none")
    summary_path = Path(correction.manifest.get("correction_summary_parquet")
                        or out / "correction_summary.parquet")
    summary_tbl = pq.read_table(summary_path)
    summary = summary_tbl.to_pandas()
    spos = pd.Index(cells_df["cell_id"].astype(str)).get_indexer(summary["cell_id"].astype(str))
    after = np.where(spos >= 0, counts[np.maximum(spos, 0)], 0)
    before = summary["n_molecules_before"].to_numpy()
    summary["n_molecules_after"] = after.astype(np.int32)
    summary["n_removed"] = (before - after).astype(np.int32)
    summary["fraction_removed"] = (before - after) / np.maximum(before, 1)
    summary["modified"] = ((before - after) > 0).astype(np.int32)
    pq.write_table(pa.Table.from_pandas(summary, schema=summary_tbl.schema, preserve_index=False),
                   summary_path, row_group_size=row_group, compression="none")
    manifest_path = out / "run.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["n_transcripts"] = int(keep_new.sum())
    manifest["cap"] = dict(amount_molecules=float(amount.sum()), cap_molecules=int(cap_count.sum()),
                           restored=int(restore.sum()), seed=int(seed), scale=float(scale))
    manifest_path.write_text(json.dumps(manifest, indent=1))
    correction.manifest["n_transcripts"] = int(keep_new.sum())
    correction.manifest["n_removed"] = int(n - keep_new.sum())
    correction.manifest["cap"] = manifest["cap"]
    return correction
