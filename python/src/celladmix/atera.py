"""Read 10x Atera output bundles by writing their Xenium-layout tables.

An Atera output directory carries an ``experiment.spatial`` manifest and
stores molecules, cells and the cell-feature matrix as zipped Zarr stores
(``transcripts.zarr.zip``, ``cells.zarr.zip``, ``cell_feature_matrix.zarr.zip``).
The cellAdmix store builder reads the Xenium layout, so this module writes the
same records into Xenium-named Parquet tables plus an ``experiment.xenium``
manifest, once per bundle, and the rest of the pipeline reads that directory
as an ordinary Xenium bundle:

- ``transcripts.parquet``: transcript_id, cell_id, overlaps_nucleus,
  feature_name, x_location, y_location, z_location, qv, codeword_index,
  codeword_category, is_gene
- ``cells.parquet``: cell_id, x_centroid, y_centroid, transcript_counts,
  cell_area, nucleus_area, nucleus_count
- ``cell_boundaries.parquet`` and ``nucleus_boundaries.parquet``: cell_id,
  vertex_x, vertex_y

Cell ids take the Xenium string form: the eight hexadecimal digits of the
integer id prefix written with the letters a-p (0 -> a, ..., 15 -> p), then
``-`` and the dataset suffix, so (25465, 1) becomes ``aaaagdhj-1``. Molecules
that the Zarr store marks with the id pair (4294967295, 4294967295) are
unassigned and get ``UNASSIGNED``. The feature type of each gene in the
cell-feature matrix ("Gene Expression", "Negative Control Probe", ...) becomes
the molecule's ``codeword_category``, and ``is_gene`` is true for
"Gene Expression" features only.

Reading the Zarr stores requires the ``zarr`` package (version 3 or later).
The conversion also runs as a command, which is how the R package reaches it::

    python -m celladmix.atera [--force] <atera_outs_dir> <out_dir>
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ATERA_MANIFEST = "experiment.spatial"
XENIUM_MANIFEST = "experiment.xenium"
ZARR_STORES = ("transcripts.zarr.zip", "cells.zarr.zip", "cell_feature_matrix.zarr.zip")
XENIUM_TABLES = ("transcripts.parquet", "cells.parquet", "cell_boundaries.parquet", "nucleus_boundaries.parquet")
# Written last, so its presence marks a complete conversion of the recorded source files.
SOURCE_RECORD = "atera_source.json"
# Subdirectory of a cellAdmix output directory that holds the converted tables.
CACHE_DIRNAME = "atera_xenium"

_LETTERS = np.frombuffer(b"abcdefghijklmnop", dtype=np.uint8)
_UNASSIGNED = np.uint32(0xFFFFFFFF)
_MANIFEST_FIELDS = ("run_name", "region_name", "pixel_size", "z_step_size", "num_cells",
                    "panel_wta", "chemistry_version", "analysis_sw_version")


def atera_bundle_dir(path) -> Path | None:
    """The bundle directory when ``path`` is an Atera bundle or its manifest, else None.

    A directory counts as an Atera bundle when it has ``experiment.spatial``
    and no ``experiment.xenium``; a directory with both is read as Xenium.
    """
    if path is None:
        return None
    path = Path(path)
    if path.is_file() and path.name == ATERA_MANIFEST:
        path = path.parent
    if path.is_dir() and (path / ATERA_MANIFEST).exists() and not (path / XENIUM_MANIFEST).exists():
        return path
    return None


def is_atera_bundle(path) -> bool:
    """Whether ``path`` is an Atera output directory or its ``experiment.spatial``."""
    return atera_bundle_dir(path) is not None


def xenium_cell_ids(pairs: np.ndarray) -> np.ndarray:
    """Convert (prefix, suffix) uint32 id pairs to Xenium cell id strings."""
    pairs = np.asarray(pairs)
    prefix = pairs[:, 0].astype(np.uint64)
    digits = np.stack([(prefix >> np.uint64(4 * (7 - k))) & np.uint64(15) for k in range(8)], axis=1)
    body = np.ascontiguousarray(_LETTERS[digits.astype(np.int64)]).view("S8").ravel().astype(str)
    out = np.char.add(np.char.add(body, "-"), pairs[:, 1].astype(str)).astype(object)
    out[(pairs[:, 0] == _UNASSIGNED) | (pairs[:, 1] == _UNASSIGNED)] = "UNASSIGNED"
    return out


def _source_record(src: Path) -> dict:
    """Path, sizes and modification times of the bundle files the conversion reads."""
    files = {}
    for name in (ATERA_MANIFEST,) + ZARR_STORES:
        stat = (src / name).stat()
        files[name] = {"size": int(stat.st_size), "mtime": int(stat.st_mtime)}
    return {"source_dir": str(src.resolve()), "files": files}


def xenium_tables_current(src, out) -> bool:
    """Whether ``out`` holds a complete conversion of the current files of bundle ``src``."""
    src, out = Path(src), Path(out)
    record_path = out / SOURCE_RECORD
    if not record_path.exists() or not all((out / name).exists() for name in XENIUM_TABLES + (XENIUM_MANIFEST,)):
        return False
    try:
        recorded = json.loads(record_path.read_text())
        return recorded.get("files") == _source_record(src)["files"]
    except (OSError, ValueError):
        return False


def _open_zip(path: Path):
    try:
        import zarr
    except ImportError as exc:
        raise ImportError("Reading Atera Zarr stores requires the zarr package (version 3 or later)") from exc
    return zarr.open_group(zarr.storage.ZipStore(str(path), mode="r"), mode="r")


def _strings(array) -> np.ndarray:
    """A Zarr string array (NumPy StringDType) as an object array of Python strings."""
    return np.array(array[:].tolist(), dtype=object)


def _dictionary_strings(codes, names):
    """Integer codes into ``names`` as a plain Arrow string column."""
    import pyarrow as pa

    return pa.DictionaryArray.from_arrays(
        pa.array(np.asarray(codes).astype(np.int32)), pa.array(list(names), pa.string())
    ).cast(pa.string())


def _feature_types(src: Path, gene_names: list) -> np.ndarray:
    """Feature type of each gene, checked against the transcript store's gene order."""
    matrix = _open_zip(src / "cell_feature_matrix.zarr.zip")
    var_names = list(_strings(matrix["var"]["feature_name"]))
    if var_names != list(gene_names):
        raise ValueError("Atera bundle: feature order differs between cell_feature_matrix and transcripts stores")
    feature_type = matrix["var"]["feature_type"]
    return _strings(feature_type["categories"])[feature_type["codes"][:]]


def _transcript_schema():
    import pyarrow as pa

    return pa.schema([
        ("transcript_id", pa.uint64()), ("cell_id", pa.string()), ("overlaps_nucleus", pa.uint8()),
        ("feature_name", pa.string()), ("x_location", pa.float32()), ("y_location", pa.float32()),
        ("z_location", pa.float32()), ("qv", pa.float32()), ("codeword_index", pa.int32()),
        ("codeword_category", pa.string()), ("is_gene", pa.bool_()),
    ])


def _tile_table(tile, first_id: int, codeword_gene: np.ndarray, gene_names: list, gene_types: np.ndarray):
    """One grid tile of the transcript store as a Xenium-layout Arrow table."""
    import pyarrow as pa

    n = tile["location"].shape[0]
    location = tile["location"][:]
    codeword = tile["codeword_identity"][:].ravel()
    gene = codeword_gene[codeword]
    # Convert each distinct cell id once rather than once per molecule.
    ids, inverse = np.unique(tile["cell_id"][:], axis=0, return_inverse=True)
    columns = [
        pa.array(np.arange(first_id, first_id + n, dtype=np.uint64)),
        _dictionary_strings(inverse.ravel(), xenium_cell_ids(ids)),
        pa.array(tile["overlaps_nucleus"][:].ravel().astype(np.uint8)),
        _dictionary_strings(gene, gene_names),
        pa.array(location[:, 0].astype(np.float32)),
        pa.array(location[:, 1].astype(np.float32)),
        pa.array(location[:, 2].astype(np.float32)),
        pa.array(tile["quality_score"][:].ravel().astype(np.float32)),
        pa.array(codeword.astype(np.int32)),
        _dictionary_strings(gene, gene_types),
        pa.array(gene_types[gene] == "Gene Expression"),
    ]
    return pa.Table.from_arrays(columns, schema=_transcript_schema())


def _write_transcripts(src: Path, out: Path, verbose: bool) -> int:
    """Write transcripts.parquet tile by tile; returns the number of molecules."""
    import pyarrow.parquet as pq

    store = _open_zip(src / "transcripts.zarr.zip")
    attrs = dict(store.attrs)
    gene_names = list(attrs["gene_names"])
    codeword_gene = np.asarray(attrs["codeword_gene_mapping"], dtype=np.int64)
    gene_types = _feature_types(src, gene_names)
    grid = store["grid"]
    writer = pq.ParquetWriter(out / "transcripts.parquet", _transcript_schema(), compression="zstd")
    n_written, start = 0, time.time()
    try:
        for key in dict(grid.attrs)["grid_keys"]:
            if key not in grid or grid[key]["location"].shape[0] == 0:
                continue
            table = _tile_table(grid[key], n_written, codeword_gene, gene_names, gene_types)
            writer.write_table(table, row_group_size=8_000_000)
            n_written += table.num_rows
            if verbose:
                print(f"  tile {key}: {table.num_rows:,} molecules "
                      f"({n_written:,} so far, {time.time() - start:.0f} s)", flush=True)
    finally:
        writer.close()
    return n_written


def _write_cells(src: Path, out: Path) -> int:
    """Write cells.parquet and the cell and nucleus boundary tables; returns the number of cells."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    store = _open_zip(src / "cells.zarr.zip")
    ids = xenium_cell_ids(store["cell_id"][:])
    summary = store["cell_summary"][:]
    names = list(dict(store["cell_summary"].attrs)["column_names"])
    column = {name: summary[:, i] for i, name in enumerate(names)}
    matrix = _open_zip(src / "cell_feature_matrix.zarr.zip")
    counts = dict(zip(_strings(matrix["obs"]["barcode"]), matrix["obs"]["transcript_counts"][:]))
    pq.write_table(pa.table({
        "cell_id": pa.array(ids, pa.string()),
        "x_centroid": column["cell_centroid_x"],
        "y_centroid": column["cell_centroid_y"],
        "transcript_counts": np.array([counts.get(i, 0) for i in ids], dtype=np.int64),
        "cell_area": column["cell_area"],
        "nucleus_area": column["nucleus_area"],
        "nucleus_count": column["nucleus_count"],
    }), out / "cells.parquet", compression="zstd")
    # Polygon set 1 holds cell boundaries and set 0 nucleus boundaries; each
    # row of vertices is a padded run of interleaved x, y coordinates.
    for set_index, name in ((1, "cell_boundaries.parquet"), (0, "nucleus_boundaries.parquet")):
        polygons = store[f"polygon_sets/{set_index}"]
        vertices = polygons["vertices"][:]
        n_vertices = polygons["num_vertices"][:].astype(np.int64).ravel()
        cell_index = polygons["cell_index"][:].astype(np.int64).ravel()
        rows = np.repeat(np.arange(len(n_vertices)), n_vertices)
        k = np.concatenate([np.arange(v) for v in n_vertices]) if len(n_vertices) else np.array([], dtype=np.int64)
        pq.write_table(pa.table({
            "cell_id": pa.array(ids[cell_index[rows]], pa.string()),
            "vertex_x": vertices[rows, 2 * k],
            "vertex_y": vertices[rows, 2 * k + 1],
        }), out / name, compression="zstd")
    return len(ids)


def _write_manifest(src: Path, out: Path) -> None:
    spatial = json.loads((src / ATERA_MANIFEST).read_text())
    manifest = {key: spatial[key] for key in _MANIFEST_FIELDS if key in spatial}
    manifest.update(
        transcripts_parquet_filepath="transcripts.parquet",
        source_format="Atera Zarr (experiment.spatial)",
        source_path=str(src.resolve()),
    )
    (out / XENIUM_MANIFEST).write_text(json.dumps(manifest, indent=4))


def write_xenium_tables(src, out, *, verbose: bool = True) -> dict:
    """Write the Xenium-layout tables and manifest for Atera bundle ``src`` into ``out``.

    Returns the numbers of cells and molecules written. The source record is
    written last, so an interrupted conversion is never taken as complete.
    """
    bundle = atera_bundle_dir(src)
    if bundle is None:
        raise ValueError(f"Not an Atera output directory (no {ATERA_MANIFEST} without {XENIUM_MANIFEST}): {src}")
    missing = [name for name in ZARR_STORES if not (bundle / name).exists()]
    if missing:
        raise FileNotFoundError(f"Atera bundle {bundle} is missing {', '.join(missing)}")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for name in (SOURCE_RECORD, XENIUM_MANIFEST):
        (out / name).unlink(missing_ok=True)
    if verbose:
        print(f"Writing Xenium-layout tables for Atera bundle {bundle} into {out}", flush=True)
    record = _source_record(bundle)
    n_cells = _write_cells(bundle, out)
    if verbose:
        print(f"  cells: {n_cells:,}", flush=True)
    n_molecules = _write_transcripts(bundle, out, verbose)
    _write_manifest(bundle, out)
    record.update(n_cells=n_cells, n_molecules=n_molecules)
    (out / SOURCE_RECORD).write_text(json.dumps(record, indent=2))
    return {"n_cells": n_cells, "n_molecules": n_molecules}


def ensure_xenium_tables(src, out, *, force: bool = False, verbose: bool = True) -> Path:
    """Return a Xenium bundle directory for Atera bundle ``src``, converting into ``out`` when needed.

    The tables in ``out`` are reused while the bundle's manifest and Zarr
    stores keep the sizes and modification times recorded at conversion.
    """
    out = Path(out)
    bundle = atera_bundle_dir(src)
    if bundle is None:
        raise ValueError(f"Not an Atera output directory: {src}")
    if force or not xenium_tables_current(bundle, out):
        write_xenium_tables(bundle, out, verbose=verbose)
    elif verbose:
        print(f"Reusing Xenium-layout tables for Atera bundle {bundle} in {out}", flush=True)
    return out


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    force = "--force" in argv
    paths = [arg for arg in argv if arg != "--force"]
    if len(paths) != 2:
        print("Usage: python -m celladmix.atera [--force] <atera_outs_dir> <out_dir>", file=sys.stderr)
        return 2
    ensure_xenium_tables(paths[0], paths[1], force=force)
    return 0


if __name__ == "__main__":
    sys.exit(main())
