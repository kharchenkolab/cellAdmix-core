"""Tests for reading Atera output bundles through Xenium-layout tables."""
import json
import os
import tempfile
import time
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import zarr
except ImportError:  # pragma: no cover - the conversion needs zarr
    zarr = None

UNASSIGNED = 4294967295
GENES = ["GA", "GB", "NegControlProbe_1"]
FEATURE_TYPES = ["Gene Expression", "Gene Expression", "Negative Control Probe"]
# Codeword -> gene index; codeword 3 is a second codeword of GA.
CODEWORD_GENE = [0, 1, 2, 0]
CELL_PAIRS = [(0, 1), (25465, 1), (15, 1)]
CELL_IDS = ["aaaaaaaa-1", "aaaagdhj-1", "aaaaaaap-1"]
CENTROIDS = [(10.0, 12.0), (30.0, 14.0), (20.0, 40.0)]


def _string_array(values):
    return np.array(values, dtype=np.dtypes.StringDType())


def _zip_group(path: Path, attributes=None):
    # Attributes are given at creation: a zip entry cannot be rewritten in place.
    store = zarr.storage.ZipStore(str(path), mode="w")
    return store, zarr.open_group(store, mode="w", attributes=attributes or {})


def _tile_molecules(tile_index):
    """Molecules of one grid tile as (cell pair, codeword, x, y, z, qv, overlaps_nucleus) rows."""
    rows = []
    for cell, (cx, cy) in enumerate(CENTROIDS):
        for k, codeword in enumerate((0, 1, 3, 2, 0)):
            rows.append((CELL_PAIRS[cell], codeword, cx + k + tile_index, cy - k, 1.5 + tile_index,
                         30.0 + k, int(k % 2 == 0)))
    for k in range(3):
        rows.append(((UNASSIGNED, UNASSIGNED), k, 50.0 + k, 50.0 + tile_index, 2.0, 20.0, 0))
    return rows


def make_atera_bundle(root: Path) -> dict:
    """Write a small Atera bundle: three cells, two grid tiles, unassigned molecules, one control probe."""
    root.mkdir(parents=True, exist_ok=True)
    (root / "experiment.spatial").write_text(json.dumps({
        "run_name": "Mock Atera", "region_name": "region", "pixel_size": 0.2125,
        "z_step_size": 3.0, "num_cells": len(CELL_IDS), "analysis_sw_version": "atera-1.0.4",
    }))

    tiles = {"0,0": _tile_molecules(0), "1,0": _tile_molecules(1)}
    store, group = _zip_group(root / "transcripts.zarr.zip", {
        "gene_names": GENES, "codeword_gene_mapping": CODEWORD_GENE,
        "codeword_gene_names": [GENES[g] for g in CODEWORD_GENE]})
    # A listed key without a stored tile must be skipped.
    grid = group.create_group("grid", attributes={"grid_keys": list(tiles) + ["2,0"]})
    for key, rows in tiles.items():
        tile = grid.create_group(key)
        tile.create_array("location", data=np.array([r[2:5] for r in rows], dtype=np.float32))
        tile.create_array("codeword_identity", data=np.array([[r[1]] for r in rows], dtype=np.uint32))
        tile.create_array("cell_id", data=np.array([r[0] for r in rows], dtype=np.uint32))
        tile.create_array("overlaps_nucleus", data=np.array([[r[6]] for r in rows], dtype=np.uint8))
        tile.create_array("quality_score", data=np.array([[r[5]] for r in rows], dtype=np.float16))
        tile.create_array("gene_offset", data=np.zeros((len(GENES), 2), dtype=np.uint32))
    store.close()

    store, group = _zip_group(root / "cells.zarr.zip")
    group.create_array("cell_id", data=np.array(CELL_PAIRS, dtype=np.uint32))
    columns = ["cell_centroid_x", "cell_centroid_y", "cell_area", "nucleus_centroid_x",
               "nucleus_centroid_y", "nucleus_area", "z_level", "nucleus_count"]
    summary = np.array([[cx, cy, 100.0, cx, cy, 30.0, 0.0, 1.0] for cx, cy in CENTROIDS])
    group.create_array("cell_summary", data=summary, attributes={"column_names": columns})
    for set_index, half_width in ((0, 2.0), (1, 5.0)):
        polygons = group.create_group(f"polygon_sets/{set_index}")
        n_vertices = np.array([4, 4, 5], dtype=np.uint32)
        vertices = np.zeros((len(CELL_IDS), 50), dtype=np.float32)
        for cell, (cx, cy) in enumerate(CENTROIDS):
            for k in range(n_vertices[cell]):
                angle = 2 * np.pi * k / n_vertices[cell]
                vertices[cell, 2 * k] = cx + half_width * np.cos(angle)
                vertices[cell, 2 * k + 1] = cy + half_width * np.sin(angle)
        polygons.create_array("vertices", data=vertices)
        polygons.create_array("num_vertices", data=n_vertices)
        polygons.create_array("cell_index", data=np.arange(len(CELL_IDS), dtype=np.uint32))
    store.close()

    store, group = _zip_group(root / "cell_feature_matrix.zarr.zip")
    group.create_array("var/feature_name", data=_string_array(GENES))
    categories = ["Gene Expression", "Negative Control Probe"]
    group.create_array("var/feature_type/categories", data=_string_array(categories))
    group.create_array("var/feature_type/codes",
                       data=np.array([categories.index(t) for t in FEATURE_TYPES], dtype=np.int8))
    group.create_array("obs/barcode", data=_string_array(CELL_IDS))
    group.create_array("obs/transcript_counts", data=np.array([8, 8, 8], dtype=np.int64))
    store.close()
    n_molecules = sum(len(rows) for rows in tiles.values())
    return {"n_molecules": n_molecules, "n_unassigned": 6}


@unittest.skipIf(zarr is None, "zarr is not installed")
class AteraConversionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from celladmix import atera

        cls.atera = atera
        cls.tmp = tempfile.TemporaryDirectory()
        root = Path(cls.tmp.name)
        cls.bundle = root / "outs"
        cls.expected = make_atera_bundle(cls.bundle)
        cls.out = root / "converted"
        cls.summary = atera.write_xenium_tables(cls.bundle, cls.out, verbose=False)
        cls.tx = pd.read_parquet(cls.out / "transcripts.parquet")
        cls.cells = pd.read_parquet(cls.out / "cells.parquet")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_bundle_detection(self):
        self.assertTrue(self.atera.is_atera_bundle(self.bundle))
        self.assertTrue(self.atera.is_atera_bundle(self.bundle / "experiment.spatial"))
        # The converted directory has experiment.xenium and is read as Xenium.
        self.assertFalse(self.atera.is_atera_bundle(self.out))

    def test_cell_ids(self):
        ids = self.atera.xenium_cell_ids(np.array(CELL_PAIRS + [(UNASSIGNED, UNASSIGNED)], dtype=np.uint32))
        self.assertEqual(list(ids), CELL_IDS + ["UNASSIGNED"])
        self.assertEqual(list(self.cells["cell_id"]), CELL_IDS)

    def test_transcripts(self):
        tx = self.tx
        self.assertEqual(len(tx), self.expected["n_molecules"])
        self.assertEqual(self.summary["n_molecules"], self.expected["n_molecules"])
        self.assertEqual(list(tx["transcript_id"]), list(range(len(tx))))
        self.assertEqual((tx["cell_id"] == "UNASSIGNED").sum(), self.expected["n_unassigned"])
        self.assertEqual(set(tx["cell_id"]), set(CELL_IDS) | {"UNASSIGNED"})
        self.assertFalse(tx["cell_id"].str.startswith("pppppppp").any())
        self.assertEqual(set(tx["feature_name"]), set(GENES))
        control = tx["feature_name"] == "NegControlProbe_1"
        self.assertTrue((tx.loc[control, "codeword_category"] == "Negative Control Probe").all())
        self.assertFalse(tx.loc[control, "is_gene"].any())
        self.assertTrue((tx.loc[~control, "codeword_category"] == "Gene Expression").all())
        self.assertTrue(tx.loc[~control, "is_gene"].all())
        # Codeword 3 maps to GA.
        self.assertTrue((tx.loc[tx["codeword_index"] == 3, "feature_name"] == "GA").all())
        first = tx[tx["cell_id"] == "aaaaaaaa-1"].iloc[0]
        self.assertAlmostEqual(first["x_location"], 10.0)
        self.assertAlmostEqual(first["y_location"], 12.0)
        self.assertAlmostEqual(first["z_location"], 1.5)
        self.assertAlmostEqual(first["qv"], 30.0)
        self.assertEqual(first["overlaps_nucleus"], 1)

    def test_cells_and_boundaries(self):
        cells = self.cells.set_index("cell_id")
        for cell_id, (cx, cy) in zip(CELL_IDS, CENTROIDS):
            self.assertAlmostEqual(cells.loc[cell_id, "x_centroid"], cx)
            self.assertAlmostEqual(cells.loc[cell_id, "y_centroid"], cy)
        self.assertEqual(list(cells["transcript_counts"]), [8, 8, 8])
        for name, half_width in (("cell_boundaries.parquet", 5.0), ("nucleus_boundaries.parquet", 2.0)):
            bounds = pd.read_parquet(self.out / name)
            self.assertEqual(bounds.groupby("cell_id").size().to_dict(),
                             {"aaaaaaaa-1": 4, "aaaagdhj-1": 4, "aaaaaaap-1": 5})
            first = bounds[bounds["cell_id"] == "aaaagdhj-1"].iloc[0]
            self.assertAlmostEqual(first["vertex_x"], 30.0 + half_width, places=5)
            self.assertAlmostEqual(first["vertex_y"], 14.0, places=5)

    def test_manifest_and_reuse(self):
        manifest = json.loads((self.out / "experiment.xenium").read_text())
        self.assertEqual(manifest["transcripts_parquet_filepath"], "transcripts.parquet")
        self.assertEqual(manifest["run_name"], "Mock Atera")
        self.assertTrue(self.atera.xenium_tables_current(self.bundle, self.out))
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "cache"
            self.atera.ensure_xenium_tables(self.bundle, out, verbose=False)
            written = (out / "transcripts.parquet").stat().st_mtime_ns
            self.atera.ensure_xenium_tables(self.bundle, out, verbose=False)
            self.assertEqual((out / "transcripts.parquet").stat().st_mtime_ns, written)
            # A changed source file invalidates the conversion.
            spatial = self.bundle / "experiment.spatial"
            stamp = time.time() + 5
            os.utime(spatial, (stamp, stamp))
            try:
                self.assertFalse(self.atera.xenium_tables_current(self.bundle, out))
            finally:
                os.utime(spatial, None)


@unittest.skipIf(zarr is None, "zarr is not installed")
class AteraStoreTests(unittest.TestCase):
    def test_dataset_builds_store_from_atera_bundle(self):
        from celladmix import CellAdmix, _core

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            expected = make_atera_bundle(root / "outs")
            ds = CellAdmix(str(root / "outs"), output_dir=str(root / "out"), num_threads=1)
            ds.ensure_store(verbose=False)
            self.assertTrue((root / "out" / "atera_xenium" / "experiment.xenium").exists())
            manifest = _core.read_input_store_manifest(str(ds.store_dir))
            molecules = pd.read_parquet(ds.store_dir / "molecules.parquet")
            # Unassigned molecules and the control probe are dropped by default.
            n_kept = expected["n_molecules"] - expected["n_unassigned"] - 2 * len(CELL_IDS)
            self.assertEqual(len(molecules), n_kept)
            self.assertEqual(int(manifest["n_molecules"]), n_kept)
            self.assertEqual(int(manifest["n_cells"]), len(CELL_IDS))
            # Boundary discovery follows the dataset to the converted tables.
            self.assertEqual(ds.xenium_source(verbose=False), root / "out" / "atera_xenium")
