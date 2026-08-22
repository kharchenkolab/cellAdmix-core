"""End-to-end test of the generative admixture correction."""
import gzip
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


def make_split_bundle(root: Path, seed: int = 11):
    """A split grid: type A on the left, type B on the right, with A-marker
    molecules planted into B cells in proportion to boundary proximity."""
    rng = np.random.default_rng(seed)
    genes_a = [f"a{i}" for i in range(3)]
    genes_b = [f"b{i}" for i in range(3)]
    side = 18
    cells = []
    molecules = []
    idx = 0
    for gy in range(side):
        for gx in range(side):
            cell_type = "A" if gx < side // 2 else "B"
            cell_id = f"c{idx:03d}"
            cx, cy = gx * 30.0 + 15.0, gy * 30.0 + 15.0
            cells.append((cell_id, cx, cy, 0.0, cell_type))
            own = genes_a if cell_type == "A" else genes_b
            pool = list(rng.choice(own, 60)) + ["h1"] * 10
            if cell_type == "B":
                n_adm = max(0, 12 - 6 * (gx - side // 2))
                pool += list(rng.choice(genes_a, n_adm))
            for gene in pool:
                molecules.append((
                    f"t{len(molecules)}", cell_id, gene,
                    cx + rng.uniform(-10, 10), cy + rng.uniform(-10, 10),
                    rng.uniform(0, 2)))
            idx += 1
    root.mkdir(parents=True, exist_ok=True)
    (root / "experiment.xenium").write_text(json.dumps({
        "run_name": "Mock Xenium", "pixel_size": 0.2125, "z_step_size": 3.0,
    }))
    tx = pd.DataFrame(molecules, columns=[
        "transcript_id", "cell_id", "feature_name", "x_location", "y_location",
        "z_location"])
    tx["qv"] = 30.0
    tx["overlaps_nucleus"] = (np.arange(len(tx)) % 3 == 0).astype(int)
    tx["nucleus_distance"] = (np.arange(len(tx)) % 11).astype(float)
    with gzip.open(root / "transcripts.csv.gz", "wt") as handle:
        tx.to_csv(handle, index=False)
    cell_df = pd.DataFrame(cells, columns=[
        "cell_id", "x_centroid", "y_centroid", "z_centroid", "cell_type"])
    with gzip.open(root / "cells.csv.gz", "wt") as handle:
        cell_df.to_csv(handle, index=False)
    planted = int(sum(1 for m in molecules
                      if m[1].startswith("c") and m[2].startswith("a")
                      and cell_df.set_index("cell_id").cell_type[m[1]] == "B"))
    return cell_df[["cell_id", "cell_type"]], planted


class GenerativeCorrectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from celladmix import CellAdmix

        cls.tmp = tempfile.TemporaryDirectory()
        root = Path(cls.tmp.name)
        annotation, cls.planted = make_split_bundle(root / "bundle")
        cls.ds = CellAdmix(
            str(root / "bundle"),
            output_dir=str(root / "out"),
            annotation=annotation.set_index("cell_id")["cell_type"],
            num_threads=2,
        )
        cls.fit = cls.ds.fit(
            rank=2, ncv_k=6, graph_k=4, nmf_iterations=20, nmf_n_runs=1,
            nmf_train_max_rows=200, seed=11, verbose=False,
        )
        cls.audit = cls.fit.audit_admixture(
            neighbor_k=6, min_target_cells=50, min_reference_cells=20,
            min_excess=50)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_correct_generative_removes_planted_admixture(self):
        correction = self.audit.correct_generative(num_threads=2)
        before, genes, cells = self.audit.fit.counts()
        after, genes_a, cells_a = correction.counts()
        self.assertEqual(list(genes), list(genes_a))
        gi = {g: i for i, g in enumerate(genes)}
        b_cells = np.array([self.audit._ctypes.get(c) == "B" for c in cells_a])
        a_genes = [gi[g] for g in genes if str(g).startswith("a")]
        b_genes = [gi[g] for g in genes if str(g).startswith("b")]
        removed_a = (before[a_genes][:, b_cells].sum()
                     - after[a_genes][:, b_cells].sum())
        removed_b = (before[b_genes][:, b_cells].sum()
                     - after[b_genes][:, b_cells].sum())
        self.assertGreater(removed_a, 0.5 * self.planted)
        # target-owned genes are structurally untouchable
        self.assertAlmostEqual(removed_b, 0.0, delta=1e-6)
        # composition and evaluation
        comp = correction.composition(source="A", target="B")
        self.assertGreater(comp["contamination"].max(), 0.0)
        report = self.audit.evaluate(correction)
        pairs = report.pairs()
        row = pairs[(pairs["source"] == "A") & (pairs["target"] == "B")]
        self.assertGreater(float(row["sensitivity"].iloc[0]), 0.5)


if __name__ == "__main__":
    unittest.main()
