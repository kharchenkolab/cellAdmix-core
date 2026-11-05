"""The cap on the factor-based cleanup: the removal in every cell and gene
is limited to the count-level amount of a fitted generative model (or any
removal matrix), the vote's choice of molecules is kept under the cap, the
rewritten corrected store stays consistent, and a cap above every removal
reproduces the uncapped store."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_generative import make_split_bundle  # noqa: E402


class CapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from celladmix import CellAdmix

        cls.tmp = tempfile.TemporaryDirectory()
        root = Path(cls.tmp.name)
        annotation, _ = make_split_bundle(root / "bundle")
        cls.ds = CellAdmix(
            str(root / "bundle"), output_dir=str(root / "out"),
            annotation=annotation.set_index("cell_id")["cell_type"], num_threads=2)
        cls.fit = cls.ds.fit(
            rank=2, ncv_k=6, graph_k=4, nmf_iterations=20, nmf_n_runs=3,
            nmf_train_max_rows=200, seed=11, verbose=False)
        cls.audit = cls.fit.audit_admixture(
            neighbor_k=6, min_target_cells=50, min_reference_cells=20, min_excess=50)
        cls.model = cls.audit.fit_generative(num_threads=2)
        cls.score = cls.fit.score_bridge(
            candidate_k=5, crossing_k=5, min_type_pair_contacts=1,
            min_factor_molecules=1, min_pairs=1, null_iterations=1,
            null_max_iterations=2, compute_null=True, verbose=False)
        # both factors removed from B cells: without the cap the correction
        # empties the B cells, so the cap has to do all the work
        cls.rules = pd.DataFrame({
            "factor": [1, 2], "target_cell_type": ["B", "B"],
            "source_cell_type": ["A", "A"]})
        cls.uncapped = cls.score.correct(cls.rules, name="cap_none")
        cls.capped = cls.score.correct(cls.rules, name="cap_model", cap=cls.model)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def _removed(self, correction):
        before, genes, cells = self.fit.counts()
        after, g2, c2 = correction.counts()
        self.assertEqual(list(genes), list(g2))
        self.assertEqual(list(cells), list(c2))
        return np.asarray((before - after).todense()), genes, cells

    def test_removal_never_exceeds_the_amount(self):
        removed, genes, cells = self._removed(self.capped)
        amount = sparse.csc_matrix(
            (self.model._removed, self.model._matrix.indices, self.model._matrix.indptr),
            shape=self.model._matrix.shape)
        amount = pd.DataFrame(np.asarray(amount.todense()), index=self.model._genes,
                              columns=self.model._cells).reindex(index=genes, columns=cells).fillna(0.0).to_numpy()
        self.assertTrue(np.all(removed <= np.floor(amount) + 1 + 1e-9))
        self.assertLessEqual(removed.sum(), amount.sum() + (removed > 0).sum())
        uncapped, _, _ = self._removed(self.uncapped)
        self.assertGreater(uncapped.sum(), removed.sum())
        self.assertGreater(removed.sum(), 0.0)
        self.assertTrue(np.all(removed <= uncapped + 1e-9))

    def test_capped_store_is_consistent(self):
        out = Path(self.capped.run_path)
        mol = pq.read_table(out / "molecules.parquet",
                            columns=["obs_id", "cell_idx", "factor_label"]).to_pandas()
        kept_uncapped = pq.read_table(Path(self.uncapped.run_path) / "molecules.parquet",
                                      columns=["obs_id"]).column(0).to_numpy()
        self.assertTrue(set(kept_uncapped) <= set(mol.obs_id))
        cells = pq.read_table(out / "cells.parquet").to_pandas()
        counts = mol.groupby("cell_idx").size().reindex(cells.cell_idx).fillna(0).to_numpy()
        self.assertTrue(np.array_equal(counts, cells.transcript_count.to_numpy()))
        nz = cells.transcript_count.to_numpy() > 0
        fr = cells[[c for c in cells.columns if c.endswith("_fraction") and c.startswith("factor_")]].to_numpy()
        # fractions are labeled molecules per factor over all the cell's
        # molecules (unlabeled ones count in the total only), as the native
        # writer computes them
        labeled = mol[mol.factor_label >= 0]
        expected = (pd.crosstab(labeled.cell_idx, labeled.factor_label)
                    .reindex(index=cells.cell_idx, columns=range(fr.shape[1])).fillna(0).to_numpy()
                    / np.maximum(counts, 1)[:, None])
        self.assertTrue(np.allclose(fr, expected))
        self.assertTrue(np.array_equal(cells.dominant_factor.to_numpy()[nz], fr[nz].argmax(axis=1) + 1))
        self.assertTrue(np.all(cells.dominant_factor.to_numpy()[~nz] == -1))
        summary = pd.read_parquet(self.capped.manifest["correction_summary_parquet"])
        self.assertEqual(int(summary.n_removed.sum()), int(self.capped.manifest["n_removed"]))
        self.assertEqual(int(summary.n_molecules_after.sum()), len(mol))
        self.assertEqual(int(self.capped.manifest["n_transcripts"]), len(mol))
        self.assertIn("cap", self.capped.manifest)

    def test_a_cap_above_every_removal_reproduces_the_uncapped_store(self):
        before, genes, cells = self.fit.counts()
        huge = (sparse.csc_matrix(np.full(before.shape, 1e6)), genes, cells)
        same = self.score.correct(self.rules, name="cap_huge", cap=huge)
        a = pq.read_table(Path(same.run_path) / "molecules.parquet", columns=["obs_id"]).column(0).to_numpy()
        b = pq.read_table(Path(self.uncapped.run_path) / "molecules.parquet", columns=["obs_id"]).column(0).to_numpy()
        self.assertTrue(np.array_equal(a, b))
        self.assertEqual(same.manifest["n_removed"], self.uncapped.manifest["n_removed"])
        self.assertEqual(same.manifest["cap"]["restored"], 0)

    def test_a_scale_relaxes_the_cap(self):
        relaxed = self.score.correct(self.rules, name="cap_scale2", cap=self.model, cap_scale=2.0)
        self.assertGreaterEqual(relaxed.manifest["n_removed"], self.capped.manifest["n_removed"])
        self.assertLessEqual(relaxed.manifest["n_removed"], self.uncapped.manifest["n_removed"])
        self.assertEqual(relaxed.manifest["cap"]["scale"], 2.0)

    def test_single_fit_correction_takes_the_cap(self):
        single = self.score.correct(self.rules, name="cap_single", cap=self.model, ensemble=1)
        removed, genes, cells = self._removed(single)
        self.assertGreater(removed.sum(), 0.0)
        self.assertIn("cap", single.manifest)


if __name__ == "__main__":
    unittest.main()
