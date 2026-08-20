"""Unit tests for the audit's ambient reference estimation and the
induced-gene screening of marker panels."""
import unittest

import numpy as np
from scipy import sparse


class ReferenceRateTests(unittest.TestCase):
    def test_reference_uses_largest_populated_neighborhood(self):
        from celladmix.audit import _reference_rate

        rng = np.random.default_rng(5)
        n = 4000
        totals = np.full(n, 100.0)
        idx = np.arange(n)
        zero_by_K = {"k15": np.full(n, True), "k60": idx > 1000,
                     "k240": idx > 3000}
        rate = np.where(idx <= 1000, 0.03, np.where(idx <= 3000, 0.012, 0.005))
        markers = rng.poisson(rate * totals).astype(float)
        ref, kind, pooled = _reference_rate(markers, totals, zero_by_K)
        self.assertEqual(kind, "k240")
        self.assertLess(ref, 0.008)
        self.assertGreater(pooled, ref)

    def test_underpopulated_neighborhood_falls_back(self):
        from celladmix.audit import _reference_rate

        rng = np.random.default_rng(6)
        n = 4000
        totals = np.full(n, 100.0)
        idx = np.arange(n)
        zero_by_K = {"k15": np.full(n, True), "k60": idx > 1000,
                     "k240": idx > 3990}
        markers = rng.poisson(0.01 * totals).astype(float)
        _, kind, _ = _reference_rate(markers, totals, zero_by_K)
        self.assertEqual(kind, "k60")

    def test_flat_profile_keeps_base_rate(self):
        from celladmix.audit import _reference_rate

        rng = np.random.default_rng(7)
        n = 4000
        totals = np.full(n, 100.0)
        idx = np.arange(n)
        zero_by_K = {"k15": np.full(n, True), "k60": idx > 1000,
                     "k240": idx > 3000}
        markers = rng.poisson(0.01 * totals).astype(float)
        ref, _, pooled = _reference_rate(markers, totals, zero_by_K)
        self.assertAlmostEqual(ref, pooled, delta=0.15 * pooled)


class ScreenPanelTests(unittest.TestCase):
    def test_induced_gene_excluded_and_replaced(self):
        from celladmix.audit import _screen_panel

        rng = np.random.default_rng(9)
        n_cells = 400
        expo = np.tile([0, 1, 2, 3], n_cells // 4)
        # Genes 0..5: transferred material proportional to the source
        # profile; gene 6: induced (large excess, tiny profile share).
        psi = np.array([600, 500, 400, 300, 200, 100, 30.0])
        rows = [rng.poisson(0.02 * psi[g] * expo + 1) for g in range(6)]
        rows.append(rng.poisson(40 * expo + 1))
        # Background expression independent of exposure, as carried by the
        # rest of the transcriptome in real cells.
        rows.append(rng.poisson(np.full(n_cells, 500.0)))
        matrix = sparse.csc_matrix(np.vstack(rows).astype(float))
        pool, induced = _screen_panel(matrix, [6, 0, 1, 2, 3, 4, 5],
            np.flatnonzero(expo > 0), np.flatnonzero(expo == 0), psi,
            n_pool=6)
        self.assertIn(6, induced.tolist())
        self.assertNotIn(6, pool.tolist())
        self.assertEqual(sorted(pool.tolist()), [0, 1, 2, 3, 4, 5])


if __name__ == "__main__":
    unittest.main()
