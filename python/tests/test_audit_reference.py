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
        ref, kind, pooled, _ = _reference_rate(markers, totals, zero_by_K)
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
        _, kind, _, _ = _reference_rate(markers, totals, zero_by_K)
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
        ref, _, pooled, _ = _reference_rate(markers, totals, zero_by_K)
        self.assertAlmostEqual(ref, pooled, delta=0.15 * pooled)


class ScreenPanelTests(unittest.TestCase):
    def _screen_case(self, profile_cv=None):
        from celladmix.audit import _screen_panel

        rng = np.random.default_rng(9)
        n_cells = 400
        expo = np.tile([0, 1, 2, 3], n_cells // 4)
        # Genes 0..9: transferred material proportional to the source
        # profile; gene 10: induced (large excess, tiny profile share);
        # gene 11: dominant transfer channel deviating ~15% from
        # proportionality - within profile uncertainty, although its
        # counting-noise significance is large.
        psi = np.array(list(range(600, 100, -50)) + [30.0, 5000.0])
        rows = [rng.poisson(0.02 * psi[g] * expo + 1) for g in range(10)]
        rows.append(rng.poisson(40 * expo + 1))
        rows.append(rng.poisson(0.02 * 5000 * expo * 1.15 + 5))
        # Background expression independent of exposure, as carried by the
        # rest of the transcriptome in real cells.
        rows.append(rng.poisson(np.full(n_cells, 500.0)))
        matrix = sparse.csc_matrix(np.vstack(rows).astype(float))
        kwargs = {} if profile_cv is None else dict(profile_cv=profile_cv)
        return _screen_panel(matrix, [11, 10] + list(range(10)),
            np.flatnonzero(expo > 0), np.flatnonzero(expo == 0), psi,
            n_pool=11, **kwargs)

    def test_induced_gene_excluded_and_replaced(self):
        pool, induced, stats = self._screen_case()
        self.assertIn(10, induced.tolist())
        self.assertNotIn(10, pool.tolist())
        self.assertNotIn(11, induced.tolist())
        self.assertEqual(sorted(pool.tolist()), sorted([11] + list(range(10))))
        fold = stats.loc[stats["gene"] == 10, "fold"].iloc[0]
        self.assertGreater(fold, 4)

    def test_profile_uncertainty_shields_large_channels(self):
        # Without the profile-uncertainty term the large channel's small
        # relative deviation becomes formally significant.
        _, induced, _ = self._screen_case(profile_cv=0.0)
        self.assertIn(11, induced.tolist())

    def test_interface_profile_uses_bordering_cells(self):
        from celladmix.audit import _interface_profile

        matrix = sparse.csc_matrix(
            np.array([[30000.0, 0.0], [10000.0, 40000.0]]))
        prof = _interface_profile(matrix, [0, 1], [10.0, 500.0])
        self.assertAlmostEqual(prof[0], 0.75)
        # A sparse near subset widens until enough molecules are available.
        matrix2 = sparse.csc_matrix(np.array([[75.0, 0.0], [25.0, 40000.0]]))
        prof2 = _interface_profile(matrix2, [0, 1], [10.0, 500.0])
        self.assertLess(prof2[0], 0.01)


if __name__ == "__main__":
    unittest.main()
