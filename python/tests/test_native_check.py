"""Tests for the native-factor false-positive check on correction rules."""
import unittest

import numpy as np
import pandas as pd


def make_fit():
    """Source cluster at the origin; target cells split into an exposed group
    near the source and a distant group far away. Factor 1 is native to the
    target (expressed everywhere); factor 2 is true admixture (expressed only
    in exposed target cells)."""
    rng = np.random.default_rng(0)
    rows = []
    for i in range(20):
        rows.append(("S%d" % i, rng.uniform(0, 10), rng.uniform(0, 10), "Source", 0.05, 0.0))
    for i in range(20):
        rows.append(("TE%d" % i, rng.uniform(5, 15), rng.uniform(0, 10), "Target", 0.30, 0.40))
    for i in range(20):
        rows.append(("TD%d" % i, rng.uniform(100, 110), rng.uniform(0, 10), "Target", 0.30, 0.0))
    cells = pd.DataFrame(rows, columns=["cell_id", "x", "y", "cell_type",
                                        "factor_1_fraction", "factor_2_fraction"])
    annotation = pd.Series(cells["cell_type"].to_numpy(), index=cells["cell_id"])

    class Dataset:
        pass

    class Fit:
        pass

    fit = Fit()
    fit.dataset = Dataset()
    fit.dataset.annotation = annotation
    fit.cell_factors = lambda: cells.drop(columns=["cell_type"])
    return fit


def make_rules():
    return pd.DataFrame(
        {
            "factor": [1, 2, 1],
            "source_cell_type": ["Source", "Source", "Missing"],
            "target_cell_type": ["Target", "Target", "Target"],
            "p_value": [0.01, 0.01, 0.01],
            "neg_log10_p": [2.0, 2.0, 2.0],
            "rule_id": ["1_Target", "2_Target", "1_Target"],
        }
    )


class NativeCheckTests(unittest.TestCase):
    def test_native_factor_vetoed_and_admixture_kept(self):
        from celladmix._score_utils import apply_native_check

        out = apply_native_check(make_rules(), make_fit(), neighbor_k=10)
        self.assertEqual(out["native_check"].tolist()[:2], ["native_median", "pass"])
        self.assertEqual(out["keep"].tolist(), [False, True, False])
        # Vetoed native factor: distant target cells still carry ~0.3.
        self.assertGreater(out.loc[0, "native_distant_median"], 0.1)
        # Kept admixture rule: strong positive exposure gradient.
        self.assertGreater(out.loc[1, "native_exposure_gradient"], 0.2)
        # Unresolvable source type is excluded to be safe.
        self.assertEqual(out.loc[2, "native_check"], "source_not_in_annotation")

    def test_empty_rules_gain_columns(self):
        from celladmix._score_utils import apply_native_check

        out = apply_native_check(make_rules().iloc[:0], make_fit())
        self.assertIn("keep", out.columns)
        self.assertEqual(len(out), 0)


if __name__ == "__main__":
    unittest.main()
