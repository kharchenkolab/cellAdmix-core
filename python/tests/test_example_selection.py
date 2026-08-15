"""Tests for explicit-cell example selection and cell-type resolution."""
import unittest

import numpy as np
import pandas as pd


class StubFitDataset:
    def __init__(self, annotation):
        self.annotation = annotation
        self.source = None


class StubFit:
    def __init__(self, cell_factors, annotation):
        self._cell_factors = cell_factors
        self.dataset = StubFitDataset(annotation)

    def cell_factors(self):
        return self._cell_factors


class StubScore:
    def __init__(self, pairs, cell_factors, annotation):
        self.pairs = pairs
        self.fit = StubFit(cell_factors, annotation)

    def annotation(self, **kwargs):
        return {"source_calls": {}}

    def rules(self, **kwargs):
        return pd.DataFrame()


def make_score():
    pairs = pd.DataFrame(
        {
            "target_cell": ["c1", "c1", "c2"],
            "target_cell_type": ["T", "T", "B"],
            "factor": [1, 2, 1],
            "mean_score": [0.5, 0.2, -0.1],
            "factor_count": [2.0, 1.0, 1.0],
            "used_in_summary": [True, False, False],
        }
    )
    cell_factors = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3"],
            "transcript_count": [10.0, 20.0, 30.0],
            "dominant_factor": [3, 3, 3],
        }
    )
    annotation = pd.Series({"c1": "T", "c2": "B", "c3": "NK"})
    return StubScore(pairs, cell_factors, annotation)


class ExplicitCellSelectionTests(unittest.TestCase):
    def test_requested_cells_return_in_order_with_fallback(self):
        from celladmix.examples import select_example_cells

        out = select_example_cells(make_score(), cells=["c3", "c1"])
        self.assertEqual(out["target_cell"].tolist(), ["c3", "c1"])
        # c3 has no score evidence: NA factor columns, type from annotation.
        self.assertEqual(out.loc[0, "target_cell_type"], "NK")
        self.assertTrue(pd.isna(out.loc[0, "top_admix_factor"]))
        self.assertEqual(float(out.loc[0, "transcript_count"]), 30.0)
        # c1 keeps its best-evidence factor even though filters would normally
        # drop it (low molecule counts, used_in_summary False rows allowed).
        self.assertEqual(int(out.loc[1, "top_admix_factor"]), 1)

    def test_negative_score_cell_still_returns(self):
        from celladmix.examples import select_example_cells

        out = select_example_cells(make_score(), cells="c2")
        self.assertEqual(len(out), 1)
        self.assertEqual(int(out.loc[0, "top_admix_factor"]), 1)
        self.assertEqual(float(out.loc[0, "top_admix_score"]), -0.1)

    def test_automatic_selection_unchanged(self):
        from celladmix.examples import select_example_cells

        out = select_example_cells(make_score())
        # Only c1/F1 passes used_in_summary + positive-score gating; molecule
        # minimums then drop it (transcript_count 10 < 50).
        self.assertTrue(out.empty)
        out = select_example_cells(make_score(), min_molecules=0, min_factor_molecules=0)
        self.assertEqual(out["target_cell"].tolist(), ["c1"])


class CellTypeResolutionTests(unittest.TestCase):
    def test_resolve_cell_types(self):
        from celladmix.examples import resolve_example_cell_types

        annotation = pd.Series({"c1": "T"})
        fit = StubFit(pd.DataFrame(), annotation)
        self.assertIs(resolve_example_cell_types(fit, "auto"), annotation)
        self.assertIsNone(resolve_example_cell_types(fit, None))
        frame = pd.DataFrame({"cell_id": ["c1"], "cell_type": ["T"]})
        resolved = resolve_example_cell_types(fit, frame)
        self.assertEqual(resolved.get("c1"), "T")
        with self.assertRaises(ValueError):
            resolve_example_cell_types(fit, "bogus")


if __name__ == "__main__":
    unittest.main()
