"""Native-check vetoes propagate to heatmaps and example selection."""
import unittest

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd


def make_rules(keep):
    return pd.DataFrame(
        {
            "factor": [1, 1, 2, 3],
            "source_cell_type": ["A", "B", "A", "A"],
            "target_cell_type": ["T", "T", "T", "U"],
            "keep": keep,
        }
    )


class VetoedCallsTests(unittest.TestCase):
    def test_group_level_veto(self):
        from celladmix._score_utils import vetoed_calls_from_rules

        # Factor 1 on T keeps one surviving rule; factor 2 on T loses all of
        # its rules; the NA on factor 3 counts as kept (matching correction).
        out = vetoed_calls_from_rules(make_rules([False, True, False, np.nan]))
        self.assertEqual(out, {(2, "T")})

        out = vetoed_calls_from_rules(make_rules([False] * 4))
        self.assertEqual(out, {(1, "T"), (2, "T"), (3, "U")})

        self.assertIsNone(vetoed_calls_from_rules(make_rules([True] * 4)))
        self.assertIsNone(vetoed_calls_from_rules(None))
        self.assertIsNone(vetoed_calls_from_rules(make_rules([True] * 4).iloc[:0]))
        self.assertIsNone(
            vetoed_calls_from_rules(make_rules([True] * 4).drop(columns="keep")))


class HeatmapVetoTests(unittest.TestCase):
    def make_summary(self):
        return pd.DataFrame(
            {
                "target_cell_type": ["T", "T"],
                "source_cell_type": ["S", "S"],
                "factor": [1, 2],
                "p_value": [1e-4, 1e-4],
                "neg_log10_p": [4.0, 4.0],
            }
        )

    def labels(self, ax):
        return [t.get_text() for t in ax.texts]

    def test_vetoed_removal_marked(self):
        from celladmix.plotting import plot_score_heatmap

        ax = plot_score_heatmap(self.make_summary(), vetoed_calls={(2, "T")})
        labels = self.labels(ax)
        self.assertIn("*", labels)
        self.assertIn("(*)", labels)
        self.assertTrue(any("native-factor check" in t for t in labels))

        ax = plot_score_heatmap(self.make_summary())
        labels = self.labels(ax)
        self.assertIn("*", labels)
        self.assertNotIn("(*)", labels)
        self.assertFalse(any("native-factor check" in t for t in labels))


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


class ExampleSelectionVetoTests(unittest.TestCase):
    def make_score(self):
        pairs = pd.DataFrame(
            {
                "target_cell": ["c1", "c2"],
                "target_cell_type": ["T", "T"],
                "factor": [1, 2],
                "mean_score": [0.6, 0.5],
                "factor_count": [5.0, 5.0],
            }
        )
        cell_factors = pd.DataFrame(
            {
                "cell_id": ["c1", "c2"],
                "transcript_count": [100.0, 100.0],
                "dominant_factor": [9, 9],
            }
        )
        annotation = pd.Series({"c1": "T", "c2": "T"})
        return StubScore(pairs, cell_factors, annotation)

    def select(self, keep):
        from celladmix.examples import select_example_cells

        rules = pd.DataFrame(
            {
                "factor": [1, 2],
                "source_cell_type": ["S", "S"],
                "target_cell_type": ["T", "T"],
                "keep": keep,
                "native_check": ["pass", "pass"],
            }
        )
        return select_example_cells(
            self.make_score(),
            rules=rules,
            score_annotation={"source_calls": {1: "S", 2: "S"}},
            min_molecules=0,
            min_factor_molecules=0,
        )

    def test_vetoed_combination_not_presented(self):
        out = self.select([True, False])
        # Factor 2's rule is vetoed: c2, whose only evidence is factor 2,
        # is not presented as an admixture example; c1 survives on factor 1.
        self.assertEqual(out["target_cell"].tolist(), ["c1"])
        self.assertEqual(int(out.loc[0, "top_admix_factor"]), 1)

        out = self.select([True, True])
        self.assertEqual(sorted(out["target_cell"]), ["c1", "c2"])


if __name__ == "__main__":
    unittest.main()
