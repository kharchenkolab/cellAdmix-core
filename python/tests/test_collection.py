"""Consensus factorization across two synthetic samples."""
import gzip
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


def make_bundle(root: Path, seed: int):
    """Two spatially mixed cell types on a grid, with planted cross-type
    molecules; the same cell-id sequence in every sample."""
    rng = np.random.default_rng(seed)
    genes = {"A": [f"GA{i}" for i in range(5)], "B": [f"GB{i}" for i in range(5)]}
    cells, molecules, idx = [], [], 0
    for gx in range(8):
        for gy in range(8):
            cell_type = "A" if (gx + gy) % 2 == 0 else "B"
            cell_id = f"c{idx}"
            cx, cy = gx * 20.0 + 10.0, gy * 20.0 + 10.0
            cells.append((cell_id, cx, cy, 0.0, cell_type))
            pool = genes[cell_type] * 8 + genes["A" if cell_type == "B" else "B"][:2]
            for gene in pool:
                molecules.append((f"t{len(molecules)}", cell_id, gene,
                                  cx + rng.uniform(-8, 8), cy + rng.uniform(-8, 8), rng.uniform(0, 2)))
            idx += 1
    root.mkdir(parents=True, exist_ok=True)
    (root / "experiment.xenium").write_text(json.dumps({
        "run_name": "Mock Xenium", "pixel_size": 0.2125, "z_step_size": 3.0}))
    tx = pd.DataFrame(molecules, columns=["transcript_id", "cell_id", "feature_name",
                                          "x_location", "y_location", "z_location"])
    tx["qv"] = 30.0
    tx["overlaps_nucleus"] = (np.arange(len(tx)) % 3 == 0).astype(int)
    tx["nucleus_distance"] = (np.arange(len(tx)) % 11).astype(float)
    with gzip.open(root / "transcripts.csv.gz", "wt") as handle:
        tx.to_csv(handle, index=False)
    cell_df = pd.DataFrame(cells, columns=["cell_id", "x_centroid", "y_centroid", "z_centroid", "cell_type"])
    with gzip.open(root / "cells.csv.gz", "wt") as handle:
        cell_df.to_csv(handle, index=False)
    return cell_df[["cell_id", "cell_type"]]


FIT = dict(rank=2, ncv_k=6, graph_k=4, nmf_iterations=60, nmf_train_max_rows=120,
           nmf_n_runs=3, seed=5, verbose=False)


class CollectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from celladmix import CellAdmix, CellAdmixCollection

        cls.tmp = tempfile.TemporaryDirectory()
        root = Path(cls.tmp.name)
        cls.datasets = {}
        for name, seed in (("a", 11), ("b", 12)):
            annotation = make_bundle(root / name, seed)
            cls.datasets[name] = CellAdmix(root / name, output_dir=root / f"out_{name}",
                                           annotation=annotation, num_threads=2)
        cls.collection = CellAdmixCollection(cls.datasets)
        cls.collection.fit(**FIT)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_fixed_profiles_keep_order_and_carry_members(self):
        from celladmix.collection import _fit_profiles_raw

        fit = self.collection.fits["a"]
        profiles = _fit_profiles_raw(fit)
        reversed_profiles = profiles.iloc[::-1].reset_index(drop=True)
        members = [reversed_profiles * 1.1, reversed_profiles * 0.9]
        fixed = self.datasets["a"].fit(nmf_fixed_h=reversed_profiles, nmf_fixed_h_members=members,
                                        run_id="fixed_members", overwrite=True, nmf_n_runs=1,
                                        ncv_k=6, graph_k=4, seed=5, verbose=False)
        self.assertTrue(fixed.manifest["fixed_profiles"])
        self.assertEqual(fixed.manifest["fixed_profile_members"], 3)
        loadings = fixed.factor_loadings().T
        loadings.columns = [str(c) for c in loadings.columns]
        np.testing.assert_allclose(loadings[reversed_profiles.columns].to_numpy(),
                                   reversed_profiles.to_numpy(), rtol=1e-6, atol=1e-9)
        from celladmix import _core

        self.assertEqual(_core.ensemble_prepare(str(fixed.run_path), 1), 3)
        self.assertEqual(fixed.manifest["nmf_diagnostics"]["stability_comparison_runs"], 2)

    def test_matching_and_training_loss(self):
        from celladmix.collection import _fit_profiles_raw, _match_factors, _training_loss

        fit = self.collection.fits["a"]
        profiles = _fit_profiles_raw(fit)
        match = _match_factors(profiles.iloc[::-1], profiles)
        self.assertEqual(match["match"], [1, 0])
        self.assertTrue(all(c > 0.99 for c in match["correlation"]))
        own = _training_loss(str(fit.run_path), profiles)
        flat = pd.DataFrame(np.ones(profiles.shape), columns=profiles.columns)
        self.assertGreater(_training_loss(str(fit.run_path), flat), own)

    def test_consensus_end_to_end(self):
        from celladmix import (CellAdmixCollectionAudit, CellAdmixCollectionCorrection,
                               CellAdmixCollectionScore, CellAdmixConsensusFit)

        consensus = self.collection.fit_consensus()
        self.assertIsInstance(consensus, CellAdmixConsensusFit)
        self.assertEqual(consensus.samples(), ["a", "b"])
        profiles = consensus.loadings()
        self.assertTrue(2 <= len(profiles) <= 4)
        self.assertEqual(len(consensus.matching), 4)
        self.assertIn("own_stability", consensus.matching.columns)
        self.assertTrue(np.isfinite(consensus.matching["own_stability"]).all())
        self.assertEqual(len(consensus.loss), 2)
        self.assertTrue(np.isfinite(consensus.loss["consensus_loss"]).all())
        for fit in consensus.fits.values():
            self.assertTrue(fit.manifest["fixed_profiles"])
            self.assertEqual(int(fit.manifest["n_factors"]), len(profiles))
            self.assertEqual(fit.manifest["fixed_profile_members"], 7)
        stability = consensus.stability()
        self.assertEqual(list(stability["sample"]), ["a", "b"])
        usage = consensus.sample_usage()
        self.assertTrue({"sample", "cell_type", "factor", "share"} <= set(usage.columns))
        divergence = consensus.usage_divergence()
        self.assertTrue(((divergence["divergence"] >= 0) & (divergence["divergence"] <= 1)).all())

        cells = consensus.cells()
        self.assertTrue(all(c.startswith(s + "_") for c, s in zip(cells["cell_id"], cells["sample"])))
        matrix, genes, cell_ids = consensus.counts()
        self.assertEqual(matrix.shape[1], len(cells))
        self.assertTrue(all(c.startswith(("a_", "b_")) for c in cell_ids))

        score = consensus.score_bridge(candidate_pairs_per_type_pair=40, max_cells_per_type_pair=40,
                                       null_iterations=1, verbose=False)
        self.assertIsInstance(score, CellAdmixCollectionScore)
        pooled = score.summary
        self.assertTrue({"n_pairs", "mean_score", "p_value", "n_samples"} <= set(pooled.columns))
        self.assertTrue(pooled["n_samples"].between(1, 2).all())
        rules = score.rules(p_thresh=0.5)
        self.assertTrue({"factor", "source_cell_type", "target_cell_type", "sample_support",
                         "keep", "native_check"} <= set(rules.columns))
        self.assertTrue(rules["sample_support"].between(0, 1).all())
        if rules.empty:
            # The grid fixture carries no bridge evidence; an explicit shared
            # rule still exercises the per-sample ensemble correction.
            rules = pd.DataFrame({"factor": [1], "source_cell_type": ["A"],
                                  "target_cell_type": ["B"]})

        correction = score.correct(rules=rules, vote=0.3)
        self.assertIsInstance(correction, CellAdmixCollectionCorrection)
        self.assertEqual(correction.samples(), ["a", "b"])
        self.assertIn("support", correction.rules.columns)
        self.assertTrue(correction.rules["support"].between(0, 1).all())
        self.assertIn("sample", correction.summary().columns)
        self.assertIn("fraction_molecules_removed", correction.summary(by_sample=False).columns)
        corrected, _, corrected_cells = correction.counts()
        self.assertEqual(corrected.shape, matrix.shape)
        self.assertLessEqual(corrected.sum(), matrix.sum())
        self.assertEqual(correction.counts(sample="a")[0].shape, consensus.counts(sample="a")[0].shape)

        audit = consensus.audit_admixture(min_excess=1, min_target_cells=2, min_reference_cells=2)
        self.assertIsInstance(audit, CellAdmixCollectionAudit)
        self.assertIn("sample", audit.pairs(by_sample=True).columns)
        self.assertTrue({"admixed_molecules", "rate", "n_samples"} <= set(audit.pairs().columns))
        report = audit.evaluate(correction)
        self.assertEqual(report.samples(), ["a", "b"])
        self.assertIn("sample", report.summary().columns)

        again = self.collection.fit_consensus()
        self.assertEqual(again.params["run_id"], consensus.params["run_id"])

        # Projection onto given profiles: sample a's own profiles, with its
        # restart pool as the shared member pool.
        from celladmix.collection import _fit_member_profiles, _fit_profiles_raw

        own = _fit_profiles_raw(self.collection.fits["a"])
        projected = self.collection.fit_profiles(own, members=_fit_member_profiles(self.collection.fits["a"]),
                                                 name="reference_a")
        self.assertEqual(projected.params["source"], "profiles")
        self.assertEqual(len(projected.loadings()), len(own))
        self.assertEqual(len(projected.matching), 4)
        self.assertTrue(np.isfinite(projected.loss["loss_gap"]).all())
        for fit in projected.fits.values():
            self.assertEqual(fit.manifest["fixed_profile_members"], 4)


if __name__ == "__main__":
    unittest.main()
