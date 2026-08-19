"""End-to-end tests for the molecule-vote ensemble correction."""
import gzip
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


def make_bundle(root: Path, seed: int = 3):
    """Write a small synthetic Xenium bundle with two spatially mixed types."""
    rng = np.random.default_rng(seed)
    genes = {"A": [f"GA{i}" for i in range(5)], "B": [f"GB{i}" for i in range(5)]}
    cells = []
    molecules = []
    idx = 0
    for gx in range(8):
        for gy in range(8):
            cell_type = "A" if (gx + gy) % 2 == 0 else "B"
            cell_id = f"c{idx}"
            cx, cy = gx * 20.0 + 10.0, gy * 20.0 + 10.0
            cells.append((cell_id, cx, cy, 0.0, cell_type))
            pool = genes[cell_type] * 8 + genes["A" if cell_type == "B" else "B"][:2]
            for gene in pool:
                molecules.append((
                    f"t{len(molecules)}", cell_id, gene,
                    cx + rng.uniform(-8, 8), cy + rng.uniform(-8, 8),
                    rng.uniform(0, 2),
                ))
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
    return cell_df[["cell_id", "cell_type"]]


class EnsembleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from celladmix import CellAdmix

        cls.tmp = tempfile.TemporaryDirectory()
        root = Path(cls.tmp.name)
        annotation = make_bundle(root / "bundle")
        cls.ds = CellAdmix(
            str(root / "bundle"),
            output_dir=str(root / "out"),
            annotation=annotation.set_index("cell_id")["cell_type"],
            num_threads=2,
        )
        cls.fit = cls.ds.fit(
            rank=2, ncv_k=8, nmf_n_runs=4, nmf_iterations=60,
            nmf_train_max_rows=60, seed=3, verbose=False,
        )

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_member_pool_persisted_and_labeled(self):
        from celladmix import _core

        run_path = Path(self.fit.run_path)
        self.assertTrue((run_path / "ensemble_h.parquet").exists())
        count = _core.ensemble_prepare(str(run_path), 2)
        self.assertEqual(count, 4)
        for member in range(4):
            self.assertTrue((run_path / f"ensemble_labels_m{member}.parquet").exists())
        fractions = _core.ensemble_member_fractions(str(run_path), 0)
        matrix = np.asarray(fractions["fractions"]["data"], dtype=float)
        self.assertTrue(np.all(matrix >= 0.0) and np.all(matrix <= 1.0))

    def test_vote_thresholds_are_monotone(self):
        from celladmix import _core

        _core.ensemble_prepare(str(self.fit.run_path), 2)
        rule = pd.DataFrame({"factor": [1], "target_cell_type": ["B"]})
        rules4 = pd.concat([rule] * 4, ignore_index=True)

        removed = {}
        for min_votes in (1, 2, 4, 5):
            corr = self.fit.correct(
                rules4, name=f"vote{min_votes}",
                rule_members=[0, 1, 2, 3], min_votes=min_votes)
            removed[min_votes] = corr.manifest["n_removed"]
        self.assertGreaterEqual(removed[1], removed[2])
        self.assertGreaterEqual(removed[2], removed[4])
        self.assertEqual(removed[5], 0)

        selected = int(self.fit.manifest["nmf_diagnostics"]["selected_run"]) - 1
        primary = self.fit.correct(rule, name="vote_primary")
        via_member = self.fit.correct(
            rule, name="vote_selected", rule_members=[selected], min_votes=1)
        self.assertEqual(primary.manifest["n_removed"], via_member.manifest["n_removed"])

        union = self.fit.correct(
            rules4, name="vote_hist", rule_members=[0, 1, 2, 3], min_votes=1)
        self.assertEqual(sum(union.manifest["vote_histogram"]), union.manifest["n_removed"])

    def test_member_rules_disk_cache(self):
        rules = pd.DataFrame({
            "factor": [1], "target_cell_type": ["B"], "source_cell_type": ["A"]})
        kwargs = dict(candidate_k=5, crossing_k=5, min_type_pair_contacts=1,
            min_factor_molecules=1, min_pairs=1, null_iterations=1,
            null_max_iterations=2, compute_null=True, verbose=False)
        score = self.fit.score_bridge(**kwargs)
        c1 = score.correct(rules, p_thresh=0.9, name="cache_a")
        cache = sorted((Path(self.fit.run_path) / "scores").glob(
            "ensemble_rules_bridge_m*.pkl"))
        self.assertEqual(len(cache), 3)
        mtimes = [p.stat().st_mtime_ns for p in cache]

        # A fresh score object reuses the disk cache instead of re-scoring.
        score2 = self.fit.score_bridge(**kwargs)
        c2 = score2.correct(rules, p_thresh=0.9, name="cache_b")
        self.assertEqual(c1.manifest["n_removed"], c2.manifest["n_removed"])
        self.assertEqual(mtimes, [p.stat().st_mtime_ns for p in cache])

        # A changed rule threshold invalidates and rewrites the cache.
        score2.correct(rules, p_thresh=0.5, name="cache_c")
        self.assertNotEqual(mtimes, [p.stat().st_mtime_ns for p in cache])

    def test_overremoval_warnings(self):
        from celladmix import audit as audit_mod
        from celladmix import fit as fit_mod

        rules = pd.DataFrame({
            "factor": [1, 2], "target_cell_type": ["B", "B"],
            "source_cell_type": ["A", "A"]})
        old_fit = fit_mod.OVERREMOVAL_MIN_MOLECULES
        fit_mod.OVERREMOVAL_MIN_MOLECULES = 10
        try:
            with self.assertWarnsRegex(UserWarning, "erasing native expression"):
                corr = self.fit.correct(rules, name="erase_b")
        finally:
            fit_mod.OVERREMOVAL_MIN_MOLECULES = old_fit

        audit = self.fit.audit_admixture(min_target_cells=10, min_reference_cells=5)
        old_audit = audit_mod.OWN_MARKER_WARN_MIN
        audit_mod.OWN_MARKER_WARN_MIN = 1
        try:
            with self.assertWarnsRegex(UserWarning, "own-marker molecules"):
                audit.evaluate(corr, warn_uncovered=False)
        finally:
            audit_mod.OWN_MARKER_WARN_MIN = old_audit

    def test_score_correct_defaults_to_ensemble(self):
        score = self.fit.score_bridge(
            candidate_k=5, crossing_k=5, min_type_pair_contacts=1,
            min_factor_molecules=1, min_pairs=1, null_iterations=1,
            null_max_iterations=2, compute_null=True, verbose=False,
        )
        detected = score.rules(p_thresh=0.9)
        if detected.empty or not detected.get("keep", pd.Series(dtype=bool)).fillna(True).any():
            # Detection is not guaranteed on the tiny synthetic bundle;
            # explicit rules still exercise the full ensemble path.
            detected = pd.DataFrame({
                "factor": [1],
                "target_cell_type": ["B"],
                "source_cell_type": ["A"],
            })
        correction = score.correct(detected, p_thresh=0.9, name="ens_default")
        self.assertEqual(correction.manifest["min_votes"], 2)
        self.assertIn("support", correction.rules.columns)
        self.assertTrue(((correction.rules["support"] >= 0) & (correction.rules["support"] <= 1)).all())
        single = score.correct(detected, p_thresh=0.9, name="ens_single", ensemble=1)
        self.assertEqual(single.ensemble()["members"], 1)


if __name__ == "__main__":
    unittest.main()
