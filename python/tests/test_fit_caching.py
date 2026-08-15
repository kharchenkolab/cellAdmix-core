"""Tests for parameter-aware fit caching helpers."""
import unittest

import pandas as pd


def make_manifest():
    return {
        "pipeline_options": {
            "rank": 8, "nmf_variant": "ls_nmf", "nmf_init": "auto",
            "molecule_scoring": "gene_loadings", "ncv_k": 71,
            "nmf_iterations": 150, "seed": 1,
        },
        "annotation_hash": "abc",
        "package_version": "0.0.1",
    }


class FitParamDiffTests(unittest.TestCase):
    def test_matching_request_and_auto_params(self):
        from celladmix.dataset import _fit_param_diff

        requested = {"rank": 8, "nmf_variant": "ls_nmf", "nmf_init": "auto",
                     "molecule_scoring": "gene_loadings"}
        self.assertEqual(_fit_param_diff(make_manifest(), requested, "abc"), [])

    def test_explicit_parameter_change_detected(self):
        from celladmix.dataset import _fit_param_diff

        requested = {"rank": 8, "nmf_variant": "ls_nmf", "nmf_init": "auto",
                     "molecule_scoring": "gene_loadings", "ncv_k": 20}
        diff = _fit_param_diff(make_manifest(), requested, "abc")
        self.assertEqual(diff, ["ncv_k: 71 -> 20"])

    def test_rank_variant_and_annotation_changes_detected(self):
        from celladmix.dataset import _fit_param_diff

        base = {"nmf_variant": "ls_nmf", "nmf_init": "auto",
                "molecule_scoring": "gene_loadings"}
        self.assertTrue(any("rank" in d for d in _fit_param_diff(
            make_manifest(), dict(base, rank=9), "abc")))
        self.assertTrue(any("nmf_variant" in d for d in _fit_param_diff(
            make_manifest(), dict(base, rank=8, nmf_variant="invsqrt_kl"), "abc")))
        self.assertEqual(_fit_param_diff(make_manifest(), dict(base, rank=8), "zzz"),
                         ["annotation content"])
        # Old manifests without a recorded hash never trigger annotation diffs.
        manifest = make_manifest()
        manifest["annotation_hash"] = ""
        self.assertEqual(_fit_param_diff(manifest, dict(base, rank=8), "zzz"), [])

    def test_annotation_hash_is_order_independent(self):
        from celladmix.dataset import _annotation_hash

        a = pd.Series({"c1": "T", "c2": "B"})
        b = pd.Series({"c2": "B", "c1": "T"})
        self.assertEqual(_annotation_hash(a), _annotation_hash(b))
        self.assertNotEqual(_annotation_hash(a), _annotation_hash(pd.Series({"c1": "T", "c2": "NK"})))
        self.assertEqual(_annotation_hash(None), "")


if __name__ == "__main__":
    unittest.main()
