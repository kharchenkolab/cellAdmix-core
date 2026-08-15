"""Tests for Xenium stain discovery, crop reading, and auto-resolution."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np


def make_bundle(root: Path, *, split_files: bool = True) -> Path:
    """Write a minimal Xenium bundle with tiny stain images.

    split_files=True mimics the multimodal layout (one file per channel with
    shared OME-XML); False keeps only the base focus file with two channels so
    discovery must fall back to channel selection.
    """
    import tifffile

    bundle = root / "bundle"
    focus_dir = bundle / "morphology_focus"
    focus_dir.mkdir(parents=True)
    manifest = {
        "images": {"morphology_focus_filepath": "morphology_focus/morphology_focus_0000.ome.tif"},
        "pixel_size": 0.5,
    }
    (bundle / "experiment.xenium").write_text(json.dumps(manifest))

    dapi = np.full((64, 64), 10, dtype=np.uint16)
    membrane = np.full((64, 64), 200, dtype=np.uint16)
    # Shared OME-XML naming both files as channels of one image, so tifffile
    # with default settings would aggregate them into one series.
    if split_files:
        tifffile.imwrite(focus_dir / "morphology_focus_0000.ome.tif", dapi)
        tifffile.imwrite(focus_dir / "morphology_focus_0001.ome.tif", membrane)
    else:
        stack = np.stack([dapi, membrane])
        tifffile.imwrite(focus_dir / "morphology_focus_0000.ome.tif", stack)
    return bundle


class StainDiscoveryTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def test_discovery_resolves_split_channel_files(self):
        from celladmix.io import discover_xenium_stain_image

        bundle = make_bundle(self.root)
        dapi = discover_xenium_stain_image(bundle, "dapi")
        membrane = discover_xenium_stain_image(bundle, "membrane")
        self.assertTrue(dapi["image_path"].endswith("morphology_focus_0000.ome.tif"))
        self.assertTrue(membrane["image_path"].endswith("morphology_focus_0001.ome.tif"))
        self.assertIsNone(dapi["channel"])
        self.assertIsNone(membrane["channel"])

    def test_crops_read_distinct_channels(self):
        from celladmix.io import discover_xenium_stain_image
        from celladmix.examples import read_stain_crop

        bundle = make_bundle(self.root)
        bbox = (0.0, 16.0, 0.0, 16.0)
        dapi = read_stain_crop(discover_xenium_stain_image(bundle, "dapi"), bbox)
        membrane = read_stain_crop(discover_xenium_stain_image(bundle, "membrane"), bbox)
        self.assertEqual(float(dapi["values"].mean()), 10.0)
        self.assertEqual(float(membrane["values"].mean()), 200.0)

    def test_multichannel_fallback_selects_channel(self):
        from celladmix.io import discover_xenium_stain_image
        from celladmix.examples import read_stain_crop

        bundle = make_bundle(self.root, split_files=False)
        membrane = discover_xenium_stain_image(bundle, "membrane")
        self.assertTrue(membrane["image_path"].endswith("morphology_focus_0000.ome.tif"))
        self.assertEqual(membrane["channel"], 1)
        crop = read_stain_crop(membrane, (0.0, 16.0, 0.0, 16.0))
        self.assertEqual(float(crop["values"].mean()), 200.0)

    def test_discover_images_skips_missing_and_non_xenium(self):
        from celladmix.io import discover_xenium_stain_images

        bundle = make_bundle(self.root)
        out = discover_xenium_stain_images(bundle)
        self.assertEqual(sorted(out), ["dapi", "membrane"])
        self.assertEqual(discover_xenium_stain_images(self.root / "nowhere"), {})
        self.assertEqual(discover_xenium_stain_images(None), {})

    def test_resolve_example_stains_semantics(self):
        from celladmix.examples import resolve_example_stains

        bundle = make_bundle(self.root)

        class FakeDataset:
            source = bundle

        class FakeFit:
            dataset = FakeDataset()

        fit = FakeFit()
        auto = resolve_example_stains(fit, "auto")
        self.assertEqual(sorted(auto), ["dapi", "membrane"])
        self.assertEqual(resolve_example_stains(fit, None), {})
        self.assertEqual(resolve_example_stains(fit, {}), {})
        explicit = {"membrane": {"image_path": "x.tif", "pixel_size": 1.0}}
        self.assertEqual(resolve_example_stains(fit, explicit), explicit)
        single = {"image_path": "x.tif", "pixel_size": 1.0, "stain": "dapi"}
        self.assertEqual(resolve_example_stains(fit, single), {"dapi": single})
        with self.assertRaises(ValueError):
            resolve_example_stains(fit, "bogus")

        class NoSourceDataset:
            source = None

        class NoSourceFit:
            dataset = NoSourceDataset()

        self.assertEqual(resolve_example_stains(NoSourceFit(), "auto"), {})

    def test_compose_uses_palettes_not_grey(self):
        from celladmix.examples import compose_stain_background

        base = np.zeros((4, 4))
        base[2:, 2:] = 100.0
        crops = {
            "dapi": {"values": base},
            "membrane": {"values": np.zeros((4, 4))},
        }
        bg = compose_stain_background(crops)
        px = bg[3, 3]
        # DAPI-only signal keeps the blue-leaning palette (no grey collapse).
        self.assertGreater(px[2], px[0])
        self.assertGreater(px[0], px[1])


if __name__ == "__main__":
    unittest.main()
