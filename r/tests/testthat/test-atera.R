library(cellAdmixCore)

# A minimal Atera bundle layout: the manifest plus placeholder Zarr store files.
write_atera_placeholder <- function(dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  writeLines('{"run_name": "Mock Atera"}', file.path(dir, "experiment.spatial"))
  for (name in c("transcripts.zarr.zip", "cells.zarr.zip", "cell_feature_matrix.zarr.zip")) {
    writeLines("placeholder", file.path(dir, name))
  }
  dir
}

# Python interpreter that can run the conversion, or "" when none is configured.
atera_test_python <- function() {
  python <- Sys.getenv("CELLADMIX_PYTHON", "")
  if (!nzchar(python)) {
    return("")
  }
  ok <- suppressWarnings(system2(python, c("-c", shQuote("import zarr, celladmix.atera")),
    stdout = FALSE, stderr = FALSE))
  if (identical(as.integer(ok), 0L)) python else ""
}

test_that("Atera bundles are detected and read as Xenium-layout input", {
  td <- tempfile("atera_detect_")
  bundle <- write_atera_placeholder(file.path(td, "outs"))
  expect_equal(cellAdmixCore:::.celladmix_infer_format(bundle), "xenium")
  expect_equal(cellAdmixCore:::.celladmix_infer_format(file.path(bundle, "experiment.spatial")), "xenium")
  expect_false(is.null(cellAdmixCore:::.celladmix_atera_bundle_dir(bundle)))
  # A directory with both manifests is a Xenium bundle.
  writeLines("{}", file.path(td, "experiment.xenium"))
  file.copy(file.path(bundle, "experiment.spatial"), td)
  expect_null(cellAdmixCore:::.celladmix_atera_bundle_dir(td))
})

test_that("Atera conversion without a usable Python names the conversion command", {
  td <- tempfile("atera_nopython_")
  bundle <- write_atera_placeholder(file.path(td, "outs"))
  old <- options(cellAdmix.python = file.path(td, "no-such-python"))
  on.exit(options(old), add = TRUE)
  expect_error(
    cellAdmix(bundle, output_dir = file.path(td, "out")),
    "-m celladmix\\.atera"
  )
})

test_that("Atera bundles convert through the Python command and build a store", {
  python <- atera_test_python()
  skip_if(!nzchar(python), "CELLADMIX_PYTHON with zarr and celladmix.atera is not configured")
  td <- tempfile("atera_store_")
  dir.create(td)
  bundle <- file.path(td, "outs")
  # Three cells (id prefixes 0, 25465 and 15) over two grid tiles, three
  # unassigned molecules per tile, and one negative-control probe.
  script <- file.path(td, "make_bundle.py")
  writeLines(c(
    "import json, sys",
    "from pathlib import Path",
    "import numpy as np, zarr",
    "root = Path(sys.argv[1]); root.mkdir(parents=True)",
    "(root / 'experiment.spatial').write_text(json.dumps({'run_name': 'Mock Atera', 'pixel_size': 0.2125}))",
    "S = lambda v: np.array(v, dtype=np.dtypes.StringDType())",
    "def zip_group(name, attributes=None):",
    "    store = zarr.storage.ZipStore(str(root / name), mode='w')",
    "    return store, zarr.open_group(store, mode='w', attributes=attributes or {})",
    "genes = ['GA', 'GB', 'NegControlProbe_1']",
    "pairs = [(0, 1), (25465, 1), (15, 1)]",
    "centroids = [(10.0, 12.0), (30.0, 14.0), (20.0, 40.0)]",
    "store, g = zip_group('transcripts.zarr.zip', {'gene_names': genes, 'codeword_gene_mapping': [0, 1, 2]})",
    "grid = g.create_group('grid', attributes={'grid_keys': ['0,0', '1,0']})",
    "for t, key in enumerate(['0,0', '1,0']):",
    "    rows = [(pairs[c], cw, cx + k, cy - k) for c, (cx, cy) in enumerate(centroids) for k, cw in enumerate([0, 1, 0, 1, 2])]",
    "    rows += [((4294967295, 4294967295), 0, 50.0 + k, 50.0) for k in range(3)]",
    "    tile = grid.create_group(key)",
    "    tile.create_array('location', data=np.array([(r[2], r[3], 1.0 + t) for r in rows], dtype=np.float32))",
    "    tile.create_array('codeword_identity', data=np.array([[r[1]] for r in rows], dtype=np.uint32))",
    "    tile.create_array('cell_id', data=np.array([r[0] for r in rows], dtype=np.uint32))",
    "    tile.create_array('overlaps_nucleus', data=np.zeros((len(rows), 1), dtype=np.uint8))",
    "    tile.create_array('quality_score', data=np.full((len(rows), 1), 30, dtype=np.float16))",
    "store.close()",
    "store, g = zip_group('cells.zarr.zip')",
    "g.create_array('cell_id', data=np.array(pairs, dtype=np.uint32))",
    "cols = ['cell_centroid_x', 'cell_centroid_y', 'cell_area', 'nucleus_centroid_x', 'nucleus_centroid_y', 'nucleus_area', 'z_level', 'nucleus_count']",
    "g.create_array('cell_summary', data=np.array([[x, y, 100, x, y, 30, 0, 1] for x, y in centroids], dtype=float), attributes={'column_names': cols})",
    "for s in (0, 1):",
    "    p = g.create_group(f'polygon_sets/{s}')",
    "    v = np.zeros((3, 50), dtype=np.float32)",
    "    for c, (x, y) in enumerate(centroids):",
    "        v[c, :8] = [x - 3, y - 3, x + 3, y - 3, x + 3, y + 3, x - 3, y + 3]",
    "    p.create_array('vertices', data=v)",
    "    p.create_array('num_vertices', data=np.full(3, 4, dtype=np.uint32))",
    "    p.create_array('cell_index', data=np.arange(3, dtype=np.uint32))",
    "store.close()",
    "store, g = zip_group('cell_feature_matrix.zarr.zip')",
    "g.create_array('var/feature_name', data=S(genes))",
    "g.create_array('var/feature_type/categories', data=S(['Gene Expression', 'Negative Control Probe']))",
    "g.create_array('var/feature_type/codes', data=np.array([0, 0, 1], dtype=np.int8))",
    "g.create_array('obs/barcode', data=S(['aaaaaaaa-1', 'aaaagdhj-1', 'aaaaaaap-1']))",
    "g.create_array('obs/transcript_counts', data=np.array([8, 8, 8]))",
    "store.close()"
  ), script)
  expect_equal(system2(python, c(shQuote(script), shQuote(bundle))), 0L)

  old <- options(cellAdmix.python = python)
  on.exit(options(old), add = TRUE)
  ds <- suppressMessages(cellAdmix(bundle, output_dir = file.path(td, "out")))
  expect_true(file.exists(file.path(td, "out", "atera_xenium", "atera_source.json")))
  counts <- ds$counts()
  expect_setequal(colnames(counts), c("aaaaaaaa-1", "aaaagdhj-1", "aaaaaaap-1"))
  # Unassigned molecules and the control probe are left out of the store:
  # each cell keeps two GA and two GB molecules per tile.
  expect_setequal(rownames(counts), c("GA", "GB"))
  expect_equal(as.numeric(Matrix::colSums(counts)), c(8, 8, 8))
  if (requireNamespace("arrow", quietly = TRUE)) {
    cells <- as.data.frame(arrow::read_parquet(file.path(td, "out", "atera_xenium", "cells.parquet")))
    expect_equal(cells$x_centroid[cells$cell_id == "aaaagdhj-1"], 30)
    expect_equal(cells$y_centroid[cells$cell_id == "aaaagdhj-1"], 14)
  }

  # A second construction reuses the conversion without running Python.
  options(cellAdmix.python = file.path(td, "no-such-python"))
  expect_no_error(cellAdmix(bundle, output_dir = file.path(td, "out")))
})
