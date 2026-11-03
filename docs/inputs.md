# Inputs

`cellAdmix-core` supports two file-level input families through the public
`cellAdmix()` constructor:

- 10x Xenium and Atera bundles
- generic tabular molecule tables

The constructor also accepts Seurat objects and prepared-project/store
directories directly (see [r-bindings.md](r-bindings.md)).

The normal workflow is object-oriented:

```r
ds <- cellAdmix(source, output_dir = ..., annotation = ...)
fit <- ds$fit(...)
score <- fit$score("membrane", ...)
clean <- score$correct(...)
```

The `celladmix_prepare_*()`, `celladmix_fit()`, and
`celladmix_score_*()` functions are internal backend helpers.

## Xenium

Use `cellAdmix()` directly for 10x Xenium bundles.

Expected source:

- Xenium bundle directory
- transcript table in Parquet or CSV/TSV form
- optional cell metadata when available

Key constructor options:

- `source`
- `output_dir`
- `analysis_bbox`
- `keep_unassigned`
- `keep_non_gene`
- `read_cells`
- `prefer_parquet`

Notes:

- `analysis_bbox` restricts the analysis scope at load time.
- `keep_unassigned = FALSE` drops molecules without a usable cell assignment.
- `keep_non_gene = FALSE` is the default. Xenium control/codeword/non-gene
  features are filtered during input-store construction, so molecule and gene
  thresholds refer to biological gene expression. For newer Xenium parquet
  bundles, cellAdmix uses `codeword_category` or `is_gene` when available; for
  older bundles it falls back to stable Xenium feature names such as
  `NegControl*`, `UnassignedCodeword*`, `DeprecatedCodeword*`, and
  `Intergenic_Region*`. Set `keep_non_gene = TRUE` only for control-feature
  diagnostics.
- When `cell_type` is absent, later clustering is expected to provide training strata.

### Xenium Membrane Images

`fit$score("membrane", ...)` can use Xenium morphology images after a run has
been fitted.

Image selection is intentionally soft:

- `image_path` can be supplied explicitly and always takes precedence.
- If `image_path` is omitted, the internal image-discovery helper reads the
  Xenium manifest and tries the morphology-focus sibling selected by
  `focus_index`.
- If the requested focus sibling is absent, discovery falls back to the
  manifest focus image for the DAPI stain, and otherwise to the manifest
  morphology image.

For the current Xenium pancreas membrane example, `focus_index = 1` resolves to:

```r
morphology_focus/morphology_focus_0001.ome.tif
```

Current image assumptions:

- physical coordinates are mapped as `pixel = coordinate / pixel_size`
- default offsets are `x_offset = 0`, `y_offset = 0`
- explicit `pixel_size`, `x_offset`, and `y_offset` can override manifest
  defaults
- single-channel tiled OME-TIFF focus images are supported, including the
  JPEG2000-compressed Xenium focus-image variant

Because Xenium bundle manifests and channel naming can change across 10x bundle
versions, notebooks should call `fit$stain_image()` and show the resolved path
before scoring. Use explicit `image_path` when discovery picks the wrong
channel.

## Atera

10x Atera output directories are passed exactly like Xenium bundles:

```r
ds <- cellAdmix("/path/to/atera/outs", output_dir = "out", annotation = annotation)
```

```python
ds = ca.CellAdmix("/path/to/atera/outs", output_dir="out", annotation=annotation)
```

An Atera directory is recognized by its `experiment.spatial` manifest (a
directory that also has `experiment.xenium` is read as Xenium). Atera stores
molecules, cells, and the cell-feature matrix as zipped Zarr stores
(`transcripts.zarr.zip`, `cells.zarr.zip`, `cell_feature_matrix.zarr.zip`)
rather than the Parquet tables the store builder reads. On first use these
stores are written once as Xenium-layout tables under
`<output_dir>/atera_xenium` (`transcripts.parquet`, `cells.parquet`,
`cell_boundaries.parquet`, `nucleus_boundaries.parquet`, and an
`experiment.xenium` manifest), and the rest of the workflow reads that
directory as a Xenium bundle. The tables are reused as long as the Atera
manifest and Zarr stores keep the sizes and modification times recorded at
conversion. The conversion reads the transcript store one spatial grid tile
at a time, so its memory use is bounded by the largest tile.

The conversion is done by the Python module `celladmix.atera`, which needs
the `zarr` package (version 3 or later; `pip install "celladmix[atera]"`).
The R package runs the same conversion as

```bash
python -m celladmix.atera /path/to/atera/outs out/atera_xenium
```

with the interpreter given by `options(cellAdmix.python = ...)`, the
`CELLADMIX_PYTHON` environment variable, or `python3` on the `PATH`. If no
interpreter with `celladmix` is available, `cellAdmix()` stops with this
command; running it once by hand lets later R sessions reuse the tables
without Python.

Record conversion:

- Cell ids take the Xenium string form: the eight hexadecimal digits of the
  integer id prefix written with the letters `a`-`p`, then `-` and the dataset
  suffix (the id pair (25465, 1) becomes `aaaagdhj-1`).
- Molecules that the Zarr store marks with the id pair (4294967295,
  4294967295) are unassigned and get `UNASSIGNED`, so `keep_unassigned`
  applies to them as for Xenium.
- The feature type of each gene in the cell-feature matrix becomes the
  molecule's `codeword_category`, and `is_gene` is true only for
  "Gene Expression" features; negative-control probes and codewords, genomic
  controls, and unassigned or deprecated codewords are therefore excluded
  under the default `keep_non_gene = FALSE`.
- The Atera per-molecule quality score becomes `qv`, so `min_qv` applies to it.
- Cell centroids come from the cell summary table; cell and nucleus
  boundaries come from the cell store's polygon sets.

The Atera manifest lists no morphology-focus images, so membrane scoring on
Atera input needs the stain image and its pixel size passed explicitly
(`image_path` and `pixel_size`).

## Generic Tabular

Use `cellAdmix(format = "tabular", ...)` for CSV/TSV/Parquet molecule tables.
CSV and TSV inputs may be gzip-compressed.

Supported molecule file formats:

- `.csv`
- `.csv.gz`
- `.tsv`
- `.tsv.gz`
- `.parquet`
- `.pq`

Required schema fields:

- `x`
- `y`
- `gene`

Optional schema fields:

- `z`
- `qv`
- `cell_type`
- `sample_id`
- `fov_id`
- `cell`
- `segmentation_mask`
- `cell_metadata`
- `cell_metadata_cell`
- `cell_metadata_cell_type`

Tabular input requires one of:

- `cell`
- `segmentation_mask`

If neither is supplied, preparation fails.

### Explicit Cell IDs

If the molecule table already contains cell assignments, pass a schema with:

- `cell = "<column name>"`

This is the simplest tabular mode. Molecules with missing or zero-like cell ids are treated as unassigned.

### TIFF Segmentation Masks

If the molecule table does not contain cell ids, pass a schema with:

- `segmentation_mask = "<mask.tif>"`

Supported TIFF mask variants:

- Single-channel labeled masks
  - pixel values are treated as cell ids directly
  - supported bit depths: `8`, `16`, `32`
- Single-channel binary masks
  - if all nonzero pixels share the same value, the mask is treated as binary
  - connected-component labeling is run first
  - the resulting connected components become cell ids

Current TIFF assumptions:

- mask must be single-channel
- molecule coordinates are mapped to pixels by:
  - `col = round(x) - 1`
  - `row = round(y) - 1`
- molecules outside the image bounds are unassigned
- molecules landing on background pixels are unassigned

Implementation detail:

- connected components are 4-connected, matching the Baysor-style binary-mask convention for ISS/DAPI masks

### Unassigned Molecules

For tabular input:

- `keep_unassigned = FALSE` drops unassigned molecules
- `keep_unassigned = TRUE` keeps them under the internal cell id `__unassigned__`

For current cell-centric workflows, `keep_unassigned = FALSE` is usually the sensible default.

## Cropping and Filtering

Xenium, Atera, and tabular input all support restricting the analysis scope.

Common filtering options:

- `analysis_bbox`
- `min_qv`

Tabular-only filtering options:

- `qv_col`

Behavior:

- rows outside `analysis_bbox` are excluded
- if `min_qv` is set and `qv_col` is available, low-quality rows are excluded

## Current Limitations

Not currently supported:

- tabular molecule tables with no cell ids and no segmentation mask
- multi-channel segmentation masks
- polygon or boundary segmentation files in the `cellAdmix` prep path

If a dataset comes with a binary nucleus mask, use `celladmix_schema(segmentation_mask = ...)`; connected components become cell ids automatically.

## Minimal Examples

### Xenium

```r
ds <- cellAdmix(
  "/path/to/xenium_bundle",
  output_dir = "/path/to/output",
  analysis_bbox = NULL,
  annotation = annotation,
  prefer_parquet = TRUE
)
```

### Tabular With Cell IDs

```r
ds <- cellAdmix(
  format = "tabular",
  molecules = "/path/to/molecules.csv",
  output_dir = "/path/to/output",
  schema = celladmix_schema(x = "x", y = "y", gene = "gene", cell = "cell_id"),
  annotation = annotation
)
```

### Tabular With TIFF Mask

```r
ds <- cellAdmix(
  format = "tabular",
  molecules = "/path/to/molecules.csv",
  output_dir = "/path/to/output",
  schema = celladmix_schema(x = "x", y = "y", gene = "gene",
    segmentation_mask = "/path/to/segmentation_mask.tif")
)
```
