# Inputs

`cellAdmix-core` currently supports two input families through the public
`cellAdmix()` constructor:

- Xenium bundles
- generic tabular molecule tables

The normal workflow is object-oriented:

```r
ds <- cellAdmix(source, output_dir = ..., annotation = ...)
fit <- ds$fit(...)
score <- fit$score("membrane", ...)
clean <- score$correct(...)
```

The older `celladmix_prepare_*()`, `celladmix_fit()`, and
`celladmix_score_*()` functions are internal backend helpers in the initial
package API.

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
- If the requested focus sibling is absent, discovery falls back to the manifest
  focus image, then to the manifest morphology image.

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

- connected components are currently 4-connected, matching the Baysor-style binary-mask path we adopted for ISS/DAPI masks

### Unassigned Molecules

For tabular input:

- `keep_unassigned = FALSE` drops unassigned molecules
- `keep_unassigned = TRUE` keeps them under the internal cell id `__unassigned__`

For current cell-centric workflows, `keep_unassigned = FALSE` is usually the sensible default.

## Cropping and Filtering

Both Xenium and tabular input support restricting the analysis scope.

Common filtering options:

- `analysis_bbox`

Tabular-only filtering options:

- `qv_col`
- `min_qv`

Behavior:

- rows outside `analysis_bbox` are excluded
- if `min_qv` is set and `qv_col` is available, low-quality rows are excluded

## Current Limitations

Not currently supported:

- tabular molecule tables with no cell ids and no segmentation mask
- multi-channel segmentation masks
- polygon or boundary segmentation files in the `cellAdmix` prep path

If a dataset comes with a binary nucleus mask, use `segmentation_mask_path`; the loader now converts connected components into cell ids automatically.

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
