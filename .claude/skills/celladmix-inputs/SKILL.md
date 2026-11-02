---
name: celladmix-inputs
description: Supported inputs of cellAdmix-core and how a dataset is built - Xenium bundles, CosMx and other tabular molecule tables, segmentation masks, annotation formats, cell id conventions, control-feature filtering, and the on-disk input store. Use when loading data or debugging dataset construction.
---

Use this skill when constructing a dataset from a 10x Xenium bundle, a CosMx export, or any molecule table, when deciding how to pass a cell-type annotation, when a loader rejects a file or a schema, or when you need to know what the input store on disk contains. The authoritative page is `docs/inputs.md`; the loaders are `src/xenium.cpp`, `src/tabular.cpp`, `src/tabular_mask.cpp`, and `src/input_store.cpp`, wrapped by `r/R/constructor.R`, `r/R/schema.R`, `r/R/api.R` (`celladmix_prepare_xenium`, `celladmix_prepare_tabular`), and `python/src/celladmix/dataset.py`.

## Input families

| Source | R | Python |
|---|---|---|
| Xenium bundle directory (contains `experiment.xenium`) | `cellAdmix(bundle_dir, output_dir = , annotation = )` | `ca.CellAdmix(bundle_dir, output_dir=, annotation=)` |
| Tabular molecule table (`.csv`, `.csv.gz`, `.tsv`, `.tsv.gz`, `.parquet`, `.pq`) | `cellAdmix(path, format = "tabular", schema = celladmix_schema(...))` | not implemented (`ensure_store()` raises `NotImplementedError`) |
| Seurat object plus the original bundle | `cellAdmix(seu, xenium_dir = , output_dir = , annotation = "cell_type")` | - |
| SpatialData object plus the original bundle | - | `ca.from_spatialdata(sdata, xenium_dir=, output_dir=, annotation=)` |
| A `celladmix_prep` object from an earlier call | `cellAdmix(prep)` | `ca.CellAdmix.attach_existing(output_dir, source=)` |

`format = "auto"` (R) detects Xenium bundles, prep objects, and Seurat objects; pass `format = "tabular"` for tables. Constructing directly from an input-store directory (`format = "store"`) is not implemented in R and stops with a message.

## Xenium bundles

Constructor options (R `cellAdmix(...)` forwards `...` to `celladmix_prepare_xenium()`; Python names in parentheses):

- `analysis_bbox = c(xmin, xmax, ymin, ymax)` restricts the analysis scope at load time; `analysis_crops` accepts several named crops (R only).
- `cell_filter`, `gene_filter`: ids to retain while building the store.
- `min_qv`: minimum transcript quality (`min_qv=-1.0` in Python means no filter).
- `keep_unassigned = FALSE`: drop molecules without a usable cell assignment.
- `keep_non_gene = FALSE`: drop control, codeword, and other non-gene features during store construction, so molecule and gene thresholds refer to biological genes. The loader uses `codeword_category` or `is_gene` when the parquet exposes them, otherwise feature-name prefixes (`NegControl*`, `UnassignedCodeword*`, `DeprecatedCodeword*`, `Intergenic_Region*`, `GenomicControl*`, `Blank*`; see `xenium_feature_name_is_control` in `include/celladmix/xenium.hpp`). Set `TRUE` only for control-feature diagnostics.
- `read_cells = TRUE`, `prefer_parquet = TRUE` (R): read the bundle's cell metadata; prefer Parquet over CSV when both exist.

If the bundle's cell metadata has no `cell_type` column, supply an annotation or cluster first (`ds$cluster()` in R).

### Stain images

Membrane and coherence scoring can use the bundle's morphology images. Discovery is soft: an explicit `image_path` always wins; otherwise the manifest's morphology-focus sibling selected by `focus_index` (default 1) is tried, then the DAPI focus image, then the manifest morphology image. Physical coordinates map to pixels as `pixel = coordinate / pixel_size` with `x_offset = y_offset = 0` by default; single-channel tiled OME-TIFF focus images including the JPEG2000-compressed variant are read. Check what was resolved before scoring:

```r
fit$stain_image()                 # membrane stain; fit$stain("dapi") for DAPI
fit$stains()                      # both "dapi" and "membrane"
```

Python: `fit.stain("membrane")`, `fit.stains()`; `celladmix.io.discover_xenium_stain_image()`.

## Tabular molecule tables (CosMx and others)

Describe the columns with `celladmix_schema()` (defaults in parentheses):

- Required: `x` ("x"), `y` ("y"), `gene` ("gene").
- Optional: `z` ("z"), `qv` (NULL), `cell_type` (NULL, molecule-level labels), `sample_id` (NULL), `fov_id` (NULL).
- Cell assignment, exactly one of: `cell` ("cell", a column of cell ids; missing or zero-like ids are unassigned) or `segmentation_mask` (NULL, a TIFF path). Because `cell` defaults to `"cell"`, pass `cell = NULL` when using a mask, otherwise `celladmix_prepare_tabular()` stops with "Specify either cell_id_col or segmentation_mask_path, not both".
- Cell metadata sidecar: `cell_metadata` (path), `cell_metadata_cell` ("cell"), `cell_metadata_cell_type` (NULL). When the last is set and no `annotation` is given, the sidecar becomes the initial annotation.

```r
schema <- celladmix_schema(x = "x_global_px", y = "y_global_px", gene = "target",
  cell = "cell_ID", cell_metadata = "cells.csv.gz",
  cell_metadata_cell = "cell_ID", cell_metadata_cell_type = "cell_type")
ds <- cellAdmix("tx_file.csv.gz", output_dir = "out", format = "tabular", schema = schema)
```

TIFF segmentation masks must be single-channel: labeled masks (8-, 16-, or 32-bit pixel values are cell ids) or binary masks (all nonzero pixels equal; 4-connected components become cells). Molecules map to pixels by `col = round(x) - 1`, `row = round(y) - 1`; molecules outside the image or on background are unassigned. Multi-channel masks and polygon boundary files are not supported.

Filtering: `analysis_bbox`, `min_qv` (needs `qv`), `keep_unassigned` (kept molecules get the cell id `__unassigned__`). The CosMx example is `examples/cosmx_nsclc_giotto/` (`prepare_cosmx_tabular.R`, `celladmix_cosmx_minimal.Rmd`); the command-line route is `Rscript scripts/celladmix_batch.R --format tabular --input molecules.csv.gz --output out --cell-metadata cells.csv.gz --cell-metadata-cell-type cell_type`.

## Annotations

- R: a named character vector (names are cell ids) or a data frame with a cell id column (`cell_id_col`, default `"cell_id"`) and a label column (`annotation_col`). Stored under a name (`annotation_name`, default `"manual"`); `ds$set_annotation()`, `ds$use_annotation()`, `ds$annotations()`, `ds$annotation(name, as_vector = TRUE)` manage several annotations, one active.
- Python: `celladmix.io.read_annotation()` accepts a CSV/CSV.GZ/Parquet path, a `pd.Series` indexed by cell id, a dict, or a DataFrame; the label column is `annotation_col` or the first of `cell_type`, `merged_annotation`, `annotation`, `cluster`, `cluster_id`. `ds.set_annotation()` replaces it. Empty and missing labels are dropped.
- Without an annotation, `ds$cluster(name = "louvain", active = NULL, ...)` (R) produces provisional labels `cluster_<k>`; the batch scripts do this with `--auto-annotate`. Curated annotations are preferred for biological analysis. Annotate at a level the scoring can resolve; very fine labels make the auto-selected rank large and the per-type evidence sparse.

The annotation content is hashed into the run manifest, so changing labels triggers a refit.

## Cell id conventions

- Cell ids are strings taken from the source (Xenium `cell_id`, the tabular `cell` column, or mask labels). Count matrices are gene-by-cell with these ids as column names; corrected counts keep them.
- Seurat integration uses the object's cell names and reorders counts to match (`celladmix_add_corrected_assay`); SpatialData uses the table's `obs_names` unless `cell_id_key` names an obs column.
- Multi-sample collections prefix combined outputs with the sample name: `<sample>_<cell_id>`; per-sample accessors return the original ids (`docs/multisample.md`).

## The input store

`<output_dir>/input_store/` holds `store.json` (manifest: `source_type`, `source_path`, `source_fingerprint`, `filter_signature`, counts of molecules, cells, and genes, column availability flags), `genes.parquet`, `cells.parquet` (ids, centroids, sample and FOV ids when present), `cell_gene_counts.parquet`, `molecules.parquet` (encoded gene and cell indices with coordinates), and `cell_offsets.parquet` (molecule row ranges per cell for block-wise processing). The store is built on first use and reused when its fingerprint and filters match; `force = TRUE`/`force=True` rebuilds it (R `ds$ensure_store()`, Python `ds.ensure_store(force=True)`). Counts and cell positions are available without any fit: `ds$counts()`, `ds$cells()` (R); `ds.counts()`, `ds.cells()` (Python).
