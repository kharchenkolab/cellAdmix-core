# R Bindings

The R package `cellAdmixCore` exposes the C++ cellAdmix implementation through
an R6 workflow. The public interface is centered on the `cellAdmix()`
constructor, which records the input source, creates an output/cache directory,
and returns a dataset object.

The R layer keeps small cell-level state in memory, such as annotations and
metadata. Large molecule tables and fitted run artifacts are handled by the C++
core through file-backed stores.

See [install.md](install.md) for R package installation and system
dependencies.

## Basic Workflow

```r
library(cellAdmixCore)

ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)
fit <- ds$fit()
score <- fit$score_membrane()
rules <- score$rules()
correction <- score$correct(rules = rules)
corrected_counts <- correction$counts()
```

The main objects are:

| Object | Purpose |
|---|---|
| `CellAdmixDataset` | Input source, output directory, annotations, clustering results, and fitted runs. |
| `CellAdmixFit` | NMF factors, molecule assignments, smoothed labels, cell-level factor fractions, and plotting helpers. |
| `CellAdmixScore` | Source/target/factor admixture evidence from bridge, membrane, or coherence scoring. |
| `CellAdmixCorrection` | Rule-based molecule removal results and corrected sparse counts. |

## Input Sources

### Xenium Bundles

```r
ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)
```

Xenium is the primary high-performance path. The C++ backend reads the Xenium
bundle and stores encoded molecules, cells, genes, and reusable intermediate
artifacts in `output_dir`.

By default, Xenium control/codeword/non-gene features are excluded while
building the input store. Use `keep_non_gene = TRUE` only when explicitly
diagnosing control features.

See [inputs.md](inputs.md) for supported Xenium files, stain-image discovery,
and cropping/filtering options.

### Tabular Molecule Tables

```r
schema <- celladmix_schema(cell_metadata = "prepared/cell_metadata.csv.gz",
  cell_metadata_cell = "cell", cell_metadata_cell_type = "cell_type")

ds <- cellAdmix("prepared/molecules.csv.gz", output_dir = "out",
  format = "tabular", schema = schema)
```

Tabular inputs support CSV, gzip-compressed CSV, TSV, gzip-compressed TSV, and
Parquet molecule tables.

### Seurat Objects

```r
seu <- LoadXenium("data", molecule.coordinates = FALSE)
seu$cell_type <- ...

ds <- cellAdmix(seu, xenium_dir = "data", output_dir = "out",
  annotation = "cell_type")
fit <- ds$fit()
score <- fit$score_membrane()
correction <- score$correct(rules = score$rules())

seu <- celladmix_add_to_seurat(seu, fit = fit, score = score,
  correction = correction, corrected_assay = "cellAdmix")
```

The Seurat object supplies the active cells, features, annotations, and
downstream storage target. The Xenium bundle supplies the full transcript table
needed by the molecule-level cellAdmix workflow.

## Fitting

```r
fit <- ds$fit()
```

`fit()` builds molecule-centered neighborhood composition vectors, fits NMF
factors, assigns molecules to factors, applies spatial smoothing, and summarizes
factor fractions at cell level.

The default rank is derived from the active annotation. The default NMF method
is `ls_nmf`, the weighted least-squares formulation of the original cellAdmix:
across datasets it recovers cell-type-native factors far more reproducibly
than the KL variants (see [nmf_stability.md](nmf_stability.md)), and the
factor decomposition it aims for — native factors per cell type or state,
with admixture read from their minor contributions in non-native cells —
matches the scoring model directly. `invsqrt_kl` remains available for
marker-driven loadings.

The NCV neighborhood size `ncv_k` is resolved automatically from the data:
it grows with the square root of the panel size (a ~400-gene panel keeps the
historical value of 20; a 5,000-gene panel resolves to ~70) and is capped at
half the median cell's molecule count to preserve sub-cellular locality. The
resolved value is logged during the fit and recorded in the run manifest.
Passing an explicit `ncv_k` overrides resolution; a warning is raised when an
explicit value leaves fewer than one neighborhood draw per 100 panel genes,
since such neighborhoods carry too little gene co-occurrence signal for
reproducible factors.

Supported `nmf_variant` values:

| `nmf_variant` | Description |
|---|---|
| `ls_nmf` | Default weighted least-squares NMF (original cellAdmix formulation). |
| `invsqrt_kl` | Sparse KL-NMF with inverse-square-root gene-prevalence weighting. |
| `kl` | Sparse KL-NMF without gene weighting. |
| `sqrt_kl` | Sparse KL-NMF after square-root transformation of NCV counts. |

Fitted runs are cached on disk by run id and reused **only when the requested
parameters match the cached run**: on each `fit()` call the requested
parameters (and the annotation content) are compared against the values
recorded in the run manifest. A match reuses the run with a
`Reusing cached run ...` message; any difference triggers an automatic refit
with the changed parameters spelled out
(`Parameters changed for run ... (ncv_k: 20 -> 71); refitting`).
Automatically resolved values (an omitted `ncv_k` or `nmf_n_runs`) match
whatever the cached run recorded. `overwrite = TRUE` forces a refit even when
parameters match — useful after package upgrades, which are reported in the
reuse message but do not trigger refits on their own. Scores and corrections
are recomputed on every call and never go stale.

Useful fit diagnostics:

```r
fit$plot_loadings()
fit$plot_stability()
fit$cell_factors()
```

With `nmf_n_runs > 1` (default: one restart per thread), the fit records a
per-factor restart stability diagnostic. Factors below the stability
threshold (0.3) are seed artifacts rather than reproducible structure and
should not drive factor interpretation or correction rules; the number of
stable factors is also the primary rank signal. See
[NMF restart stability](nmf_stability.md) for the definition, the restart
initialization recipe, and interpretation guidance.

## Scoring

Scoring interprets NMF factors as potential source/target admixture patterns.

```r
membrane_score <- fit$score_membrane()
bridge_score <- fit$score_bridge()
coherence_score <- fit$score_coherence()
```

Common score operations:

```r
score$plot_heatmap()
score$plot_pairs()
rules <- score$rules(p_thresh = 0.1)
```

Rules pass a native-factor false-positive check by default: each rule is
tested against target cells that have no source-type cells among their
nearest neighbors, and rules whose factor persists in those source-distant
cells are flagged `keep = FALSE` with the reason in the `native_check` column
(supporting evidence columns included). Correction skips flagged rules; pass
`native_check = FALSE` to disable, or adjust `native_median_thresh` and
related thresholds. See the Native-Factor Check section in
[scoring_methods.md](scoring_methods.md).

## Auditing Admixture and Verifying Cleanup

Independently of factorization and scoring, the amount of admixture in a
dataset can be estimated from its spatial structure: source-marker content
in target cells rises with the number of source-type neighbor cells, while
unexposed target cells provide an internal negative control (see
[benchmarks.md](benchmarks.md) for the methodology). The audit measures
this for every ordered cell-type pair:

```r
audit <- fit$audit_admixture()
audit$pairs()                      # per-pair leaked-molecule estimates
audit$plot_map()                   # source x target admixture overview
audit$plot_exposure()              # cumulative exposure profile, 95% intervals
audit$plot_exposure("Exocrine epithelial", "Endothelial", correction = correction)
audit$plot_remaining(list(membrane = correction))  # admixture left per correction
```

`audit$evaluate(correction)` verifies a correction against the same
measurements: per-pair cleanup sensitivity, the own-marker false-removal
rate per cell type (removal of near-surely-genuine molecules), and a
warning for any detected pair that no removal rule covers.

```r
report <- audit$evaluate(correction)
report$summary()
report$plot_cleanup()
```

Estimates are conservative lower bounds: contamination that reaches even
unexposed cells raises the reference level and is not counted.

Example-cell overlays show molecule-level evidence around selected target
cells. For Xenium-backed fits, the DAPI/membrane stain background, cell
boundaries, and cell-type contour coloring are discovered automatically:

```r
score$plot_examples()                              # auto-selected cells
score$plot_examples(cells = c("cell-1", "cell-2")) # specific cells
score$plot_example("cell-1", stains = NULL)        # one cell, no background
```

Pass `color_cell_types = FALSE` to drop the cell-type contours, or explicit
`stains` / `boundaries` / `cell_types` arguments to override discovery.

See [scoring_methods.md](scoring_methods.md) for method definitions and tradeoffs.

## Correction

```r
correction <- score$correct(rules = rules, name = "clean")
corrected_counts <- correction$counts()
correction$summary()
correction$plot_removed_molecules()
```

Correction removes molecules assigned to rule-supported admixture factors in
their target cell types and returns sparse corrected counts.

## Clustering and Annotation

If no annotation is supplied, cellAdmix can generate provisional clustering
labels:

```r
clust <- ds$cluster(name = "cluster", active = TRUE)
plot(clust)
```

These clusters are intended as quick labels for testing or exploratory runs.
Curated or externally derived annotations should be preferred for biological
analysis.

## Batch Processing

For CLI-style processing:

```bash
Rscript scripts/celladmix_batch.R --input data --output out \
  --annotation annotations/annotation.csv.gz \
  --annotation-col merged_annotation \
  --score auto \
  --report \
  --threads 10
```

`--score auto` uses membrane scoring when a membrane stain is discovered and
falls back to bridge scoring otherwise. If no annotation is supplied, the script
stops unless `--auto-annotate` is passed.

Run:

```bash
Rscript scripts/celladmix_batch.R --help
```

for the full option list.

## Example Notebooks

- [Minimal CosMx NSCLC tutorial](../examples/cosmx_nsclc_giotto/celladmix_cosmx_minimal.ipynb): a compact tabular example mirroring the original cellAdmix [NSCLC tutorial](https://github.com/kharchenkolab/cellAdmix/blob/main/vignettes/NSCLC_tutorial_fulldata.ipynb).
- [Xenium pancreas membrane/bridge scoring tutorial](../examples/xenium_pancreas_membrane_377_full/pancreas_membrane_scoring_clean.ipynb): a modern Xenium bundle example with membrane cell staining.
- [Seurat Xenium integration tutorial](../examples/xenium_pancreas_membrane_377_full/pancreas_seurat_integration.ipynb): the same pancreas dataset, using Seurat for cell-level state and cellAdmix for molecule-complete fitting, scoring, and correction.

## Detailed Pages

- [Inputs](inputs.md)
- [Installation](install.md)
- [Scoring methods](scoring_methods.md)
- [Documentation index](README.md)
