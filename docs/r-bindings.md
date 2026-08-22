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
than the KL variants (gene-ownership correlation 0.92 versus 0.48 between
restarts; see [benchmarks.md](benchmarks.md), Figure 2a), and the
factor decomposition it aims for — native factors per cell type or state,
with admixture read from their minor contributions in non-native cells —
matches the scoring model directly. In cleanup benchmarks
([benchmarks.md](benchmarks.md)) the most effective variant follows the
scoring method: `invsqrt_kl` paired with membrane scoring on stained data,
`ls_nmf` paired with bridge scoring — the pancreas quickstart follows this
pairing.

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

With `nmf_n_runs > 1` (default: at least 10 restarts, one per thread above
that), the fit records a per-factor restart stability diagnostic. Factors
below the stability threshold (0.3) are seed-dependent rather than
reproducible structure and should be interpreted with caution — the default
ensemble correction votes across the restarts rather than trusting any
single one; the number of stable factors is also the primary rank signal. See
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
[benchmarks.md](benchmarks.md) for the methodology). Because a tissue
section shows only a slab of the tissue, cells with zero observed
source-type neighbors can still be contaminated by source cells just above
or below the section plane; the audit therefore takes its reference level
from target cells with zero source-type cells among a progressively larger
set of nearest neighbors (up to 240, spanning only a ~100 um radius, so the
comparison stays within the local tissue neighborhood), using the largest
neighborhood that retains enough reference molecules. The `reference_kind`
and `reference_inflation` columns of `pairs()` record the neighborhood used
and how much contamination this removed from the comparison group. Marker panels are also screened for likely induced
genes — genes whose exposure-linked excess is disproportionate to the
source expression profile reflect a transcriptional response to proximity
rather than transferred material. The comparison profile is measured from
the source cells bordering the target type (source cells at an interface
can genuinely express activation genes above the source average, and
material transferred from them carries that share), and the test allows
for both counting noise and profile uncertainty, so small relative
deviations on large transfer channels are not flagged. Flagged genes are
excluded from the panels, replaced by the next-ranked source-specific
genes, and listed with their measured disproportionality in
`markers()$induced` and `markers()$induced_stats`. The audit measures every ordered cell-type
pair:

```r
audit <- fit$audit_admixture()
audit$pairs()                      # per-pair admixture rates and molecule estimates
audit$plot_map()                   # admixture-rate map (% of target-type molecules)
audit$plot_exposure()              # pooled excess-exposure profile, 95% intervals
audit$plot_exposure("Exocrine epithelial", "Endothelial", correction = correction)
audit$plot_remaining(list(membrane = correction))  # admixture left per correction
```

The audit can also correct the admixture it measures, through the
generative model described in [generative.md](generative.md): each target
cell's counts are decomposed into the cell's own expression, contamination
from each detected source type, ambient background, and
neighborhood-induced expression; the contamination and ambient shares are
removed while the induced ones are retained. The fit runs in the C++ core
(seconds to a few minutes depending on dataset size) and returns a
correction object that `evaluate()` accepts, together with the model's
per-cell decomposition:

```r
correction <- audit$correct_generative(num_threads = 8)
correction$counts()                # corrected gene x cell matrix
correction$pairs                   # per-pair removed and induced molecule totals
correction$induced                 # retained induced genes with disproportionality
correction$composition("Exocrine epithelial", "Ductal/tumor epithelial")
                                   # per-cell dose, contamination, induced activity
audit$evaluate(correction)         # verified like any other correction
```

`audit$evaluate(correction)` verifies a correction against the same
measurements: per-pair cleanup sensitivity, the own-marker false-removal
rate per cell type (removal of near-surely-genuine molecules), a warning
for any detected pair whose molecules the correction measurably failed to
remove (based on the corrected counts themselves, since ensemble members
can remove molecules for pairs absent from the primary rule list and rules
can exist yet remove nothing), and a severe-over-removal warning when a
type loses more than 25% of its own-marker molecules.

```r
report <- audit$evaluate(correction)
report$summary()
report$plot_cleanup()
```

The reported rates and molecule counts include the contamination present
in zero-neighbor cells above the ambient reference, and extrapolate the
marker-panel measurement to the full transcriptome by the markers' share
of the source expression profile, kept as the `coverage` column of
`audit$pairs()`. Pair detection itself remains based on the rise of
exposed cells over unexposed ones, which no reference choice can inflate.

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
correction$ensemble()   # member count, vote threshold, vote histogram
```

Correction removes molecules assigned to rule-supported admixture factors in
their target cell types and returns sparse corrected counts.

By default the correction is a molecule-vote ensemble over the fit's NMF
restarts: each restart's factorization labels every molecule, is scored with
the same method and parameters, and is vetted by its own native-factor
check; a molecule is removed when at least `vote` (default 0.3) of the
members remove it. Single-fit corrections are highly seed-dependent, and the
vote both stabilizes them and consistently reaches the upper range of the
individual restarts on the cleanup benchmark (see
[benchmarks.md](benchmarks.md)), with the threshold acting as a
sensitivity/specificity dial:

```r
score$correct()               # ensemble over up to 10 restarts, vote = 0.3
score$correct(vote = 0.5)     # stricter: majority vote, higher specificity
score$correct(ensemble = 1)   # single-fit correction from the selected restart
```

The first ensemble correction computes and caches per-member molecule
labelings in the run directory (one projection and smoothing pass per
member); re-voting at a different threshold reuses the cached member rules.
When `rules` is passed explicitly, its source→target pairs restrict every
member, so vetoed pairs stay vetoed across the ensemble; the returned
correction's `rules` carry a `support` column giving the fraction of members
that independently kept each pair. Runs fitted with `nmf_n_runs = 1` (or
cached runs whose stores lack a restart member pool) fall back to the
single-fit correction with a message.

`correct()` warns when a rule set removes more than 60% of a cell type's
molecules — usually a sign that the factorization has no native factor for
that type; the minimum-molecule gate for the warning is
`options(celladmix.overremoval_min_molecules = )`.

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

- [Pancreas quickstart](../examples/xenium_pancreas_membrane_377_full/pancreas_quickstart.md): the recommended workflow in its shortest form - fit, audit, correct, verify.
- [Minimal CosMx NSCLC tutorial](../examples/cosmx_nsclc_giotto/celladmix_cosmx_minimal.md): a compact tabular example mirroring the original cellAdmix [NSCLC tutorial](https://github.com/kharchenkolab/cellAdmix/blob/main/vignettes/NSCLC_tutorial_fulldata.ipynb).
- [Detailed Xenium pancreas tutorial](../examples/xenium_pancreas_membrane_377_full/pancreas_membrane_scoring_clean.md): the full walkthrough on a membrane-stained Xenium bundle - audit, membrane and bridge scoring on their recommended factorizations, and the variant comparison behind that pairing.
- [Seurat Xenium integration tutorial](../examples/xenium_pancreas_membrane_377_full/pancreas_seurat_integration.md): the same pancreas dataset, using Seurat for cell-level state and cellAdmix for molecule-complete fitting, scoring, and correction.
- [Xenium breast 5K membrane scoring tutorial](../examples/xenium_breast_membrane_5k_full/breast_5k_membrane_scoring.md): the same workflow at full 5K-panel scale.
- [Generative correction in brief](../examples/xenium_pancreas_membrane_377_full/pancreas_generative_minimal.md): the generative admixture model in four calls.
- [Generative correction on Xenium pancreas](../examples/xenium_pancreas_membrane_377_full/pancreas_generative.md): the admixture pattern between cell types and the induced expression the model preserves, characterized and visualized.
- [Generative correction on CosMx NSCLC](../examples/cosmx_nsclc_giotto/nsclc_generative.md): the same characterization on a CosMx dataset.

## Detailed Pages

- [Inputs](inputs.md)
- [Installation](install.md)
- [Scoring methods](scoring_methods.md)
- [NMF restart stability](nmf_stability.md)
- [Cleanup benchmarks](benchmarks.md)
- [Documentation index](README.md)
