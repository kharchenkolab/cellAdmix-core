# cellAdmix-core

`cellAdmix-core` is a C++ implementation of the core
[cellAdmix](https://github.com/kharchenkolab/cellAdmix) workflow for detecting
and correcting transcript-level admixture in spatial transcriptomics data.

The package is designed around a high-performance C++ core with language
bindings. The current public interfaces are the R package `cellAdmixCore` and
an initial Python package, with integrations for Seurat and SpatialData-style
workflows.

The implementation is optimized for large molecule-level datasets. Xenium
bundles are handled as first-class inputs, with molecule tables streamed into an
encoded on-disk store so large experiments do not need to be kept fully in R or
Python memory. As a rough local benchmark, on Xenium 5K data (breast 5K;
96.6M kept molecules and 694K cells), a 10-thread batch run with
auto-annotation, rank-30 fit, membrane scoring, a single-fit correction
(`--ensemble 1`), and a minimal report completed in about 62 minutes with
peak RSS about 42 GB; the default ensemble correction repeats the
molecule-assignment and scoring stages once per ensemble member.

The quickest way to see the workflow end to end is the
[pancreas quickstart notebook](examples/xenium_pancreas_membrane_377_full/pancreas_quickstart.ipynb):
audit the dataset's admixture, fit, correct, and verify the cleanup in a few
short steps.

## Highlights

- C++ core for neighborhood composition construction, NMF fitting, molecule
  assignment, spatial smoothing, admixture scoring, and correction.
- R bindings with a high-level `cellAdmix()` constructor and R6 workflow.
- Python bindings with a high-level `CellAdmix` class and optional SpatialData
  adapter.
- Native Xenium bundle support, including transcript tables, cell metadata,
  boundaries, and morphology images when available.
- Generic tabular input support for `.csv`, `.csv.gz`, `.tsv`, `.tsv.gz`, and
  Parquet molecule tables.
- Seurat integration: use a Seurat object for cells, genes, annotations, and
  downstream storage while using the original Xenium bundle as the
  molecule-complete backing source.
- Optional quick cell clustering and marker diagnostics when no annotation is
  available.
- Corrected sparse cell-by-gene count collection and add-back to Seurat.
- Admixture audit: exposure-gradient estimates of leaked molecules per
  cell-type pair, independent of the factorization, with post-correction
  verification of cleanup sensitivity and false removal
  ([docs/benchmarks.md](docs/benchmarks.md)).
- Ensemble correction by default: molecule-level voting across NMF restarts
  stabilizes seed-dependent corrections, with the vote threshold as a
  sensitivity/specificity dial.

## NMF Methods

cellAdmix-core supports several factorization modes. `ls_nmf` is the
default; in cleanup benchmarks ([docs/benchmarks.md](docs/benchmarks.md))
the most effective variant follows the scoring method — `invsqrt_kl` for
membrane-scored corrections on stained data, `ls_nmf` for bridge-scored
ones.

| `nmf_variant` | Description |
|---|---|
| `ls_nmf` | Default. Weighted least-squares NMF following the original cellAdmix formulation; produces the most reproducible cell-type-native factors across restarts. |
| `invsqrt_kl` | Sparse KL-NMF with inverse-square-root gene-prevalence weighting; marker-driven loadings at the cost of lower restart stability. Best cleanup with membrane scoring. |
| `kl` | Sparse KL-NMF without gene weighting. |
| `sqrt_kl` | Sparse KL-NMF after square-root transformation of NCV counts. |

Multirun fits report a per-factor restart stability diagnostic used for
factor trust and rank selection; see
[docs/nmf_stability.md](docs/nmf_stability.md).

## Installation

See [docs/install.md](docs/install.md) for system dependencies and R/Python
binding installation. The repository includes CMake presets, a vcpkg manifest,
and GitHub Actions smoke builds for the C++ core plus R and Python bindings.

## R Bindings

The R interface is documented in [docs/r-bindings.md](docs/r-bindings.md).

Minimal R workflow:

```r
library(cellAdmixCore)

ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)
fit <- ds$fit()
score <- fit$score_membrane()
rules <- score$rules()
correction <- score$correct(rules = rules)
corrected_counts <- correction$counts()
```

The R documentation page links to more detailed notes on inputs, scoring
methods, batch processing, and example notebooks.

## Python Bindings

The Python interface is documented in
[docs/python-bindings.md](docs/python-bindings.md). It is a thin pybind11 layer
over the same C++ analysis store and run artifacts used by the R package.

Minimal Python workflow:

```python
import celladmix as ca

ds = ca.CellAdmix("data", output_dir="out", annotation=cell_annotation)
fit = ds.fit()
score = fit.score_membrane()
rules = score.rules()
correction = score.correct(rules=rules)
counts, genes, cells = correction.counts()
```

Batch processing from Python is available through
[`scripts/celladmix_batch.py`](scripts/celladmix_batch.py). The current Python
batch entrypoint supports Xenium bundle inputs and mirrors the core R batch
workflow:

```bash
python scripts/celladmix_batch.py --input data --output out \
  --annotation annotations/annotation.csv.gz \
  --annotation-col merged_annotation \
  --score auto --report --threads 10
```

See [docs/python-batch.md](docs/python-batch.md) for options and current
Python-side limitations.

## Examples

Rendered notebooks with embedded output are grouped by binding.

### R Examples

- [Pancreas quickstart](examples/xenium_pancreas_membrane_377_full/pancreas_quickstart.ipynb): the recommended workflow in its shortest form - audit, fit, correct, verify.
- [Minimal CosMx NSCLC tutorial](examples/cosmx_nsclc_giotto/celladmix_cosmx_minimal.ipynb): a compact tabular example mirroring the original cellAdmix [NSCLC tutorial](https://github.com/kharchenkolab/cellAdmix/blob/main/vignettes/NSCLC_tutorial_fulldata.ipynb).
- [Detailed Xenium pancreas tutorial](examples/xenium_pancreas_membrane_377_full/pancreas_membrane_scoring_clean.ipynb): the full walkthrough on a membrane-stained Xenium bundle - audit, membrane and bridge scoring on their recommended factorizations (`invsqrt_kl` and `ls_nmf`), and the variant comparison behind that pairing.
- [Seurat Xenium integration tutorial](examples/xenium_pancreas_membrane_377_full/pancreas_seurat_integration.ipynb): the same pancreas dataset, using Seurat for cell-level state and cellAdmix for molecule-complete fitting, scoring, and correction.
- [Xenium breast 5K membrane scoring tutorial](examples/xenium_breast_membrane_5k_full/breast_5k_membrane_scoring.ipynb): the full-scale 5K-panel workflow behind the timing benchmark above.

### Python Examples

- [Standalone Xenium pancreas tutorial](examples/xenium_pancreas_membrane_377_full/pancreas_python_standalone.ipynb): the pancreas dataset through the standalone Python API.
- [SpatialData Xenium integration tutorial](examples/xenium_pancreas_membrane_377_full/pancreas_spatialdata_integration.ipynb): the same pancreas dataset through a SpatialData-first workflow.
- [Standalone Xenium breast 5K tutorial](examples/xenium_breast_membrane_5k_full/breast_5k_python_standalone.ipynb): the same core workflow at 5K-panel scale.

## Documentation

- [docs/README.md](docs/README.md): documentation index.
- [docs/install.md](docs/install.md): installation instructions.
- [docs/r-bindings.md](docs/r-bindings.md): R interface overview.
- [docs/inputs.md](docs/inputs.md): supported input formats.
- [docs/scoring_methods.md](docs/scoring_methods.md): bridge, membrane, and coherence scoring definitions.
- [docs/nmf_stability.md](docs/nmf_stability.md): the per-factor restart stability diagnostic.
- [docs/benchmarks.md](docs/benchmarks.md): the cleanup benchmark methodology and factorization/scoring comparisons.
- [docs/python-bindings.md](docs/python-bindings.md): Python interface and SpatialData integration overview.
- [docs/python-batch.md](docs/python-batch.md): Python command-line batch workflow.

## Status

This is an initial public release. The C++ core, R bindings, and Python
bindings are usable, but interfaces may still evolve as the package matures.
