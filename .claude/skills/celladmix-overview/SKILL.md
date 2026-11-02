---
name: celladmix-overview
description: Orientation for the cellAdmix-core repository - what the package does, the pipeline stages, where code and docs live, how to build and test each layer, and the project conventions. Use first when starting any task in this repo.
---

Use this skill when you begin work in the cellAdmix-core repository and need to know what the software does, which stage of the pipeline a task touches, where the authoritative documentation and the code for each layer live, how to build and test the C++ core, the R package, and the Python package, and which writing and commit conventions apply. The other `celladmix-*` skills go deeper on inputs, the two binding APIs, scoring and correction, the audit, the C++ core, multi-sample processing, and the generative model.

## What the package does

cellAdmix-core detects and removes molecular admixture in imaging spatial transcriptomics (10x Xenium, CosMx): molecules that were transcribed in one cell but assigned by segmentation to a neighboring cell. It is a C++ core (`include/`, `src/`) with an R package (`r/`, package name `cellAdmixCore`) and a Python package (`python/`, package name `celladmix`). Both bindings drive the same store-backed C++ pipeline and read the same on-disk run artifacts, so results are identical across languages.

## Pipeline stages

1. **Input store.** The molecule table of a Xenium bundle or a tabular file is streamed into an encoded Parquet store under `<output_dir>/input_store/` (`store.json`, `genes.parquet`, `cells.parquet`, `cell_gene_counts.parquet`, `molecules.parquet`, `cell_offsets.parquet`). Xenium control, codeword, and other non-gene features are dropped at this step (`keep_non_gene = FALSE`). See `celladmix-inputs`.
2. **Neighborhood composition vectors.** For each molecule, the gene composition of its `ncv_k` nearest molecules within the same cell. `ncv_k` is resolved from the panel size (20 for a ~400-gene panel, growing with the square root of the gene count, capped at half the median cell's molecule count).
3. **NMF with restarts.** Up to `nmf_train_max_rows` (10,000) neighborhood rows, stratified by annotated cell type, are factorized at rank `ceil(1.2 * number of cell types)` (capped at 30). Variants: `ls_nmf` (default, weighted least squares), `invsqrt_kl`, `kl`, `sqrt_kl`. At least 10 random restarts run by default; the selected restart defines the factors, a per-factor restart stability is recorded, and every restart's gene profiles are kept as the ensemble member pool (`ensemble_h.parquet`).
4. **Molecule labeling.** Every molecule receives factor potentials from the gene loadings (`molecule_scoring = "gene_loadings"`) and a hard label after smoothing on a within-cell k-nearest-neighbor graph (`graph_k = 10`, `same_label_ratio = 5`). Cell-level factor fractions summarize the labels.
5. **Scoring.** Each (target cell type, source cell type, factor) triple is tested for spatial evidence of admixture by bridge scoring (molecule contacts across the cell boundary), membrane scoring (stain intensity along source-facing rays), or coherence scoring (within-cell patches), plus a factor source score that matches factors to annotated cell types by marker similarity.
6. **Rules and the native-factor check.** Removal rules (`factor`, `source_cell_type`, `target_cell_type`) are the triples passing `p_thresh` (0.1). Each rule is then tested against target cells with no source-type neighbors; a factor that persists there is native expression and the rule is kept in the table with `keep = FALSE`.
7. **Molecule-vote ensemble correction.** Each NMF restart labels and scores the molecules independently; a molecule is removed when at least `vote` (0.3) of up to `ensemble` (10) members remove it. The corrected run holds corrected molecules, sparse counts, and a per-cell removal summary.
8. **Audit.** Independently of the factorization, source-marker content in target cells is compared across source-neighbor exposure bins to estimate per-pair admixture rates and to verify a correction's sensitivity and own-marker false removal. The audit also fits a generative model that separates transferred material from neighborhood-induced expression. See `celladmix-audit-and-stability` and `celladmix-generative`.

## Repository layout

- `include/celladmix/*.hpp`, `src/*.cpp`: the C++ core (see `celladmix-core-cpp`). `include/` also vendors header-only libraries (aarand, irlba, kmeans, knncolle, subpar, umappp, hnswlib, nanoflann).
- `r/`: R package. `r/R/*.R` (R6 classes and backend wrappers), `r/man/*.Rd`, `r/src/bindings.cpp` plus a `Makevars` that compiles `../../src/*.cpp`, `r/tests/testthat/`.
- `python/`: Python package built with scikit-build-core. `python/src/celladmix/*.py`, `python/src/celladmix/_core_bindings.cpp` (pybind11 module `celladmix._core`), `python/tests/`.
- `docs/`: user documentation. `docs/development/` is gitignored design material.
- `tests/`: C++ unit tests with the minimal framework in `tests/test_framework.hpp`.
- `examples/`: rendered tutorials per dataset (`xenium_pancreas_membrane_377_full/`, `xenium_breast_membrane_5k_full/`, `cosmx_nsclc_giotto/`); R tutorials are `.Rmd` with rendered `.md`, Python tutorials are `.ipynb`. Data and output folders inside are gitignored.
- `scripts/`: `celladmix_batch.R` and `celladmix_batch.py` command-line workflows.
- `analysis/`: evaluation harnesses behind the documentation (`cleanup_benchmark/` for the audit comparisons, `audit_guided/`, `generative_model/`, `multisample/`). Not part of the packages.
- `CMakeLists.txt`, `CMakePresets.json`, `cmake/CellAdmixDependencies.cmake`, `vcpkg.json`, `.github/workflows/`.

## Authoritative documentation

- `README.md`: overview, highlights, NMF variants, minimal R and Python workflows, example index.
- `docs/README.md`: documentation index.
- `docs/install.md`: system dependencies, R and Python installation, CMake presets.
- `docs/inputs.md`: Xenium and tabular inputs, schema fields, stain-image discovery, filtering.
- `docs/r-bindings.md`: R workflow, objects, fitting, scoring, audit, correction, Seurat, batch.
- `docs/python-bindings.md`: Python workflow, run reattachment, factor indexing, audit, SpatialData.
- `docs/python-batch.md`: options and outputs of `scripts/celladmix_batch.py`.
- `docs/scoring_methods.md`: definitions of bridge, membrane, coherence, and factor source scores and the native-factor check.
- `docs/nmf_stability.md`: the per-factor restart stability diagnostic.
- `docs/audit.md`: the admixture audit, the stochasticity findings, and the ensemble correction.
- `docs/generative.md`: the generative model of admixture and induced expression.
- `docs/multisample.md`: consensus factorization and pooled rules across samples.

## Build and test

C++ core and unit tests (CMake 3.20+ for presets; dependencies in `docs/install.md`):

```bash
cmake --preset tests
cmake --build --preset tests -j8
ctest --preset tests            # or run build/tests/celladmix_core_tests directly
```

`cmake --preset user` builds the optimized library without tests; `user-vcpkg` resolves dependencies through vcpkg.

R package (compiles the core from the repository root, so install from a checkout):

```bash
MAKEFLAGS=-j8 R CMD INSTALL r
Rscript -e 'library(cellAdmixCore); packageVersion("cellAdmixCore")'
Rscript -e 'library(cellAdmixCore); testthat::test_dir("r/tests/testthat")'
```

Python package (Python 3.8+; SpatialData extras need 3.11+):

```bash
python -m pip install -e python                        # isolated build fetches build tools
python -m pip install scikit-build-core pybind11 cmake ninja
python -m pip install -e python --no-build-isolation   # faster repeated rebuilds
python -m pip install -e "python[plot,spatialdata]"    # optional extras
python -c "import celladmix as ca; print(ca.__version__)"
python -m pytest python/tests
```

Set `CMAKE_PREFIX_PATH` (CMake, Python) or `PKG_CONFIG_PATH` (R) when dependencies live in a custom prefix.

## Project conventions (from `CLAUDE.md`)

- Commits: a single short imperative subject line; no body unless essential; no Co-Authored-By or other trailers.
- Documentation is kept coherently current: when behavior, defaults, or APIs change, rework `README.md`, `docs/*.md`, roxygen, and docstrings in the same change so each page reads as one description of the present code. No change journals, no "new in this version" or "previously/now" notes.
- Plain scientific language: describe methods and quantities by what they are or do; do not coin shorthand labels. If a shorthand is needed, define it at first use.

## Local environment on this machine

- R 4.4.3 is at `/usr/local/R/R-4.4.3/bin`. The worktree uses a local R library `.Rlib`; install the package with `R_LIBS=<worktree>/.Rlib R CMD INSTALL --library=<worktree>/.Rlib r` (add `MAKEFLAGS=-j12` to speed compilation) and load it with `R_LIBS=<worktree>/.Rlib` set.
- The Python virtual environment is `.venv-sdata` in the worktree (Python 3.11); use `.venv-sdata/bin/python -m pip install -e python --no-build-isolation` and `.venv-sdata/bin/python -m pytest python/tests`.
- The C++ build directory is `build/tests`, produced by `cmake --preset tests`.
- Long-running R jobs read their script files lazily: do not edit a script while a job started from it is still running, and do not reinstall the package into a library (`.Rlib`) that running jobs have loaded; install and test development code in a second library (`.Rlib_dev`) until the jobs finish.
- Branch layout of the local checkouts: `dev` carries the audit/ensemble line, `dev_model` the generative-model work (the main checkout `cellAdmix-core`), and `dev_multisample` the multi-sample work (worktree `wt-dev_multisample`, branched from `dev`).
