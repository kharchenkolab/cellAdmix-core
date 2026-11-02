---
name: celladmix-core-cpp
description: The C++ core of cellAdmix-core - the main headers and what each does, the run store layout on disk (manifest and parquet files), the NMF variants with their objectives and feature spaces, how the R and Python bindings call the core, and how to add and run C++ tests. Use when reading or changing code under include/, src/, or tests/.
---

Use this skill when a task touches the C++ implementation: changing an algorithm, adding a core function and exposing it to R or Python, reading or writing run artifacts, or adding a unit test. The code is in `include/celladmix/` and `src/`; the bindings are `r/src/bindings.cpp` (Rcpp) and `python/src/celladmix/_core_bindings.cpp` (pybind11); tests are in `tests/`. Build instructions are in `docs/install.md` and the `celladmix-overview` skill.

## Headers and what they do

Data model and utilities:
- `types.hpp`: `TranscriptTable` (molecule rows with coordinates, genes, cell ids, optional z, qv, cell types, sample and FOV ids, nucleus flags; `finalize()` builds gene and cell dictionaries), `CellTable`, centroid helpers, `kCelladmixVersion`.
- `matrix.hpp`: `DenseMatrix` (row-major) and `SparseRowMatrix` (CSR) with small linear-algebra helpers.
- `spatial.hpp`: exact molecule k-nearest-neighbor index (nanoflann) used by neighborhood construction and smoothing. `hnsw.hpp`: approximate neighbors for cell and report UMAPs.
- `graph.hpp`: within-cell molecule kNN graph, iterative-conditional-modes label smoothing, `assign_factors_per_cell` (with prebuilt `CellGraphs` for ensemble members).

Loaders and stores:
- `xenium.hpp`: manifest parsing, `load_xenium_bundle`, crop boxes, control-feature classification (`xenium_is_biological_feature`).
- `tabular.hpp`, `tabular_mask.hpp`: `TabularSourceSpec` column mapping, CSV/Parquet loading, TIFF mask reading and cell assignment.
- `input_store.hpp`: `build_xenium_input_store`, `build_tabular_input_store`, manifest and paths, block-wise molecule access (`load_input_store_molecule_range`, `load_input_store_cell_offsets`), counts and cells loaders.
- `run_store.hpp`: persisted run layout, `RunManifest`, `write_basic_run`, `write_corrected_run`, `load_run_data`, `collect_run_counts`, ensemble member pool and labels (`write_ensemble_h_parquet`, `ensure_ensemble_labels`, `load_ensemble_member_labels`).

Factorization and labeling:
- `ncv.hpp`: neighborhood composition vectors (`NcvOptions`: `k = 20`, `within_cell`, `include_self`), the KL feature transform (`make_ncv_feature_transform`), projections of rows or gene loadings onto fixed factors.
- `nmf_euclidean.hpp`: weighted least-squares NMF (`weighted_nmf`, `WeightedNmfOptions`), multirun selection, member pool, `weighted_ls_usages_fixed_h`.
- `nmf_kl.hpp`: sparse NMF with KL or Euclidean loss (`sparse_nmf`, `SparseNmfOptions`: `init_mode`, `loss_mode`, `h_l1_penalty`, `h_diversity_penalty`), `order_nmf_factors_by_training_importance`.
- `nmf_stability.hpp`: ownership profiles, Hungarian matching, `match_factor_profiles`, `matched_ownership_stability`, `kStableFactorThreshold = 0.3`.
- `pipeline_common.hpp`: `resolve_auto_ncv_k`, `resolve_molecule_scoring`, `is_weighted_ls_variant`, `align_profiles_to_genes`, conversion of weighted LS results to the common `SparseNmfResult`.
- `pipeline.hpp`: `BasicPipelineOptions` (`ncv_k`, `rank`, `graph_k = 10`, `same_label_ratio = 5`, `nmf_iterations = 150`, `nmf_n_runs`, `nmf_init`, `nmf_variant = "ls_nmf"`, `molecule_scoring = "gene_loadings"`, `nmf_train_max_rows`, `nmf_min_molecules = 10`, `seed`, `training_cell_strata`, `training_scope_cell_types`, `nmf_fixed_h`, `nmf_fixed_h_members`, `annotation_hash`), `run_basic_pipeline` on an in-memory table, `apply_removal_rule`, `build_cell_gene_counts`.
- `pipeline_store.hpp`: `run_basic_pipeline_store` (the production path: processes the input store by cell blocks and writes the run), `training_loss_under_profiles`.

Scoring, audit model, and auxiliary analyses:
- `bridge.hpp` (`BridgeTestOptions`, `run_bridge_test`), `membrane.hpp` (`Image2D`, `MembraneImageOptions`, `MembraneTestOptions`, `run_membrane_test`, stain window reading), `coherence.hpp` (`CoherenceTestOptions`, `run_coherence_test`). Struct defaults differ from the binding defaults in two places (bridge `candidate_mode = "molecule_global"`, `max_cells_per_type_pair = 200`); the R and Python wrappers pass `"cell_center"` and 400, so change user-facing defaults in the bindings.
- `clustering.hpp`: cell clustering and UMAP for provisional annotations, `cell_neighbor_type_counts`, `cell_nearest_type_distance` (used by the native-factor check and the audit).
- `generative.hpp`: `GenerativeOptions`, `GenerativePairSpec`, `GenerativeResult`, `fit_generative` (the per-target-type Poisson mixture of `docs/generative.md`), `pava_increasing`.
- `domains.hpp` (tissue domains from cell positions and annotations), `simulation.hpp` (`simulate_nsclc_admixture` synthetic data for tests), `report_embedding.hpp` (training-molecule UMAP sidecar), `workflow.hpp` (table subsetting, `summarize_cells`).

Vendored header libraries under `include/`: aarand, irlba, kmeans, knncolle, subpar, umappp, `third_party/hnswlib`, `third_party/nanoflann.hpp` (see `THIRD_PARTY.md`).

## Run store layout

`make_run_paths(root)` defines a run directory:

```
<run>/run.json                 RunManifest (format_version "0.2", package_version, run_type "fit" or
                               corrected, annotation_hash, fixed_profiles, fixed_profile_members,
                               source, paths, pipeline_options, storage_options, genes,
                               nmf_gene_weights, nmf_diagnostics, n_transcripts, n_cells, n_factors,
                               n_training_rows, has_* column flags, parent_run for corrected runs)
<run>/factors.parquet          gene loadings per factor
<run>/cells.parquet            per-cell summary: dominant factor, factor fractions, counts, centroid
<run>/molecules.parquet        molecule rows with zero-based factor label and margin, tiled by tile_size
<run>/training_rows.parquet    ids of the sampled neighborhood training rows
<run>/ensemble_h.parquet       restart member pool (one full-gene-space H per restart)
<run>/ensemble_labels_m<k>.parquet   member k's labels aligned row-for-row with molecules.parquet
<run>/scores/                  score outputs written by the bindings
<run>/corrected/<name>/        corrected child run (same layout, plus correction_summary.parquet,
                               vote statistics in run.json)
```

The input store (`input_store.hpp`): `store.json`, `genes.parquet`, `cells.parquet`, `cell_gene_counts.parquet`, `molecules.parquet`, `cell_offsets.parquet`. Parquet is read and written through Arrow; `RunStorageOptions` sets `tile_size` (100) and `parquet_row_group_size` (65536). `RunLoadOptions.ensemble_member >= 0` swaps molecule labels for a member's when loading a run for scoring.

## NMF variants

`nmf_variant` selects both the objective and the feature space of the training rows:

| variant | objective | training features | molecule potentials |
|---|---|---|---|
| `ls_nmf` | weighted least squares (`weighted_nmf`, inverse column-frequency gene weights) | raw neighborhood counts | gene loadings (default `molecule_scoring`) |
| `kl` | generalized KL divergence (`sparse_nmf`, `loss_mode = "kl"`) | raw counts | gene loadings; `ncv_projection` solves fixed-H KL per molecule |
| `sqrt_kl` | KL | square-root transformed counts (`make_ncv_feature_transform`) | as above |
| `invsqrt_kl` | KL | counts scaled per gene by inverse square-root prevalence; weights stored in `nmf_gene_weights`, row sums rescaled to `nmf_transform_target_row_sum` | as above |

`is_weighted_ls_variant()` is true only for `ls_nmf`. Fixed-profile runs (`nmf_fixed_h` non-empty) skip factorization, label against the given profiles in order, and store `nmf_fixed_h_members` as member pool with the fixed profiles as member 0. All multirun results are converted to `SparseNmfResult`, so downstream code sees one shape (`h`, `w`, `candidate_h` pool, stability fields).

## How the bindings call the core

- R: `r/src/bindings.cpp` exports `.Call` entry points named `_cellAdmixCore_celladmix_*` (store building, `celladmix_fit_store_run`, `celladmix_bridge_scores_run`, `celladmix_membrane_scores_run`, `celladmix_coherence_scores_run`, `celladmix_correct_run`, `celladmix_ensemble_prepare`, `celladmix_ensemble_member_fractions`, `celladmix_fit_generative`, `celladmix_match_factors`, `celladmix_training_loss`, collectors for cells, counts, transcripts, and summaries). `r/R/native.R` wraps them; `r/R/api.R` adds defaults and R-side tables; the R6 classes in `r/R/*.R` sit on top. `r/src/Makevars` compiles `../../src/*.cpp` with `-I../../include`, so a core change needs `R CMD INSTALL r`.
- Python: `_core_bindings.cpp` defines the `celladmix._core` module (`build_xenium_store`, `fit_store`, `score_bridge`, `score_membrane`, `correct_run`, `ensemble_prepare`, `ensemble_h`, `ensemble_member_fractions`, `fit_generative`, `match_factors`, `training_loss`, `cluster_store_counts`, `cluster_counts_matrix`, `cluster_run_counts`, `cluster_factors`, `cell_neighbor_type_counts`, `cell_nearest_type_distance`, `collect_counts`, `collect_input_store_counts`, `read_run_manifest`, `read_input_store_manifest`, `load_run_cells_path`, `training_strata_for_cells`, `core_version`). Keyword defaults live in the `py::arg(...) = ...` lists. scikit-build-core builds it via the top-level CMake with `CELLADMIX_BUILD_PYTHON=ON` (`python/pyproject.toml`), so a core change needs `pip install -e python --no-build-isolation`.
- Both bindings pass the annotation as parallel `cell_ids`/`labels` vectors at scoring time, so a stored run can be scored against different labels without refitting. Rule-level logic (rule derivation, native-factor check, ensemble member rules, audit statistics) is implemented in R and Python separately; keep the two in step when changing either.

Adding a core function: declare it in the matching header, implement it in `src/`, add it to `add_library(celladmix_core ...)` in `CMakeLists.txt` and to `CORE_SOURCES` in `r/src/Makevars` if it is a new file, export it in both bindings, and document it.

## Tests

`tests/test_framework.hpp` is a minimal registry: `TEST_CASE("name") { ... }` registers a function; assertions are `REQUIRE(cond)`, `REQUIRE_NEAR(a, b, tol)`, `REQUIRE_EQ/GT/GE/LT/LE`. `tests/test_main.cpp` runs every registered case, prints `[PASS]`/`[FAIL]`, and returns nonzero on failure. Test files: `test_bridge`, `test_clustering`, `test_generative`, `test_graph`, `test_membrane`, `test_ncv`, `test_nmf`, `test_nmf_stability`, `test_run_store`, `test_simulation`, `test_tabular`, `test_workflow`, `test_xenium`. Most build a small dataset with `simulate_nsclc_admixture()` and write to a temporary directory.

To add a test, create or extend a `tests/test_<topic>.cpp`, include `"test_framework.hpp"` and the headers under test, and list a new file in `add_executable(celladmix_core_tests ...)` in `CMakeLists.txt`. Build and run:

```bash
cmake --preset tests && cmake --build --preset tests -j8
ctest --preset tests                 # or: build/tests/celladmix_core_tests
```

The binary has no case filter; it runs all registered cases. CI (`.github/workflows/`) builds the core on Ubuntu, macOS, and Windows, and runs the tests plus R and Python install smoke jobs on Ubuntu.
