---
name: celladmix-r-api
description: The R user API of cellAdmixCore - the cellAdmix() constructor, R6 dataset/fit/score/correction/audit objects, their methods and default parameters, plotting helpers, Seurat export, caching, and the output directory layout. Use when writing or reading R code that drives the pipeline.
---

Use this skill when writing R code against `cellAdmixCore`, when you need a method name or a default, or when you must explain what a run directory contains. Sources: `r/NAMESPACE`, `r/R/constructor.R`, `dataset.R`, `fit.R`, `score.R`, `correction.R`, `audit.R`, `seurat.R`, `api.R` (backend defaults), `r/man/*.Rd`, and `docs/r-bindings.md`. Scoring semantics are in `celladmix-scoring-and-correction`; audit details in `celladmix-audit-and-stability`.

## Minimal workflow

```r
library(cellAdmixCore)
ds <- cellAdmix("xenium_bundle", output_dir = "out", annotation = annotation)
fit <- ds$fit()                              # NMF, molecule labels, restart pool
audit <- fit$audit_admixture()               # factorization-independent admixture estimate
score <- fit$score_membrane()                # or fit$score_bridge(), fit$score_coherence()
rules <- score$rules()                       # p_thresh = 0.1, native-factor check applied
correction <- score$correct(rules = rules)   # molecule-vote ensemble, vote = 0.3
counts <- correction$counts()                # sparse gene x cell dgCMatrix
report <- audit$evaluate(correction); report$summary()
```

## Constructor

`cellAdmix(source, output_dir, annotation, format = c("auto","xenium","tabular","prep","store","seurat"), schema, annotation_name = "manual", annotation_col, cell_id_col = "cell_id", num_threads, molecules, xenium_dir, assay, image, ...)` returns a `CellAdmixDataset`. `...` goes to `celladmix_prepare_xenium()` (`analysis_bbox`, `analysis_crops`, `cell_filter`, `gene_filter`, `min_qv`, `keep_unassigned = FALSE`, `keep_non_gene = FALSE`, `read_cells = TRUE`, `prefer_parquet = TRUE`) or `celladmix_prepare_tabular()`. `num_threads` defaults to `CELLADMIX_NUM_THREADS` or the core count capped at 10. `output_dir = NULL` creates a temporary directory. Inputs and schemas: `celladmix-inputs`.

## CellAdmixDataset

Fields: `prep`, `format`, `output_dir`, `source_info`, `active_annotation`, `active_domain`, `fits`, `num_threads`.

- Annotations: `annotations()`, `annotation(name = NULL, as_vector = FALSE)`, `set_annotation(annotation, name, active = TRUE, annotation_col, cell_id_col)`, `use_annotation(name)`, `drop_annotation(name)`.
- Store and cell-level data: `ensure_store(required = c("molecules","counts"), ...)`, `counts()`, `cells()` (cell id and centroid), `cell_state_umap(...)`.
- Clustering and domains: `cluster(name = "louvain", active = NULL, ...)` (forwards to `celladmix_cluster_cells()`: `min_molecules = 10`, `min_genes = 5`, `n_variable_genes = 1000`, `pca_dims = 30`, `graph_k = 15`, `resolution = 1`, `compute_umap = TRUE`); `identify_domains(name, annotation, active = TRUE, ...)`, `load_domains(name)`, `domains(name = NULL)`.
- Fits: `fit(annotation = NULL, rank = "auto", rank_multiplier = 1.2, rank_cap = 30L, run_id = NULL, overwrite = FALSE, ...)`, `list_fits()`, `get_fit(run_id)`, `read_fit(run_dir, annotation = NULL, name = NULL)` to reattach a persisted run.
- Audit without a fit: `audit_admixture(annotation = NULL, ...)`.

`fit()` `...` reaches `celladmix_fit()`: `ncv_k = NULL` (auto), `graph_k = 10`, `same_label_ratio = 5`, `nmf_iterations = 150`, `nmf_init = "auto"|"random"|"cluster"`, `nmf_variant = "ls_nmf"|"invsqrt_kl"|"kl"|"sqrt_kl"`, `molecule_scoring = "gene_loadings"|"ncv_projection"|"auto"`, `nmf_n_runs = NA` (auto: at least 10, one per thread above that), `nmf_train_max_rows = 10000`, `nmf_min_molecules = 10`, `seed = 1`, `verbose = FALSE`, `nmf_fixed_h`, `nmf_fixed_h_members` (fixed-profile projection; see `docs/multisample.md`). Auto rank is `min(rank_cap, max(2, ceiling(rank_multiplier * n_labels)))`.

## CellAdmixFit

Fields: `dataset`, `run`, `run_id`, `run_dir`, `annotation_name`, `annotation_hash`, `rank`, `scores_registry`, `corrections_registry`, `params`.

- Data: `loadings()` (factor x gene matrix), `cell_factors()` (per cell: `dominant_factor`, `factor_<k>_fraction`, `transcript_count`), `counts(...)`, `molecules(...)`, `training_molecules()`, `region(bbox, ...)`, `stability(importance = "molecule_fraction")`.
- Images and boundaries: `stain(stain = "membrane", ...)`, `stain_image()`, `stains(names = c("dapi","membrane"))`, `stain_crop(image, bbox)`, `cell_boundaries(boundary_path, cells, bbox)`.
- Scoring: `score(method = c("membrane","bridge","coherence"), name, annotation, ...)`, `score_membrane(name = "membrane", annotation, ...)`, `score_bridge(name = "bridge", ...)`, `score_coherence(name = "coherence", ...)`, `score_factor_sources(name = "factor_sources", annotation, counts, ...)`, `scores()`, `get_score(name)`.
- Correction: `correct(score, name = NULL, p_thresh = 0.1, adjust_p = FALSE, targets = NULL, rules = NULL, ensemble = NULL, vote = 0.3, ...)`; `corrections()`, `get_correction(name)`.
- Audit: `audit_admixture(neighbor_k = 15, n_pool = 20, q_thresh = 0.01, min_excess = 200, min_target_cells = 200, min_reference_cells = 100)`.
- Plots: `plot_loadings(n_genes = 10, ncol = 5)`, `plot_stability()`, `prepare_cell_example()`, `plot_cell_example()`.

## CellAdmixScore

Fields: `name`, `method`, `fit`, `result`, `annotation_name`, `paths`, `params`.

- `summary()` (target/source/factor table with `mean_score`, `p_value`, `neg_log10_p`, method-specific columns), `pairs()` (per cell pair or cell), `annotation(p_thresh = 0.1, adjust_p = FALSE, ...)` (source calls and removal calls).
- `rules(p_thresh = 0.1, adjust_p = FALSE, targets = NULL, native_check = TRUE, native_median_thresh = 0.1, native_expr_thresh = 0.05, native_outlier_min_frac = 0.1, neighbor_k = 15L, ...)` returns `factor`, `source_cell_type`, `target_cell_type`, `p_value`, `neg_log10_p`, `rule_id`, and the check columns `keep`, `native_check`, `native_distant_n`, `native_exposed_n`, `native_distant_median`, `native_exposure_gradient`, `native_distant_expr_frac`.
- `correct(name, p_thresh, adjust_p, targets, rules, ensemble = NULL, vote = 0.3, ...)` delegates to `fit$correct()`; the correction is named `<score>_clean` by default.
- Plots: `plot_heatmap(..., p_thresh = 0.1, adjust_p = FALSE, native_check = TRUE)` (removal calls marked `*`, vetoed ones `(*)`), `plot_pairs(...)`, `examples(...)`, `plot_example(example, ...)`, `plot_examples(examples = NULL, targets, n_per_target = 4, cells, ncol = 2, ...)` with automatic stain backgrounds, boundaries, and cell-type contours for Xenium runs (`stains`, `boundaries`, `cell_types`, `color_cell_types` override discovery).

`CellAdmixFactorSourceScore` (from `score_factor_sources()`): `scores()`, `marker_weights()`, `top_markers(n_genes = 10)`, `top_factor_genes(n_genes = 8, sources)`, `annotation(min_score = 0.15, min_margin = 0.03)`, `plot_heatmap()`, `plot_top_genes()`.

## CellAdmixCorrection

Fields: `name`, `fit`, `score_name`, `run`, `params`, `rules` (the applied rules; with a `support` column in ensemble corrections).

`counts(...)` (sparse `dgCMatrix`, genes x cells), `cells()`, `molecules(...)`, `summary()` (per cell type: `n_cells`, `n_modified_cells`, `molecules_before`, `molecules_after`, `molecules_removed`, `fraction_molecules_removed`, medians per modified cell; first row `all`), `cell_summary()` (per cell), `ensemble()` (`members`, `vote`, `min_votes`, `n_removed`, `vote_histogram`), `plot_removed_molecules()`, `cell_state_umap()`. `correct()` warns when more than 60% of a cell type's molecules are removed (gate `options(celladmix.overremoval_min_molecules = )`).

## Audit objects

`CellAdmixAudit`: `pairs(detected_only = FALSE)`, `markers(source, target)`, `plot_map(value = c("rate","molecules"), detected_only = TRUE)`, `plot_exposure(source = NULL, target = NULL, correction = NULL, strict = TRUE)`, `plot_remaining(corrections = list())`, `plot_reference(source, target)`, `evaluate(correction, warn_uncovered = TRUE)`, `fit_generative(init = NULL, n_programs = 4L, seed = 1L, num_threads = NULL, ...)`, `correct_generative(name = "generative", ...)`. `CellAdmixCleanupReport`: `pairs()`, `false_removal()`, `summary()`, `plot_cleanup()`. See `celladmix-audit-and-stability` and `celladmix-generative`.

## Multi-sample objects

`cellAdmixCollection(datasets, sample_allocation = c("proportional","equal"))` and the `CellAdmixCollection`, `CellAdmixConsensusFit`, `CellAdmixCollectionScore`, `CellAdmixCollectionCorrection`, `CellAdmixCollectionAudit`, `CellAdmixCollectionCleanupReport` classes: see `celladmix-multisample` and `docs/multisample.md`.

## Seurat and other exported helpers

```r
seu <- celladmix_add_to_seurat(seu, fit = fit, score = score, correction = correction,
  corrected_assay = "cellAdmix", prefix = "celladmix_")
```

`celladmix_add_factors_to_seurat(object, fit, prefix)` adds `dominant_factor` and `factor_<k>_fraction` metadata; `celladmix_add_corrected_assay(object, correction, assay = "cellAdmix")` adds the corrected counts, reordered to the object's cells. Differential expression and impact: `celladmix_de()`, `celladmix_de_by_group()`, `celladmix_correction_impact()`, `celladmix_plot_volcano()`, `celladmix_plot_marker_expression()`, `celladmix_plot_de_shift()`. Score views: `celladmix_plot_score_pairs()`, `celladmix_plot_score_agreement()`. Spatial and domain plots: `celladmix_plot_spatial(data, color_by, x = "x", y = "y", ...)`, `celladmix_plot_domains_spatial()`, `celladmix_plot_domain_composition()`, `celladmix_plot_domain_umap()`, `celladmix_domain_composition()`, `celladmix_domain_diagnostics()`. Cluster reports: `celladmix_score_cluster_markers()`, `celladmix_top_cluster_markers()`, `celladmix_plot_cluster_marker_heatmap()`, `_dotplot()`, `_umap()`. Example cells: `celladmix_select_example_cells()`, `celladmix_prepare_cell_example()`, `celladmix_plot_cell_example()`, `celladmix_discover_cell_boundaries()`, `celladmix_read_cell_boundaries()`, `celladmix_discover_stain_image(run, stain = "membrane", image_path, focus_index, pixel_size, x_offset = 0, y_offset = 0)`. `celladmix_core_version()` reports the compiled core version. Functions named `celladmix_prepare_*`, `celladmix_fit`, `celladmix_score_*`, `celladmix_correct` are internal backends.

## Caching and output directory layout

```
out/
  prep.json                       # prepared-project description
  annotations/r6/<name>/          # annotation sidecars
  input_store/                    # encoded molecules, cells, genes, counts
  runs/<run_id>/                  # run_id = fit_<annotation>_rank<K>_<variant>
    run.json                      # manifest: options, genes, nmf_diagnostics, counts
    factors.parquet  cells.parquet  molecules.parquet  training_rows.parquet
    ensemble_h.parquet            # restart member pool (gene profiles)
    ensemble_labels_m<k>.parquet  # per-member molecule labels (ensemble correction)
    fit_metadata.json
    scores/<score_name>/          # score tables and score_metadata.json
    corrected/<correction_name>/  # corrected run: run.json, molecules.parquet,
                                  # cells.parquet, correction_summary.parquet,
                                  # correction_metadata.json
```

`fit()` compares the requested parameters and the annotation hash with `run.json`; a match prints `Reusing cached run ...`, a difference prints the changed parameters and refits, `overwrite = TRUE` always refits. Omitted auto-resolved values (`ncv_k`, `nmf_n_runs`) match whatever the cached run recorded. Scores and corrections are recomputed on every call. Per-member labelings for the ensemble are computed once per run and reused.

## Batch script

`Rscript scripts/celladmix_batch.R --input data --output out --annotation annotation.csv.gz --annotation-col merged_annotation --score auto --report --threads 10`; `--help` lists all options (`--score auto` uses membrane scoring when a stain is found, otherwise bridge; `--auto-annotate` clusters when no annotation is given; `--ensemble`, `--vote`, `--nmf-variant` default `ls_nmf`).
