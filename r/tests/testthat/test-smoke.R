library(cellAdmixCore)

celladmix_simulate_nsclc <- cellAdmixCore:::celladmix_simulate_nsclc
celladmix_resolve_rank <- cellAdmixCore:::.celladmix_resolve_rank
celladmix_weighted_nmf_matrix <- cellAdmixCore:::.celladmix_weighted_nmf_matrix
celladmix_build_reference_ncv_matrix <- cellAdmixCore:::.celladmix_build_reference_ncv_matrix
celladmix_build_pipeline_ncv_matrix <- cellAdmixCore:::.celladmix_build_pipeline_ncv_matrix
celladmix_bridge_annotation <- cellAdmixCore:::celladmix_bridge_annotation
celladmix_cluster_cells <- cellAdmixCore:::celladmix_cluster_cells
celladmix_collect_cells <- cellAdmixCore:::celladmix_collect_cells
celladmix_collect_counts <- cellAdmixCore:::celladmix_collect_counts
celladmix_collect_region <- cellAdmixCore:::celladmix_collect_region
celladmix_collect_training_molecules <- cellAdmixCore:::celladmix_collect_training_molecules
celladmix_collect_transcripts <- cellAdmixCore:::celladmix_collect_transcripts
celladmix_correct <- cellAdmixCore:::celladmix_correct
celladmix_fit <- cellAdmixCore:::celladmix_fit
celladmix_prepare_tabular <- cellAdmixCore:::celladmix_prepare_tabular
celladmix_prepare_xenium <- cellAdmixCore:::celladmix_prepare_xenium
celladmix_read_run <- cellAdmixCore:::celladmix_read_run
celladmix_report_data <- cellAdmixCore:::celladmix_report_data
celladmix_score_bridge <- cellAdmixCore:::celladmix_score_bridge
celladmix_score_coherence <- cellAdmixCore:::celladmix_score_coherence
celladmix_select_factor_by_markers <- cellAdmixCore:::celladmix_select_factor_by_markers
celladmix_store <- cellAdmixCore:::celladmix_store

test_that("core version and simulation helpers work", {
  expect_match(celladmix_core_version(), "^0\\.0\\.1$")
  expect_setequal(
    getNamespaceExports("cellAdmixCore"),
    c("CellAdmixCorrection", "CellAdmixDataset", "CellAdmixFactorSourceScore",
      "CellAdmixFit",
      "CellAdmixScore", "cellAdmix", "celladmix_core_version",
      "celladmix_correction_impact", "celladmix_de",
      "celladmix_de_by_group", "celladmix_domain_composition",
      "celladmix_domain_diagnostics", "celladmix_discover_cell_boundaries",
      "celladmix_discover_stain_image", "celladmix_add_corrected_assay",
      "celladmix_add_factors_to_seurat", "celladmix_add_to_seurat",
      "celladmix_plot_de_shift", "celladmix_plot_cell_example",
      "celladmix_plot_cluster_marker_dotplot",
      "celladmix_plot_cluster_marker_heatmap",
      "celladmix_plot_cluster_marker_umap",
      "celladmix_plot_domain_composition", "celladmix_plot_domain_umap",
      "celladmix_plot_domains_spatial",
      "celladmix_plot_marker_expression",
      "celladmix_plot_molecule_removal", "celladmix_plot_score_pairs",
      "celladmix_plot_score_agreement", "celladmix_plot_spatial",
      "celladmix_plot_volcano", "celladmix_prepare_cell_example",
      "celladmix_read_cell_boundaries", "celladmix_schema",
      "celladmix_score_cluster_markers", "celladmix_select_example_cells",
      "celladmix_top_cluster_markers")
  )

  sim <- celladmix_simulate_nsclc(
    transcripts_per_cell = 30L,
    admixture_per_target_cell = 6L,
    seed = 11L
  )

  expect_true(is.data.frame(sim$transcripts))
  expect_true(is.data.frame(sim$cells))
  expect_true(nrow(sim$transcripts) > 0)
  expect_true(nrow(sim$cells) == 4)
})

test_that("diagnostic helpers work on sparse count matrices", {
  testthat::skip_if_not_installed("ggplot2")
  counts <- Matrix::Matrix(c(
    10, 8, 1, 1,
    0, 0, 5, 6,
    2, 2, 9, 8
  ), nrow = 3, byrow = TRUE, sparse = TRUE)
  rownames(counts) <- c("g1", "g2", "g3")
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  groups <- c(c1 = "a", c2 = "a", c3 = "b", c4 = "b")

  de <- celladmix_de(counts, groups, contrast = c("a", "b"))
  expect_setequal(c("gene", "logFC", "p_adj"), intersect(c("gene", "logFC", "p_adj"), names(de)))
  expect_equal(nrow(de), nrow(counts))

  expect_s3_class(celladmix_plot_volcano(de, markers = "g1", label = "g1"), "ggplot")
  expect_s3_class(celladmix_plot_volcano(de,
    markers = list(source_a = "g1", source_b = "g2")), "ggplot")
  expect_s3_class(celladmix_plot_marker_expression(counts, counts,
    cells = colnames(counts), markers = list(test = c("g1", "g2"))), "ggplot")
  expect_s3_class(celladmix_plot_de_shift(de, de, markers = "g1", label = "g1"), "ggplot")
  removal <- data.frame(
    cell_type = c("a", "a", "b"),
    n_molecules_before = c(10, 12, 9),
    n_molecules_after = c(8, 11, 7),
    n_removed = c(2, 1, 2),
    fraction_removed = c(0.2, 1 / 12, 2 / 9)
  )
  expect_s3_class(celladmix_plot_molecule_removal(removal), "ggplot")
  score_pairs <- data.frame(
    target_cell_type = c("A", "A", "B"),
    source_cell_type = c("B", "C", "A"),
    factor = c(1, 2, 1),
    mean_score = c(0.2, -0.1, 0.5),
    p_value = c(0.01, 0.2, 0.03),
    neg_log10_p = -log10(c(0.01, 0.2, 0.03)),
    used_in_summary = TRUE,
    stringsAsFactors = FALSE
  )
  expect_s3_class(celladmix_plot_score_pairs(score_pairs, method_label = "Smoke",
    show_top_pairs = TRUE), "ggplot")
  bridge_pairs <- score_pairs
  bridge_pairs$score <- bridge_pairs$mean_score
  bridge_pairs$mean_score <- NULL
  expect_s3_class(celladmix_plot_score_pairs(bridge_pairs, method_label = "Bridge",
    show_top_pairs = TRUE), "ggplot")
  testthat::skip_if_not_installed("cowplot")
  expect_s3_class(celladmix_plot_score_pairs(score_pairs, method_label = "Smoke"),
    "ggplot")
  score_summary <- data.frame(
    target_cell_type = c("A", "A", "B", "B"),
    source_cell_type = c("A", "B", "A", "B"),
    factor = c(1, 1, 1, 1),
    mean_score = c(0.1, 0.4, -0.2, 0.05),
    p_value = c(0.5, 0.01, 0.05, 0.8),
    neg_log10_p = -log10(c(0.5, 0.01, 0.05, 0.8)),
    stringsAsFactors = FALSE
  )
  expect_s3_class(celladmix_plot_score_pairs(score_summary,
    method_label = "Smoke", factor = "F1"), "ggplot")
  expect_s3_class(celladmix_plot_score_pairs(score_summary,
    method_label = "Smoke", factor = c("F1", "F2")), "ggplot")
  ann_a <- list(
    score_matrix = matrix(c(2, 0.5, 1, 3), nrow = 2,
      dimnames = list(c("A", "B"), c("F1", "F2"))),
    source_calls = list(f_1 = "B", f_2 = "A")
  )
  ann_b <- list(
    score_matrix = matrix(c(1.5, 0.25, 2, 2.5), nrow = 2,
      dimnames = list(c("A", "B"), c("F1", "F2"))),
    source_calls = list(f_1 = "B", f_2 = "B")
  )
  expect_s3_class(celladmix_plot_score_agreement(list(MethodA = ann_a, MethodB = ann_b)),
    "ggplot")
})

test_that("generic spatial plot handles labels and continuous values", {
  testthat::skip_if_not_installed("ggplot2")
  cell_df <- data.frame(
    cell_id = paste0("c", 1:4),
    x = c(0, 1, 0, 1),
    y = c(0, 0, 1, 1),
    cell_type = c("A", "A", "B", "B"),
    dominant_factor = c(1L, 2L, 1L, 2L),
    factor_1_fraction = c(0.1, 0.4, 0.7, 0.9),
    stringsAsFactors = FALSE
  )
  expect_s3_class(celladmix_plot_spatial(cell_df, color_by = "cell_type"), "ggplot")
  expect_s3_class(celladmix_plot_spatial(cell_df, color_by = "dominant_factor"), "ggplot")
  expect_s3_class(celladmix_plot_spatial(cell_df, color_by = "factor_1_fraction",
    reverse_y = TRUE), "ggplot")
  expect_error(celladmix_plot_spatial(cell_df, color_by = "missing"), "missing required columns")
})

test_that("internal weighted_nmf_matrix helper returns consistent matrix shapes", {
  x <- matrix(c(
    4, 1, 0,
    3, 1, 0,
    0, 1, 4,
    0, 1, 3
  ), nrow = 4, byrow = TRUE)
  rownames(x) <- paste0("row_", seq_len(nrow(x)))
  colnames(x) <- paste0("gene_", seq_len(ncol(x)))
  beta <- 1 / colSums(x)

  fit <- celladmix_weighted_nmf_matrix(
    x,
    column_weights = beta,
    rank = 2L,
    max_iterations = 50L,
    n_runs = 3L,
    num_threads = 2L,
    seed = 7L
  )

  expect_equal(dim(fit$w), c(4L, 2L))
  expect_equal(dim(fit$h), c(2L, 3L))
  expect_equal(rownames(fit$w), rownames(x))
  expect_equal(colnames(fit$h), colnames(x))
  expect_true(length(fit$losses) >= 1L)
  expect_true(is.finite(fit$final_loss))
  expect_length(fit$candidate_final_losses, 3L)
  expect_true(fit$selected_run >= 1)

  fit_cluster <- celladmix_weighted_nmf_matrix(
    x,
    column_weights = beta,
    rank = 2L,
    max_iterations = 50L,
    init_mode = "cluster",
    row_groups = c(1L, 1L, 2L, 2L),
    n_runs = 1L,
    num_threads = 1L,
    seed = 7L
  )
  expect_equal(dim(fit_cluster$w), c(4L, 2L))
  expect_equal(dim(fit_cluster$h), c(2L, 3L))

  fit_variant <- celladmix_weighted_nmf_matrix(
    x,
    column_weights = beta,
    rank = 2L,
    max_iterations = 50L,
    n_runs = 1L,
    num_threads = 1L,
    tolerance = 1e-8,
    update_epsilon = 1e-9,
    renormalize_each_iteration = FALSE,
    lee_style_epsilon = TRUE,
    random_init = "rnmf",
    seed = 7L
  )
  expect_equal(dim(fit_variant$w), c(4L, 2L))
  expect_equal(dim(fit_variant$h), c(2L, 3L))
})

test_that("internal reference NCV builder matches reference neighborhood cardinality convention", {
  df <- data.frame(
    mol_id = paste0("m", 1:6),
    gene = c("g1", "g2", "g1", "g2", "g3", "g3"),
    cell = c("c1", "c1", "c1", "c2", "c2", "c2"),
    x = c(0, 1, 2, 0, 1, 2),
    y = c(0, 0, 0, 1, 1, 1),
    z = c(0, 0, 0, 0, 0, 0)
  )

  x <- celladmix_build_reference_ncv_matrix(df, k = 2L)
  expect_equal(dim(x$x), c(6L, 3L))
  expect_true(all(rowSums(x$x) == 3))
  expect_equal(colnames(x$x), c("g1", "g2", "g3"))
})

test_that("internal pipeline NCV builder matches current pipeline cardinality convention", {
  df <- data.frame(
    mol_id = paste0("m", 1:6),
    gene = c("g1", "g2", "g1", "g2", "g3", "g3"),
    cell = c("c1", "c1", "c1", "c2", "c2", "c2"),
    x = c(0, 1, 2, 0, 1, 2),
    y = c(0, 0, 0, 1, 1, 1),
    z = c(0, 0, 0, 0, 0, 0)
  )

  x <- celladmix_build_pipeline_ncv_matrix(df, k = 2L)
  expect_equal(dim(x$x), c(6L, 3L))
  expect_true(all(rowSums(x$x) == 3))
})

write_mock_xenium_bundle <- function(td, sim, include_cell_type = TRUE) {
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      "{",
      '  "run_name": "Mock Xenium",',
      '  "pixel_size": 0.2125,',
      '  "z_step_size": 3.0',
      "}"
    ),
    file.path(td, "experiment.xenium")
  )

  tx <- data.frame(
    transcript_id = sim$transcripts$transcript_id,
    cell_id = sim$transcripts$cell,
    overlaps_nucleus = seq_len(nrow(sim$transcripts)) %% 3L == 0L,
    feature_name = sim$transcripts$gene,
    x_location = sim$transcripts$x,
    y_location = sim$transcripts$y,
    z_location = sim$transcripts$z,
    qv = rep(30, nrow(sim$transcripts))
  )
  tx$nucleus_distance <- as.numeric(seq_len(nrow(tx)) %% 11L)
  utils::write.csv(tx, gzfile(file.path(td, "transcripts.csv.gz")), row.names = FALSE)

  cells <- data.frame(
    cell_id = sim$cells$cell_id,
    x_centroid = sim$cells$x,
    y_centroid = sim$cells$y,
    z_centroid = sim$cells$z
  )
  if (include_cell_type) {
    cells$cell_type <- sim$cells$cell_type
  }
  utils::write.csv(cells, gzfile(file.path(td, "cells.csv.gz")), row.names = FALSE)
}

write_mock_tabular_bundle <- function(path, sim, include_cell_type = TRUE) {
  tx <- data.frame(
    x = sim$transcripts$x,
    y = sim$transcripts$y,
    z = sim$transcripts$z,
    gene = sim$transcripts$gene,
    cell_id = sim$transcripts$cell
  )
  if (include_cell_type) {
    tx$cell_type <- sim$transcripts$cell_type
  }
  utils::write.csv(tx, path, row.names = FALSE)
}

test_that("R6 front door supports a compact tabular fit workflow", {
  csv_path <- tempfile("tabular_r6_input_", fileext = ".csv")
  output_dir <- tempfile("celladmix_r6_output_")
  sim <- celladmix_simulate_nsclc(
    transcripts_per_cell = 24L,
    admixture_per_target_cell = 4L,
    seed = 19L
  )
  write_mock_tabular_bundle(csv_path, sim, include_cell_type = TRUE)

  annotation <- stats::setNames(sim$cells$cell_type, sim$cells$cell_id)
  ds <- cellAdmix(
    csv_path,
    output_dir = output_dir,
    schema = celladmix_schema(x = "x", y = "y", z = "z", gene = "gene",
      cell = "cell_id", cell_type = "cell_type"),
    annotation = annotation,
    annotation_name = "manual"
  )
  expect_true(inherits(ds, "CellAdmixDataset"))
  expect_equal(nrow(ds$annotations()), 1L)
  expect_true(file.exists(ds$prep$paths$project_json))

  fit <- ds$fit(
    rank = 2L,
    run_id = "r6_smoke",
    overwrite = TRUE,
    ncv_k = 6L,
    graph_k = 4L,
    nmf_iterations = 30L,
    nmf_n_runs = 2L,
    nmf_train_max_rows = 24L,
    report_ncv_umap = TRUE,
    num_threads = 2L,
    seed = 19L
  )
  expect_true(inherits(fit, "CellAdmixFit"))
  expect_true(file.exists(file.path(fit$run_dir, "run.json")))
  expect_true(is.data.frame(fit$cell_factors()))
  expect_true(is.data.frame(fit$training_molecules()))
  expect_equal(nrow(fit$stability()), 2L)
})

test_that("prepare_xenium does cheap preflight checks and reports missing cell_type metadata", {
  td <- tempfile("xenium_prep_mock_")
  sim <- celladmix_simulate_nsclc(
    transcripts_per_cell = 20L,
    admixture_per_target_cell = 4L,
    seed = 13L
  )
  write_mock_xenium_bundle(td, sim, include_cell_type = FALSE)

  prep <- celladmix_prepare_xenium(td)

  expect_s3_class(prep, "celladmix_prep")
  expect_true(file.exists(prep$paths$project_json))
  expect_identical(prep$source$type, "xenium")
  expect_true(any(grepl("cell_type", prep$warnings, fixed = TRUE)))
  expect_true("feature_name" %in% prep$probe$transcripts$columns)
})

test_that("input store modes expose consistent capabilities", {
  td <- tempfile("xenium_store_modes_mock_")
  sim <- celladmix_simulate_nsclc(
    transcripts_per_cell = 20L,
    admixture_per_target_cell = 4L,
    seed = 17L
  )
  write_mock_xenium_bundle(td, sim, include_cell_type = TRUE)

  prep <- celladmix_prepare_xenium(td)
  counts_store <- celladmix_store(prep, store_mode = "counts")
  expect_s3_class(counts_store, "celladmix_store")
  expect_identical(counts_store$store_mode, "counts")
  expect_true(isTRUE(counts_store$capabilities$has_cell_gene_counts))
  expect_false(isTRUE(counts_store$capabilities$has_molecule_rows))

  full_store <- celladmix_store(prep, store_mode = "full")
  expect_s3_class(full_store, "celladmix_store")
  expect_identical(full_store$store_mode, "full")
  expect_true(isTRUE(full_store$capabilities$has_cell_gene_counts))
  expect_true(isTRUE(full_store$capabilities$has_molecule_rows))
  expect_true(isTRUE(full_store$capabilities$has_cell_offsets))

  reread <- cellAdmixCore:::.celladmix_read_input_store(full_store$path)
  expect_identical(reread$store_mode, "full")
  expect_true(isTRUE(reread$capabilities$has_molecule_rows))
})

test_that("cluster_cells writes persisted annotations and fit can use cluster labels", {
  td <- tempfile("xenium_cluster_mock_")
  sim <- celladmix_simulate_nsclc(
    transcripts_per_cell = 30L,
    admixture_per_target_cell = 6L,
    seed = 23L
  )
  write_mock_xenium_bundle(td, sim, include_cell_type = FALSE)

  prep <- celladmix_prepare_xenium(td)
  clust <- celladmix_cluster_cells(
    prep,
    cells_max = 2L,
    n_variable_genes = 20L,
    pca_dims = 5L,
    graph_k = 2L,
    resolution = 1,
    umap_neighbors = 2L,
    umap_epochs = 20L,
    seed = 23L
  )

  expect_s3_class(clust, "celladmix_clusters")
  expect_true(file.exists(clust$paths$clusters_path))
  expect_true(file.exists(clust$paths$embedding_path))
  expect_equal(nrow(clust$clusters), 2L)
  expect_true(all(c("umap_1", "umap_2") %in% colnames(clust$embedding)))

  run <- celladmix_fit(
    prep,
    training_labels = clust,
    out_dir = tempfile("celladmix_run_clustered_"),
    ncv_k = 6L,
    rank = 2L,
    graph_k = 4L,
    nmf_iterations = 50L,
    nmf_init = "cluster",
    nmf_n_runs = 3L,
    nmf_train_max_rows = 20L,
    seed = 23L
  )

  expect_s3_class(run, "celladmix_run")
  expect_true(file.exists(file.path(run$path, "run.json")))
  expect_true(run$n_training_rows <= sum(clust$clusters$transcript_count))
  expect_equal(run$pipeline_options$nmf_n_runs, 3L)
  expect_identical(run$pipeline_options$nmf_init, "cluster")
  expect_length(run$nmf_candidate_final_objectives, 3L)
  expect_length(run$nmf_candidate_best_match_correlations, 3L)
  expect_true(is.finite(run$nmf_candidate_final_objective_sd))
  expect_true(is.finite(run$nmf_candidate_best_match_correlation_mean))

  run_auto <- celladmix_fit(
    prep,
    training_labels = "cluster",
    out_dir = tempfile("celladmix_run_clustered_auto_"),
    ncv_k = 6L,
    rank = 2L,
    graph_k = 4L,
    nmf_iterations = 50L,
    nmf_train_max_rows = 20L,
    report_ncv_umap = TRUE,
    seed = 23L
  )
  expect_s3_class(run_auto, "celladmix_run")
  expect_true(run_auto$n_training_rows <= sum(clust$clusters$transcript_count))
  expect_true(file.exists(run_auto$paths$training_molecule_umap_parquet))
  expect_equal(nrow(celladmix_collect_training_molecules(run_auto)), run_auto$n_training_rows)

  report_data <- celladmix_report_data(run, clust = clust, what = c("cell_umap", "ncv_umap"))
  expect_true(is.data.frame(report_data$cell_umap))
  expect_true(is.data.frame(report_data$ncv_umap))
  expect_true(all(c("umap_1", "umap_2", "dominant_factor") %in% colnames(report_data$cell_umap)))
  expect_true(all(c("umap_x", "umap_y", "factor_label") %in% colnames(report_data$ncv_umap)))
})

test_that("analysis_crop fit stays on store-backed cell-block path", {
  td <- tempfile("xenium_crop_store_mock_")
  sim <- celladmix_simulate_nsclc(
    transcripts_per_cell = 32L,
    admixture_per_target_cell = 6L,
    seed = 29L
  )
  write_mock_xenium_bundle(td, sim, include_cell_type = TRUE)

  prep <- celladmix_prepare_xenium(td, analysis_crops = data.frame(
    crop_id = c("lower", "upper"),
    xmin = c(0, 0),
    xmax = c(100, 100),
    ymin = c(0, 52),
    ymax = c(52, 110)
  ))

  run <- celladmix_fit(
    prep,
    training_labels = "none",
    out_dir = tempfile("celladmix_run_crop_store_"),
    analysis_crop = "lower",
    ncv_k = 6L,
    rank = 2L,
    graph_k = 4L,
    nmf_iterations = 30L,
    nmf_train_max_rows = 20L,
    seed = 29L
  )

  expect_s3_class(run, "celladmix_run")
  expect_identical(run$analysis_crop, "lower")
  expect_true(run$n_transcripts < nrow(sim$transcripts))

  cells <- celladmix_collect_cells(run)
  tx <- celladmix_collect_transcripts(run)
  expect_equal(nrow(cells), 2L)
  expect_true(all(tx$analysis_crop == "lower"))
})

test_that("prepare_tabular supports explicit cell-id input through clustering and fit", {
  csv_path <- tempfile("tabular_input_", fileext = ".csv")
  sim <- celladmix_simulate_nsclc(
    transcripts_per_cell = 30L,
    admixture_per_target_cell = 6L,
    seed = 31L
  )
  write_mock_tabular_bundle(csv_path, sim, include_cell_type = TRUE)

  prep <- celladmix_prepare_tabular(
    csv_path,
    cell_id_col = "cell_id",
    cell_type_col = "cell_type"
  )
  expect_s3_class(prep, "celladmix_prep")
  expect_identical(prep$source$type, "tabular")
  expect_true(file.exists(prep$paths$project_json))
  expect_true("gene" %in% prep$probe$molecules$columns)

  counts_store <- celladmix_store(prep, store_mode = "counts", force = TRUE)
  expect_identical(counts_store$store_mode, "counts")
  expect_true(isTRUE(counts_store$capabilities$has_cell_gene_counts))
  expect_false(isTRUE(counts_store$capabilities$has_molecule_rows))

  full_store <- celladmix_store(prep, store_mode = "full")
  expect_identical(full_store$store_mode, "full")
  expect_true(isTRUE(full_store$capabilities$has_cell_gene_counts))
  expect_true(isTRUE(full_store$capabilities$has_molecule_rows))
  expect_true(isTRUE(full_store$capabilities$has_cell_offsets))

  clust <- celladmix_cluster_cells(
    prep,
    cells_max = 4L,
    n_variable_genes = 20L,
    pca_dims = 5L,
    graph_k = 2L,
    resolution = 1,
    umap_neighbors = 2L,
    umap_epochs = 20L,
    seed = 31L
  )
  expect_s3_class(clust, "celladmix_clusters")
  expect_true(file.exists(clust$paths$clusters_path))

  run <- celladmix_fit(
    prep,
    training_labels = clust,
    out_dir = tempfile("celladmix_tabular_run_"),
    ncv_k = 6L,
    rank = 2L,
    graph_k = 4L,
    nmf_iterations = 50L,
    nmf_train_max_rows = 40L,
    seed = 31L
  )
  expect_s3_class(run, "celladmix_run")
  expect_identical(run$source$type, "tabular")
  expect_true(file.exists(file.path(run$path, "run.json")))
})

test_that("tabular store can attach cell types from a cell metadata sidecar", {
  csv_path <- tempfile("tabular_input_", fileext = ".csv")
  meta_path <- tempfile("tabular_cell_metadata_", fileext = ".csv")
  sim <- celladmix_simulate_nsclc(
    transcripts_per_cell = 30L,
    admixture_per_target_cell = 6L,
    seed = 41L
  )
  write_mock_tabular_bundle(csv_path, sim, include_cell_type = FALSE)
  label_by_cell <- stats::setNames(sim$cells$cell_type, sim$cells$cell_id)
  cell_meta <- unique(data.frame(
    cell_id = sim$transcripts$cell,
    label = label_by_cell[sim$transcripts$cell],
    stringsAsFactors = FALSE
  ))
  utils::write.csv(cell_meta, meta_path, row.names = FALSE)

  prep <- celladmix_prepare_tabular(
    csv_path,
    cell_id_col = "cell_id",
    cell_metadata_path = meta_path,
    cell_metadata_cell_id_col = "cell_id",
    cell_metadata_cell_type_col = "label"
  )
  expect_s3_class(prep, "celladmix_prep")
  expect_false(any(grepl("No cell_type_col", prep$warnings, fixed = TRUE)))

  full_store <- celladmix_store(prep, store_mode = "full", force = TRUE)
  expect_true(isTRUE(full_store$capabilities$has_cell_gene_counts))
  expect_true(isTRUE(full_store$capabilities$has_molecule_rows))
  expect_true(isTRUE(full_store$has_cell_type))

  run <- celladmix_fit(
    prep,
    training_labels = "cell_type",
    out_dir = tempfile("celladmix_tabular_sidecar_run_"),
    ncv_k = 6L,
    rank = 2L,
    graph_k = 4L,
    nmf_iterations = 40L,
    nmf_train_max_rows = 40L,
    seed = 41L
  )
  expect_s3_class(run, "celladmix_run")
  expect_true(isTRUE(run$has_cell_type))

  scoped_run <- celladmix_fit(
    prep,
    training_labels = "cell_type",
    scope = "fibroblast",
    out_dir = tempfile("celladmix_tabular_scoped_run_"),
    ncv_k = 6L,
    rank = 2L,
    graph_k = 4L,
    nmf_iterations = 40L,
    nmf_train_max_rows = 40L,
    report_ncv_umap = TRUE,
    seed = 41L
  )
  expect_s3_class(scoped_run, "celladmix_run")
  expect_equal(scoped_run$n_cells, run$n_cells)
  expect_identical(scoped_run$pipeline_options$training_scope_cell_types, "fibroblast")
  scoped_training <- celladmix_collect_training_molecules(scoped_run)
  expect_true(all(label_by_cell[unique(scoped_training$cell_id)] == "fibroblast"))
})

test_that("fit persists supported KL-NMF matrix variants", {
  csv_path <- tempfile("tabular_variant_input_", fileext = ".csv")
  sim <- celladmix_simulate_nsclc(
    transcripts_per_cell = 24L,
    admixture_per_target_cell = 4L,
    seed = 43L
  )
  write_mock_tabular_bundle(csv_path, sim, include_cell_type = TRUE)

  prep <- celladmix_prepare_tabular(
    csv_path,
    cell_id_col = "cell_id",
    cell_type_col = "cell_type"
  )

  sqrt_run <- celladmix_fit(
    prep,
    training_labels = "none",
    out_dir = tempfile("celladmix_sqrt_variant_run_"),
    ncv_k = 6L,
    rank = 2L,
    graph_k = 4L,
    nmf_iterations = 25L,
    nmf_train_max_rows = 24L,
    nmf_variant = "sqrt_kl",
    seed = 43L
  )
  expect_identical(sqrt_run$pipeline_options$nmf_variant, "sqrt_kl")
  expect_equal(length(sqrt_run$nmf_gene_weights), 0L)
  expect_gt(sqrt_run$nmf_transform_target_row_sum, 0)

  invsqrt_run <- celladmix_fit(
    prep,
    training_labels = "none",
    out_dir = tempfile("celladmix_invsqrt_variant_run_"),
    ncv_k = 6L,
    rank = 2L,
    graph_k = 4L,
    nmf_iterations = 25L,
    nmf_train_max_rows = 24L,
    nmf_variant = "invsqrt_kl",
    seed = 44L
  )
  expect_identical(invsqrt_run$pipeline_options$nmf_variant, "invsqrt_kl")
  expect_equal(length(invsqrt_run$nmf_gene_weights), length(invsqrt_run$genes))
  expect_gt(invsqrt_run$nmf_transform_target_row_sum, 0)

  reread <- celladmix_read_run(invsqrt_run$path)
  expect_identical(reread$pipeline_options$nmf_variant, "invsqrt_kl")
  expect_equal(reread$nmf_gene_weights, invsqrt_run$nmf_gene_weights)
})

test_that("automatic factor rank follows clustering size and caps large clusterings", {
  prep <- structure(
    list(
      project_dir = normalizePath(tempdir(), winslash = "/", mustWork = FALSE),
      paths = list(annotations_dir = tempdir())
    ),
    class = "celladmix_prep"
  )

  clust_small <- structure(
    list(
      output_dir = prep$project_dir,
      summary = list(n_clusters = 8L),
      paths = list(clusters_path = tempfile("clusters_small_", fileext = ".parquet"))
    ),
    class = "celladmix_clusters"
  )
  expect_equal(celladmix_resolve_rank(prep, clust_small, rank = NULL), 10L)
  expect_equal(celladmix_resolve_rank(prep, clust_small, rank = 7L), 7L)

  clust_large <- structure(
    list(
      output_dir = prep$project_dir,
      summary = list(n_clusters = 24L),
      paths = list(clusters_path = tempfile("clusters_large_", fileext = ".parquet"))
    ),
    class = "celladmix_clusters"
  )
  expect_warning(
    expect_equal(celladmix_resolve_rank(prep, clust_large, rank = NULL), 24L),
    "capping"
  )
})

test_that("persisted runs support read, collect, bridge scoring, and correction", {
  td <- tempfile("xenium_run_mock_")
  sim <- celladmix_simulate_nsclc(
    transcripts_per_cell = 40L,
    admixture_per_target_cell = 10L,
    seed = 31L
  )
  write_mock_xenium_bundle(td, sim, include_cell_type = TRUE)

  annotation <- data.frame(
    cell_id = as.character(sim$cells$cell_id),
    cell_type = as.character(sim$cells$cell_type),
    stringsAsFactors = FALSE
  )
  ds <- cellAdmix(td, output_dir = tempfile("celladmix_output_"),
    annotation = annotation, annotation_name = "cell_type",
    annotation_col = "cell_type", cell_id_col = "cell_id")
  fit <- ds$fit(
    annotation = "cell_type",
    run_id = "celladmix_run_public",
    overwrite = TRUE,
    ncv_k = 8L,
    rank = 2L,
    graph_k = 4L,
    nmf_iterations = 60L,
    nmf_train_max_rows = 40L,
    seed = 31L
  )
  run <- fit$run

  expect_s3_class(run, "celladmix_run")
  expect_true(file.exists(file.path(run$path, "run.json")))
  expect_true(isTRUE(run$has_overlaps_nucleus))
  expect_true(isTRUE(run$has_nucleus_distance))

  reread <- celladmix_read_run(run$path)
  expect_s3_class(reread, "celladmix_run")
  expect_equal(reread$n_transcripts, run$n_transcripts)
  expect_true(isTRUE(reread$has_overlaps_nucleus))
  expect_true(isTRUE(reread$has_nucleus_distance))

  tx_view <- celladmix_collect_transcripts(run, sample_n = 10L, seed = 1L)
  expect_true(nrow(tx_view) <= 10L)
  expect_true(all(c("gene", "cell", "x", "y", "z", "factor_label") %in% colnames(tx_view)))

  tx_region <- celladmix_collect_region(run, bbox = c(0, 1000, 0, 1000), sample_n = 10L, seed = 2L)
  expect_true(nrow(tx_region) <= 10L)

  cells <- celladmix_collect_cells(run)
  report_data <- celladmix_report_data(run, what = "ncv_umap")
  training_molecules <- report_data$ncv_umap
  counts <- celladmix_collect_counts(run)
  expect_true(is.data.frame(cells))
  expect_true(is.data.frame(training_molecules))
  expect_equal(nrow(training_molecules), run$n_training_rows)
  expect_true(all(c("obs_id", "cell_id", "factor_label", "umap_x", "umap_y", "component_1") %in%
    colnames(training_molecules)))
  expect_true(is.matrix(counts))
  expect_equal(ncol(counts), nrow(cells))

  factor_id <- celladmix_select_factor_by_markers(run, c("KRT19", "KRT8", "KRT17"))
  expect_true(factor_id %in% c(1L, 2L))

  bridge <- celladmix_score_bridge(
    run,
    factor = factor_id,
    candidate_k = 5L,
    crossing_k = 5L,
    min_type_pair_contacts = 1L,
    min_factor_molecules = 1L,
    min_pairs = 1L,
    null_iterations = 1L,
    null_max_iterations = 2L,
    compute_null = TRUE
  )
  expect_s3_class(bridge, "celladmix_bridge_result")
  expect_true(is.data.frame(bridge$scores))
  expect_true(is.data.frame(bridge$summary))
  expect_true(all(c("target_cell", "source_cell", "target_cell_type", "source_cell_type", "factor", "score") %in%
    colnames(bridge$scores)))
  expect_true(all(c("target_cell_type", "source_cell_type", "factor", "n_pairs", "p_value") %in%
    colnames(bridge$summary)))
  bridge_ann <- celladmix_bridge_annotation(bridge, p_thresh = 0.5)
  expect_true(is.matrix(bridge_ann$score_matrix))
  expect_true(is.list(bridge_ann$source_calls))

  coherence <- celladmix_score_coherence(
    run,
    factor = factor_id,
    stain_mode = "none",
    k_neighbors = 4L,
    min_factor_molecules = 1L,
    min_cells = 1L,
    null_iterations = 2L,
    compute_null = TRUE
  )
  expect_s3_class(coherence, "celladmix_coherence_result")
  expect_true(is.data.frame(coherence$scores))
  expect_true(is.data.frame(coherence$summary))
  expect_identical(coherence$null_matching$method, "label_permutation")
  expect_true(isTRUE(coherence$null_matching$nucleus_overlap))

  partial_annotation <- stats::setNames(annotation$cell_type, annotation$cell_id)
  partial_annotation <- partial_annotation[-1]
  bridge_partial <- celladmix_score_bridge(
    run,
    annotation = partial_annotation,
    factor = factor_id,
    candidate_k = 5L,
    crossing_k = 5L,
    min_type_pair_contacts = 1L,
    min_factor_molecules = 1L,
    min_pairs = 1L,
    compute_null = FALSE
  )
  expect_false("__unannotated__" %in% bridge_partial$cell_types)
  expect_false(any(bridge_partial$scores$target_cell_type == "__unannotated__"))
  expect_false(any(bridge_partial$scores$source_cell_type == "__unannotated__"))

  corrected <- celladmix_correct(
    run,
    rules = data.frame(factor = factor_id, target_cell_type = "fibroblast"),
    out_dir = tempfile("celladmix_corrected_run_")
  )
  expect_s3_class(corrected, "celladmix_correction_run")
  expect_true(file.exists(file.path(corrected$path, "run.json")))

  corrected_counts <- celladmix_collect_counts(corrected)
  expect_true(sum(corrected_counts) <= sum(counts))
})

test_that("example stain auto-discovery resolves quietly", {
  fake_run <- structure(list(source = list(type = "tabular", path = "")),
    class = "celladmix_run")
  expect_identical(cellAdmixCore:::.celladmix_discover_stain_images(fake_run), list())

  fit_stub <- list(stains = function(...) list())
  resolve <- cellAdmixCore:::.celladmix_resolve_example_stains
  expect_identical(resolve(fit_stub, "auto"), list())
  expect_identical(resolve(fit_stub, NULL), list())

  descriptor <- list(path = "x.tif", pixel_size = 0.5, stain = "membrane")
  expect_identical(resolve(fit_stub, descriptor), list(membrane = descriptor))
  named <- list(dapi = descriptor)
  expect_identical(resolve(fit_stub, named), named)

  fit_with_stains <- list(stains = function(...) list(dapi = descriptor))
  expect_identical(resolve(fit_with_stains, "auto"), list(dapi = descriptor))
})

test_that("explicit cell selection returns requested cells with relaxed filters", {
  pairs <- data.frame(
    target_cell = c("c1", "c1", "c2"),
    target_cell_type = c("T", "T", "B"),
    factor = c(1, 2, 1),
    mean_score = c(0.5, 0.2, -0.1),
    factor_count = c(2, 1, 1),
    used_in_summary = c(TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  cell_factors <- data.frame(
    cell_id = c("c1", "c2", "c3"),
    transcript_count = c(10, 20, 30),
    dominant_factor = c(3, 3, 3),
    stringsAsFactors = FALSE
  )
  stub_score <- structure(list(
    pairs = function() pairs,
    annotation = function(...) list(source_calls = list()),
    rules = function(...) data.frame(),
    fit = list(
      cell_factors = function() cell_factors,
      annotation_name = "manual",
      dataset = list(annotation = function(...) c(c1 = "T", c2 = "B", c3 = "NK"))
    )
  ), class = "CellAdmixScore")

  out <- celladmix_select_example_cells(stub_score, cells = c("c3", "c1"))
  expect_equal(as.character(out$target_cell), c("c3", "c1"))
  expect_equal(as.character(out$target_cell_type[[1]]), "NK")
  expect_true(is.na(out$top_admix_factor[[1]]))
  expect_equal(out$transcript_count[[1]], 30)
  expect_equal(as.integer(out$top_admix_factor[[2]]), 1L)

  neg <- celladmix_select_example_cells(stub_score, cells = "c2")
  expect_equal(nrow(neg), 1L)
  expect_equal(neg$top_admix_score[[1]], -0.1)

  auto <- celladmix_select_example_cells(stub_score,
    min_molecules = 0, min_factor_molecules = 0)
  expect_equal(as.character(auto$target_cell), "c1")
})

test_that("example cell-type resolution handles auto, NULL, and frames", {
  resolve <- cellAdmixCore:::.celladmix_resolve_example_cell_types
  fit_stub <- list(
    annotation_name = "manual",
    dataset = list(annotation = function(...) c(c1 = "T"))
  )
  expect_identical(resolve(fit_stub, "auto"), c(c1 = "T"))
  expect_null(resolve(fit_stub, NULL))
  frame <- data.frame(cell_id = "c1", cell_type = "T", stringsAsFactors = FALSE)
  expect_identical(resolve(fit_stub, frame), c(c1 = "T"))
})

test_that("native-factor check vetoes persistent factors and keeps true admixture", {
  set.seed(0)
  n <- 20L
  cells <- data.frame(
    cell_id = c(paste0("S", 1:n), paste0("TE", 1:n), paste0("TD", 1:n)),
    x = c(runif(n, 0, 10), runif(n, 5, 15), runif(n, 100, 110)),
    y = runif(3L * n, 0, 10),
    factor_1_fraction = c(rep(0.05, n), rep(0.30, 2L * n)),
    factor_2_fraction = c(rep(0, n), rep(0.40, n), rep(0, n)),
    stringsAsFactors = FALSE
  )
  annotation <- setNames(rep(c("Source", "Target"), c(n, 2L * n)), cells$cell_id)
  fit_stub <- list(
    cell_factors = function() cells,
    annotation_name = "manual",
    dataset = list(annotation = function(...) annotation)
  )
  rules <- data.frame(
    factor = c(1L, 2L, 1L),
    source_cell_type = c("Source", "Source", "Missing"),
    target_cell_type = "Target",
    p_value = 0.01,
    neg_log10_p = 2,
    rule_id = c("1_Target", "2_Target", "1_Target"),
    stringsAsFactors = FALSE
  )

  out <- cellAdmixCore:::.celladmix_apply_native_check(rules, fit_stub, neighbor_k = 10L)
  expect_equal(out$native_check[1:2], c("native_median", "pass"))
  expect_equal(out$keep, c(FALSE, TRUE, FALSE))
  expect_gt(out$native_distant_median[[1]], 0.1)
  expect_gt(out$native_exposure_gradient[[2]], 0.2)
  expect_equal(out$native_check[[3]], "source_not_in_annotation")

  empty <- cellAdmixCore:::.celladmix_apply_native_check(rules[0, ], fit_stub)
  expect_true("keep" %in% names(empty))
  expect_equal(nrow(empty), 0L)
})
