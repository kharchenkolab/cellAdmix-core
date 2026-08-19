celladmix_simulate_nsclc <- cellAdmixCore:::celladmix_simulate_nsclc
celladmix_correct <- cellAdmixCore:::celladmix_correct
celladmix_select_factor_by_markers <- cellAdmixCore:::celladmix_select_factor_by_markers

make_ensemble_fit <- function(seed = 7L, nmf_n_runs = 4L) {
  td <- tempfile("xenium_ens_mock_")
  sim <- celladmix_simulate_nsclc(
    transcripts_per_cell = 40L,
    admixture_per_target_cell = 10L,
    seed = seed
  )
  write_mock_xenium_bundle(td, sim, include_cell_type = TRUE)
  annotation <- data.frame(
    cell_id = as.character(sim$cells$cell_id),
    cell_type = as.character(sim$cells$cell_type),
    stringsAsFactors = FALSE
  )
  ds <- cellAdmix(td, output_dir = tempfile("celladmix_ens_out_"),
    annotation = annotation, annotation_name = "cell_type",
    annotation_col = "cell_type", cell_id_col = "cell_id")
  fit <- ds$fit(
    annotation = "cell_type",
    run_id = "ens_run",
    overwrite = TRUE,
    ncv_k = 8L,
    rank = 2L,
    graph_k = 4L,
    nmf_iterations = 60L,
    nmf_train_max_rows = 40L,
    nmf_n_runs = nmf_n_runs,
    seed = seed
  )
  fit
}

test_that("multirun fits persist a restart member pool with per-member labels", {
  fit <- make_ensemble_fit(nmf_n_runs = 4L)
  run <- fit$run

  expect_true(file.exists(file.path(run$path, "ensemble_h.parquet")))
  n_members <- cellAdmixCore:::.celladmix_ensemble_prepare(run$path, 2L)
  expect_identical(n_members, 4L)
  expect_true(all(file.exists(
    file.path(run$path, sprintf("ensemble_labels_m%d.parquet", 0:3)))))

  # The selected member reproduces the run's own labeling: its cell factor
  # fractions must match the fit's cell factors exactly.
  selected <- as.integer(run$nmf_diagnostics$selected_run) - 1L
  member_fr <- cellAdmixCore:::.celladmix_ensemble_member_fractions(run$path, selected)
  cells <- fit$cell_factors()
  idx <- match(as.character(cells$cell_id), as.character(member_fr$cell_id))
  expect_false(anyNA(idx))
  for (k in seq_len(ncol(member_fr$fractions))) {
    expect_lt(max(abs(member_fr$fractions[idx, k] -
      cells[[paste0("factor_", k, "_fraction")]])), 1e-9)
  }

  # Every member yields valid fractions.
  for (m in 0:3) {
    fr <- cellAdmixCore:::.celladmix_ensemble_member_fractions(run$path, m)
    expect_true(all(fr$fractions >= 0 & fr$fractions <= 1))
  }
})

test_that("vote-threshold correction is monotone and matches the primary path", {
  fit <- make_ensemble_fit(nmf_n_runs = 4L)
  run <- fit$run
  invisible(cellAdmixCore:::.celladmix_ensemble_prepare(run$path, 2L))
  selected <- as.integer(run$nmf_diagnostics$selected_run) - 1L
  factor_id <- celladmix_select_factor_by_markers(run, c("KRT19", "KRT8", "KRT17"))

  rule <- data.frame(factor = factor_id, target_cell_type = "fibroblast",
    stringsAsFactors = FALSE)
  rules4 <- rule[rep(1L, 4L), , drop = FALSE]

  corr_primary <- celladmix_correct(run, rule,
    out_dir = tempfile("ens_corr_primary_"))
  corr_selected <- celladmix_correct(run, rule,
    out_dir = tempfile("ens_corr_selected_"),
    rule_member = selected, min_votes = 1L)
  # The selected member's labeling is the run's own labeling.
  expect_identical(corr_selected$n_removed, corr_primary$n_removed)

  corr_union <- celladmix_correct(run, rules4,
    out_dir = tempfile("ens_corr_union_"),
    rule_member = 0:3, min_votes = 1L)
  corr_majority <- celladmix_correct(run, rules4,
    out_dir = tempfile("ens_corr_majority_"),
    rule_member = 0:3, min_votes = 2L)
  corr_unanimous <- celladmix_correct(run, rules4,
    out_dir = tempfile("ens_corr_all_"),
    rule_member = 0:3, min_votes = 4L)

  expect_true(corr_union$n_removed >= corr_majority$n_removed)
  expect_true(corr_majority$n_removed >= corr_unanimous$n_removed)
  # The union removes at least what any single member removes.
  expect_true(corr_union$n_removed >= corr_primary$n_removed)

  # Vote histogram accounts for every molecule the union removes.
  expect_identical(sum(corr_union$vote_histogram), as.integer(corr_union$n_removed))
  # A threshold above the member count removes nothing.
  corr_none <- celladmix_correct(run, rules4,
    out_dir = tempfile("ens_corr_none_"),
    rule_member = 0:3, min_votes = 5L)
  expect_identical(corr_none$n_removed, 0)
})

test_that("fit$correct defaults to the ensemble and falls back cleanly", {
  fit <- make_ensemble_fit(nmf_n_runs = 4L)
  score <- fit$score_bridge(
    candidate_k = 5L,
    crossing_k = 5L,
    min_type_pair_contacts = 1L,
    min_factor_molecules = 1L,
    min_pairs = 1L,
    null_iterations = 1L,
    null_max_iterations = 2L,
    compute_null = TRUE
  )

  expect_message(
    corr <- fit$correct(score, name = "ens_default", p_thresh = 0.9),
    "Ensemble correction over 4 members")
  expect_s3_class(corr, "CellAdmixCorrection")
  info <- corr$ensemble()
  expect_identical(info$members, 4L)
  expect_identical(info$min_votes, 2L)
  if (!is.null(corr$rules) && nrow(corr$rules)) {
    expect_true("support" %in% names(corr$rules))
    expect_true(all(corr$rules$support >= 0 & corr$rules$support <= 1))
  }

  corr_single <- fit$correct(score, name = "ens_single", p_thresh = 0.9,
    ensemble = 1)
  expect_identical(corr_single$ensemble()$members, 1L)
  # Re-correcting at a different vote threshold reuses cached member rules.
  corr_union <- fit$correct(score, name = "ens_union", p_thresh = 0.9,
    vote = 1e-6)
  expect_true(corr_union$run$n_removed >= corr$run$n_removed)
})

test_that("member rules are cached on disk and reused across score objects", {
  fit <- make_ensemble_fit(seed = 13L, nmf_n_runs = 4L)
  score_args <- list(candidate_k = 5L, crossing_k = 5L,
    min_type_pair_contacts = 1L, min_factor_molecules = 1L, min_pairs = 1L,
    null_iterations = 1L, null_max_iterations = 2L, compute_null = TRUE)
  score <- do.call(fit$score_bridge, score_args)
  corr1 <- fit$correct(score, name = "cache_a", p_thresh = 0.9)

  cache_files <- list.files(fit$run$paths$scores_dir,
    pattern = "^ensemble_rules_bridge_m", full.names = TRUE)
  expect_length(cache_files, 3L)
  mtimes <- file.info(cache_files)$mtime
  keys <- vapply(cache_files, function(f) readRDS(f)$key, character(1))

  # A fresh score object (a new session, effectively) reuses the disk cache
  # instead of re-scoring the members.
  Sys.sleep(1.1)
  score2 <- do.call(fit$score_bridge, score_args)
  corr2 <- fit$correct(score2, name = "cache_b", p_thresh = 0.9)
  expect_identical(corr2$run$n_removed, corr1$run$n_removed)
  expect_equal(file.info(cache_files)$mtime, mtimes)

  # A changed rule threshold invalidates and rewrites the cache.
  corr3 <- fit$correct(score2, name = "cache_c", p_thresh = 0.5)
  keys_after <- vapply(cache_files, function(f) readRDS(f)$key, character(1))
  expect_false(any(keys_after == keys))
})

test_that("member labels recompute deterministically after cache removal", {
  fit <- make_ensemble_fit(seed = 17L, nmf_n_runs = 3L)
  n_members <- cellAdmixCore:::.celladmix_ensemble_prepare(fit$run$path, 2L)
  selected <- as.integer(fit$run$nmf_diagnostics$selected_run) - 1L
  member <- setdiff(seq_len(n_members) - 1L, selected)[[1]]
  before <- cellAdmixCore:::.celladmix_ensemble_member_fractions(fit$run$path, member)
  file.remove(file.path(fit$run$path,
    sprintf("ensemble_labels_m%d.parquet", member)))
  cellAdmixCore:::.celladmix_ensemble_prepare(fit$run$path, 2L)
  after <- cellAdmixCore:::.celladmix_ensemble_member_fractions(fit$run$path, member)
  expect_equal(before$fractions, after$fractions)
})

test_that("erasing a cell type's whole content triggers the over-removal warning", {
  fit <- make_ensemble_fit(seed = 19L, nmf_n_runs = 1L)
  score <- fit$score_bridge(
    candidate_k = 5L, crossing_k = 5L, min_type_pair_contacts = 1L,
    min_factor_molecules = 1L, min_pairs = 1L, null_iterations = 1L,
    null_max_iterations = 2L, compute_null = TRUE)
  rules <- data.frame(factor = c(1L, 2L),
    target_cell_type = "fibroblast", source_cell_type = "malignant",
    stringsAsFactors = FALSE)
  old <- options(celladmix.overremoval_min_molecules = 10)
  on.exit(options(old), add = TRUE)
  expect_warning(
    fit$correct(score, rules = rules, ensemble = 1, name = "erase_all"),
    "erasing native expression")
})

test_that("single-restart fits fall back to the single-fit correction", {
  fit <- make_ensemble_fit(seed = 11L, nmf_n_runs = 1L)
  expect_false(file.exists(file.path(fit$run$path, "ensemble_h.parquet")))
  score <- fit$score_bridge(
    candidate_k = 5L,
    crossing_k = 5L,
    min_type_pair_contacts = 1L,
    min_factor_molecules = 1L,
    min_pairs = 1L,
    null_iterations = 1L,
    null_max_iterations = 2L,
    compute_null = TRUE
  )
  expect_message(
    corr <- fit$correct(score, name = "single_fallback", p_thresh = 0.9),
    "no ensemble member pool")
  expect_identical(corr$ensemble()$members, 1L)
})
