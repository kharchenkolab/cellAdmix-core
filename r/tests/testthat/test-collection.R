celladmix_simulate_nsclc <- cellAdmixCore:::celladmix_simulate_nsclc

make_sample_dataset <- function(seed, label) {
  td <- tempfile(paste0("xenium_collection_", label, "_"))
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
  cellAdmix(td, output_dir = tempfile(paste0("celladmix_collection_out_", label, "_")),
    annotation = annotation, annotation_name = "cell_type",
    annotation_col = "cell_type", cell_id_col = "cell_id", num_threads = 2L)
}

fit_args <- list(ncv_k = 8L, rank = 2L, graph_k = 4L, nmf_iterations = 60L,
  nmf_train_max_rows = 60L, nmf_n_runs = 3L, seed = 7L)

test_that("fixed-profile fits keep the given factor order and carry a member pool", {
  ds <- make_sample_dataset(31L, "fixed")
  base <- do.call(ds$fit, c(list(annotation = "cell_type", run_id = "base", overwrite = TRUE), fit_args))
  H <- cellAdmixCore:::.celladmix_fit_profiles_raw(base)
  # Reverse the factor order: the fixed-profile run must report the
  # reversed profiles in the given order, not re-sorted by importance.
  H_rev <- H[rev(seq_len(nrow(H))), , drop = FALSE]
  members <- list(H_rev * 1.1, H_rev * 0.9)
  fit <- ds$fit(annotation = "cell_type", run_id = "fixed_members", overwrite = TRUE,
    nmf_fixed_h = H_rev, nmf_fixed_h_members = members, nmf_n_runs = 1L,
    ncv_k = 8L, graph_k = 4L, seed = 7L)
  expect_true(isTRUE(fit$run$fixed_profiles))
  expect_equal(fit$run$fixed_profile_members, 3L)
  loadings <- fit$loadings()
  colnames(loadings) <- fit$run$genes
  expect_equal(unname(loadings[, colnames(H_rev)]), unname(H_rev), tolerance = 1e-8)
  expect_equal(cellAdmixCore:::.celladmix_ensemble_prepare(fit$run$path, 1L), 3L)
  pool <- cellAdmixCore:::.celladmix_ensemble_h(fit$run$path)
  expect_length(pool, 3L)
  expect_equal(unname(pool[[1]][, colnames(H_rev)]), unname(H_rev), tolerance = 1e-8)
  expect_equal(unname(pool[[2]][, colnames(H_rev)]), unname(H_rev * 1.1), tolerance = 1e-8)
  # Agreement of the members with the fixed profiles is reported in the
  # stability slots (scaled copies own exactly the same genes).
  expect_equal(fit$run$nmf_diagnostics$stability_comparison_runs, 2L)
  expect_true(all(fit$run$nmf_diagnostics$selected_factor_stability > 0.99))
})

test_that("factor matching and training loss are reachable from R", {
  ds <- make_sample_dataset(41L, "loss")
  fit <- do.call(ds$fit, c(list(annotation = "cell_type", run_id = "base", overwrite = TRUE), fit_args))
  H <- cellAdmixCore:::.celladmix_fit_profiles_raw(fit)
  match <- cellAdmixCore:::.celladmix_match_factors(H[2:1, ], H)
  expect_equal(match$match, c(2L, 1L))
  expect_true(all(match$correlation > 0.99))
  own <- cellAdmixCore:::.celladmix_training_loss(fit$run$path, H)
  expect_true(is.finite(own) && own > 0)
  # A profile set that ignores the gene structure fits the rows worse.
  flat <- matrix(1, nrow = 2, ncol = ncol(H), dimnames = dimnames(H))
  expect_gt(cellAdmixCore:::.celladmix_training_loss(fit$run$path, flat), own)
})

test_that("a consensus factorization pools two samples end to end", {
  datasets <- list(a = make_sample_dataset(51L, "a"), b = make_sample_dataset(52L, "b"))
  collection <- cellAdmixCollection(datasets)
  expect_equal(collection$samples()$sample, c("a", "b"))
  do.call(collection$fit, c(list(annotation = "cell_type", overwrite = TRUE), fit_args))
  expect_equal(names(collection$fits), c("a", "b"))

  consensus <- collection$fit_consensus()
  expect_s3_class(consensus, "CellAdmixConsensusFit")
  expect_equal(consensus$samples(), c("a", "b"))
  profiles <- consensus$loadings()
  expect_true(nrow(profiles) >= 2L && nrow(profiles) <= 4L)
  expect_equal(nrow(consensus$matching), 4L)
  expect_true("own_stability" %in% names(consensus$matching))
  expect_true(all(is.finite(consensus$matching$own_stability)))
  factors <- consensus$factors
  expect_equal(nrow(factors), nrow(profiles))
  expect_true(all(c("consensus_factor", "samples", "coverage", "restarts", "components",
    "agreement") %in% names(factors)))
  # one factorization contributes at most one component to a shared factor
  expect_true(all(factors$components == factors$restarts))
  expect_true(all(factors$restarts >= 2L))
  expect_true(all(factors$coverage > 0 & factors$coverage <= 1))
  recurrence <- consensus$recurrence()
  expect_equal(nrow(recurrence), nrow(profiles) * 2L)
  expect_true(all(recurrence$recurrence >= 0 & recurrence$recurrence <= 1))
  expect_true(all(c("sample", "restart", "factor", "cluster") %in% names(consensus$unreproduced)))
  # every pooled component is either in a shared factor or unreproduced
  expect_equal(sum(!is.na(consensus$components$consensus_factor)) + nrow(consensus$unreproduced),
    nrow(consensus$components))
  expect_equal(nrow(consensus$loss), 2L)
  expect_true(all(is.finite(consensus$loss$own_loss)))
  expect_true(all(is.finite(consensus$loss$consensus_loss)))
  for (sample in c("a", "b")) {
    fit <- consensus$fits[[sample]]
    expect_true(isTRUE(fit$run$fixed_profiles))
    expect_equal(as.integer(fit$run$n_factors), nrow(profiles))
    # consensus itself plus both samples' three aligned restarts, shared
    expect_equal(fit$run$fixed_profile_members, 7L)
  }
  stability <- consensus$stability()
  expect_equal(stability$sample, c("a", "b"))
  usage <- consensus$sample_usage()
  expect_true(all(c("sample", "cell_type", "factor", "share") %in% names(usage)))
  divergence <- consensus$usage_divergence()
  expect_true(all(divergence$divergence >= 0 & divergence$divergence <= 1))

  cells <- consensus$cells()
  expect_true(all(startsWith(cells$cell_id, paste0(cells$sample, "_"))))
  counts <- consensus$counts()
  expect_equal(ncol(counts), nrow(cells))
  expect_true(all(startsWith(colnames(counts), "a_") | startsWith(colnames(counts), "b_")))

  score <- consensus$score_bridge(candidate_pairs_per_type_pair = 40L,
    max_cells_per_type_pair = 40L, null_iterations = 1L)
  expect_s3_class(score, "CellAdmixCollectionScore")
  pooled <- score$summary()
  expect_true(all(c("n_pairs", "mean_score", "p_value", "n_samples") %in% names(pooled)))
  expect_true(all(pooled$n_samples >= 1 & pooled$n_samples <= 2))
  per_sample <- score$sample_summaries()
  expect_true("sample" %in% names(per_sample))
  expect_true(all(per_sample$sample %in% c("a", "b")))

  rules <- score$rules(p_thresh = 0.5)
  expect_true(all(c("factor", "source_cell_type", "target_cell_type",
    "sample_support", "keep", "native_check") %in% names(rules)))
  expect_true(all(rules$sample_support >= 0 & rules$sample_support <= 1))
  if (!nrow(rules)) {
    # The simulated fixture carries no bridge evidence; an explicit shared
    # rule still exercises the per-sample ensemble correction.
    rules <- data.frame(factor = 1L, source_cell_type = "malignant",
      target_cell_type = "fibroblast", stringsAsFactors = FALSE)
  }

  correction <- score$correct(rules = rules, vote = 0.3)
  expect_s3_class(correction, "CellAdmixCollectionCorrection")
  expect_equal(correction$samples(), c("a", "b"))
  expect_true("support" %in% names(correction$rules))
  expect_true(all(correction$rules$support >= 0 & correction$rules$support <= 1))
  diagnostics <- consensus$diagnostics()
  expect_equal(diagnostics$sample, c("a", "b"))
  expect_true(all(c("loss_gap", "restart_spread", "representation_flag", "shared_fraction",
    "shared_factors_flag", "within_sample_stability", "consensus_agreement",
    "max_usage_divergence", "divergent_cell_type") %in% names(diagnostics)))
  expect_true(all(diagnostics$shared_fraction >= 0 & diagnostics$shared_fraction <= 1))
  by_sample <- consensus$usage_divergence(by_sample = TRUE)
  expect_true(all(c("cell_type", "sample", "divergence") %in% names(by_sample)))
  correction_diagnostics <- correction$diagnostics()
  expect_equal(correction_diagnostics$sample, c("a", "b"))
  expect_true(all(correction_diagnostics$rules_without_local_evidence >= 0))
  expect_true(all(is.na(correction_diagnostics$removed_without_local_evidence) |
    (correction_diagnostics$removed_without_local_evidence >= 0 &
      correction_diagnostics$removed_without_local_evidence <= 1)))
  per_sample_summary <- correction$summary()
  expect_true("sample" %in% names(per_sample_summary))
  pooled_summary <- correction$summary(by_sample = FALSE)
  expect_true("fraction_molecules_removed" %in% names(pooled_summary))
  corrected <- correction$counts()
  expect_equal(nrow(corrected), nrow(counts))
  expect_equal(ncol(corrected), ncol(counts))
  expect_true(sum(corrected) <= sum(counts))
  expect_equal(dim(correction$counts(sample = "a")), dim(consensus$counts(sample = "a")))

  audit <- consensus$audit_admixture(min_excess = 1, min_target_cells = 2L,
    min_reference_cells = 2L)
  expect_s3_class(audit, "CellAdmixCollectionAudit")
  by_sample <- audit$pairs(by_sample = TRUE)
  expect_true("sample" %in% names(by_sample))
  pooled_pairs <- audit$pairs()
  expect_true(all(c("admixed_molecules", "rate", "n_samples") %in% names(pooled_pairs)))
  expect_true(all(c("source", "target") %in% names(by_sample)))
  report <- audit$evaluate(correction)
  expect_equal(report$samples(), c("a", "b"))
  expect_true("sample" %in% names(report$summary()))

  # Rerunning the consensus reuses the cached fixed-profile runs.
  again <- collection$fit_consensus()
  expect_equal(again$params$run_id, consensus$params$run_id)

  # Projection onto given profiles: sample a's own profiles with its restart
  # pool as the shared member pool.
  own <- cellAdmixCore:::.celladmix_fit_profiles_raw(collection$fits$a)
  projected <- collection$fit_profiles(own,
    members = cellAdmixCore:::.celladmix_fit_member_profiles(collection$fits$a),
    name = "reference_a")
  expect_equal(projected$params$source, "profiles")
  expect_equal(nrow(projected$loadings()), nrow(own))
  expect_equal(nrow(projected$matching), 4L)
  expect_true(all(is.finite(projected$loss$loss_gap)))
  for (fit in projected$fits) {
    expect_equal(fit$run$fixed_profile_members, 4L)
  }
})
