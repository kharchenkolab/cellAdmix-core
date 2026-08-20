celladmix_simulate_nsclc <- cellAdmixCore:::celladmix_simulate_nsclc

test_that("fixed-H fits with more factors than the derived rank keep all fractions", {
  td <- tempfile("xenium_fixedh_mock_")
  sim <- celladmix_simulate_nsclc(
    transcripts_per_cell = 40L,
    admixture_per_target_cell = 10L,
    seed = 23L
  )
  write_mock_xenium_bundle(td, sim, include_cell_type = TRUE)
  annotation <- data.frame(
    cell_id = as.character(sim$cells$cell_id),
    cell_type = as.character(sim$cells$cell_type),
    stringsAsFactors = FALSE
  )
  ds <- cellAdmix(td, output_dir = tempfile("celladmix_fixedh_out_"),
    annotation = annotation, annotation_name = "cell_type",
    annotation_col = "cell_type", cell_id_col = "cell_id")
  fit <- ds$fit(annotation = "cell_type", run_id = "base", overwrite = TRUE,
    ncv_k = 8L, rank = 2L, graph_k = 4L, nmf_iterations = 60L,
    nmf_train_max_rows = 40L, nmf_n_runs = 1L, seed = 23L)

  # Extend the fitted loadings by one anchor row and refit with fixed H:
  # the run must carry three factors end to end, including cell fractions.
  H0 <- fit$loadings()
  anchor <- rep(1 / ncol(H0), ncol(H0)) * stats::median(rowSums(H0))
  H3 <- rbind(H0, anchor)
  colnames(H3) <- fit$run$genes

  fit_a <- ds$fit(annotation = "cell_type", run_id = "anchored",
    overwrite = TRUE, nmf_fixed_h = H3, nmf_variant = "kl", nmf_n_runs = 1L,
    ncv_k = 8L, graph_k = 4L, seed = 23L)
  expect_equal(as.integer(fit_a$rank), 3L)
  cells <- fit_a$cell_factors()
  expect_true(all(paste0("factor_", 1:3, "_fraction") %in% names(cells)))
  fr <- as.matrix(cells[, paste0("factor_", 1:3, "_fraction")])
  expect_true(all(rowSums(fr) <= 1 + 1e-8))

  # The native check must be able to evaluate rules on the extra factor.
  rules <- data.frame(factor = 3L, target_cell_type = "fibroblast",
    source_cell_type = "malignant", stringsAsFactors = FALSE)
  corr <- cellAdmixCore:::celladmix_correct(fit_a$run, rules,
    out_dir = tempfile("fixedh_corr_"))
  expect_true(is.numeric(corr$n_removed))
})
