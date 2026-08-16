test_that("audit metric internals behave on planted contamination", {
  genes <- c(paste0("smk", 1:4), paste0("tmk", 1:4), "shared1")
  n_T <- 200; n_S <- 100
  cells <- c(paste0("T", seq_len(n_T)), paste0("S", seq_len(n_S)))
  types <- setNames(c(rep("T", n_T), rep("S", n_S)), cells)
  exposure <- setNames(c(rep(0:4, length.out = n_T), rep(0, n_S)), cells)
  counts <- Matrix::Matrix(0, length(genes), length(cells), sparse = TRUE,
    dimnames = list(genes, cells))
  counts[paste0("tmk", 1:4), types == "T"] <- c(50, 40, 30, 20)
  counts[paste0("smk", 1:4), types == "T"] <- 2
  counts["shared1", ] <- 30
  counts[paste0("smk", 1:4), types == "S"] <- 60
  planted <- counts * 0
  for (g in paste0("smk", 1:4)) {
    planted[g, types == "T"] <- 8 * exposure[types == "T"]
  }
  before <- counts + planted

  T_cells <- names(types)[types == "T"]
  totals <- Matrix::colSums(before)
  bins <- cellAdmixCore:::.celladmix_audit_bins(exposure[T_cells])
  profiles <- cellAdmixCore:::.celladmix_audit_profiles(before, types)
  baseline <- cellAdmixCore:::.celladmix_audit_pseudobulk(before,
    T_cells[exposure[T_cells] == 0])

  pool <- cellAdmixCore:::.celladmix_audit_marker_pool(profiles, "S", baseline,
    n_pool = 4)
  expect_setequal(pool, paste0("smk", 1:4))
  strict <- cellAdmixCore:::.celladmix_audit_pool_strict(
    c(pool, "shared1"), profiles, "S", baseline)
  expect_false("shared1" %in% strict)

  mcount <- cellAdmixCore:::.celladmix_audit_mcount
  rb <- cellAdmixCore:::.celladmix_audit_bin_rates(
    mcount(before, pool, T_cells), totals[T_cells], bins)
  expect_true(all(rb$rate_lo <= rb$rate & rb$rate <= rb$rate_hi))
  det <- cellAdmixCore:::.celladmix_audit_excess(rb)
  expect_lt(det$p, 1e-10)
  expect_equal(det$excess, sum(planted), tolerance = 0.15)

  for (case in list(list(before - planted, 1), list(before, 0),
                    list(before - planted / 2, 0.5))) {
    ra <- cellAdmixCore:::.celladmix_audit_bin_rates(
      mcount(case[[1]], pool, T_cells), totals[T_cells], bins)
    expect_equal(cellAdmixCore:::.celladmix_audit_power(rb, ra), case[[2]],
      tolerance = 0.06)
  }
})

test_that("audit end-to-end detects planted admixture on a grid dataset", {
  set.seed(11)
  n_side <- 18
  spacing <- 30
  cells <- expand.grid(ix = seq_len(n_side), iy = seq_len(n_side))
  cells$cell_id <- sprintf("c%03d", seq_len(nrow(cells)))
  cells$cell_type <- ifelse(cells$ix <= n_side / 2, "A", "B")
  cells$x <- cells$ix * spacing
  cells$y <- cells$iy * spacing

  genes_a <- paste0("a", 1:3); genes_b <- paste0("b", 1:3)
  tx <- list()
  for (i in seq_len(nrow(cells))) {
    ci <- cells[i, ]
    own <- if (ci$cell_type == "A") genes_a else genes_b
    n_own <- 60L
    g <- c(sample(own, n_own, replace = TRUE), rep("h1", 10L))
    # plant A-marker admixture into B cells by proximity to the boundary
    if (ci$cell_type == "B") {
      dist_cols <- ci$ix - n_side / 2
      n_adm <- max(0L, 12L - 6L * (as.integer(dist_cols) - 1L))
      if (n_adm > 0) g <- c(g, sample(genes_a, n_adm, replace = TRUE))
    }
    tx[[i]] <- data.frame(gene = g,
      x = ci$x + stats::runif(length(g), -10, 10),
      y = ci$y + stats::runif(length(g), -10, 10),
      z = 0, cell_id = ci$cell_id, cell_type = ci$cell_type,
      stringsAsFactors = FALSE)
  }
  tx <- do.call(rbind, tx)
  csv_path <- tempfile("audit_grid_", fileext = ".csv")
  utils::write.csv(tx, csv_path, row.names = FALSE)

  ds <- cellAdmix(csv_path, output_dir = tempfile("audit_grid_out_"),
    schema = celladmix_schema(x = "x", y = "y", z = "z", gene = "gene",
      cell = "cell_id", cell_type = "cell_type"),
    annotation = stats::setNames(cells$cell_type, cells$cell_id),
    annotation_name = "manual")
  fit <- ds$fit(rank = 2L, run_id = "audit_grid", ncv_k = 6L, graph_k = 4L,
    nmf_iterations = 20L, nmf_n_runs = 1L, nmf_train_max_rows = 200L,
    num_threads = 2L, seed = 11L)
  audit <- fit$audit_admixture(neighbor_k = 6L, min_target_cells = 50L,
    min_reference_cells = 20L, min_excess = 50)
  pairs <- audit$pairs()
  expect_true(all(c("source", "target", "excess", "q_value", "detected")
    %in% names(pairs)))
  ab <- pairs[pairs$source == "A" & pairs$target == "B", ]
  expect_equal(nrow(ab), 1L)
  expect_true(ab$detected)
  planted <- sum(tx$cell_type == "B" & tx$gene %in% genes_a)
  expect_gt(ab$excess, 0.5 * planted)
})
