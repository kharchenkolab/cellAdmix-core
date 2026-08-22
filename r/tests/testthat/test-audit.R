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

test_that("the reference rate uses the largest sufficiently populated neighborhood", {
  set.seed(5)
  n <- 4000
  totals <- rep(100, n)
  # Cells that pass the base zero-neighbor definition can still sit in
  # source-rich surroundings and carry contamination; requiring zero source
  # cells among progressively more neighbors isolates the clean ones.
  z15 <- rep(TRUE, n)
  z60 <- seq_len(n) > 1000
  z240 <- seq_len(n) > 3000
  rate <- ifelse(seq_len(n) <= 1000, 0.03,
    ifelse(seq_len(n) <= 3000, 0.012, 0.005))
  markers <- rpois(n, rate * totals)
  ref <- cellAdmixCore:::.celladmix_audit_reference_rate(markers, totals,
    list(k15 = z15, k60 = z60, k240 = z240))
  expect_identical(ref$kind, "k240")
  expect_lt(ref$rate, 0.008)
  expect_gt(ref$rate_unexposed, ref$rate)

  # Deepest neighborhood underpopulated: fall back to the next level.
  ref2 <- cellAdmixCore:::.celladmix_audit_reference_rate(markers, totals,
    list(k15 = z15, k60 = z60, k240 = seq_len(n) > 3990))
  expect_identical(ref2$kind, "k60")

  # Flat contamination-free profile: reference matches the base rate.
  markers_flat <- rpois(n, 0.01 * totals)
  ref_flat <- cellAdmixCore:::.celladmix_audit_reference_rate(
    markers_flat, totals, list(k15 = z15, k60 = z60, k240 = z240))
  expect_equal(ref_flat$rate, ref_flat$rate_unexposed, tolerance = 0.15)
})

test_that("panel screening excludes disproportionate induced genes", {
  set.seed(7)
  smk <- paste0("smk", 1:10)
  genes <- c(smk, "induced1", "bigch")
  n_T <- 400
  cells <- paste0("T", seq_len(n_T))
  expo <- rep(0:3, length.out = n_T)
  counts <- Matrix::Matrix(0, length(genes), n_T, sparse = TRUE,
    dimnames = list(genes, cells))
  # Transferred material: proportional to the source profile across genes.
  psi <- setNames(c(seq(600, 150, by = -50), 30, 5000), genes)
  for (g in smk) {
    counts[g, ] <- rpois(n_T, 0.02 * psi[[g]] * expo + 1)
  }
  # Induced gene: large exposure-linked excess despite a tiny profile share.
  counts["induced1", ] <- rpois(n_T, 40 * expo + 1)
  # Dominant transfer channel deviating ~15% from proportionality: within
  # profile uncertainty, although its counting-noise significance is large.
  counts["bigch", ] <- rpois(n_T, 0.02 * 5000 * expo * 1.15 + 5)
  totals <- Matrix::colSums(counts) + 500
  screen <- function(...) cellAdmixCore:::.celladmix_audit_screen_panel(
    counts, candidates = c("bigch", "induced1", smk),
    T_cells = cells, expo = expo, totals = setNames(totals, cells),
    source_profile = psi, n_pool = 11L, ...)
  screened <- screen()
  expect_true("induced1" %in% screened$induced)
  expect_false("induced1" %in% screened$pool)
  expect_false("bigch" %in% screened$induced)
  expect_setequal(screened$pool, c("bigch", smk))
  st <- screened$induced_stats
  expect_gt(st$fold[st$gene == "induced1"], 4)
  # Without the profile-uncertainty term the large channel's small relative
  # deviation becomes formally significant; the term is what prevents that.
  expect_true("bigch" %in% screen(profile_cv = 0)$induced)
})

test_that("the interface profile reflects the bordering source cells", {
  counts <- Matrix::Matrix(c(30000, 10000, 0, 40000), 2, 2, sparse = TRUE,
    dimnames = list(c("g1", "g2"), c("s_near", "s_far")))
  dist_T <- c(s_near = 10, s_far = 500)
  prof <- cellAdmixCore:::.celladmix_audit_interface_profile(
    counts, c("s_near", "s_far"), dist_T)
  expect_equal(prof[["g1"]], 0.75)
  # A sparse near subset widens until enough molecules are available.
  counts2 <- Matrix::Matrix(c(75, 25, 0, 40000), 2, 2, sparse = TRUE,
    dimnames = list(c("g1", "g2"), c("s_near", "s_far")))
  prof2 <- cellAdmixCore:::.celladmix_audit_interface_profile(
    counts2, c("s_near", "s_far"), dist_T)
  expect_lt(prof2[["g1"]], 0.01)
})

make_audit_fit <- function(seed = 11L) {
  set.seed(seed)
  n_side <- 18; spacing <- 30
  cells <- expand.grid(ix = seq_len(n_side), iy = seq_len(n_side))
  cells$cell_id <- sprintf("c%03d", seq_len(nrow(cells)))
  cells$cell_type <- ifelse(cells$ix <= n_side / 2, "A", "B")
  cells$x <- cells$ix * spacing; cells$y <- cells$iy * spacing
  genes_a <- paste0("a", 1:3); genes_b <- paste0("b", 1:3)
  tx <- list()
  for (i in seq_len(nrow(cells))) {
    ci <- cells[i, ]
    own <- if (ci$cell_type == "A") genes_a else genes_b
    g <- c(sample(own, 60L, replace = TRUE), rep("h1", 10L))
    if (ci$cell_type == "B") {
      n_adm <- max(0L, 12L - 6L * (as.integer(ci$ix - n_side / 2) - 1L))
      if (n_adm > 0) g <- c(g, sample(genes_a, n_adm, replace = TRUE))
    }
    tx[[i]] <- data.frame(gene = g,
      x = ci$x + stats::runif(length(g), -10, 10),
      y = ci$y + stats::runif(length(g), -10, 10),
      z = 0, cell_id = ci$cell_id, cell_type = ci$cell_type,
      stringsAsFactors = FALSE)
  }
  tx <- do.call(rbind, tx)
  csv_path <- tempfile("audit_fx_", fileext = ".csv")
  utils::write.csv(tx, csv_path, row.names = FALSE)
  ds <- cellAdmix(csv_path, output_dir = tempfile("audit_fx_out_"),
    schema = celladmix_schema(x = "x", y = "y", z = "z", gene = "gene",
      cell = "cell_id", cell_type = "cell_type"),
    annotation = stats::setNames(cells$cell_type, cells$cell_id),
    annotation_name = "manual")
  ds$fit(rank = 2L, run_id = "audit_fx", ncv_k = 6L, graph_k = 4L,
    nmf_iterations = 20L, nmf_n_runs = 1L, nmf_train_max_rows = 200L,
    num_threads = 2L, seed = seed)
}

test_that("evaluate warns from measured removal, not the rule list", {
  fit <- make_audit_fit()
  audit <- fit$audit_admixture(neighbor_k = 6L, min_target_cells = 50L,
    min_reference_cells = 20L, min_excess = 5)
  p <- audit$pairs()
  expect_true(all(c("reference_kind", "reference_inflation", "reference_trend", "n_induced")
    %in% names(p)))
  mk <- audit$markers(p$source[[1]], p$target[[1]])
  expect_true(all(c("pool", "strict", "induced") %in% names(mk)))
  if (any(p$detected)) {
    # A correction that removes nothing must trigger the low-removal
    # warning for every sufficiently large detected pair, regardless of
    # what its rule list claims.
    identity_correction <- list(counts = function() fit$counts(),
      rules = data.frame(source_cell_type = p$source[p$detected][[1]],
        target_cell_type = p$target[p$detected][[1]]))
    class(identity_correction) <- "list"
    big <- p[p$detected & p$excess >= audit$params$min_excess * 5, ,
      drop = FALSE]
    if (nrow(big)) {
      expect_warning(audit$evaluate(identity_correction),
        "Correction removed only")
    }
  }
})

test_that("the generative correction removes planted admixture", {
  fit <- make_audit_fit()
  audit <- fit$audit_admixture(neighbor_k = 6L, min_target_cells = 50L,
    min_reference_cells = 20L, min_excess = 50)
  correction <- audit$correct_generative(num_threads = 2L)
  before <- fit$counts()
  after <- correction$counts()
  genes_a <- paste0("a", 1:3); genes_b <- paste0("b", 1:3)
  ann <- fit$dataset$annotation(fit$annotation_name, as_vector = TRUE)
  b_cells <- colnames(before)[!is.na(ann[colnames(before)]) &
    ann[colnames(before)] == "B"]
  planted <- sum(before[genes_a, b_cells])
  removed_a <- sum(before[genes_a, b_cells]) - sum(after[genes_a, b_cells])
  removed_b <- sum(before[genes_b, b_cells]) - sum(after[genes_b, b_cells])
  expect_gt(removed_a, 0.4 * planted)
  # target-owned genes are structurally untouchable
  expect_equal(removed_b, 0, tolerance = 1e-9)
  comp <- correction$composition(source = "A", target = "B")
  expect_gt(max(comp$contamination), 0)
  report <- audit$evaluate(correction)
  p <- report$pairs()
  ab <- p[p$source == "A" & p$target == "B", ]
  expect_gt(ab$sensitivity, 0.5)
})
