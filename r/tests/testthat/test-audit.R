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
  # Dominant transfer channel deviating ~30% from the other genes' proportion:
  # within profile uncertainty once it carries the regression, although its
  # counting-noise significance is large.
  counts["bigch", ] <- rpois(n_T, 0.02 * 5000 * expo * 1.3 + 5)
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

test_that("genes below the proportion do not pull the screen's slope down", {
  # Three large channels and many small genes the profile overstates five-fold:
  # a fit weighted by the observed counts would take the small genes' level as
  # the proportion and call the large channels induced.
  set.seed(3)
  n_small <- 100L
  genes <- c(paste0("big", 1:3), paste0("small", seq_len(n_small)))
  psi <- setNames(c(rep(0.2, 3), rep(0.004, n_small)), genes)
  n_cells <- 2000L
  cells <- paste0("c", seq_len(n_cells))
  expo <- setNames(as.numeric(seq_len(n_cells) <= n_cells / 2), cells)
  scale <- 2 * 50000 / n_cells
  rows <- lapply(seq_along(genes), function(g) {
    rpois(n_cells, scale * psi[[g]] * expo * (if (g <= 3) 1 else 0.1) + 0.01)
  })
  counts <- Matrix::Matrix(do.call(rbind, rows), sparse = TRUE,
    dimnames = list(genes, cells))
  totals <- Matrix::colSums(counts) + 100
  screened <- cellAdmixCore:::.celladmix_audit_screen_panel(
    counts, candidates = genes, T_cells = cells, expo = expo,
    totals = setNames(totals, cells), source_profile = psi, n_pool = 20L)
  expect_length(intersect(screened$induced, paste0("big", 1:3)), 0)
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

# A split grid, A left and B right, with A material planted into B cells
# near the boundary. With shared both types express s0 (half of A's
# molecules, a little more of B's) and twice the A material is planted,
# carrying s0 at A's rate; with response the B cells of the two columns next
# to A express r0, which neither type expresses otherwise.
make_audit_fit <- function(seed = 11L, shared = FALSE, response = FALSE) {
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
    if (shared) {
      own_p <- if (ci$cell_type == "A") c(rep(0.5 / 3, 3), 0.5) else c(rep(0.15, 3), 0.55)
      g <- c(sample(c(own, "s0"), 60L, replace = TRUE, prob = own_p), rep("h1", 10L))
    } else {
      g <- c(sample(own, 60L, replace = TRUE), rep("h1", 10L))
    }
    if (ci$cell_type == "B") {
      d <- as.integer(ci$ix - n_side / 2) - 1L
      n_adm <- max(0L, 12L - 6L * d)
      if (n_adm > 0) {
        g <- c(g, if (shared) sample(c(genes_a, "s0"), 2L * n_adm, replace = TRUE,
          prob = c(rep(0.5 / 3, 3), 0.5)) else sample(genes_a, n_adm, replace = TRUE))
      }
      if (response) g <- c(g, rep("r0", max(0L, 8L - 4L * d)))
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
  model <- audit$fit_generative(num_threads = 2L)
  comp <- model$composition(source = "A", target = "B")
  expect_gt(max(comp$contamination), 0)
  # the fitted expectation adds up to the cells' molecules, and the planted
  # A molecules in B cells are admixture on the admixed expectation, not a
  # response on the contact term
  E <- model$expectation("B")
  expect_equal(dim(E$own), c(length(E$cells), nrow(before)))
  total <- E$own + E$ambient + E$floor + Reduce(`+`, E$admixed)
  for (I in E$induced) total <- total + I
  expect_lt(abs(sum(total) - sum(E$totals)) / sum(E$totals), 0.2)
  expect_error(model$expectation("no such type"), "no fitted pairs")
  det <- model$detect_responses(source = "A", target = "B", verbose = FALSE)
  expect_equal(nrow(det), nrow(before))
  expect_true(all(c("owner", "joint_z", "joint_called") %in% names(det)))
  a_rows <- det[det$gene %in% genes_a, ]
  expect_gt(sum(a_rows$admixed_fit), 0)
  expect_false(any(a_rows$called))
  # the joint test leaves out the owner's terms: A's genes are not tested against A
  expect_true(all(a_rows$owner == "A"))
  expect_false(any(a_rows$joint_called))
  det2 <- model$detect_responses(source = "A", target = "B", kernels = "both",
    range_um = c(20, 60), verbose = FALSE)
  expect_true(all(c("range_rate", "range_z", "range_share") %in% names(det2)))
  # retention of what was called changes nothing here; a response declared
  # on a1 reduces a1's removal in B cells and nothing else
  same <- model$retain_responses(det)
  expect_equal(sum(abs(same$correct()$counts() - after)), 0, tolerance = 1e-9)
  forced <- det[det$gene == "a1", ]
  forced$called <- TRUE; forced$contact_rate <- 0.5
  kept <- model$retain_responses(forced)
  after_kept <- kept$correct()$counts()
  expect_lt(sum(before["a1", b_cells] - after_kept["a1", b_cells]),
    sum(before["a1", b_cells] - after["a1", b_cells]))
  other <- setdiff(rownames(before), "a1")
  expect_equal(sum(abs(after_kept[other, ] - after[other, ])), 0, tolerance = 1e-9)
  c2 <- model$correct()
  expect_equal(sum(abs(c2$counts() - after)), 0, tolerance = 1e-6)
  removed_default <- sum(before) - sum(after)
  for (ini in c("pseudobulk", "clusters")) {
    mi <- audit$fit_generative(init = ini, num_threads = 2L)$correct()$counts()
    removed_i <- sum(before) - sum(mi)
    expect_gt(removed_i, 0.5 * removed_default)
    expect_lt(removed_i, 2 * removed_default)
  }
  report <- audit$evaluate(correction, yardstick_min_molecules = 50)
  p <- report$pairs()
  ab <- p[p$source == "A" & p$target == "B", ]
  expect_gt(ab$sensitivity, 0.5)
  # the excess yardstick: the A-marker molecules planted into B cells near
  # the boundary lie entirely above their level in B cells far from A, and
  # the correction removes most of them
  yg <- report$yardstick_genes()
  a_in_b <- yg[yg$cell_type == "B" & yg$gene %in% genes_a, ]
  expect_equal(nrow(a_in_b), 3)
  expect_true(all(a_in_b$shedding_type == "A"))
  expect_true(all(a_in_b$excess_share > 0.5))
  expect_true(all(a_in_b$removed_share > 0.5))
  y <- report$yardstick()
  expect_true("all" %in% y$cell_type)
  s <- report$summary()
  expect_true(all(c("yardstick_excess_share", "yardstick_removed_share",
    "yardstick_deviation") %in% names(s)))
  identity_correction <- list(counts = function() before)
  y0 <- audit$evaluate(identity_correction, warn_uncovered = FALSE,
    yardstick_min_molecules = 50)$yardstick()
  expect_equal(y0$removed_share[y0$cell_type == "all"], 0, tolerance = 1e-9)
})

test_that("detection calls a response and not admixture on a gene both types express", {
  fit <- make_audit_fit(shared = TRUE, response = TRUE)
  audit <- fit$audit_admixture(neighbor_k = 6L, min_target_cells = 50L,
    min_reference_cells = 20L, min_excess = 50)
  model <- audit$fit_generative(num_threads = 2L)
  det <- model$detect_responses(target = "B", verbose = FALSE)
  r0 <- det[det$source == "A" & det$gene == "r0", ]
  expect_true(r0$called)
  expect_true(r0$joint_called)
  expect_gt(r0$response_share, 0.5)
  # s0 is B's gene, carried into B cells by A's material: the admixed term
  # extended to genes both types express takes it, so it is no response
  s0 <- det[det$source == "A" & det$gene == "s0", ]
  expect_equal(s0$owner, "B")
  expect_gt(s0$admixed_fit, 0)
  expect_false(s0$called)
  expect_false(s0$joint_called)
})

test_that("the joint test finds an association shared by co-occurring terms", {
  set.seed(3)
  n <- 4000
  region <- stats::runif(n)
  c1 <- stats::rpois(n, 3 * region); c2 <- stats::rpois(n, 3 * region)
  own <- matrix(2, n, 1)
  y <- matrix(as.numeric(stats::rpois(n, 2 + 0.15 * (c1 + c2))), n, 1)
  fit <- .celladmix_poisson_terms(y, matrix(0, n, 1), list(own), cbind(c1, c2))
  b <- fit$coef[2:3, 1]
  S <- fit$cov[2:3, 2:3, 1]
  z <- b / sqrt(diag(S))
  joint <- as.numeric(t(b) %*% solve(S) %*% b)
  expect_gt(stats::qnorm(stats::pchisq(joint, 2, lower.tail = FALSE), lower.tail = FALSE), max(z) + 1)
  expect_equal(sum(b), 0.3, tolerance = 0.08 / 0.3)
})

test_that("the generative fit reports progress only when verbose", {
  fit <- make_audit_fit()
  audit <- fit$audit_admixture(neighbor_k = 6L, min_target_cells = 50L,
    min_reference_cells = 20L, min_excess = 50)
  quiet <- capture.output(audit$fit_generative(num_threads = 2L))
  expect_false(any(grepl("[INFO", quiet, fixed = TRUE)))
  log <- capture.output(m <- audit$fit_generative(num_threads = 2L,
    verbose = TRUE))
  info <- grep("^\\[INFO ", log, value = TRUE)
  expect_match(info[[1]], "Generative fit: .*2 types")
  expect_true(any(grepl("Fitting round 1/2, target B (1/2)", info,
    fixed = TRUE)))
  expect_true(any(grepl("expectation-maximization: .*log-likelihood", info)))
  expect_match(info[[length(info)]], "Generative fit finished")
  every <- capture.output(audit$fit_generative(num_threads = 2L,
    verbose = TRUE, progress_interval_sec = 0))
  expect_true(any(grepl("expectation-maximization: iteration 2/", every,
    fixed = TRUE)))
})

test_that("dose options run and the covariate control collapses", {
  fit <- make_audit_fit()
  audit <- fit$audit_admixture(neighbor_k = 6L, min_target_cells = 50L,
    min_reference_cells = 20L, min_excess = 50)
  before <- fit$counts()
  ann <- fit$dataset$annotation(fit$annotation_name, as_vector = TRUE)
  b_cells <- colnames(before)[!is.na(ann[colnames(before)]) &
    ann[colnames(before)] == "B"]
  genes_a <- paste0("a", 1:3)
  removed_a <- function(model) {
    after <- model$correct()$counts()
    sum(before[genes_a, b_cells]) - sum(after[genes_a, b_cells])
  }
  plain <- audit$fit_generative(num_threads = 2L)
  joint <- audit$fit_generative(num_threads = 2L, dose_mode = "joint",
    dose_ladder = c(12L, 24L), dose_reference_min_molecules = 500)
  expect_true(all(c("prior_unexposed_molecules", "cross_attributed_molecules")
    %in% names(joint$pairs)))
  expect_gt(removed_a(joint), 0.5 * removed_a(plain))
  control <- audit$fit_generative(num_threads = 2L, dose_mode = "joint",
    dose_ladder = c(12L, 24L), dose_reference_min_molecules = 500,
    use_ambient = FALSE, permute_covariates_seed = 3L)
  expect_lt(removed_a(control), 0.1 * removed_a(plain))
  expect_error(audit$fit_generative(num_threads = 2L, dose_mode = "bogus"))
})

test_that("the dataset-level audit matches the fit-level audit and runs NMF-free", {
  fit <- make_audit_fit()
  audit_fit <- fit$audit_admixture(neighbor_k = 6L, min_target_cells = 50L,
    min_reference_cells = 20L, min_excess = 50)
  ds <- fit$dataset
  audit_ds <- ds$audit_admixture(neighbor_k = 6L, min_target_cells = 50L,
    min_reference_cells = 20L, min_excess = 50)
  pf <- audit_fit$pairs()
  pd <- audit_ds$pairs()
  expect_equal(nrow(pd), nrow(pf))
  key <- function(p) paste(p$source, p$target)
  pd <- pd[match(key(pf), key(pd)), ]
  expect_equal(pd$excess, pf$excess, tolerance = 1e-6)
  expect_equal(pd$detected, pf$detected)

  # The full NMF-free path: dataset -> audit -> generative model ->
  # correction -> evaluation, with no factorization involved anywhere.
  model <- audit_ds$fit_generative(num_threads = 2L)
  expect_identical(model$init, "clusters")
  correction <- model$correct()
  report <- audit_ds$evaluate(correction)
  p <- report$pairs()
  ab <- p[p$source == "A" & p$target == "B", ]
  expect_gt(ab$sensitivity, 0.5)
})

test_that("the molecule-level realization and the cap follow the count-level removal", {
  fit <- make_audit_fit()
  audit <- fit$audit_admixture(neighbor_k = 6L, min_target_cells = 50L,
    min_reference_cells = 20L, min_excess = 50)
  model <- audit$fit_generative(num_threads = 2L)
  before <- fit$counts()
  expected <- model$correct()$counts()
  # realization: per cell and gene the realized removal is the expected one
  # up to the rounding of a fraction
  realized <- model$realize(fit, name = "generative_molecules", seed = 3L)
  got <- realized$counts()
  got <- got[rownames(before), colnames(before)]
  diff <- abs(as.matrix(before - got) - as.matrix(before - expected))
  expect_lte(max(diff), 1 + 1e-9)
  removed_expected <- sum(before - expected)
  expect_lt(abs(sum(before - got) - removed_expected), 0.1 * removed_expected + 5)
  expect_equal(realized$run$n_removed, sum(before - got))
  expect_true(!is.null(realized$run$realization))
  # the cap: both factors removed from B cells empties them without it;
  # with it the removal per cell and gene stays within the amount + 1
  rules <- data.frame(factor = c(1L, 2L), target_cell_type = c("B", "B"),
    stringsAsFactors = FALSE)
  uncapped <- celladmix_correct(fit$run, rules, out_dir = tempfile("cap_none_"))
  amounts <- cellAdmixCore:::.celladmix_removal_amounts(model)
  capped <- celladmix_correct(fit$run, rules, out_dir = tempfile("cap_model_"),
    limit = c(amounts, list(scale = 1, seed = 0)))
  removed_capped <- as.matrix(before - celladmix_collect_counts_sparse(capped)[rownames(before), colnames(before)])
  removed_uncapped <- as.matrix(before - celladmix_collect_counts_sparse(uncapped)[rownames(before), colnames(before)])
  amount <- as.matrix(model$removal())[rownames(before), colnames(before)]
  expect_true(all(removed_capped <= floor(amount) + 1 + 1e-9))
  expect_true(all(removed_capped <= removed_uncapped + 1e-9))
  expect_gt(sum(removed_uncapped), sum(removed_capped))
  expect_gt(sum(removed_capped), 0)
  expect_true(!is.null(capped$cap))
  expect_equal(capped$cap$restored, uncapped$n_removed - capped$n_removed)
  # a cap above every removal reproduces the uncapped correction
  huge <- list(removed = Matrix::Matrix(1e6, nrow(before), ncol(before),
    dimnames = dimnames(before), sparse = TRUE))
  same <- celladmix_correct(fit$run, rules, out_dir = tempfile("cap_huge_"),
    limit = c(cellAdmixCore:::.celladmix_removal_amounts(huge), list(scale = 1, seed = 0)))
  expect_equal(same$n_removed, uncapped$n_removed)
  expect_equal(same$cap$restored, 0)
})

test_that("the prior strength of the per-cell fraction can be estimated per pair", {
  fit <- make_audit_fit()
  audit <- fit$audit_admixture(neighbor_k = 6L, min_target_cells = 50L,
    min_reference_cells = 20L, min_excess = 50)
  est <- audit$fit_generative(num_threads = 2L, alpha_prior_strength = 0)
  expect_true("prior_strength" %in% names(est$pairs))
  expect_true(all(est$pairs$prior_strength >= 5 & est$pairs$prior_strength <= 1000))
  fixed <- audit$fit_generative(num_threads = 2L, alpha_prior_strength = 300)
  expect_true(all(fixed$pairs$prior_strength == 300))
  default <- audit$fit_generative(num_threads = 2L)
  expect_equal(default$pairs$prior_strength, est$pairs$prior_strength)
  expect_error(audit$fit_generative(num_threads = 2L, alpha_prior_strength = -1))
})
