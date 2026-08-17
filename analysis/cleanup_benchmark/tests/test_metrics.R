# Synthetic-ground-truth tests for the cleanup benchmark metrics.
# A planted-contamination scenario with known exposure structure checks that
# each metric reads 1/0/intermediate where it must.

suppressMessages(library(Matrix))
suppressMessages(library(testthat))
source(file.path(Sys.getenv("BENCH_DIR", "benchmarks/cleanup_gradient"), "metrics.R"))

make_scenario <- function() {
  genes <- c(paste0("smk", 1:4), paste0("tmk", 1:4), paste0("bg", 1:4), "shared1")
  n_T <- 200; n_S <- 100; n_O <- 100
  cells <- c(paste0("T", seq_len(n_T)), paste0("S", seq_len(n_S)), paste0("O", seq_len(n_O)))
  types <- setNames(c(rep("T", n_T), rep("S", n_S), rep("O", n_O)), cells)
  exposure <- setNames(c(rep(0:4, length.out = n_T), rep(0, n_S + n_O)), cells)

  counts <- Matrix(0, nrow = length(genes), ncol = length(cells), sparse = TRUE,
    dimnames = list(genes, cells))
  # T cells: native tmk high, bg mid, shared mid, smk small native baseline
  counts[paste0("tmk", 1:4), types == "T"] <- c(50, 40, 30, 20)
  counts[paste0("bg", 1:4), types == "T"] <- 20
  counts["shared1", types == "T"] <- 30
  counts[paste0("smk", 1:4), types == "T"] <- 2
  # planted contamination: 8 molecules per marker per exposure unit
  planted <- Matrix(0, nrow = length(genes), ncol = length(cells), sparse = TRUE,
    dimnames = list(genes, cells))
  for (g in paste0("smk", 1:4)) {
    planted[g, types == "T"] <- 8 * exposure[types == "T"]
  }
  # S cells: smk high, shared mid; O cells: bg only
  counts[paste0("smk", 1:4), types == "S"] <- 60
  counts["shared1", types == "S"] <- 30
  counts[paste0("bg", 1:4), types == "O"] <- 40
  before <- counts + planted
  list(genes = genes, cells = cells, types = types, exposure = exposure,
    base = counts, planted = planted, before = before)
}

sc <- make_scenario()
totals <- Matrix::colSums(sc$before)
T_cells <- names(sc$types)[sc$types == "T"]
bins <- bench_exposure_bins(sc$exposure[T_cells])
profiles <- bench_type_profiles(sc$before, sc$types)
baseline_T0 <- bench_pseudobulk(sc$before, T_cells[sc$exposure[T_cells] == 0])

marker_counts <- function(m, genes, cells) Matrix::colSums(m[genes, cells, drop = FALSE])

test_that("marker pool ranks source-specific genes above shared ones", {
  pool <- bench_marker_pool(profiles, "S", "T", baseline_T0, n_pool = 4)
  expect_setequal(pool, paste0("smk", 1:4))
  full <- bench_marker_pool(profiles, "S", "T", baseline_T0, n_pool = 20)
  # shared1 is top-expressed in S so it may enter a wide pool, but only
  # after every truly specific gene
  expect_true(all(match(paste0("smk", 1:4), full) < match("shared1", full)))
})

test_that("detection fires on the contaminated pair and not on a clean one", {
  pool <- paste0("smk", 1:4)
  rates <- bench_bin_rates(marker_counts(sc$before, pool, T_cells), totals[T_cells], bins)
  det <- bench_detect(rates)
  expect_lt(det$p, 1e-10)
  # O markers in T cells: no contamination planted
  rates_o <- bench_bin_rates(marker_counts(sc$before, paste0("bg", 1:4), T_cells),
    totals[T_cells], bins)
  expect_gt(bench_detect(rates_o)$p, 0.05)
})

test_that("power reads 1 / 0 / 0.5 for perfect / no / half cleanup", {
  pool <- paste0("smk", 1:4)
  rb <- bench_bin_rates(marker_counts(sc$before, pool, T_cells), totals[T_cells], bins)
  perfect <- sc$before - sc$planted
  none <- sc$before
  half <- sc$before - sc$planted / 2
  for (case in list(list(perfect, 1), list(none, 0), list(half, 0.5))) {
    ra <- bench_bin_rates(marker_counts(case[[1]], pool, T_cells), totals[T_cells], bins)
    expect_equal(bench_power(rb, ra), case[[2]], tolerance = 0.02)
  }
})

test_that("baseline erosion shows in safety, not power", {
  pool <- paste0("smk", 1:4)
  # overzealous: removes planted AND half the native tmk in exposed cells
  over <- sc$before - sc$planted
  exposed <- T_cells[sc$exposure[T_cells] > 0]
  over[paste0("tmk", 1:4), exposed] <- over[paste0("tmk", 1:4), exposed] / 2
  rb <- bench_bin_rates(marker_counts(sc$before, pool, T_cells), totals[T_cells], bins)
  ra <- bench_bin_rates(marker_counts(over, pool, T_cells), totals[T_cells], bins)
  expect_equal(bench_power(rb, ra), 1, tolerance = 0.02)
  expect_equal(bench_removal_depth(rb, ra), 0, tolerance = 0.02)
  ret <- bench_safety_native(
    marker_counts(sc$before, paste0("tmk", 1:4), T_cells),
    marker_counts(over, paste0("tmk", 1:4), T_cells),
    totals[T_cells], bins)
  expect_lt(ret, 0.6)
  # gene guillotine (all marker molecules removed everywhere): full power,
  # full removal depth
  over2 <- sc$before
  over2[paste0("smk", 1:4), T_cells] <- 0
  ra2 <- bench_bin_rates(marker_counts(over2, pool, T_cells), totals[T_cells], bins)
  expect_equal(bench_power(rb, ra2), 1, tolerance = 1e-6)
  expect_gt(bench_removal_depth(rb, ra2), 0.99)
  ret_perfect <- bench_safety_native(
    marker_counts(sc$before, paste0("tmk", 1:4), T_cells),
    marker_counts(sc$before - sc$planted, paste0("tmk", 1:4), T_cells),
    totals[T_cells], bins)
  expect_equal(ret_perfect, 1, tolerance = 0.02)
})

test_that("strict subpool drops genes with real native baseline", {
  pool <- c(paste0("smk", 1:4), "shared1")
  strict <- bench_pool_strict(pool, profiles, "S", baseline_T0)
  expect_false("shared1" %in% strict)
  expect_setequal(strict, paste0("smk", 1:4))
})

test_that("per-gene power exposes a single escaping gene", {
  pool <- paste0("smk", 1:4)
  partial <- sc$before - sc$planted
  partial["smk1", ] <- sc$before["smk1", ]  # smk1 escapes cleanup entirely
  pg <- bench_power_per_gene(sc$before, partial, pool, T_cells,
    totals[T_cells], bins)
  expect_equal(pg$power[pg$gene == "smk1"], 0, tolerance = 0.02)
  expect_true(all(pg$power[pg$gene != "smk1"] > 0.98))
})

test_that("native-restricted integrity ignores foreign-gene removal", {
  e0 <- T_cells[sc$exposure[T_cells] == 0]
  expect_equal(bench_profile_integrity(sc$before, sc$before - sc$planted, e0), 1,
    tolerance = 1e-6)
  # deep cleanup that also strips the S-gene baseline from e0 cells:
  # T-native-restricted integrity stays 1
  deep <- sc$before - sc$planted
  deep[paste0("smk", 1:4), e0] <- 0
  nat <- bench_profile_integrity(sc$before, deep, e0, profiles = profiles, T_type = "T")
  expect_equal(nat, 1, tolerance = 1e-6)
  # but eroding a native gene lowers it
  worse <- deep
  worse["tmk1", e0] <- worse["tmk1", e0] / 4
  expect_lt(bench_profile_integrity(sc$before, worse, e0,
    profiles = profiles, T_type = "T"), 0.995)
})

test_that("nnls recovers a planted mixture", {
  X <- cbind(a = c(10, 0, 5), b = c(0, 8, 5))
  y <- 0.7 * X[, "a"] + 0.3 * X[, "b"]
  w <- bench_nnls(y, X)
  expect_equal(unname(w / sum(w)), c(0.7, 0.3), tolerance = 0.02)
})

test_that("nnls contamination index drops after perfect cleanup", {
  profiles_e0 <- bench_type_profiles(sc$base, sc$types)  # uncontaminated refs
  exposed <- T_cells[sc$exposure[T_cells] > 0]
  before_idx <- bench_nnls_contamination(sc$before, profiles_e0, exposed, "S", "T")
  after_idx <- bench_nnls_contamination(sc$before - sc$planted, profiles_e0, exposed, "S", "T")
  expect_gt(before_idx, after_idx + 0.02)
  expect_lt(after_idx, 0.05)
})

cat("all metric tests defined\n")
