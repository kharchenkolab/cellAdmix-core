# Cleanup-benchmark metrics: spatial-exposure self-grading of admixture removal.
#
# The evaluation unit is an ordered cell-type pair (source S -> target T).
# Source-cell exposure of each target cell (S-cells among its k nearest
# cells) stratifies target cells; source-specific marker excess in exposed
# cells over the exposure-0 baseline is the contamination signal. Cleanup
# power is the fraction of that excess removed; safety metrics check that
# baseline and native expression survive.
#
# All "after" rates share the "before" per-cell totals so that removal reads
# as removal instead of being hidden by renormalization.

suppressMessages(library(Matrix))

bench_exposure_bins <- function(exposure) {
  cut(exposure, breaks = c(-0.5, 0.5, 1.5, 2.5, Inf), labels = c("0", "1", "2", "3+"))
}

# Column-normalized cpm pseudobulk profile for one set of cells.
bench_pseudobulk <- function(counts, cells, scale = 1e6) {
  cells <- intersect(cells, colnames(counts))
  if (!length(cells)) {
    return(setNames(rep(NA_real_, nrow(counts)), rownames(counts)))
  }
  sums <- Matrix::rowSums(counts[, cells, drop = FALSE])
  sums / max(sum(sums), 1) * scale
}

# Per-type cpm profiles, ordered gene x type matrix.
bench_type_profiles <- function(counts, cell_types) {
  types <- sort(unique(cell_types[!is.na(cell_types)]))
  out <- sapply(types, function(t) {
    bench_pseudobulk(counts, names(cell_types)[!is.na(cell_types) & cell_types == t])
  })
  rownames(out) <- rownames(counts)
  out
}

# Source-specific marker pool for an ordered pair: genes strongly specific to
# S among all types, with a near-zero native baseline in exposure-0 T cells.
bench_marker_pool <- function(profiles, S, T_type, baseline_cpm_T0,
                              n_pool = 20L, baseline_frac = 0.05,
                              min_source_cpm = 50) {
  others <- setdiff(colnames(profiles), S)
  max_other <- apply(profiles[, others, drop = FALSE], 1, max)
  spec <- profiles[, S] / pmax(profiles[, S] + max_other, 1e-9)
  ok <- profiles[, S] >= min_source_cpm &
    baseline_cpm_T0 < baseline_frac * profiles[, S]
  cand <- rownames(profiles)[ok]
  cand[order(spec[cand], profiles[cand, S], decreasing = TRUE)][
    seq_len(min(n_pool, length(cand)))]
}

# Abundance-matched control genes: similar overall cpm to the markers but not
# S-specific. Used to separate admixture removal from microenvironment shifts.
bench_control_genes <- function(profiles, S, markers, n_per_marker = 2L) {
  others <- setdiff(colnames(profiles), S)
  max_other <- apply(profiles[, others, drop = FALSE], 1, max)
  spec <- profiles[, S] / pmax(profiles[, S] + max_other, 1e-9)
  total <- rowMeans(profiles)
  pool <- setdiff(rownames(profiles)[spec < stats::median(spec, na.rm = TRUE)], markers)
  unique(unlist(lapply(markers, function(g) {
    pool[order(abs(log1p(total[pool]) - log1p(total[[g]])))][seq_len(n_per_marker)]
  })))
}

# Exposure-binned pooled rates: marker molecules / cell totals per bin.
bench_bin_rates <- function(marker_counts, totals, bins) {
  agg_m <- tapply(marker_counts, bins, sum, default = 0)
  agg_t <- tapply(totals, bins, sum, default = 0)
  n <- tapply(rep(1, length(bins)), bins, sum, default = 0)
  data.frame(bin = names(agg_m), markers = as.numeric(agg_m),
    totals = as.numeric(agg_t), n_cells = as.numeric(n),
    rate = as.numeric(agg_m) / pmax(as.numeric(agg_t), 1))
}

# Detection test: is there a significant exposure-linked marker excess?
# Poisson test of pooled exposed counts against the exposure-0 rate.
bench_detect <- function(rates) {
  r0 <- rates$rate[rates$bin == "0"]
  exposed <- rates[rates$bin != "0", , drop = FALSE]
  m <- sum(exposed$markers)
  expect <- r0 * sum(exposed$totals)
  if (!length(r0) || is.na(r0) || expect <= 0 || m <= expect) {
    return(list(p = 1, excess_molecules = 0))
  }
  p <- stats::ppois(m - 1, expect, lower.tail = FALSE)
  list(p = p, excess_molecules = m - expect)
}

# Core power metric: fraction of the exposure-linked excess eliminated.
# excess is measured per bin against that state's own exposure-0 rate, so
# baseline erosion registers in safety, not power.
bench_power <- function(rates_before, rates_after) {
  eb <- merge(rates_before, rates_after, by = "bin", suffixes = c("_b", "_a"))
  eb <- eb[eb$bin != "0", , drop = FALSE]
  r0b <- rates_before$rate[rates_before$bin == "0"]
  r0a <- rates_after$rate[rates_after$bin == "0"]
  excess_b <- pmax(eb$rate_b - r0b, 0) * eb$totals_b
  excess_a <- pmax(eb$rate_a - r0a, 0) * eb$totals_b
  if (sum(excess_b) <= 0) {
    return(NA_real_)
  }
  1 - sum(excess_a) / sum(excess_b)
}

# Dose-response power: same idea on a (near-)continuous exposure dose,
# weighting each dose level's excess by its molecule mass.
bench_dose_power <- function(marker_before, marker_after, totals, dose, max_dose = 10L) {
  dose <- pmin(dose, max_dose)
  rb <- bench_bin_rates(marker_before, totals, factor(dose))
  ra <- bench_bin_rates(marker_after, totals, factor(dose))
  names(rb)[1] <- names(ra)[1] <- "bin"
  r0b <- rb$rate[rb$bin == "0"]
  r0a <- ra$rate[ra$bin == "0"]
  keep <- rb$bin != "0"
  excess_b <- pmax(rb$rate[keep] - r0b, 0) * rb$totals[keep]
  excess_a <- pmax(ra$rate[keep] - r0a, 0) * rb$totals[keep]
  if (!length(excess_b) || sum(excess_b) <= 0) {
    return(NA_real_)
  }
  1 - sum(excess_a) / sum(excess_b)
}

# Safety 1: marker retention in exposure-0 target cells (native floor).
bench_safety_baseline <- function(rates_before, rates_after) {
  r0b <- rates_before$rate[rates_before$bin == "0"]
  r0a <- rates_after$rate[rates_after$bin == "0"]
  if (!length(r0b) || r0b <= 0) {
    return(NA_real_)
  }
  r0a / r0b
}

# Safety 2: retention of the target's own markers, worst bin. Catches
# overcorrection that strips native expression from exposed cells.
bench_safety_native <- function(native_before, native_after, totals, bins) {
  rb <- bench_bin_rates(native_before, totals, bins)
  ra <- bench_bin_rates(native_after, totals, bins)
  ret <- ifelse(rb$rate > 0, ra$rate / rb$rate, NA)
  min(ret, na.rm = TRUE)
}

# Safety 3: whole-profile integrity of exposure-0 target cells.
bench_profile_integrity <- function(counts_before, counts_after, cells_e0) {
  pb <- bench_pseudobulk(counts_before, cells_e0)
  pa <- bench_pseudobulk(counts_after, cells_e0)
  ok <- is.finite(pb) & is.finite(pa)
  stats::cor(log1p(pb[ok]), log1p(pa[ok]))
}

# Small non-negative least squares via multiplicative coordinate descent.
bench_nnls <- function(y, X, n_iter = 2000L) {
  w <- rep(mean(y) / max(colMeans(X), 1e-9) / ncol(X), ncol(X))
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  for (i in seq_len(n_iter)) {
    grad <- XtX %*% w - Xty
    step <- 1 / max(diag(XtX))
    w <- pmax(w - step * as.numeric(grad), 0)
  }
  setNames(as.numeric(w), colnames(X))
}

# Marker-free contamination estimate: decompose the exposed-target pseudobulk
# over all type profiles; the S coefficient share is the contamination index.
# Linear cpm space: mixtures are linear in rates, not in log space.
bench_nnls_contamination <- function(counts, profiles_e0, cells_exposed, S, T_type) {
  y <- bench_pseudobulk(counts, cells_exposed)
  ok <- is.finite(y) & apply(profiles_e0, 1, function(r) all(is.finite(r)))
  w <- bench_nnls(y[ok], profiles_e0[ok, , drop = FALSE])
  sum_w <- sum(w)
  if (sum_w <= 0) {
    return(NA_real_)
  }
  unname(w[S] / sum_w)
}

# Per-cell doublet-axis shift: project target cells on the S-minus-T profile
# axis; contamination pulls exposed cells toward S. Unlike the rate metrics,
# this is a compositional measure, so each state uses its own cell totals.
bench_doublet_shift <- function(counts, cells_T, exposure, profile_S, profile_T,
                                n_genes = 40L) {
  spec_S <- profile_S / pmax(profile_S + profile_T, 1e-9)
  genes <- unique(c(
    names(sort(spec_S, decreasing = TRUE))[seq_len(n_genes)],
    names(sort(spec_S, decreasing = FALSE))[seq_len(n_genes)]))
  genes <- intersect(genes, rownames(counts))
  axis <- log1p(profile_S[genes]) - log1p(profile_T[genes])
  axis <- axis / sqrt(sum(axis^2))
  cells_T <- intersect(cells_T, colnames(counts))
  totals <- Matrix::colSums(counts[, cells_T, drop = FALSE])
  cpm <- t(t(as.matrix(counts[genes, cells_T, drop = FALSE])) /
    pmax(totals, 1)) * 1e6
  scores <- as.numeric(t(log1p(cpm)) %*% axis)
  e <- exposure[cells_T]
  mean(scores[e > 0]) - mean(scores[e == 0])
}
