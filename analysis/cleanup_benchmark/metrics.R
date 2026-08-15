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

# Source-specific marker pool for an ordered pair: among genes whose
# top-expressing type is S, the top n by contrast against the exposure-0
# target baseline. Rank-based on purpose - an absolute baseline cutoff would
# be skewed by the very contamination being measured, and the power metric
# subtracts the baseline anyway.
bench_marker_pool <- function(profiles, S, T_type, baseline_cpm_T0,
                              n_pool = 20L, min_source_cpm = 50,
                              baseline_floor = 20) {
  top_type <- colnames(profiles)[max.col(profiles, ties.method = "first")]
  cand <- rownames(profiles)[top_type == S & profiles[, S] >= min_source_cpm]
  contrast <- profiles[cand, S] / pmax(baseline_cpm_T0[cand], baseline_floor)
  cand[order(contrast, profiles[cand, S], decreasing = TRUE)][
    seq_len(min(n_pool, length(cand)))]
}

# Abundance-matched control genes: similar overall cpm to the markers but not
# clearly owned by ANY cell type, so that leakage from other colocated types
# does not masquerade as a confound signal. What remains measures
# microenvironment/density gradients.
bench_control_genes <- function(profiles, S, markers, n_per_marker = 2L) {
  norm <- t(t(profiles) / pmax(colSums(profiles), 1))
  ownership <- apply(norm, 1, max) / pmax(rowSums(norm), 1e-12)
  total <- rowMeans(profiles)
  pool <- setdiff(rownames(profiles)[ownership < stats::median(ownership, na.rm = TRUE)],
    markers)
  unique(unlist(lapply(markers, function(g) {
    pool[order(abs(log1p(total[pool]) - log1p(total[[g]])))][seq_len(n_per_marker)]
  })))
}

# Strict subpool: pool genes with near-zero native baseline in exposure-0
# target cells. Their excess is unambiguous contamination; the broad pool adds
# shared genes whose exposure gradient may include real biology.
bench_pool_strict <- function(pool, profiles, S, baseline_cpm_T0, baseline_frac = 0.05) {
  pool[baseline_cpm_T0[pool] < baseline_frac * profiles[pool, S]]
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

# Removal depth: how much of the exposure-0 marker signal was removed. For
# highly source-specific markers that signal is mostly ambient/dispersed
# contamination, so this measures cleanup reach, not damage.
bench_removal_depth <- function(rates_before, rates_after) {
  r0b <- rates_before$rate[rates_before$bin == "0"]
  r0a <- rates_after$rate[rates_after$bin == "0"]
  if (!length(r0b) || r0b <= 0) {
    return(NA_real_)
  }
  1 - r0a / r0b
}

# Per-gene power: the excess-removal fraction gene by gene, exposing genes
# that escape cleanup inside an otherwise well-cleaned pool.
bench_power_per_gene <- function(counts_before, counts_after, pool, cells_T,
                                 totals, bins) {
  do.call(rbind, lapply(pool, function(g) {
    gb <- as.numeric(counts_before[g, cells_T])
    ga <- as.numeric(counts_after[g, cells_T])
    rb <- bench_bin_rates(gb, totals, bins)
    ra <- bench_bin_rates(ga, totals, bins)
    det <- bench_detect(rb)
    data.frame(gene = g, excess_before = round(det$excess_molecules, 1),
      power = round(bench_power(rb, ra), 3))
  }))
}

# Safety 2: retention of the target's own markers, worst bin. Catches
# overcorrection that strips native expression from exposed cells.
bench_safety_native <- function(native_before, native_after, totals, bins) {
  rb <- bench_bin_rates(native_before, totals, bins)
  ra <- bench_bin_rates(native_after, totals, bins)
  ret <- ifelse(rb$rate > 0, ra$rate / rb$rate, NA)
  min(ret, na.rm = TRUE)
}

# Safety 3: profile integrity of exposure-0 target cells, restricted to the
# target's own genes (top-expressing type is T). Removal of drifted foreign
# material should not count against integrity; erosion of native genes does.
bench_profile_integrity <- function(counts_before, counts_after, cells_e0,
                                    profiles = NULL, T_type = NULL) {
  pb <- bench_pseudobulk(counts_before, cells_e0)
  pa <- bench_pseudobulk(counts_after, cells_e0)
  genes <- rownames(counts_before)
  if (!is.null(profiles) && !is.null(T_type)) {
    top_type <- colnames(profiles)[max.col(profiles, ties.method = "first")]
    genes <- rownames(profiles)[top_type == T_type]
  }
  ok <- genes[is.finite(pb[genes]) & is.finite(pa[genes])]
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

