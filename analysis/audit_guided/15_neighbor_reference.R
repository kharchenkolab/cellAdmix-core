# Can the ambient reference be recovered from the in-plane neighbor-count
# dependency alone, without comparing to distant tissue regions?
#
# Three measurements per pair:
#  1. Reference versus neighborhood size: the pooled marker rate of cells
#     with zero source neighbors among their K nearest, for growing K.
#     Hidden exposure shrinks as K grows (a cell with 0 of 120 neighbors is
#     far cleaner than one with 0 of 15), while the lateral radius spanned
#     by K neighbors grows slowly - if the curve converges before the
#     radius reaches different tissue compartments, a neighbor-only
#     reference is viable.
#  2. Beta-binomial extrapolation at K = 15: the neighbor count k samples
#     the cell's true source surroundings, so k follows a beta-binomial
#     whose alpha parameter measures the hidden exposure of zero-count
#     cells; under that model the dose-response line rate(k) reaches the
#     true zero-exposure level at k = -alpha, an extrapolation that never
#     leaves the local neighborhood mix.
#  3. Compartment drift: correlation of the target cells' own expression
#     profile (marker panel excluded) between zero-neighbor cells near
#     (<60 um) and far (>250 um) from source regions. Low correlation
#     supports the concern that distant same-type cells are biologically
#     different, favoring the neighbor-based reference.
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(cellAdmixCore)); suppressMessages(library(Matrix))
ag_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/audit_guided"

setwd("/home/pkharchenko/cellAdmix/cellAdmix-core/examples/xenium_pancreas_membrane_377_full")
annotation <- read.csv("annotations/annotation.csv.gz", stringsAsFactors = FALSE)
cell_annotation <- setNames(annotation$merged_annotation, annotation$cell_id)
ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)
fit <- ds$fit(nmf_variant = "invsqrt_kl", verbose = FALSE)

cells <- fit$cell_factors()
counts <- fit$counts()
ann <- cell_annotation
audit <- fit$audit_admixture()
pairs <- audit$pairs(detected_only = TRUE)
totals <- Matrix::colSums(counts)
nearest <- cellAdmixCore:::.celladmix_source_nearest_distance(cells, ann)

Ks <- c(15L, 30L, 60L, 120L, 240L)
expo_K <- lapply(Ks, function(K)
  cellAdmixCore:::.celladmix_source_exposure_counts(cells, ann, K))
for (i in seq_along(Ks)) rownames(expo_K[[i]]$counts) <- as.character(cells$cell_id)

# Median lateral radius spanned by K neighbors (sample of cells).
set.seed(3)
sample_cells <- sample(nrow(cells), 3000)
radius_K <- vapply(Ks, function(K) {
  d <- FNN_dist <- NULL
  # distance to the K-th nearest cell via the count helper is unavailable;
  # approximate from density: r_K = sqrt(K / (pi * local_density)).
  area <- diff(range(cells$x)) * diff(range(cells$y))
  sqrt(K / (pi * nrow(cells) / area))
}, 0)

bb_fit <- function(k, N = 15L) {
  nll <- function(par) {
    a <- exp(par[1]); b <- exp(par[2])
    -sum(lchoose(N, k) + lbeta(k + a, N - k + b) - lbeta(a, b))
  }
  o <- optim(c(log(0.2), log(3)), nll, method = "Nelder-Mead")
  list(alpha = exp(o$par[1]), beta = exp(o$par[2]))
}

rows <- list()
for (pi_ in seq_len(nrow(pairs))) {
  S <- pairs$source[[pi_]]; T_type <- pairs$target[[pi_]]
  mk <- audit$markers(S, T_type)
  T_cells <- names(ann)[!is.na(ann) & ann == T_type]
  T_cells <- intersect(T_cells, colnames(counts))
  T_cells <- intersect(T_cells, rownames(expo_K[[1]]$counts))
  pool_counts <- cellAdmixCore:::.celladmix_audit_mcount(counts, mk$pool, T_cells)
  tt <- totals[T_cells]
  dist_S <- nearest[T_cells, S]

  # (1) zero-neighbor reference versus K
  ref_K <- vapply(seq_along(Ks), function(i) {
    e <- expo_K[[i]]$counts[T_cells, S]
    z <- e == 0
    if (sum(tt[z]) < 2e4) return(NA_real_)
    sum(pool_counts[z]) / sum(tt[z])
  }, 0)

  # (2) beta-binomial extrapolation at K = 15
  e15 <- expo_K[[1]]$counts[T_cells, S]
  bb <- tryCatch(bb_fit(as.integer(round(e15))), error = function(e) NULL)
  ref_bb <- NA_real_
  if (!is.null(bb)) {
    ks <- sort(unique(as.integer(round(e15))))
    rate_k <- vapply(ks, function(k) {
      z <- as.integer(round(e15)) == k
      if (sum(z) < 30) return(NA_real_)
      sum(pool_counts[z]) / max(sum(tt[z]), 1)
    }, 0)
    wt_k <- vapply(ks, function(k) sum(tt[as.integer(round(e15)) == k]), 0)
    okk <- is.finite(rate_k) & ks <= 8
    if (sum(okk) >= 3) {
      fit_lm <- stats::lm(rate_k[okk] ~ ks[okk], weights = wt_k[okk])
      ref_bb <- max(0, unname(stats::coef(fit_lm)[1] -
        stats::coef(fit_lm)[2] * bb$alpha))
    }
  }

  # (3) compartment drift among zero-neighbor cells
  z15 <- e15 == 0
  near <- T_cells[z15 & dist_S < 60]
  far <- T_cells[z15 & dist_S > 250]
  drift_cor <- NA_real_
  if (length(near) >= 100 && length(far) >= 100) {
    pb <- function(cs) {
      v <- Matrix::rowSums(counts[, cs, drop = FALSE])
      v[setdiff(rownames(counts), mk$pool)]
    }
    a <- log1p(pb(near) / max(sum(tt[near]), 1) * 1e4)
    b <- log1p(pb(far) / max(sum(tt[far]), 1) * 1e4)
    drift_cor <- stats::cor(a, b)
  }

  rows[[length(rows) + 1]] <- data.frame(pair = paste(S, "->", T_type),
    ref_k0_15 = ref_K[1], ref_k0_30 = ref_K[2], ref_k0_60 = ref_K[3],
    ref_k0_120 = ref_K[4], ref_k0_240 = ref_K[5],
    ref_bb = ref_bb, bb_alpha = if (is.null(bb)) NA else bb$alpha,
    ref_distance = pairs$rate[[pi_]] * 0 +  # placeholder replaced below
      NA_real_,
    profile_drift_cor = drift_cor, stringsAsFactors = FALSE)
  # distance-based reference recomputed directly for comparison
  ref_d <- cellAdmixCore:::.celladmix_audit_reference_rate(
    pool_counts[z15], tt[z15], dist_S[z15])
  rows[[length(rows)]]$ref_distance <- ref_d$rate
}
res <- do.call(rbind, rows)
res_scaled <- res
num <- vapply(res, is.numeric, TRUE)
res_scaled[num] <- lapply(res[num], function(v) round(1e3 * v, 3))
res_scaled$bb_alpha <- round(res$bb_alpha, 3)
res_scaled$profile_drift_cor <- round(res$profile_drift_cor, 3)
write.csv(res_scaled, file.path(ag_dir, "results", "neighbor_reference.csv"),
  row.names = FALSE)

cat("lateral radius spanned by K neighbors (approx, um):",
  paste(Ks, "->", round(radius_K), collapse = ", "), "\n\n")
ok <- is.finite(res$ref_k0_240) & is.finite(res$ref_distance)
cat("pairs with usable K=240 reference:", sum(ok), "of", nrow(res), "\n")
cat("median ratio K=15 zero-neighbor reference / distance reference:",
  round(median(res$ref_k0_15 / res$ref_distance, na.rm = TRUE), 2), "\n")
cat("median ratio K=240 zero-neighbor reference / distance reference:",
  round(median(res$ref_k0_240[ok] / res$ref_distance[ok]), 2), "\n")
okb <- is.finite(res$ref_bb) & is.finite(res$ref_distance) & res$ref_distance > 0
cat("median ratio beta-binomial extrapolated / distance reference:",
  round(median(res$ref_bb[okb] / res$ref_distance[okb]), 2), "\n")
cat("median target-profile correlation near vs far (compartment check):",
  round(median(res$profile_drift_cor, na.rm = TRUE), 3), "\n")
print(head(res_scaled[order(-res$ref_k0_15), c("pair", "ref_k0_15", "ref_k0_60",
  "ref_k0_240", "ref_bb", "bb_alpha", "ref_distance", "profile_drift_cor")], 12),
  row.names = FALSE)
cat("NEIGHBOR REFERENCE DONE\n")
