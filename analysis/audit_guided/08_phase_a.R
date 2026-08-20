# Phase A: exposure-regression corrector v3 (option 3, regression tier).
#
# Upgrades over the E4 prototype, each fixing a measured shortfall:
#  - isotonic (monotone) dose-response on raw exposure counts, replacing the
#    0/1/2/3+ bins whose "3+" pooling carried 65-95% of the excess;
#  - budgets delivered exactly at the exposure-stratum level (waterfilling
#    across genes within each stratum), immune to per-cell sparsity;
#  - ambient tier: strict genes (near-zero native baseline) are removed
#    entirely from target cells - their baseline is contamination too;
#  - top-up pass: re-measure each pair's gradient after applying, spend the
#    residual once more; stop when flat.
#
# Arms:
#  validation: budgets from A halves, ambient off  -> honest power_B
#  production: budgets from full pools, ambient on -> gate metrics
#  shuffled:   production machinery, permuted exposure -> control
# Plus a downstream cell-state kNN-purity comparison against the default
# membrane ensemble and the original counts.
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(cellAdmixCore))
ag_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/audit_guided"
source(file.path(ag_dir, "01_split_eval.R"))

setwd("/home/pkharchenko/cellAdmix/cellAdmix-core/examples/xenium_pancreas_membrane_377_full")
annotation <- read.csv("annotations/annotation.csv.gz", stringsAsFactors = FALSE)
cell_annotation <- setNames(annotation$merged_annotation, annotation$cell_id)
ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)
fit <- ds$fit(nmf_variant = "invsqrt_kl", verbose = FALSE)

defs <- build_split_defs(fit, cell_annotation)
counts0 <- defs$counts_before
profiles <- defs$profiles
top_type <- colnames(profiles)[max.col(profiles, ties.method = "first")]
totals0 <- defs$totals_before

# Weighted pool-adjacent-violators: monotone-increasing fit of y with weights w.
pava <- function(y, w) {
  n <- length(y)
  val <- y; wt <- w; idx <- as.list(seq_len(n))
  i <- 1
  while (i < length(val)) {
    if (val[i] > val[i + 1] + 1e-15) {
      merged_w <- wt[i] + wt[i + 1]
      merged_v <- (val[i] * wt[i] + val[i + 1] * wt[i + 1]) / merged_w
      val[i] <- merged_v; wt[i] <- merged_w; idx[[i]] <- c(idx[[i]], idx[[i + 1]])
      val <- val[-(i + 1)]; wt <- wt[-(i + 1)]; idx[[i + 1]] <- NULL
      i <- max(1, i - 1)
    } else i <- i + 1
  }
  out <- numeric(n)
  for (j in seq_along(val)) out[idx[[j]]] <- val[j]
  out
}

# Per-pair stratum budgets from an isotonic dose-response on guide genes,
# measured on the CURRENT counts (per-cell totals stay at their original
# values throughout, matching the harness convention).
stratum_budgets <- function(cur, d, guide, e) {
  strata <- sort(unique(e))
  mk <- ag_mcount(cur, guide, d$T_cells)
  tt <- totals0[d$T_cells]
  m_s <- tapply(mk, e, sum)[as.character(strata)]
  t_s <- tapply(tt, e, sum)[as.character(strata)]
  rate <- pava(as.numeric(m_s) / pmax(as.numeric(t_s), 1), as.numeric(t_s))
  r0 <- rate[strata == 0]
  if (!length(r0)) return(NULL)
  excess <- pmax(rate - r0, 0)
  data.frame(stratum = strata, excess_rate = excess, totals = as.numeric(t_s))
}

exposure_correct_v3 <- function(guide_half = c("A", "pool"), ambient = FALSE,
                                shuffle = FALSE, passes = 2L, seed = 1L) {
  guide_half <- match.arg(guide_half)
  set.seed(seed)
  cur <- counts0
  removed_ambient <- 0
  if (ambient) {
    # Strict genes: baseline in target cells is contamination (near-zero
    # native expression). Remove their content from target cells outright.
    for (nm in names(defs$pairs)) {
      d <- defs$pairs[[nm]]
      gs <- intersect(d$strict, rownames(cur))
      if (!length(gs)) next
      removed_ambient <- removed_ambient + sum(cur[gs, d$T_cells])
      cur[gs, d$T_cells] <- 0
    }
  }
  for (pass in seq_len(passes)) {
    total_budget <- 0; total_delivered <- 0
    for (nm in names(defs$pairs)) {
      d <- defs$pairs[[nm]]
      e <- d$e[d$T_cells]
      if (shuffle) e <- setNames(sample(as.numeric(e)), names(e))
      guide <- if (guide_half == "A") d$A else d$pool
      if (ambient) guide <- setdiff(guide, d$strict)
      if (length(guide) < 3) next
      sb <- stratum_budgets(cur, d, guide, e)
      if (is.null(sb) || sum(sb$excess_rate) <= 0) next
      gset <- rownames(profiles)[top_type == d$S]
      if (ambient) gset <- setdiff(gset, d$strict)
      gset <- intersect(gset, rownames(cur))
      if (!length(gset)) next
      share_guide <- sum(profiles[intersect(guide, rownames(profiles)), d$S]) /
        max(sum(profiles[, d$S]), 1)
      share_gset <- sum(profiles[gset, d$S]) / max(sum(profiles[, d$S]), 1)
      baseline_T0 <- bench_pseudobulk(counts0, d$T_cells[d$e[d$T_cells] == 0])
      psi <- profiles[gset, d$S]
      w <- psi * psi / (psi + pmax(baseline_T0[gset], 0))
      w <- w / max(sum(w), 1e-12)
      pair_mat <- as.matrix(cur[gset, d$T_cells, drop = FALSE])
      for (si in which(sb$excess_rate > 0)) {
        s <- sb$stratum[si]
        in_s <- e[d$T_cells] == s
        budget <- sb$excess_rate[si] * sb$totals[si] * share_gset / max(share_guide, 1e-6)
        total_budget <- total_budget + budget
        sub <- pair_mat[, in_s, drop = FALSE]
        avail_g <- rowSums(sub)
        rem_g <- numeric(length(gset)); resid <- budget
        for (round in 1:4) {
          if (resid <= 0.01 * budget) break
          hw <- w * (avail_g - rem_g > 0.5)
          if (sum(hw) <= 0) break
          step <- pmin(avail_g - rem_g, hw / sum(hw) * resid)
          rem_g <- rem_g + step
          resid <- resid - sum(step)
        }
        total_delivered <- total_delivered + sum(rem_g)
        frac <- ifelse(avail_g > 1e-12, pmin(1, rem_g / pmax(avail_g, 1e-12)), 0)
        pair_mat[, in_s] <- sub * (1 - frac)
      }
      cur[gset, d$T_cells] <- as(pair_mat, "CsparseMatrix")
    }
    message(sprintf("  pass %d: budget %.0f, delivered %.0f (%.0f%%)",
      pass, total_budget, total_delivered,
      100 * total_delivered / max(total_budget, 1)))
    if (total_delivered < 0.05 * sum(totals0) / 1e3) break
  }
  message(sprintf("  ambient tier removed: %.0f molecules", removed_ambient))
  cur
}

# Downstream diagnostic: same-label kNN purity in PCA space on a 5,000-cell
# subset, normalized against the label-frequency baseline.
knn_purity <- function(m, n_cells = 5000L, n_genes = 1000L, k = 15L, seed = 7L) {
  set.seed(seed)
  ct <- defs$cell_types[colnames(m)]
  keep <- names(ct)[!is.na(ct) & totals0[names(ct)] >= 30]
  cells <- sample(keep, min(n_cells, length(keep)))
  sub <- m[, cells, drop = FALSE]
  norm <- log1p(t(t(as.matrix(sub)) / pmax(Matrix::colSums(sub), 1)) * 1e4)
  vars <- apply(norm, 1, stats::var)
  norm <- norm[order(-vars)[seq_len(min(n_genes, nrow(norm)))], ]
  pc <- stats::prcomp(t(norm), rank. = 30, center = TRUE, scale. = FALSE)$x
  dm <- as.matrix(stats::dist(pc))
  diag(dm) <- Inf
  lab <- as.character(ct[cells])
  same <- vapply(seq_along(cells), function(i) {
    nn <- order(dm[i, ])[seq_len(k)]
    mean(lab[nn] == lab[i])
  }, 0)
  base <- sum((table(lab) / length(lab))^2)
  (mean(same) - base) / (1 - base)
}

message("== validation arm (A-guided, no ambient) ==")
cv <- exposure_correct_v3("A", ambient = FALSE)
rv <- eval_split(defs, cv, "v3_validation")
print(rv$summary, row.names = FALSE)
write.csv(rv$pairs, file.path(ag_dir, "results", "phaseA_validation_pairs.csv"),
  row.names = FALSE)

message("\n== production arm (pool-guided, ambient on) ==")
cp <- exposure_correct_v3("pool", ambient = TRUE)
rp <- eval_split(defs, cp, "v3_production")
print(rp$summary, row.names = FALSE)
write.csv(rp$pairs, file.path(ag_dir, "results", "phaseA_production_pairs.csv"),
  row.names = FALSE)
write.csv(rp$types, file.path(ag_dir, "results", "phaseA_production_types.csv"),
  row.names = FALSE)

message("\n== shuffled-exposure control ==")
cs <- exposure_correct_v3("pool", ambient = FALSE, shuffle = TRUE, passes = 1L)
rs <- eval_split(defs, cs, "v3_shuffled")
print(rs$summary, row.names = FALSE)

message("\n== coverage gate: pairs below 0.8 power_B (production) ==")
low <- rp$pairs[!is.na(rp$pairs$power_B) & rp$pairs$power_B < 0.8, ]
print(low[order(-low$excess_B), c("pair", "excess_B", "power_B", "power_strictB")],
  row.names = FALSE)

message("\n== downstream kNN purity ==")
score_m <- fit$score_membrane(verbose = FALSE)
corr_ens <- suppressWarnings(fit$correct(score_m, name = "ag_phaseA_ens"))
p0 <- knn_purity(counts0)
pe <- knn_purity(corr_ens$counts())
pp <- knn_purity(cp)
message(sprintf("original %.3f | membrane ensemble %.3f | exposure v3 %.3f",
  p0, pe, pp))
writeLines(sprintf("original,%.4f\nensemble,%.4f\nexposure_v3,%.4f", p0, pe, pp),
  file.path(ag_dir, "results", "phaseA_purity.csv"))
message("PHASE A DONE")
