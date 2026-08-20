# E4: minimal exposure-regression corrector (option 3's cheapest instance).
#
# Per pair (S -> T): the A-half dose-response sets a per-cell leaked-molecule
# budget (excess A-rate at the cell's exposure bin, extrapolated by the
# A-half's share of the S transcriptome). The budget is spent across S-owned
# genes weighted by source abundance times a baseline-aware likelihood ratio,
# capped at observed counts. Zero-exposure cells get budget zero, so the
# baseline is preserved by construction (depth ~ 0).
#
# Controls: shuffled exposure (budgets should collapse) is run alongside.
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
counts <- defs$counts_before
profiles <- defs$profiles
top_type <- colnames(profiles)[max.col(profiles, ties.method = "first")]

exposure_correct <- function(defs, shuffle = FALSE, seed = 1L) {
  set.seed(seed)
  removal <- counts * 0
  for (nm in names(defs$pairs)) {
    d <- defs$pairs[[nm]]
    bins <- d$bins
    if (shuffle) bins <- sample(bins)
    # A-half excess rate per bin (guide side only).
    rA <- bench_bin_rates(ag_mcount(counts, d$A, d$T_cells),
      defs$totals_before[d$T_cells], bins)
    r0 <- rA$rate[rA$bin == "0"]
    excess_rate <- setNames(pmax(rA$rate - r0, 0), rA$bin)
    # Extrapolation: A-half share of the source transcriptome.
    share_A <- sum(profiles[intersect(d$A, rownames(profiles)), d$S]) /
      max(sum(profiles[, d$S]), 1)
    # Removal gene set: S-owned genes; weights = source abundance x
    # baseline-aware likelihood ratio against the unexposed target profile.
    baseline_T0 <- bench_pseudobulk(counts, d$T_cells[d$bins == "0"])
    gset <- rownames(profiles)[top_type == d$S]
    psi <- profiles[gset, d$S]
    lr <- psi / (psi + pmax(baseline_T0[gset], 0))
    w <- psi * lr
    share_gset <- sum(psi) / max(sum(profiles[, d$S]), 1)
    w <- w / max(sum(w), 1e-12)
    # Per-cell budget on the removal set (in molecules).
    budget <- excess_rate[as.character(bins)] *
      defs$totals_before[d$T_cells] / max(share_A, 1e-6) * share_gset
    names(budget) <- d$T_cells
    budget <- budget[budget > 0]
    if (!length(budget)) next
    gset <- intersect(gset, rownames(counts))
    sub <- as.matrix(counts[gset, names(budget), drop = FALSE])
    # Expected removal per gene-cell, capped at observed counts; sparse
    # counts leave much of the budget undelivered, so redistribute the
    # residual budget over genes with remaining molecules (waterfilling).
    wv <- w[gset]
    rem <- matrix(0, nrow(sub), ncol(sub), dimnames = dimnames(sub))
    resid <- budget
    for (round in 1:4) {
      if (all(resid <= 0.01 * budget)) break
      head_w <- wv * (sub - rem > 0.5)
      head_w <- sweep(head_w, 2, pmax(colSums(head_w), 1e-12), "/")
      step <- pmin(sub - rem, sweep(head_w, 2, resid, "*"))
      rem <- rem + step
      resid <- pmax(resid - colSums(step), 0)
    }
    removal[gset, names(budget)] <- pmax(removal[gset, names(budget)],
      as(rem, "CsparseMatrix"))
  }
  out <- counts - removal
  out@x[out@x < 0] <- 0
  out
}

message("== exposure-regression corrector ==")
after <- exposure_correct(defs)
res <- eval_split(defs, after, "exposure_reg")
write.csv(res$pairs, file.path(ag_dir, "results", "e4_exposure_pairs.csv"),
  row.names = FALSE)
write.csv(res$types, file.path(ag_dir, "results", "e4_exposure_types.csv"),
  row.names = FALSE)
print(res$summary, row.names = FALSE)

message("\n== shuffled-exposure control ==")
after_sh <- exposure_correct(defs, shuffle = TRUE)
res_sh <- eval_split(defs, after_sh, "exposure_shuffled")
write.csv(res_sh$pairs, file.path(ag_dir, "results", "e4_shuffled_pairs.csv"),
  row.names = FALSE)
print(res_sh$summary, row.names = FALSE)

message("\n-- key pairs (exposure_reg): power vs depth --")
key <- res$pairs[order(-res$pairs$excess_B), ][1:12, ]
print(key[, c("pair", "excess_B", "power_B", "power_strictB", "depth_B")],
  row.names = FALSE)
message("E4 DONE")
