# Rank and neighborhood ladder on the 5,101-gene breast crop (same crop as the
# cleanup benchmark): both variants at ranks 8 and 16 with 20, 35, 71 (the
# automatic value) and 140 molecules per neighborhood row.
.libPaths(c("/home/pkharchenko/cellAdmix/wt-dev_multisample/.Rlib", .libPaths()))
suppressPackageStartupMessages(library(cellAdmixCore))
ns <- asNamespace("cellAdmixCore")
out_csv <- commandArgs(TRUE)[[1]]
ann_df <- read.csv("/home/pkharchenko/cellAdmix/cellAdmix-core/examples/xenium_breast_membrane_5k_full/annotations/annotation.csv.gz",
  stringsAsFactors = FALSE)
ann <- setNames(ann_df$merged_annotation, ann_df$cell_id)
ann <- ann[ann != "Ambiguous / low-quality"]
ds <- cellAdmix("/home/pkharchenko/baysor/rewrite/Baysor/examples/Xenium_breast_memrane_5k/data",
  output_dir = "/home/pkharchenko/cellAdmix/multisample_runs/breast_5k_crop", annotation = ann,
  num_threads = 8, analysis_bbox = c(5000, 7000, 4500, 6500))
norm_rows <- function(h) h / pmax(rowSums(h), .Machine$double.eps)
project_w <- function(X, H, iters = 400L) {
  W <- matrix(1 / nrow(H), nrow(X), nrow(H)); XHt <- X %*% t(H); HHt <- H %*% t(H)
  for (i in seq_len(iters)) W <- W * XHt / pmax(W %*% HHt, 1e-12)
  W
}
unexplained <- function(B, A) { W <- project_w(B, A); rowSums((B - W %*% A)^2) / rowSums(B^2) }
importance <- function(fit) {
  cells <- fit$cell_factors()
  cols <- grep("^factor_[0-9]+_fraction$", names(cells), value = TRUE)
  w <- cells$transcript_count
  vapply(cols, function(col) sum(cells[[col]] * w) / sum(w), numeric(1))
}
describe <- function(fit) {
  run <- fit$run
  obj <- as.numeric(run$nmf_candidate_final_objectives)
  sel <- as.integer(run$nmf_diagnostics$selected_run)
  pool <- lapply(ns$.celladmix_fit_member_profiles(fit), norm_rows)
  others <- setdiff(seq_along(pool), sel)
  u <- apply(sapply(others, function(a) unexplained(pool[[sel]], pool[[a]])), 1, median)
  imp <- importance(fit); st <- as.numeric(run$nmf_factor_stability)
  data.frame(rank = nrow(pool[[sel]]), restarts = length(pool),
    objective_spread = round(max(obj) / min(obj) - 1, 4),
    stability_mean = round(mean(st), 3), stable_factors = sum(st >= 0.5),
    unexplained_mass_weighted = round(sum(imp * u) / sum(imp), 3),
    unexplained_mean = round(mean(u), 3),
    noise_factors = sum(u > 0.3), noise_mass = round(sum(imp[u > 0.3]), 3),
    stringsAsFactors = FALSE)
}
rows <- list()
for (variant in c("ls_nmf", "invsqrt_kl")) for (k in c(20L, 35L, 71L, 140L)) for (rank in c(8L, 16L)) {
  fit <- ds$fit(nmf_variant = variant, rank = rank, ncv_k = k,
    run_id = sprintf("fit_rank%d_%s_k%d", rank, variant, k))
  rows[[length(rows) + 1]] <- cbind(sample = "breast 5k crop", variant = variant, ncv_k = k,
    training_rows = as.integer(fit$run$n_training_rows), describe(fit))
  out <- do.call(rbind, rows); print(out, row.names = FALSE); write.csv(out, out_csv, row.names = FALSE)
}
