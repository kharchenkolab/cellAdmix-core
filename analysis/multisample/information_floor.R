# How many molecules a neighborhood row needs to tell two programs apart:
# for the selected restart's unit-sum gene profiles, the symmetric
# Bhattacharyya distance per molecule between every pair of programs (the
# error exponent of telling them apart), and the row size at which the
# error bound reaches exp(-3) (k = 3 / distance), summarized over pairs (median and the 90th
# percentile, the hardest pairs). Molecule-weighted pairs count each program
# by its molecule share.
.libPaths(c("/home/pkharchenko/cellAdmix/wt-dev_multisample/.Rlib", .libPaths()))
suppressPackageStartupMessages(library(cellAdmixCore))
read_labels <- function(path) {
  labels <- read.csv(path, stringsAsFactors = FALSE)
  labels <- labels[!is.na(labels$cell_type) & labels$cell_type != "Unlabeled", ]
  setNames(as.character(labels$cell_type), as.character(labels$cell_id))
}
D <- "/home/pkharchenko/cellAdmix/datasets"; M <- "/home/pkharchenko/cellAdmix/multisample_runs"
A <- "/home/pkharchenko/cellAdmix/wt-dev_multisample/examples/xenium_breast_replicates/annotations"
ann5k <- read.csv("/home/pkharchenko/cellAdmix/cellAdmix-core/examples/xenium_breast_membrane_5k_full/annotations/annotation.csv.gz", stringsAsFactors = FALSE)
ann5k <- setNames(ann5k$merged_annotation, ann5k$cell_id); ann5k <- ann5k[ann5k != "Ambiguous / low-quality"]
specs <- list(
  list(name = "breast rep1 (321 genes), ls_nmf rank 23", ds = function() cellAdmix(file.path(D, "xenium_breast_rep1/outs"), output_dir = file.path(M, "breast/rep1"), annotation = read_labels(file.path(A, "rep1_annotation.csv.gz")), num_threads = 2), fit = function(ds) ds$fit(nmf_variant = "ls_nmf", rank = 23)),
  list(name = "breast rep1 (321 genes), invsqrt_kl rank 23, k = 50", ds = function() cellAdmix(file.path(D, "xenium_breast_rep1/outs"), output_dir = file.path(M, "breast_membrane/rep1"), annotation = read_labels(file.path(A, "rep1_annotation.csv.gz")), num_threads = 2), fit = function(ds) ds$fit(nmf_variant = "invsqrt_kl", rank = 23, ncv_k = 50L, run_id = "fit_manual_rank23_invsqrt_kl_k50")),
  list(name = "pancreas left (377 genes), ls_nmf rank 9", ds = function() cellAdmix(file.path(D, "pancreas_halves/left"), output_dir = file.path(M, "pancreas/left"), annotation = read_labels(file.path(D, "pancreas_halves/left_annotation.csv.gz")), num_threads = 2), fit = function(ds) ds$fit(nmf_variant = "ls_nmf", rank = 9)),
  list(name = "breast 5k crop (5101 genes), ls_nmf rank 8, k = 71", ds = function() cellAdmix("/home/pkharchenko/baysor/rewrite/Baysor/examples/Xenium_breast_memrane_5k/data", output_dir = file.path(M, "breast_5k_crop"), annotation = ann5k, num_threads = 2, analysis_bbox = c(5000, 7000, 4500, 6500)), fit = function(ds) ds$fit(nmf_variant = "ls_nmf", rank = 8, ncv_k = 71L, run_id = "fit_rank8_ls_nmf_k71")),
  list(name = "breast 5k crop (5101 genes), ls_nmf rank 16, k = 71", ds = function() cellAdmix("/home/pkharchenko/baysor/rewrite/Baysor/examples/Xenium_breast_memrane_5k/data", output_dir = file.path(M, "breast_5k_crop"), annotation = ann5k, num_threads = 2, analysis_bbox = c(5000, 7000, 4500, 6500)), fit = function(ds) ds$fit(nmf_variant = "ls_nmf", rank = 16, ncv_k = 71L, run_id = "fit_rank16_ls_nmf_k71")))
for (s in specs) {
  ds <- s$ds(); fit <- suppressMessages(s$fit(ds))
  h <- fit$loadings(); p <- h / rowSums(h)
  cells <- fit$cell_factors(); cols <- grep("^factor_[0-9]+_fraction$", names(cells), value = TRUE)
  share <- vapply(cols, function(col) sum(cells[[col]] * cells$transcript_count), numeric(1)); share <- share / sum(share)
  n <- nrow(p); div <- matrix(NA_real_, n, n); w <- matrix(NA_real_, n, n)
  for (a in seq_len(n)) for (b in seq_len(n)) if (a < b) {
    div[a, b] <- -log(sum(sqrt(p[a, ] * p[b, ]))); w[a, b] <- share[a] * share[b]
  }
  d <- div[upper.tri(div)]; ww <- w[upper.tri(w)]
  k_needed <- 3 / d
  wq <- function(x, wts, q) { o <- order(x); x <- x[o]; wts <- wts[o]; x[which(cumsum(wts) / sum(wts) >= q)[1]] }
  cat(sprintf("%s: Bhattacharyya distance per molecule median %.2f (90th pct of hardest pairs %.2f); molecules for 3 nats: median %.0f, hardest 10%% of pairs %.0f, molecule-weighted median %.0f\n",
    s$name, median(d), quantile(d, 0.1), median(k_needed), quantile(k_needed, 0.9), wq(k_needed, ww, 0.5)))
}
