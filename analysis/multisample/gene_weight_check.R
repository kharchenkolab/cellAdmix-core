# Usage: Rscript gene_weight_check.R <out_csv> <rank> <name>=<variant>=<bundle>=<annotation>=<output_dir> ...
# For every factor of each cached fit: the loading-weighted mean of
# log10 gene total (how high-count the genes a factor loads on are), the
# share of its loading on the bottom-quartile genes by total, next to its
# molecule share, one-to-one stability and unrepresentable share.
args <- commandArgs(TRUE)
out_csv <- args[[1]]; rank <- as.integer(args[[2]]); specs <- args[-(1:2)]
suppressPackageStartupMessages({library(cellAdmixCore); library(Matrix)})
ns <- asNamespace("cellAdmixCore")
read_labels <- function(path) {
  labels <- read.csv(path, stringsAsFactors = FALSE)
  labels <- labels[!is.na(labels$cell_type) & labels$cell_type != "Unlabeled", ]
  setNames(as.character(labels$cell_type), as.character(labels$cell_id))
}
norm_rows <- function(h) h / pmax(rowSums(h), .Machine$double.eps)
project_w <- function(X, H, iters = 400L) {
  W <- matrix(1 / nrow(H), nrow(X), nrow(H)); XHt <- X %*% t(H); HHt <- H %*% t(H)
  for (i in seq_len(iters)) W <- W * XHt / pmax(W %*% HHt, 1e-12)
  W
}
unexplained <- function(B, A) { W <- project_w(B, A); rowSums((B - W %*% A)^2) / rowSums(B^2) }
rows <- list()
for (spec in specs) {
  p <- strsplit(spec, "=", fixed = TRUE)[[1]]
  ds <- cellAdmix(p[[3]], output_dir = p[[5]], annotation = read_labels(p[[4]]), num_threads = 2)
  fit <- suppressMessages(ds$fit(nmf_variant = p[[2]], rank = rank))
  run <- fit$run
  counts <- fit$counts()
  totals <- Matrix::rowSums(counts); names(totals) <- rownames(counts)
  h <- fit$loadings(); colnames(h) <- run$genes
  totals <- totals[colnames(h)]
  pr <- norm_rows(h)
  low <- totals <= quantile(totals, 0.25)
  cells <- fit$cell_factors()
  cols <- grep("^factor_[0-9]+_fraction$", names(cells), value = TRUE)
  share <- vapply(cols, function(col) sum(cells[[col]] * cells$transcript_count), numeric(1)); share <- share / sum(share)
  sel <- as.integer(run$nmf_diagnostics$selected_run)
  pool <- lapply(ns$.celladmix_fit_member_profiles(fit), norm_rows)
  others <- setdiff(seq_along(pool), sel)
  u <- apply(sapply(others, function(a) unexplained(pool[[sel]], pool[[a]])), 1, median)
  rows[[length(rows) + 1]] <- data.frame(sample = p[[1]], variant = p[[2]], rank = rank, factor = seq_len(nrow(h)),
    molecule_share = round(share, 3), stability = round(as.numeric(run$nmf_factor_stability), 2),
    unrepresentable = round(u, 2), mean_log10_gene_total = round(as.numeric(pr %*% log10(totals)), 2),
    low_count_loading = round(as.numeric(pr %*% as.numeric(low)), 2), stringsAsFactors = FALSE)
}
out <- do.call(rbind, rows)
write.csv(out, out_csv, row.names = FALSE)
for (key in unique(paste(out$sample, out$variant))) {
  d <- out[paste(out$sample, out$variant) == key, ]
  noise <- d$unrepresentable > 0.3
  cat(sprintf("%s: genes median log10 total %.2f; program factors mean log10 total %.2f (n=%d), noise factors %.2f (n=%d); low-count loading program %.2f, noise %.2f; cor(unrep, mean log total) %.2f\n",
    key, median(log10(totals)), mean(d$mean_log10_gene_total[!noise]), sum(!noise), mean(d$mean_log10_gene_total[noise]), sum(noise),
    mean(d$low_count_loading[!noise]), mean(d$low_count_loading[noise]), cor(d$unrepresentable, d$mean_log10_gene_total)))
}
