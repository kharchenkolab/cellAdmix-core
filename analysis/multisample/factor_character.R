# Usage: Rscript factor_character.R <out_prefix> <rank,rank,...> <name>=<variant>=<bundle>=<annotation>=<output_dir> ...
# Characterizes the selected restart of each cached fit: gene concentration
# of every factor (top-gene share and effective gene number of the unit-sum
# loading row), the factor's molecule share, its one-to-one stability, the
# share of it that other restarts cannot represent, and its cell-type
# composition (molecule-weighted purity). Writes <out_prefix>_factors.csv and
# <out_prefix>_fits.csv.
args <- commandArgs(TRUE)
out_prefix <- args[[1]]; ranks <- as.integer(strsplit(args[[2]], ",")[[1]]); specs <- args[-(1:2)]
`%||%` <- function(a, b) if (is.null(a)) b else a
overrides <- list(); tag <- ""
iterations <- as.integer(Sys.getenv("CELLADMIX_NMF_ITERATIONS", NA))
if (!is.na(iterations)) { overrides$nmf_iterations <- iterations; tag <- paste0(tag, "_iter", iterations) }
ncv_k <- as.integer(Sys.getenv("CELLADMIX_NCV_K", NA))
if (!is.na(ncv_k)) { overrides$ncv_k <- ncv_k; tag <- paste0(tag, "_k", ncv_k) }
train_rows <- as.integer(Sys.getenv("CELLADMIX_TRAIN_ROWS", NA))
if (!is.na(train_rows)) { overrides$nmf_train_max_rows <- train_rows; tag <- paste0(tag, "_rows", train_rows) }
suppressPackageStartupMessages(library(cellAdmixCore))
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
effective_number <- function(p) { p <- p[p > 0]; exp(-sum(p * log(p))) }
factor_rows <- list(); fit_rows <- list()
for (spec in specs) {
  p <- strsplit(spec, "=", fixed = TRUE)[[1]]
  r <- p[[1]]; variant <- p[[2]]
  ann <- read_labels(p[[4]])
  ds <- cellAdmix(p[[3]], output_dir = p[[5]], annotation = ann, num_threads = 2)
  for (k in ranks) {
    fit <- tryCatch(suppressMessages(if (nzchar(tag)) {
      do.call(ds$fit, c(list(nmf_variant = variant, rank = k,
        run_id = sprintf("fit_manual_rank%d_%s%s", k, variant, tag)), overrides))
    } else {
      ds$fit(nmf_variant = variant, rank = k)
    }), error = function(e) NULL)
    if (is.null(fit)) { cat("missing:", r, variant, k, "\n"); next }
    run <- fit$run
    h <- fit$loadings(); colnames(h) <- run$genes
    pr <- norm_rows(h)
    top_gene <- colnames(pr)[apply(pr, 1, which.max)]
    top_share <- apply(pr, 1, max)
    eff_genes <- apply(pr, 1, effective_number)
    ownership <- t(t(h) / pmax(colSums(h), .Machine$double.eps))
    top_owned <- apply(ownership, 1, function(x) paste(names(sort(x, decreasing = TRUE))[1:5], collapse = " "))
    cells <- fit$cell_factors()
    cells$cell_type <- ann[as.character(cells$cell_id)]
    cells <- cells[!is.na(cells$cell_type), ]
    cols <- grep("^factor_[0-9]+_fraction$", names(cells), value = TRUE)
    mass <- as.matrix(cells[, cols]) * cells$transcript_count
    by_type <- rowsum(mass, cells$cell_type)
    share <- colSums(by_type) / sum(by_type)
    purity <- apply(by_type, 2, function(x) max(x) / sum(x))
    top_type <- rownames(by_type)[apply(by_type, 2, which.max)]
    eff_types <- apply(by_type, 2, function(x) effective_number(x / sum(x)))
    home <- apply(by_type, 1, function(x) max(x) / sum(x))
    sel <- as.integer(run$nmf_diagnostics$selected_run)
    pool <- lapply(ns$.celladmix_fit_member_profiles(fit), norm_rows)
    others <- setdiff(seq_along(pool), sel)
    u <- apply(sapply(others, function(a) unexplained(pool[[sel]], pool[[a]])), 1, median)
    st <- as.numeric(run$nmf_factor_stability)
    factor_rows[[length(factor_rows) + 1]] <- data.frame(sample = r, variant = paste0(variant, tag), rank = k,
      factor = seq_len(nrow(h)), molecule_share = round(share, 3), stability = round(st, 2),
      unrepresentable = round(u, 2), top_gene = top_gene, top_gene_share = round(top_share, 2),
      effective_genes = round(eff_genes, 1), top_cell_type = top_type, purity = round(purity, 2),
      effective_cell_types = round(eff_types, 1), top_owned_genes = top_owned, stringsAsFactors = FALSE)
    fit_rows[[length(fit_rows) + 1]] <- data.frame(sample = r, variant = paste0(variant, tag), rank = k,
      single_gene_factors = sum(top_share > 0.5), single_gene_mass = round(sum(share[top_share > 0.5]), 3),
      median_top_gene_share = round(median(top_share), 2), median_effective_genes = round(median(eff_genes), 1),
      mass_weighted_effective_genes = round(sum(share * eff_genes), 1),
      median_purity = round(median(purity), 2), mass_weighted_purity = round(sum(share * purity), 2),
      cell_types_with_home = sum(home >= 0.5), cell_types = length(home),
      median_home_share = round(median(home), 2),
      noise_factors = sum(u > 0.3), noise_mass = round(sum(share[u > 0.3]), 3),
      stringsAsFactors = FALSE)
  }
}
factors <- do.call(rbind, factor_rows); fits <- do.call(rbind, fit_rows)
write.csv(factors, paste0(out_prefix, "_factors.csv"), row.names = FALSE)
write.csv(fits, paste0(out_prefix, "_fits.csv"), row.names = FALSE)
print(fits, row.names = FALSE)
