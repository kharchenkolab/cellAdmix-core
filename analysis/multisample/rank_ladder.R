# Usage: Rscript rank_ladder.R <variant> <rank,rank,...> <out_csv> <name>=<bundle>=<annotation>=<output_dir> ...
# Fits every sample at every rank (cached fits reused) and summarizes the
# restart structure of each fit: objective spread, mean one-to-one stability,
# stable factors (>= 0.5), the share of the primary factors' molecules that
# another restart's factor cone cannot represent (nonnegative least squares
# on unit-sum profiles, median over the other restarts), and the factors
# more than 30% unrepresentable ("noise factors") with their molecule share.
# CELLADMIX_NMF_ITERATIONS overrides the NMF iteration count under its own run id.
args <- commandArgs(TRUE)
variant <- args[[1]]; ranks <- as.integer(strsplit(args[[2]], ",")[[1]]); out_csv <- args[[3]]; specs <- args[-(1:3)]
iterations <- as.integer(Sys.getenv("CELLADMIX_NMF_ITERATIONS", NA))
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
for (spec in specs) {
  p <- strsplit(spec, "=", fixed = TRUE)[[1]]
  ds <- cellAdmix(p[[2]], output_dir = p[[4]], annotation = read_labels(p[[3]]), num_threads = 6)
  for (k in ranks) {
    fit <- if (is.na(iterations)) {
      ds$fit(nmf_variant = variant, rank = k)
    } else {
      ds$fit(nmf_variant = variant, rank = k, nmf_iterations = iterations,
        run_id = sprintf("fit_manual_rank%d_%s_iter%d", k, variant, iterations))
    }
    rows[[length(rows) + 1]] <- cbind(sample = p[[1]], variant = variant,
      iterations = if (is.na(iterations)) as.integer(fit$run$pipeline_options$nmf_iterations) else iterations,
      describe(fit))
  }
}
out <- do.call(rbind, rows)
print(out, row.names = FALSE)
write.csv(out, out_csv, row.names = FALSE)
