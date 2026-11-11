# Usage: Rscript restart_structure.R <variant> <rank> <name>=<bundle>=<annotation>=<output_dir> ...
# For each sample: the spread of the restart objectives, how much of each
# restart's factors another restart's factor cone can represent (nonnegative
# least squares on unit-sum profiles), the one-to-one ownership correlation,
# and the primary factors least representable in the other restarts with
# their molecule share, stability, and top genes. Fits are reused from the
# output directories (fitted if absent).
args <- commandArgs(TRUE)
variant <- args[[1]]; rank <- as.integer(args[[2]]); specs <- args[-(1:2)]
suppressPackageStartupMessages(library(cellAdmixCore))
read_labels <- function(path) {
  labels <- read.csv(path, stringsAsFactors = FALSE)
  labels <- labels[!is.na(labels$cell_type) & labels$cell_type != "Unlabeled", ]
  setNames(as.character(labels$cell_type), as.character(labels$cell_id))
}
datasets <- list()
for (spec in specs) {
  p <- strsplit(spec, "=", fixed = TRUE)[[1]]
  datasets[[p[[1]]]] <- cellAdmix(p[[2]], output_dir = p[[4]], annotation = read_labels(p[[3]]), num_threads = 2)
}
collection <- cellAdmixCollection(datasets)
collection$fit(nmf_variant = variant, rank = rank)
ns <- asNamespace("cellAdmixCore")
norm_rows <- function(h) h / pmax(rowSums(h), .Machine$double.eps)
project_w <- function(X, H, iters = 400L) {
  W <- matrix(1 / nrow(H), nrow(X), nrow(H)); XHt <- X %*% t(H); HHt <- H %*% t(H)
  for (i in seq_len(iters)) W <- W * XHt / pmax(W %*% HHt, 1e-12)
  W
}
unexplained <- function(B, A) {  # per row of B: residual share after NNLS onto rows of A
  W <- project_w(B, A); rowSums((B - W %*% A)^2) / rowSums(B^2)
}
importance <- function(fit) {
  cells <- fit$cell_factors()
  cols <- grep("^factor_[0-9]+_fraction$", names(cells), value = TRUE)
  w <- cells$transcript_count
  vapply(cols, function(col) sum(cells[[col]] * w) / sum(w), numeric(1))
}
cat("variant:", variant, " rank:", rank, "\n")
first <- TRUE
for (s in names(collection$fits)) {
  fit <- collection$fits[[s]]; run <- fit$run
  if (first) {
    cat("run fields:", paste(grep("nmf|objective|iter|converg|stab", names(run), value = TRUE), collapse = ", "), "\n")
    cat("nmf_diagnostics fields:", paste(names(run$nmf_diagnostics), collapse = ", "), "\n")
    opts <- run$pipeline_options
    cat("options:", paste(names(opts)[grepl("nmf|iter|tol", names(opts))], unlist(opts[grepl("nmf|iter|tol", names(opts))]), sep = "=", collapse = ", "), "\n")
    first <- FALSE
  }
  obj <- as.numeric(run$nmf_candidate_final_objectives)
  cat(sprintf("\n[%s] restart objectives: min %.6g, spread (max/min - 1) %.4f, values rel. to min: %s\n", s, min(obj),
    max(obj) / min(obj) - 1, paste(sprintf("%.4f", obj / min(obj) - 1), collapse = " ")))
  if (!is.null(run$nmf_diagnostics$candidate_iterations)) cat("  iterations:", paste(run$nmf_diagnostics$candidate_iterations, collapse = " "), "\n")
  sel <- as.integer(run$nmf_diagnostics$selected_run)
  pool <- lapply(ns$.celladmix_fit_member_profiles(fit), norm_rows)
  n <- length(pool)
  # pairwise representability of every restart's factors within another restart's cone
  pair <- matrix(NA_real_, n, n)
  for (a in seq_len(n)) for (b in seq_len(n)) if (a != b) pair[b, a] <- mean(unexplained(pool[[b]], pool[[a]]))
  cat(sprintf("  mean unexplained share of a restart's factors within another restart's cone: median %.3f, range %.3f to %.3f\n",
    median(pair, na.rm = TRUE), min(pair, na.rm = TRUE), max(pair, na.rm = TRUE)))
  # one-to-one matched ownership correlation for the same pairs (the stability metric's ingredient)
  corr_pairs <- sapply(setdiff(seq_len(n), sel), function(a) {
    m <- ns$.celladmix_match_factors(pool[[sel]], pool[[a]]); mean(m$correlation, na.rm = TRUE)
  })
  cat(sprintf("  primary restart vs others: mean one-to-one ownership correlation %.3f (recorded stability mean %.3f)\n",
    mean(corr_pairs), mean(as.numeric(run$nmf_factor_stability))))
  # primary restart factors: unexplained within each other restart, mass share, stability
  u <- sapply(setdiff(seq_len(n), sel), function(a) unexplained(pool[[sel]], pool[[a]]))
  imp <- importance(fit)
  st <- as.numeric(run$nmf_factor_stability)
  df <- data.frame(factor = seq_len(rank), mass = round(imp, 3), stability = round(st, 2),
    unexplained_median = round(apply(u, 1, median), 3), unexplained_max = round(apply(u, 1, max), 3))
  df$top_genes <- apply(pool[[sel]], 1, function(r) paste(names(sort(r, decreasing = TRUE))[1:4], collapse = "/"))
  df <- df[order(-df$unexplained_median), ]
  cat(sprintf("  mass-weighted unexplained share of the primary factors (median over other restarts): %.3f; unweighted %.3f\n",
    sum(imp * apply(u, 1, median)) / sum(imp), mean(apply(u, 1, median))))
  cat(sprintf("  mass share of primary factors with stability < 0.5: %.3f (%d factors); with unexplained > 0.3: %.3f (%d factors)\n",
    sum(imp[st < 0.5]) / sum(imp), sum(st < 0.5), sum(imp[apply(u, 1, median) > 0.3]) / sum(imp), sum(apply(u, 1, median) > 0.3)))
  print(head(df, 8), row.names = FALSE)
}
