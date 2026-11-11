# Usage: Rscript cleanup_vs_k.R <out_csv> <variant> <k,k,...>
# Breast rep1, rank 23 of the given variant: for each neighborhood size,
# fit (cached), bridge scoring, the default 10-member vote correction, and
# the audit's evaluation of the correction (sensitivity per detected pair,
# own-marker false removal). The audit is computed once from the variant's
# default fit. The least-squares runs live under breast/rep1, the KL runs
# under breast_membrane/rep1.
args <- commandArgs(TRUE)
out_csv <- args[[1]]; variant <- args[[2]]; ks <- as.integer(strsplit(args[[3]], ",")[[1]])
out_dir <- if (variant == "ls_nmf") "breast" else "breast_membrane"
.libPaths(c("/home/pkharchenko/cellAdmix/wt-dev_multisample/.Rlib", .libPaths()))
suppressPackageStartupMessages(library(cellAdmixCore))
read_labels <- function(path) {
  labels <- read.csv(path, stringsAsFactors = FALSE)
  labels <- labels[!is.na(labels$cell_type) & labels$cell_type != "Unlabeled", ]
  setNames(as.character(labels$cell_type), as.character(labels$cell_id))
}
ann <- read_labels("/home/pkharchenko/cellAdmix/wt-dev_multisample/examples/xenium_breast_replicates/annotations/rep1_annotation.csv.gz")
ds <- cellAdmix("/home/pkharchenko/cellAdmix/datasets/xenium_breast_rep1/outs",
  output_dir = file.path("/home/pkharchenko/cellAdmix/multisample_runs", out_dir, "rep1"), annotation = ann, num_threads = 8)
fit0 <- ds$fit(nmf_variant = variant, rank = 23)
t0 <- Sys.time(); audit <- fit0$audit_admixture(); message("audit: ", format(Sys.time() - t0))
rows <- list()
for (k in ks) {
  fit <- if (k == 20L) fit0 else ds$fit(nmf_variant = variant, rank = 23, ncv_k = k,
    run_id = sprintf("fit_manual_rank23_%s_k%d", variant, k))
  t0 <- Sys.time(); score <- fit$score_bridge(); message("k=", k, " scoring: ", format(Sys.time() - t0))
  rules <- score$rules(p_thresh = 0.1)
  t0 <- Sys.time(); correction <- score$correct(name = "bridge_clean"); message("k=", k, " correction: ", format(Sys.time() - t0))
  t0 <- Sys.time(); report <- audit$evaluate(correction, warn_uncovered = FALSE); message("k=", k, " evaluation: ", format(Sys.time() - t0))
  s <- report$summary(); p <- report$pairs()
  removed <- sum(fit$counts()) - sum(correction$counts())
  rows[[length(rows) + 1]] <- data.frame(variant = variant, ncv_k = as.integer(fit$run$pipeline_options$ncv_k),
    stability_mean = round(mean(as.numeric(fit$run$nmf_factor_stability)), 3),
    rules_kept = sum(rules$keep), rules_candidate = nrow(rules),
    molecules_removed = removed,
    detected_pairs = s$detected_pairs, estimated_admixed = round(s$estimated_admixed_molecules),
    leakage_removed_overall = round(s$leakage_removed_overall, 4),
    median_pair_sensitivity = round(s$median_pair_sensitivity, 4),
    pairs_below_0.2 = sum(p$sensitivity < 0.2, na.rm = TRUE),
    own_marker_false_removal = round(s$own_marker_false_removal, 4),
    worst_false_removal = s$worst_false_removal, stringsAsFactors = FALSE)
  out <- do.call(rbind, rows); print(out, row.names = FALSE); write.csv(out, out_csv, row.names = FALSE)
  write.csv(p, sub("\\.csv$", sprintf("_pairs_k%d.csv", k), out_csv), row.names = FALSE)
}
