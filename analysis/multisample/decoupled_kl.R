# Usage: Rscript decoupled_kl.R <out_csv> <train_k> <label_k>
# Breast rep1, invsqrt_kl rank 23: the programs come from the fit with
# <train_k>-molecule rows (its selected restart as the fixed profiles, all
# its restarts as the member pool); molecules are labeled with
# <label_k>-molecule rows projected onto those profiles; then bridge
# scoring, the 10-member vote correction and the audit's evaluation, in the
# columns of cleanup_vs_k.R.
args <- commandArgs(TRUE)
out_csv <- args[[1]]; train_k <- as.integer(args[[2]]); label_k <- as.integer(args[[3]])
.libPaths(c("/home/pkharchenko/cellAdmix/wt-dev_multisample/.Rlib", .libPaths()))
suppressPackageStartupMessages(library(cellAdmixCore))
ns <- asNamespace("cellAdmixCore")
read_labels <- function(path) {
  labels <- read.csv(path, stringsAsFactors = FALSE)
  labels <- labels[!is.na(labels$cell_type) & labels$cell_type != "Unlabeled", ]
  setNames(as.character(labels$cell_type), as.character(labels$cell_id))
}
ann <- read_labels("/home/pkharchenko/cellAdmix/wt-dev_multisample/examples/xenium_breast_replicates/annotations/rep1_annotation.csv.gz")
ds <- cellAdmix("/home/pkharchenko/cellAdmix/datasets/xenium_breast_rep1/outs",
  output_dir = "/home/pkharchenko/cellAdmix/multisample_runs/breast_membrane/rep1", annotation = ann, num_threads = 8)
trained <- if (train_k == 20L) ds$fit(nmf_variant = "invsqrt_kl", rank = 23) else
  ds$fit(nmf_variant = "invsqrt_kl", rank = 23, ncv_k = train_k, run_id = sprintf("fit_manual_rank23_invsqrt_kl_k%d", train_k))
pool <- ns$.celladmix_fit_member_profiles(trained)
sel <- as.integer(trained$run$nmf_diagnostics$selected_run)
t0 <- Sys.time()
fit <- ds$fit(nmf_fixed_h = pool[[sel]], nmf_fixed_h_members = pool, nmf_variant = "invsqrt_kl", nmf_n_runs = 1L,
  ncv_k = label_k, run_id = sprintf("fit_rank23_invsqrt_kl_train%d_label%d", train_k, label_k))
message("fixed-profile labeling: ", format(Sys.time() - t0))
t0 <- Sys.time(); audit <- trained$audit_admixture(); message("audit: ", format(Sys.time() - t0))
t0 <- Sys.time(); score <- fit$score_bridge(); message("scoring: ", format(Sys.time() - t0))
rules <- score$rules(p_thresh = 0.1)
t0 <- Sys.time(); correction <- score$correct(name = "bridge_clean"); message("correction: ", format(Sys.time() - t0))
t0 <- Sys.time(); report <- audit$evaluate(correction, warn_uncovered = FALSE); message("evaluation: ", format(Sys.time() - t0))
s <- report$summary(); p <- report$pairs()
removed <- sum(fit$counts()) - sum(correction$counts())
out <- data.frame(variant = sprintf("invsqrt_kl, programs at %d, labels at %d", train_k, label_k), ncv_k = label_k,
  stability_mean = round(mean(as.numeric(trained$run$nmf_factor_stability)), 3),
  rules_kept = sum(rules$keep), rules_candidate = nrow(rules), molecules_removed = removed,
  detected_pairs = s$detected_pairs, estimated_admixed = round(s$estimated_admixed_molecules),
  leakage_removed_overall = round(s$leakage_removed_overall, 4),
  median_pair_sensitivity = round(s$median_pair_sensitivity, 4),
  pairs_below_0.2 = sum(p$sensitivity < 0.2, na.rm = TRUE),
  own_marker_false_removal = round(s$own_marker_false_removal, 4),
  worst_false_removal = s$worst_false_removal, stringsAsFactors = FALSE)
print(out, row.names = FALSE); write.csv(out, out_csv, row.names = FALSE)
write.csv(p, sub("\\.csv$", "_pairs.csv", out_csv), row.names = FALSE)
