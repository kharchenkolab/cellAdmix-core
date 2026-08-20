# E3: audit-calibrated vote threshold. The sweep is scored on the A halves
# (selection side) and the chosen threshold is reported on the B halves,
# with own-marker false removal as the selection constraint.
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(cellAdmixCore))
ag_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/audit_guided"
source(file.path(ag_dir, "01_split_eval.R"))

setwd("/home/pkharchenko/cellAdmix/cellAdmix-core/examples/xenium_pancreas_membrane_377_full")
annotation <- read.csv("annotations/annotation.csv.gz", stringsAsFactors = FALSE)
cell_annotation <- setNames(annotation$merged_annotation, annotation$cell_id)
ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)
fit <- ds$fit(nmf_variant = "invsqrt_kl", verbose = FALSE)
score <- fit$score_membrane(verbose = FALSE)
defs <- build_split_defs(fit, cell_annotation)

budget_false_removal <- 0.05
votes <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0)
sweep <- lapply(votes, function(v) {
  corr <- suppressWarnings(fit$correct(score, vote = v,
    name = sprintf("ag_e3_v%02.0f", 100 * v)))
  res <- eval_split(defs, corr$counts(), sprintf("vote=%.2f", v))
  wb <- function(p, ex) { ok <- !is.na(p); sum(p[ok] * ex[ok]) / sum(ex[ok]) }
  data.frame(vote = v,
    power_A = wb(res$pairs$power_A, res$pairs$excess_pool),
    power_B = wb(res$pairs$power_B, res$pairs$excess_B),
    depth_B = stats::median(res$pairs$depth_B, na.rm = TRUE),
    false_removal = res$summary$false_removal_pooled)
})
sweep <- do.call(rbind, sweep)
print(sweep, row.names = FALSE)
write.csv(sweep, file.path(ag_dir, "results", "e3_vote_sweep.csv"), row.names = FALSE)

ok <- sweep$false_removal <= budget_false_removal
chosen <- if (any(ok)) sweep$vote[ok][which.max(sweep$power_A[ok])] else
  sweep$vote[which.min(sweep$false_removal)]
message(sprintf("selected on A under %.0f%% budget: vote=%.2f -> reported power_B=%.3f (false removal %.3f)",
  100 * budget_false_removal, chosen, sweep$power_B[sweep$vote == chosen],
  sweep$false_removal[sweep$vote == chosen]))
message("E3 DONE")
