# Evaluate the generative-model correction arms against the audit's
# gene-split harness (eval_split from analysis/audit_guided/01_split_eval.R).
# Reads the removed-molecule matrices written by 01_model.py, applies them to
# the exported "before" counts, and writes per-pair and per-type CSVs plus a
# gate summary.
#
# Usage: Rscript 02_eval_gates.R [arm1 arm2 ...]
# Default arms: validation production production_noind shuffled
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(Matrix))
ag_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/audit_guided"
gm_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
source(file.path(ag_dir, "01_split_eval.R"))

arms <- commandArgs(trailingOnly = TRUE)
if (!length(arms)) arms <- c("validation", "production", "production_noind",
  "shuffled")
dataset <- Sys.getenv("GM_DATASET", "pancreas")
data_dir <- file.path(gm_dir,
  if (dataset == "pancreas") "data" else paste0("data_", dataset))
prefix <- if (dataset == "pancreas") "" else paste0(dataset, "_")

defs <- readRDS(file.path(data_dir, "defs.rds"))
counts0 <- defs$counts_before

summaries <- list()
for (arm in arms) {
  f <- file.path(data_dir, sprintf("gm_removed_%s.mtx", arm))
  if (!file.exists(f)) { message("missing: ", f, " - skipping"); next }
  message("== arm: ", arm, " ==")
  rem <- as(Matrix::readMM(f), "CsparseMatrix")
  dimnames(rem) <- dimnames(counts0)
  counts_after <- counts0 - rem
  counts_after@x <- pmax(counts_after@x, 0)
  message(sprintf("  removed %.0f of %.0f molecules (%.2f%%)",
    sum(rem), sum(counts0), 100 * sum(rem) / sum(counts0)))
  res <- eval_split(defs, counts_after, paste0("gm_", arm))
  print(res$summary, row.names = FALSE)
  write.csv(res$pairs,
    file.path(gm_dir, "results", sprintf("%sgm_%s_pairs.csv", prefix, arm)),
    row.names = FALSE)
  write.csv(res$types,
    file.path(gm_dir, "results", sprintf("%sgm_%s_types.csv", prefix, arm)),
    row.names = FALSE)
  summaries[[arm]] <- res$summary
  if (arm %in% c("production", "production_noind")) {
    low <- res$pairs[!is.na(res$pairs$power_B) & res$pairs$power_B < 0.8, ]
    message("  pairs below 0.8 power_B: ", nrow(low))
    if (nrow(low)) print(low[order(-low$excess_B),
      c("pair", "excess_B", "power_B", "power_strictB")], row.names = FALSE)
  }
}
if (length(summaries)) {
  write.csv(do.call(rbind, summaries),
    file.path(gm_dir, "results", paste0(prefix, "gm_gate_summaries.csv")), row.names = FALSE)
}
message("EVAL DONE")
