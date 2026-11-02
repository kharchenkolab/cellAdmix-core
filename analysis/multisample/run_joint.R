# Joint stand-in: all samples in one molecule table (built by
# make_joint_table.py with a large coordinate offset between samples and a
# sample_id column) through the ordinary single-dataset pipeline: one
# factorization on rows pooled across samples, pooled scoring, one rule set,
# one ensemble correction, one audit. Per-sample views are then cut out by
# the sample prefix of the cell ids: per-sample correction summaries, and a
# per-sample audit and evaluation on the sample's own cells and counts.
#
# Usage: Rscript run_joint.R <molecules.parquet> <annotation.csv.gz> <out_dir> [threads]
.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/home/pkharchenko/cellAdmix/wt-dev_multisample/.Rlib"), .libPaths()))
suppressMessages(library(cellAdmixCore))
suppressMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: run_joint.R <molecules.parquet> <annotation.csv.gz> <out_dir> [threads]")
molecules_path <- args[[1]]
annotation_path <- args[[2]]
out_dir <- args[[3]]
threads <- if (length(args) >= 4) as.integer(args[[4]]) else 10L
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

ann <- read.csv(annotation_path, stringsAsFactors = FALSE)
annotation <- stats::setNames(as.character(ann$cell_type), as.character(ann$cell_id))
sample_of_cell <- stats::setNames(as.character(ann$sample), as.character(ann$cell_id))
samples <- sort(unique(ann$sample))

ds <- cellAdmix(molecules_path, output_dir = out_dir, format = "tabular",
  schema = celladmix_schema(x = "x", y = "y", z = "z", gene = "gene", cell = "cell",
    sample_id = "sample_id", qv = "qv"),
  annotation = annotation, num_threads = threads)
t0 <- Sys.time()
fit <- ds$fit(nmf_variant = "ls_nmf", nmf_n_runs = 10L, verbose = TRUE)
log_line("joint fit done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min; rank", fit$rank)
write.csv(fit$summary(), file.path(out_dir, "fit_summary.csv"), row.names = FALSE)
h <- fit$loadings()
colnames(h) <- fit$run$genes
write.csv(data.frame(factor = paste0("F", seq_len(nrow(h))), h, check.names = FALSE),
  file.path(out_dir, "profiles.csv"), row.names = FALSE)

t0 <- Sys.time()
score <- fit$score_bridge()
log_line("scoring done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")
rules <- score$rules()
write.csv(rules, file.path(out_dir, "rules.csv"), row.names = FALSE)
write.csv(score$summary(), file.path(out_dir, "score_summary.csv"), row.names = FALSE)
log_line("rules:", nrow(rules), "kept", sum(rules$keep))

t0 <- Sys.time()
correction <- score$correct(rules = rules)
log_line("correction done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")
write.csv(correction$rules, file.path(out_dir, "correction_rules.csv"), row.names = FALSE)

# Per-sample views by cell-id prefix.
cell_summary <- correction$cell_summary()
cell_summary$sample <- sub("_.*$", "", as.character(cell_summary$cell_id))
per_sample <- do.call(rbind, lapply(samples, function(s) {
  compact <- cellAdmixCore:::.celladmix_compact_correction_summary(
    cell_summary[cell_summary$sample == s, , drop = FALSE])
  cbind(sample = s, compact, stringsAsFactors = FALSE)
}))
write.csv(per_sample, file.path(out_dir, "correction_summary.csv"), row.names = FALSE)

t0 <- Sys.time()
counts_before <- fit$counts()
counts_after <- correction$counts()
cells <- fit$cell_factors()
cells$sample <- sub("_.*$", "", as.character(cells$cell_id))
evaluations <- list()
for (s in samples) {
  cells_s <- cells[cells$sample == s, , drop = FALSE]
  ids <- as.character(cells_s$cell_id)
  ann_s <- annotation[names(annotation) %in% ids]
  audit_s <- cellAdmixCore:::.celladmix_audit_from_data(
    counts = counts_before[, ids, drop = FALSE],
    cells = data.frame(cell_id = ids, x = cells_s$x, y = cells_s$y, stringsAsFactors = FALSE),
    annotation = ann_s)
  pairs_s <- audit_s$pairs()
  write.csv(cbind(sample = s, pairs_s, stringsAsFactors = FALSE),
    file.path(out_dir, paste0("audit_pairs_", s, ".csv")), row.names = FALSE)
  corrected_s <- counts_after[, ids, drop = FALSE]
  report_s <- audit_s$evaluate(list(counts = function() corrected_s))
  evaluations[[s]] <- cbind(sample = s, as.data.frame(report_s$summary()), stringsAsFactors = FALSE)
  write.csv(cbind(sample = s, report_s$pairs(), stringsAsFactors = FALSE),
    file.path(out_dir, paste0("audit_evaluation_pairs_", s, ".csv")), row.names = FALSE)
}
write.csv(do.call(rbind, evaluations), file.path(out_dir, "audit_evaluation.csv"), row.names = FALSE)
log_line("per-sample audits done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")
saveRDS(list(run_dir = fit$run_dir, correction_dir = correction$run$path, samples = samples),
  file.path(out_dir, "joint_pointers.rds"))
log_line("done")
