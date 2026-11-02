# Separate per-sample processing (today's pipeline), one sample per call.
#
# Builds the input store, fits (ls_nmf, rank from the annotation, 10
# restarts), scores with the bridge method, derives rules with the
# native-factor check, applies the default ensemble correction, and
# audits before and after. Writes compact CSV summaries next to the run.
#
# Usage: Rscript run_separate.R <sample> <bundle_dir> <annotation.csv[.gz]>
#          <output_dir> [threads] [nmf_variant] [score_method]
#
# The annotation file has columns cell_id and cell_type; rows labeled
# "Unlabeled" are dropped (those cells stay in the store as unannotated).
.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/home/pkharchenko/cellAdmix/wt-dev_multisample/.Rlib"), .libPaths()))
suppressMessages(library(cellAdmixCore))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: run_separate.R <sample> <bundle_dir> <annotation> <output_dir> [threads] [nmf_variant] [score_method]")
sample <- args[[1]]
bundle_dir <- args[[2]]
annotation_path <- args[[3]]
output_dir <- args[[4]]
threads <- if (length(args) >= 5) as.integer(args[[5]]) else 10L
nmf_variant <- if (length(args) >= 6) args[[6]] else "ls_nmf"
score_method <- if (length(args) >= 7) args[[7]] else "bridge"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_annotation <- function(path) {
  ann <- read.csv(path, stringsAsFactors = FALSE)
  ann <- ann[!is.na(ann$cell_type) & ann$cell_type != "Unlabeled", ]
  stats::setNames(as.character(ann$cell_type), as.character(ann$cell_id))
}

log_line <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sample, ...,
  "\n", sep = " ")

annotation <- read_annotation(annotation_path)
log_line("annotation:", length(annotation), "cells,", length(unique(annotation)), "types")

ds <- cellAdmix(bundle_dir, output_dir = output_dir, annotation = annotation,
  num_threads = threads)
t0 <- Sys.time()
fit <- ds$fit(nmf_variant = nmf_variant, nmf_n_runs = 10L, verbose = TRUE)
log_line("fit done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min; rank", fit$rank)
write.csv(fit$summary(), file.path(output_dir, "fit_summary.csv"), row.names = FALSE)

t0 <- Sys.time()
score <- if (identical(score_method, "membrane")) fit$score_membrane() else fit$score_bridge()
log_line("scoring done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")
rules <- score$rules()
write.csv(rules, file.path(output_dir, "rules.csv"), row.names = FALSE)
write.csv(score$summary(), file.path(output_dir, "score_summary.csv"), row.names = FALSE)
log_line("rules:", nrow(rules), "of which kept", sum(rules$keep))

t0 <- Sys.time()
correction <- score$correct(rules = rules)
log_line("correction done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")
write.csv(correction$summary(), file.path(output_dir, "correction_summary.csv"), row.names = FALSE)
write.csv(correction$rules, file.path(output_dir, "correction_rules.csv"), row.names = FALSE)

t0 <- Sys.time()
audit <- fit$audit_admixture()
write.csv(audit$pairs(), file.path(output_dir, "audit_pairs.csv"), row.names = FALSE)
report <- audit$evaluate(correction)
write.csv(report$summary(), file.path(output_dir, "audit_evaluation.csv"), row.names = FALSE)
write.csv(report$pairs(), file.path(output_dir, "audit_evaluation_pairs.csv"), row.names = FALSE)
log_line("audit done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")
saveRDS(list(run_dir = fit$run_dir, correction_dir = correction$run$path,
  rank = fit$rank, nmf_variant = nmf_variant, score_method = score_method),
  file.path(output_dir, "separate_pointers.rds"))
log_line("done")
