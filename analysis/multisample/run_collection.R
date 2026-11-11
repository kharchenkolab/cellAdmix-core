# Multi-sample arms on the per-sample datasets prepared by run_separate.R.
#
# Rebuilds each sample's dataset on the same output directory (so the cached
# per-sample fits are reused), then runs one of:
#   consensus              consensus factorization (matched, averaged profiles)
#   reference <sample>     reference-sample projection: one sample's profiles
#                          and restart pool projected onto every sample
#   profiles <run_dir>     projection of a persisted run's profiles (e.g. the
#                          joint stand-in's) with that run's restart pool
# followed by scoring (bridge by default; CELLADMIX_SCORE=membrane selects
# membrane scoring, CELLADMIX_VARIANT the per-sample NMF variant), pooled
# rules, per-sample ensemble correction (CELLADMIX_APPLY_RULES=where_present
# skips, in each sample, the rules on programs its own restarts do not
# carry), and per-sample audits. Writes CSV summaries under <out_dir>.
#
# Usage: Rscript run_collection.R <out_dir> <threads> <mode> [<mode arg>]
#          <sample>=<bundle>=<annotation>=<separate_dir> ...
.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/home/pkharchenko/cellAdmix/wt-dev_multisample/.Rlib"), .libPaths()))
suppressMessages(library(cellAdmixCore))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: run_collection.R <out_dir> <threads> <mode> [<mode arg>] <sample>=<bundle>=<annotation>=<separate_dir> ...")
out_dir <- args[[1]]
threads <- as.integer(args[[2]])
mode <- args[[3]]
rest <- args[-(1:3)]
mode_arg <- NULL
if (mode %in% c("reference", "profiles")) {
  mode_arg <- rest[[1]]
  rest <- rest[-1]
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

read_annotation <- function(path) {
  ann <- read.csv(path, stringsAsFactors = FALSE)
  ann <- ann[!is.na(ann$cell_type) & ann$cell_type != "Unlabeled", ]
  stats::setNames(as.character(ann$cell_type), as.character(ann$cell_id))
}

specs <- lapply(strsplit(rest, "=", fixed = TRUE), function(parts) {
  if (length(parts) != 4) stop("sample spec must be <sample>=<bundle>=<annotation>=<separate_dir>")
  list(sample = parts[[1]], bundle = parts[[2]], annotation = parts[[3]], dir = parts[[4]])
})
datasets <- lapply(specs, function(spec) {
  cellAdmix(spec$bundle, output_dir = spec$dir, annotation = read_annotation(spec$annotation),
    num_threads = threads)
})
names(datasets) <- vapply(specs, `[[`, character(1), "sample")

collection <- cellAdmixCollection(datasets)
t0 <- Sys.time()
nmf_variant <- Sys.getenv("CELLADMIX_VARIANT", "ls_nmf")
score_method <- Sys.getenv("CELLADMIX_SCORE", "bridge")
collection$fit(nmf_variant = nmf_variant, nmf_n_runs = 10L)
log_line("per-sample fits ready in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")

t0 <- Sys.time()
min_stability <- as.numeric(Sys.getenv("CELLADMIX_MIN_STABILITY", "0"))
cfit <- if (identical(mode, "consensus")) {
  collection$fit_consensus(min_stability = min_stability, num_threads = threads)
} else if (identical(mode, "reference")) {
  own <- collection$fits[[mode_arg]]
  collection$fit_profiles(
    cellAdmixCore:::.celladmix_fit_profiles_raw(own),
    members = cellAdmixCore:::.celladmix_fit_member_profiles(own),
    name = paste0("reference_", mode_arg), num_threads = threads)
} else if (identical(mode, "profiles")) {
  run <- cellAdmixCore:::celladmix_read_run(mode_arg)
  h <- run$h
  colnames(h) <- run$genes
  members <- cellAdmixCore:::.celladmix_ensemble_h(run$path)
  collection$fit_profiles(h, members = members, name = "joint_profiles", num_threads = threads)
} else {
  stop("unknown mode: ", mode)
}
log_line(mode, "projection done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
  "min; factors", nrow(cfit$loadings()))
write.csv(cfit$matching, file.path(out_dir, "matching.csv"), row.names = FALSE)
write.csv(cfit$loss, file.path(out_dir, "loss.csv"), row.names = FALSE)
write.csv(cfit$stability(), file.path(out_dir, "stability.csv"), row.names = FALSE)
write.csv(cfit$factors, file.path(out_dir, "factors.csv"), row.names = FALSE)
write.csv(cfit$factor_stability(), file.path(out_dir, "factor_stability.csv"), row.names = FALSE)
write.csv(cfit$unreproduced, file.path(out_dir, "unreproduced.csv"), row.names = FALSE)
write.csv(cfit$sample_usage(), file.path(out_dir, "sample_usage.csv"), row.names = FALSE)
write.csv(cfit$usage_divergence(), file.path(out_dir, "usage_divergence.csv"), row.names = FALSE)
profiles <- cfit$loadings()
write.csv(data.frame(factor = rownames(profiles), profiles, check.names = FALSE),
  file.path(out_dir, "profiles.csv"), row.names = FALSE)

t0 <- Sys.time()
score <- if (identical(score_method, "membrane")) cfit$score_membrane() else cfit$score_bridge()
log_line("scoring done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")
write.csv(score$summary(), file.path(out_dir, "pooled_summary.csv"), row.names = FALSE)
write.csv(score$sample_summaries(), file.path(out_dir, "sample_summaries.csv"), row.names = FALSE)
rules <- score$rules()
write.csv(rules, file.path(out_dir, "rules.csv"), row.names = FALSE)
log_line("rules:", nrow(rules), "kept", sum(rules$keep),
  "fully supported", sum(rules$keep & rules$sample_support >= 1))

t0 <- Sys.time()
apply_rules <- Sys.getenv("CELLADMIX_APPLY_RULES", "everywhere")
correction <- score$correct(rules = rules, apply_rules = apply_rules,
  name = paste0(score_method, "_clean", if (identical(apply_rules, "where_present")) "_where_present" else ""))
log_line("correction done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")
write.csv(correction$summary(), file.path(out_dir, "correction_summary.csv"), row.names = FALSE)
write.csv(correction$rules, file.path(out_dir, "correction_rules.csv"), row.names = FALSE)

t0 <- Sys.time()
audit <- cfit$audit_admixture()
write.csv(audit$pairs(by_sample = TRUE), file.path(out_dir, "audit_pairs.csv"), row.names = FALSE)
report <- audit$evaluate(correction)
write.csv(report$summary(), file.path(out_dir, "audit_evaluation.csv"), row.names = FALSE)
write.csv(report$pairs(), file.path(out_dir, "audit_evaluation_pairs.csv"), row.names = FALSE)
log_line("audit done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")

pointers <- data.frame(
  sample = cfit$samples(),
  run_dir = vapply(cfit$fits, function(f) f$run_dir, character(1)),
  correction_dir = vapply(correction$corrections, function(c) c$run$path, character(1)),
  own_run_dir = vapply(cfit$own_fits, function(f) f$run_dir, character(1)),
  stringsAsFactors = FALSE)
write.csv(pointers, file.path(out_dir, "pointers.csv"), row.names = FALSE)
log_line("done")
