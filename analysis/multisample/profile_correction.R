# Times the stages of a collection correction (member-label preparation,
# member scoring, pooled member rules, per-sample votes) on cached per-sample
# fits. Usage: Rscript profile_correction.R <threads_per_dataset> <budget>
#          <sample>=<bundle>=<annotation>=<separate_dir> ...
.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/home/pkharchenko/cellAdmix/wt-dev_multisample/.Rlib"), .libPaths()))
suppressMessages(library(cellAdmixCore))
args <- commandArgs(trailingOnly = TRUE)
threads <- as.integer(args[[1]]); budget <- as.integer(args[[2]]); rest <- args[-(1:2)]
log_line <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
timed <- function(label, expr) {
  t0 <- Sys.time(); value <- force(expr)
  log_line(label, "in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min")
  invisible(value)
}
read_annotation <- function(path) {
  ann <- read.csv(path, stringsAsFactors = FALSE)
  ann <- ann[!is.na(ann$cell_type) & ann$cell_type != "Unlabeled", ]
  stats::setNames(as.character(ann$cell_type), as.character(ann$cell_id))
}
specs <- lapply(strsplit(rest, "=", fixed = TRUE), function(p)
  list(sample = p[[1]], bundle = p[[2]], annotation = p[[3]], dir = p[[4]]))
datasets <- lapply(specs, function(s) cellAdmix(s$bundle, output_dir = s$dir,
  annotation = read_annotation(s$annotation), num_threads = threads))
names(datasets) <- vapply(specs, `[[`, character(1), "sample")
collection <- cellAdmixCollection(datasets)
timed("fits", collection$fit(nmf_variant = "ls_nmf", nmf_n_runs = 10L))
cfit <- timed("consensus", collection$fit_consensus(min_stability = 0.5, num_threads = budget))
score <- timed("scoring", cfit$score_bridge())
rules <- timed("rules", score$rules())
kept <- rules[is.na(rules$keep) | rules$keep, , drop = FALSE]
fits <- cfit$fits
for (s in names(fits)) {
  timed(paste("prepare labels", s),
    cellAdmixCore:::.celladmix_ensemble_prepare(fits[[s]]$run$path, threads))
}
priv <- score$.__enclos_env__$private
selected <- as.integer(fits[[1]]$run$nmf_diagnostics$selected_run) - 1L
members <- unique(c(selected, setdiff(0:9, selected)))
timed("score members (concurrent)", priv$score_members(setdiff(members, selected), budget))
member_rules <- timed("pooled member rules", priv$pooled_member_rules(members, selected, kept,
  p_thresh = 0.1, adjust_p = FALSE, targets = NULL, num_threads = budget))
min_votes <- 3L
for (s in names(fits)) {
  timed(paste("vote", s), {
    out_dir <- file.path(fits[[s]]$run$paths$corrected_dir, "profile_votes")
    cellAdmixCore:::celladmix_correct(fits[[s]]$run, rules = member_rules$rules, out_dir = out_dir,
      annotation = score$sample_scores[[s]]$annotation_vector(),
      rule_member = member_rules$member, min_votes = min_votes)
  })
}
log_line("done")
