# Times a complete collection correction (shared ensemble, concurrent
# per-sample steps) on cached per-sample fits and writes the removal per
# sample, to compare against a run of the same correction done sequentially.
# Usage: Rscript time_collection_correction.R <out_csv> <threads_per_dataset>
#          <budget> <correction_name> <sample>=<bundle>=<annotation>=<separate_dir> ...
.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/home/pkharchenko/cellAdmix/wt-dev_multisample/.Rlib"), .libPaths()))
suppressMessages(library(cellAdmixCore))
args <- commandArgs(trailingOnly = TRUE)
out_csv <- args[[1]]; threads <- as.integer(args[[2]]); budget <- as.integer(args[[3]])
correction_name <- args[[4]]; rest <- args[-(1:4)]
log_line <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
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
collection$fit(nmf_variant = "ls_nmf", nmf_n_runs = 10L)
cfit <- collection$fit_consensus(min_stability = 0.5, num_threads = budget)
score <- cfit$score_bridge()
rules <- score$rules()
t0 <- Sys.time()
correction <- score$correct(rules = rules, name = correction_name, num_threads = budget)
log_line("correction done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")
removed <- data.frame(sample = names(correction$corrections),
  n_removed = vapply(correction$corrections, function(x) as.numeric(x$run$n_removed), numeric(1)))
write.csv(removed, out_csv, row.names = FALSE)
print(removed)
