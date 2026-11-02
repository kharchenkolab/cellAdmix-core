# Positive control for the per-sample diagnostics: project a deliberately
# poor profile set (one sample's rank-10 fit) onto the rank-23 breast
# replicates through the reference-projection path and print the per-sample
# diagnostics; the representation and shared-factor flags should fire.
#
# Usage: Rscript diagnostics_check.R <out_dir> <threads>
#          <sample>=<bundle>=<annotation>=<separate_dir> ...
.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/home/pkharchenko/cellAdmix/wt-dev_multisample/.Rlib_dev"), .libPaths()))
suppressMessages(library(cellAdmixCore))

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[[1]]
threads <- as.integer(args[[2]])
rest <- args[-(1:2)]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

read_annotation <- function(path) {
  ann <- read.csv(path, stringsAsFactors = FALSE)
  ann <- ann[!is.na(ann$cell_type) & ann$cell_type != "Unlabeled", ]
  stats::setNames(as.character(ann$cell_type), as.character(ann$cell_id))
}
specs <- lapply(strsplit(rest, "=", fixed = TRUE), function(parts) {
  list(sample = parts[[1]], bundle = parts[[2]], annotation = parts[[3]], dir = parts[[4]])
})
datasets <- lapply(specs, function(spec) {
  cellAdmix(spec$bundle, output_dir = spec$dir, annotation = read_annotation(spec$annotation),
    num_threads = threads)
})
names(datasets) <- vapply(specs, `[[`, character(1), "sample")

collection <- cellAdmixCollection(datasets)
collection$fit(nmf_variant = "ls_nmf", nmf_n_runs = 10L)
log_line("per-sample fits ready")

t0 <- Sys.time()
rank_arg <- intersect(c("rank", "nmf_rank"), names(formals(datasets[[1]]$fit)))[[1]]
poor <- do.call(datasets[[1]]$fit, c(list(nmf_variant = "ls_nmf", nmf_n_runs = 10L),
  stats::setNames(list(10L), rank_arg)))
log_line("rank-10 fit of", names(datasets)[[1]], "done in",
  round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")

t0 <- Sys.time()
control <- collection$fit_profiles(
  cellAdmixCore:::.celladmix_fit_profiles_raw(poor),
  members = cellAdmixCore:::.celladmix_fit_member_profiles(poor),
  name = "control_rank10", num_threads = threads)
log_line("projection done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")
diagnostics <- control$diagnostics()
write.csv(diagnostics, file.path(out_dir, "control_rank10_diagnostics.csv"), row.names = FALSE)
write.csv(control$usage_divergence(by_sample = TRUE),
  file.path(out_dir, "control_rank10_usage_divergence.csv"), row.names = FALSE)
print(diagnostics)
log_line("done")
