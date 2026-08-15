# Cleanup benchmark: pancreas Xenium dataset.
# Usage: Rscript analysis/cleanup_benchmark/run_pancreas.R [seeds]
# Passing "seeds" runs the 3-seed reproducibility arm (membrane only).

.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/tmp/celladmix_r_lib"), .libPaths()))
suppressMessages(library(cellAdmixCore))
bench_dir <- normalizePath(Sys.getenv("BENCH_DIR", "analysis/cleanup_benchmark"))
source(file.path(bench_dir, "metrics.R"))
source(file.path(bench_dir, "run_dataset.R"))
Sys.setenv(BENCH_DIR = bench_dir)

annotation_df <- read.csv("examples/xenium_pancreas_membrane_377_full/annotations/annotation.csv.gz",
  stringsAsFactors = FALSE)
cell_annotation <- setNames(annotation_df$merged_annotation, annotation_df$cell_id)
setwd("examples/xenium_pancreas_membrane_377_full")
ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)

if (identical(commandArgs(trailingOnly = TRUE), "seeds")) {
  bench_run(ds, cell_annotation, label = "pancreas_seeds",
    methods = "membrane", seeds = c(1L, 2L, 3L))
} else {
  bench_run(ds, cell_annotation, label = "pancreas")
}
message("DONE")
