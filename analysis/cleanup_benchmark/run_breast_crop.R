# Cleanup benchmark: breast 5K medium-crop Xenium dataset.
# Usage: Rscript analysis/cleanup_benchmark/run_breast_crop.R

.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/tmp/celladmix_r_lib"), .libPaths()))
suppressMessages(library(cellAdmixCore))
bench_dir <- normalizePath(Sys.getenv("BENCH_DIR", "analysis/cleanup_benchmark"))
source(file.path(bench_dir, "metrics.R"))
source(file.path(bench_dir, "run_dataset.R"))
Sys.setenv(BENCH_DIR = bench_dir)

annotation_df <- read.csv("examples/xenium_breast_membrane_5k_full/annotations/annotation.csv.gz",
  stringsAsFactors = FALSE)
cell_annotation <- setNames(annotation_df$merged_annotation, annotation_df$cell_id)
cell_annotation <- cell_annotation[cell_annotation != "Ambiguous / low-quality"]
setwd("examples/xenium_breast_membrane_5k_full")
ds <- cellAdmix("data", output_dir = "out_medium_crop", annotation = cell_annotation,
  num_threads = 10, analysis_bbox = c(5000, 7000, 4500, 6500))

bench_run(ds, cell_annotation, label = "breast_crop")
message("DONE")
