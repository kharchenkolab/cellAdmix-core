# Cleanup benchmark: CosMx NSCLC dataset (tabular; no stains, bridge only).
# Usage: Rscript analysis/cleanup_benchmark/run_nsclc.R

.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/tmp/celladmix_r_lib"), .libPaths()))
suppressMessages(library(cellAdmixCore))
bench_dir <- normalizePath(Sys.getenv("BENCH_DIR", "analysis/cleanup_benchmark"))
source(file.path(bench_dir, "metrics.R"))
source(file.path(bench_dir, "run_dataset.R"))
Sys.setenv(BENCH_DIR = bench_dir)

setwd("examples/cosmx_nsclc_giotto")
cell_meta <- read.csv(file.path("prepared", "cell_metadata_all.csv.gz"),
  stringsAsFactors = FALSE)
cell_annotation <- setNames(cell_meta$cell_type_coarse, cell_meta$cell)
ds <- cellAdmix(file.path("prepared", "molecules_all.csv.gz"),
  output_dir = "out", annotation = cell_annotation)

bench_run(ds, cell_annotation, label = "nsclc", methods = "bridge")
message("DONE")
