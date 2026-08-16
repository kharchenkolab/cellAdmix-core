# Cleanup benchmark: CosMx NSCLC, bridge only (seeds + consensus + anchor arms).
.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/tmp/celladmix_r_lib"), .libPaths()))
suppressMessages(library(cellAdmixCore))
bench_dir <- normalizePath(Sys.getenv("BENCH_DIR", "analysis/cleanup_benchmark"))
source(file.path(bench_dir, "metrics.R"))
source(file.path(bench_dir, "run_dataset.R"))
Sys.setenv(BENCH_DIR = bench_dir)

LABEL <- "nsclc"; METHODS <- "bridge"
setwd("examples/cosmx_nsclc_giotto")
cell_meta <- read.csv(file.path("prepared", "cell_metadata_all.csv.gz"), stringsAsFactors = FALSE)
cell_annotation <- setNames(cell_meta$cell_type_coarse, cell_meta$cell)
ds <- cellAdmix(file.path("prepared", "molecules_all.csv.gz"),
  output_dir = "out", annotation = cell_annotation)

seeds <- 1:3
variants <- c("ls_nmf", "invsqrt_kl")
runs_dir <- ds$prep$paths$runs_dir
python_bin <- "/home/pkharchenko/cellAdmix/cellAdmix-core/.venv-sdata/bin/python"
fit0 <- ds$fit(verbose = FALSE)
rank <- fit0$rank
extra <- list()
for (s in seeds) {
  for (v in variants) {
    ds$fit(nmf_variant = v, seed = s, verbose = FALSE,
      run_id = sprintf("bench_seed%d_%s", s, v))
  }
  seed_run <- sprintf("bench_seed%d_invsqrt_kl", s)
  h_csv <- file.path(bench_dir, "results", sprintf("%s_anchor_h_s%d.csv", LABEL, s))
  if (!file.exists(h_csv)) {
    status <- system2(python_bin, c(file.path(bench_dir, "anchor_h.py"),
      file.path(runs_dir, seed_run), rank, h_csv))
    stopifnot(status == 0)
  }
  H <- as.matrix(read.csv(h_csv, check.names = FALSE))
  anchor_id <- sprintf("bench_anchor_s%d", s)
  ds$fit(nmf_fixed_h = H, nmf_variant = "kl", nmf_n_runs = 1L, seed = s,
    run_id = anchor_id, verbose = FALSE)
  extra[[sprintf("anchor_s%d", s)]] <- anchor_id
}
bench_run(ds, cell_annotation, label = LABEL, methods = METHODS,
  seeds = seeds, consensus = TRUE, extra_fits = extra)
message("DONE")
