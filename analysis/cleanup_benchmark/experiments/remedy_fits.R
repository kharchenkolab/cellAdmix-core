.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(cellAdmixCore))
setwd("/home/pkharchenko/cellAdmix/cellAdmix-core/examples/xenium_pancreas_membrane_377_full")
annotation <- read.csv("annotations/annotation.csv.gz", stringsAsFactors = FALSE)
cell_annotation <- setNames(annotation$merged_annotation, annotation$cell_id)
ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)

for (seed in 1:3) {
  fit <- ds$fit(nmf_variant = "invsqrt_kl", seed = seed, nmf_n_runs = 30L,
    run_id = sprintf("bench_nr30_seed%d_invsqrt_kl", seed), verbose = TRUE)
  fit <- ds$fit(nmf_variant = "invsqrt_kl", seed = seed, rank = 12L,
    run_id = sprintf("bench_r12_seed%d_invsqrt_kl", seed), verbose = TRUE)
}
cat("REMEDY FITS DONE\n")
