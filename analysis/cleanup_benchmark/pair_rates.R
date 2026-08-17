# Per-pair extrapolated admixture rates for the pancreas dataset: the
# demonstrated pool leakage L, the pool's source-profile coverage s, and
# the extrapolated admixed-molecule count A = L / s and admixture rate
# r = L / (s * M_T). Writes results/pancreas_pair_rates.csv for the
# Figure 1 heatmap panel.
.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/tmp/celladmix_r_lib"), .libPaths()))
suppressMessages(library(cellAdmixCore))
bench_dir <- normalizePath(Sys.getenv("BENCH_DIR", "analysis/cleanup_benchmark"))

setwd("examples/xenium_pancreas_membrane_377_full")
annotation <- read.csv(file.path("annotations", "annotation.csv.gz"))
cell_annotation <- setNames(annotation$merged_annotation, annotation$cell_id)
ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)
fit <- ds$fit(nmf_variant = "invsqrt_kl", verbose = FALSE)

audit <- fit$audit_admixture()
pairs <- audit$pairs(detected_only = TRUE)
counts <- fit$counts()
cells <- fit$cell_factors()
cell_ids <- as.character(cells$cell_id)
types <- as.character(cell_annotation[cell_ids])
totals <- Matrix::colSums(counts)

rows <- lapply(seq_len(nrow(pairs)), function(i) {
  source_type <- pairs$source[[i]]
  target_type <- pairs$target[[i]]
  pool <- audit$markers(source_type, target_type)$pool
  source_cells <- intersect(cell_ids[!is.na(types) & types == source_type],
    colnames(counts))
  target_cells <- intersect(cell_ids[!is.na(types) & types == target_type],
    colnames(counts))
  coverage <- sum(counts[intersect(pool, rownames(counts)), source_cells]) /
    sum(counts[, source_cells])
  target_molecules <- sum(totals[target_cells])
  data.frame(
    source = source_type, target = target_type,
    excess = pairs$excess[[i]],
    coverage = coverage,
    target_molecules = target_molecules,
    admixed_molecules = pairs$excess[[i]] / coverage,
    rate = pairs$excess[[i]] / (coverage * target_molecules))
})
out <- do.call(rbind, rows)
write.csv(out, file.path(bench_dir, "results", "pancreas_pair_rates.csv"),
  row.names = FALSE)
message("PAIR RATES DONE: ", nrow(out), " pairs")
