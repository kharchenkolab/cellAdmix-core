# Figure data for docs/benchmarks.md: exposure-bin marker rates before/after
# a standard (bare ls_nmf, seed 1) membrane cleanup, for every detected pair
# on the pancreas dataset.
.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/tmp/celladmix_r_lib"), .libPaths()))
suppressMessages(library(cellAdmixCore))
suppressMessages(library(Matrix))
bench_dir <- normalizePath(Sys.getenv("BENCH_DIR", "analysis/cleanup_benchmark"))
source(file.path(bench_dir, "metrics.R"))

annotation_df <- read.csv("examples/xenium_pancreas_membrane_377_full/annotations/annotation.csv.gz",
  stringsAsFactors = FALSE)
cell_annotation <- setNames(annotation_df$merged_annotation, annotation_df$cell_id)
setwd("examples/xenium_pancreas_membrane_377_full")
ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)

fit <- ds$fit(nmf_variant = "ls_nmf", seed = 1L, verbose = FALSE,
  run_id = "bench_seed1_ls_nmf")
counts_before <- fit$counts()
score <- fit$score_membrane()
rules <- score$rules(p_thresh = 0.1)
counts_after <- score$correct(rules = rules, name = "cmp_membrane_ls_nmf_s1")$counts()

cells <- fit$cell_factors()
exp15 <- cellAdmixCore:::.celladmix_source_exposure_counts(cells, cell_annotation, 15L)
rownames(exp15$counts) <- as.character(cells$cell_id)
types <- exp15$types
# Exposure at each neighborhood size of the reference ladder, for the
# ambient reference computed on the same gene pool as the plotted rates.
ladder_K <- c(15L, 30L, 60L, 120L, 240L)
exp_ladder <- lapply(ladder_K, function(K) {
  e <- cellAdmixCore:::.celladmix_source_exposure_counts(cells, cell_annotation, K)
  rownames(e$counts) <- as.character(cells$cell_id)
  e$counts
})
names(exp_ladder) <- paste0("k", ladder_K)
cell_types <- setNames(as.character(cell_annotation[as.character(cells$cell_id)]),
  as.character(cells$cell_id))
totals_before <- Matrix::colSums(counts_before)
profiles <- bench_type_profiles(counts_before, cell_types)
mcount <- function(m, genes, cs) {
  genes <- intersect(genes, rownames(m)); cs <- intersect(cs, colnames(m))
  Matrix::colSums(m[genes, cs, drop = FALSE])
}

rows <- list()
for (S in types) for (T_type in setdiff(types, S)) {
  T_cells <- intersect(names(cell_types)[!is.na(cell_types) & cell_types == T_type],
    intersect(rownames(exp15$counts), colnames(counts_before)))
  if (length(T_cells) < 200) next
  expo <- exp15$counts[T_cells, S]
  if (sum(expo == 0) < 100) next
  baseline_T0 <- bench_pseudobulk(counts_before, T_cells[expo == 0])
  pool <- bench_marker_pool(profiles, S, T_type, baseline_T0)
  if (length(pool) < 3) next
  strict <- bench_pool_strict(pool, profiles, S, baseline_T0)
  if (length(strict) < 2) next
  bins <- bench_exposure_bins(expo)
  rb <- bench_bin_rates(mcount(counts_before, strict, T_cells), totals_before[T_cells], bins)
  ra <- bench_bin_rates(mcount(counts_after, strict, T_cells), totals_before[T_cells], bins)
  det <- bench_detect(rb)
  zero_by_K <- lapply(exp_ladder, function(ec) ec[T_cells, S] == 0)
  ref <- cellAdmixCore:::.celladmix_audit_reference_rate(
    mcount(counts_before, strict, T_cells), totals_before[T_cells], zero_by_K)
  for (i in seq_len(nrow(rb))) {
    rows[[length(rows) + 1]] <- data.frame(
      source = S, target = T_type, bin = rb$bin[[i]], n_cells = rb$n_cells[[i]],
      rate_before = rb$rate[[i]], rate_after = ra$rate[[i]],
      excess_strict = det$excess_molecules, p_detect = det$p,
      reference_rate = ref$rate, reference_kind = ref$kind,
      stringsAsFactors = FALSE)
  }
}
out <- do.call(rbind, rows)
write.csv(out, file.path(bench_dir, "results", "pancreas_fig_bins.csv"), row.names = FALSE)
message("wrote pancreas_fig_bins.csv: ", nrow(out), " rows")
