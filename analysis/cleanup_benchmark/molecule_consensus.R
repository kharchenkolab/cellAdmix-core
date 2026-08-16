# Molecule-action consensus: each seed's pipeline runs to its end product
# (the per-molecule removal set); a molecule is removed iff a majority of
# seeds remove it. This ensembles away both decision-stage and label-stage
# seed variance, unlike rule-level consensus which only pools decisions.
#
# Usage: Rscript analysis/cleanup_benchmark/molecule_consensus.R <dataset>
#   dataset in {pancreas, breast_crop, nsclc}

.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/tmp/celladmix_r_lib"), .libPaths()))
suppressMessages(library(cellAdmixCore))
suppressMessages(library(Matrix))
bench_dir <- normalizePath(Sys.getenv("BENCH_DIR", "analysis/cleanup_benchmark"))
source(file.path(bench_dir, "metrics.R"))

dataset <- commandArgs(trailingOnly = TRUE)[[1]]
cfg <- switch(dataset,
  pancreas = list(dir = "examples/xenium_pancreas_membrane_377_full",
    ann = "annotations/annotation.csv.gz", ann_col = "merged_annotation",
    ann_id = "cell_id", out = "out", methods = c("membrane", "bridge"),
    build = function(ca) cellAdmix("data", output_dir = "out", annotation = ca)),
  breast_crop = list(dir = "examples/xenium_breast_membrane_5k_full",
    ann = "annotations/annotation.csv.gz", ann_col = "merged_annotation",
    ann_id = "cell_id", out = "out_medium_crop", methods = c("membrane", "bridge"),
    build = function(ca) cellAdmix("data", output_dir = "out_medium_crop",
      annotation = ca[ca != "Ambiguous / low-quality"], num_threads = 10,
      analysis_bbox = c(5000, 7000, 4500, 6500))),
  nsclc = list(dir = "examples/cosmx_nsclc_giotto",
    ann = "prepared/cell_metadata_all.csv.gz", ann_col = "cell_type_coarse",
    ann_id = "cell", out = "out", methods = "bridge",
    build = function(ca) cellAdmix(file.path("prepared", "molecules_all.csv.gz"),
      output_dir = "out", annotation = ca)))

ann_df <- read.csv(file.path(cfg$dir, cfg$ann), stringsAsFactors = FALSE)
cell_annotation <- setNames(ann_df[[cfg$ann_col]], ann_df[[cfg$ann_id]])
setwd(cfg$dir)
ds <- cfg$build(cell_annotation)
if (dataset == "breast_crop") {
  cell_annotation <- cell_annotation[cell_annotation != "Ambiguous / low-quality"]
}

fit0 <- ds$fit(verbose = FALSE)
counts_before <- fit0$counts()
cells <- fit0$cell_factors()
exp15 <- cellAdmixCore:::.celladmix_source_exposure_counts(cells, cell_annotation, 15L)
rownames(exp15$counts) <- as.character(cells$cell_id)
types <- exp15$types
cell_types <- setNames(as.character(cell_annotation[as.character(cells$cell_id)]),
  as.character(cells$cell_id))
totals_before <- Matrix::colSums(counts_before)
profiles <- bench_type_profiles(counts_before, cell_types)
top_type <- colnames(profiles)[max.col(profiles, ties.method = "first")]
native_markers <- lapply(types, function(t) {
  cand <- rownames(profiles)[top_type == t & profiles[, t] >= 200]
  mo <- apply(profiles[cand, setdiff(types, t), drop = FALSE], 1, max)
  spec <- profiles[cand, t] / pmax(profiles[cand, t] + mo, 1e-9)
  head(cand[order(spec, decreasing = TRUE)], 10)
})
names(native_markers) <- types
mcount <- function(m, genes, cs) {
  genes <- intersect(genes, rownames(m)); cs <- intersect(cs, colnames(m))
  Matrix::colSums(m[genes, cs, drop = FALSE])
}

pair_defs <- list()
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
  bins <- bench_exposure_bins(expo)
  rates <- bench_bin_rates(mcount(counts_before, pool, T_cells), totals_before[T_cells], bins)
  det <- bench_detect(rates)
  if (det$p > 1e-4 || det$excess_molecules < 200) next
  rates_s <- bench_bin_rates(mcount(counts_before, strict, T_cells), totals_before[T_cells], bins)
  pair_defs[[paste(S, T_type, sep = " -> ")]] <- list(S = S, T_type = T_type,
    T_cells = T_cells, expo = expo, bins = bins, pool = pool, strict = strict,
    excess = det$excess_molecules, rates_before = rates,
    excess_strict = bench_detect(rates_s)$excess_molecules, rates_before_strict = rates_s)
}
message("pairs: ", length(pair_defs))

# consensus-removal deltas precomputed by consensus_delta.py
for (variant in c("ls_nmf", "invsqrt_kl")) {
  for (method in cfg$methods) {
   for (tag in c("molcons", "moluni")) {
    delta_csv <- file.path(bench_dir, "results",
      sprintf("%s_%s_%s_%s.csv.gz", dataset, tag, variant, method))
    if (!file.exists(delta_csv)) { message("missing: ", delta_csv); next }
    dd <- read.csv(gzfile(delta_csv), stringsAsFactors = FALSE)
    dd <- dd[dd$gene %in% rownames(counts_before) &
             dd$cell_id %in% colnames(counts_before), , drop = FALSE]
    delta <- sparseMatrix(
      i = match(dd$gene, rownames(counts_before)),
      j = match(dd$cell_id, colnames(counts_before)),
      x = dd$n,
      dims = dim(counts_before), dimnames = dimnames(counts_before))
    counts_after <- counts_before - delta
    counts_after@x <- pmax(counts_after@x, 0)

    pow <- vapply(pair_defs, function(d) {
      bench_power(d$rates_before, bench_bin_rates(
        mcount(counts_after, d$pool, d$T_cells), totals_before[d$T_cells], d$bins))
    }, 0)
    pow_s <- vapply(pair_defs, function(d) {
      if (length(d$strict) < 2) return(NA_real_)
      bench_power(d$rates_before_strict, bench_bin_rates(
        mcount(counts_after, d$strict, d$T_cells), totals_before[d$T_cells], d$bins))
    }, 0)
    saf <- vapply(pair_defs, function(d) {
      bench_safety_native(
        mcount(counts_before, native_markers[[d$T_type]], d$T_cells),
        mcount(counts_after, native_markers[[d$T_type]], d$T_cells),
        totals_before[d$T_cells], d$bins)
    }, 0)
    ex <- vapply(pair_defs, function(d) d$excess, 0)
    exs <- vapply(pair_defs, function(d) d$excess_strict, 0)
    oks <- !is.na(pow_s)
    message(sprintf(
      "%s %s %s/%s: overall=%.3f strict=%.3f pairs>=0.9: %d/%d identity med=%.3f min=%.3f",
      toupper(tag), dataset, variant, method,
      sum(pow * ex) / sum(ex),
      sum(pow_s[oks] * exs[oks]) / sum(exs[oks]),
      sum(pow >= 0.9, na.rm = TRUE), length(pow),
      stats::median(saf, na.rm = TRUE), min(saf, na.rm = TRUE)))
   }
  }
}
message("MOLECULE CONSENSUS DONE")
