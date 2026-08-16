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
suppressMessages(library(arrow))
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

# molecule -> (gene, cell) map from one fit run (identical across fits)
run_dir1 <- file.path(cfg$out, "runs", "bench_seed1_ls_nmf")
mol <- as.data.frame(read_parquet(file.path(run_dir1, "molecules.parquet"),
  col_select = c("obs_id", "gene_idx", "cell_idx")))
run_meta <- jsonlite::fromJSON(file.path(run_dir1, "run.json"))
genes_vec <- run_meta$genes
cells_tab <- as.data.frame(read_parquet(file.path(run_dir1, "cells.parquet"),
  col_select = c("cell_idx", "cell_id")))
cell_id_by_idx <- setNames(as.character(cells_tab$cell_id), cells_tab$cell_idx)
assigned <- mol$cell_idx >= 0
max_obs <- max(mol$obs_id) + 1L

for (variant in c("ls_nmf", "invsqrt_kl")) {
  for (method in cfg$methods) {
    votes <- integer(max_obs)
    n_seeds <- 0L
    for (s in 1:3) {
      cdir <- file.path(cfg$out, "runs", sprintf("bench_seed%d_%s", s, variant),
        "corrected", sprintf("cmp_%s_%s_s%d", method, variant, s))
      if (!dir.exists(cdir)) { message("missing: ", cdir); next }
      kept <- read_parquet(file.path(cdir, "molecules.parquet"),
        col_select = "obs_id")$obs_id
      removed <- rep(TRUE, max_obs)
      removed[kept + 1L] <- FALSE
      removed[mol$obs_id[!assigned] + 1L] <- FALSE
      votes <- votes + as.integer(removed)
      n_seeds <- n_seeds + 1L
    }
    if (n_seeds < 2) next
    consensus_removed <- votes >= 2L
    rm_mask <- consensus_removed[mol$obs_id + 1L] & assigned
    message(sprintf("%s/%s: consensus removes %d molecules (seed removals pooled from %d seeds)",
      variant, method, sum(rm_mask), n_seeds))
    delta <- sparseMatrix(
      i = mol$gene_idx[rm_mask] + 1L,
      j = mol$cell_idx[rm_mask] + 1L,
      x = 1,
      dims = c(length(genes_vec), nrow(cells_tab)),
      dimnames = list(genes_vec, {
        cn <- character(nrow(cells_tab))
        cn[cells_tab$cell_idx + 1L] <- as.character(cells_tab$cell_id)
        cn
      }))
    counts_after <- counts_before
    common_g <- intersect(rownames(counts_before), rownames(delta))
    common_c <- intersect(colnames(counts_before), colnames(delta))
    counts_after[common_g, common_c] <-
      pmax(counts_before[common_g, common_c] - delta[common_g, common_c], 0)

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
      "MOLCONS %s %s/%s: overall=%.3f strict=%.3f pairs>=0.9: %d/%d identity med=%.3f min=%.3f",
      dataset, variant, method,
      sum(pow * ex) / sum(ex),
      sum(pow_s[oks] * exs[oks]) / sum(exs[oks]),
      sum(pow >= 0.9, na.rm = TRUE), length(pow),
      stats::median(saf, na.rm = TRUE), min(saf, na.rm = TRUE)))
  }
}
message("MOLECULE CONSENSUS DONE")
