# Score the remedy fits (30-restart, rank-12) and a consensus-rule prototype
# on the pancreas benchmark's strict tier.
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(cellAdmixCore))
suppressMessages(library(Matrix))
bench_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/cleanup_benchmark"
source(file.path(bench_dir, "metrics.R"))

setwd("/home/pkharchenko/cellAdmix/cellAdmix-core/examples/xenium_pancreas_membrane_377_full")
annotation <- read.csv("annotations/annotation.csv.gz", stringsAsFactors = FALSE)
cell_annotation <- setNames(annotation$merged_annotation, annotation$cell_id)
ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)

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
  if (length(strict) < 2) next
  bins <- bench_exposure_bins(expo)
  rates_s <- bench_bin_rates(mcount(counts_before, strict, T_cells), totals_before[T_cells], bins)
  det_s <- bench_detect(rates_s)
  det <- bench_detect(bench_bin_rates(mcount(counts_before, pool, T_cells),
    totals_before[T_cells], bins))
  if (det$p > 1e-4 || det$excess_molecules < 200) next
  pair_defs[[paste(S, T_type, sep = " -> ")]] <- list(S = S, T_type = T_type,
    T_cells = T_cells, bins = bins, strict = strict,
    excess_strict = det_s$excess_molecules, rates_before_strict = rates_s)
}
message("pairs: ", length(pair_defs))

eval_correction <- function(counts_after, tag) {
  pow <- vapply(pair_defs, function(d) {
    ra <- bench_bin_rates(mcount(counts_after, d$strict, d$T_cells),
      totals_before[d$T_cells], d$bins)
    bench_power(d$rates_before_strict, ra)
  }, 0)
  ex <- vapply(pair_defs, function(d) d$excess_strict, 0)
  ok <- !is.na(pow)
  overall <- sum(pow[ok] * ex[ok]) / sum(ex[ok])
  message(sprintf("%-28s leakage_removed_strict = %.3f  (pairs>=0.9: %d/%d)",
    tag, overall, sum(pow[ok] >= 0.9), sum(ok)))
  invisible(pow)
}

score_fit <- function(run_id, tag, consensus_pairs = NULL) {
  fit <- ds$read_fit(file.path("out", "runs", run_id))
  score <- fit$score_membrane()
  rules <- score$rules(p_thresh = 0.1)
  if (!is.null(consensus_pairs)) {
    ann <- score$annotation(p_thresh = 0.1)
    src <- ann$source_calls
    have <- unique(paste(rules$source_cell_type[rules$keep],
      rules$target_cell_type[rules$keep]))
    add <- list()
    for (cp in consensus_pairs) {
      if (cp %in% have) next
      parts <- strsplit(cp, "\\|\\|")[[1]]
      S <- parts[[1]]; T_type <- parts[[2]]
      fs <- as.integer(sub("^f_", "", names(src)[vapply(src, function(v)
        identical(as.character(v), S), TRUE)]))
      for (f in fs) {
        row <- rules[1, , drop = FALSE]
        row$factor <- f; row$source_cell_type <- S; row$target_cell_type <- T_type
        row$keep <- TRUE; row$native_check <- "consensus"
        add[[length(add) + 1]] <- row
      }
    }
    if (length(add)) rules <- rbind(rules, do.call(rbind, add))
  }
  correction <- score$correct(rules = rules, name = paste0("cmp_rem_", tag))
  eval_correction(correction$counts(), paste0(tag, if (!is.null(consensus_pairs)) "+consensus"))
  list(rules = rules)
}

message("\n-- remedy fits --")
for (seed in 1:3) score_fit(sprintf("bench_nr30_seed%d_invsqrt_kl", seed), sprintf("nr30_s%d", seed))
for (seed in 1:3) score_fit(sprintf("bench_r12_seed%d_invsqrt_kl", seed), sprintf("r12_s%d", seed))

message("\n-- consensus rules over original 3 seeds --")
kept_pairs <- character()
for (seed in 1:3) {
  fit <- ds$read_fit(sprintf("out/runs/bench_seed%d_invsqrt_kl", seed))
  r <- fit$score_membrane()$rules(p_thresh = 0.1)
  kept_pairs <- union(kept_pairs, paste(r$source_cell_type[r$keep],
    r$target_cell_type[r$keep], sep = "||"))
}
message("consensus (S,T) coverage: ", length(kept_pairs), " pairs")
for (seed in 1:3) score_fit(sprintf("bench_seed%d_invsqrt_kl", seed),
  sprintf("orig_s%d", seed), consensus_pairs = kept_pairs)
message("REMEDY SCORE DONE")
