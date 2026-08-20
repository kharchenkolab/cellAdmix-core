# A/B gene-split evaluation harness for audit-guided correction experiments.
#
# Each detected pair's 20-gene marker pool is split into a guide half (A)
# and a validation half (B), alternating genes along the pool's contrast
# ranking so both halves carry similar coverage. Anything that *guides* a
# correction (anchors, priors, threshold selection) may only see A;
# performance is reported on B (and on the strict subset of B), with B's
# own re-estimated exposure-0 baseline. Own-marker retention per type is
# the specificity metric; native markers are disjoint from guide pools by
# construction check.
#
# Source with: source("01_split_eval.R"); then build_split_defs() once and
# eval_split(counts_after, tag) per correction.

suppressMessages(library(Matrix))
bench_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/cleanup_benchmark"
source(file.path(bench_dir, "metrics.R"))

ag_mcount <- function(m, genes, cs) {
  genes <- intersect(genes, rownames(m)); cs <- intersect(cs, colnames(m))
  if (!length(genes) || !length(cs)) return(setNames(numeric(length(cs)), cs))
  Matrix::colSums(m[genes, cs, drop = FALSE])
}

# Build pair definitions with split pools from a fitted dataset.
# Returns list(pairs = <defs>, native = <per-type native markers>,
#   cell_types, totals_before, counts_before).
build_split_defs <- function(fit, cell_annotation, k_neighbors = 15L,
                             n_pool = 20L, seed_offset = 0L) {
  counts_before <- fit$counts()
  cells <- fit$cell_factors()
  expo <- cellAdmixCore:::.celladmix_source_exposure_counts(cells, cell_annotation,
    k_neighbors)
  rownames(expo$counts) <- as.character(cells$cell_id)
  types <- expo$types
  cell_types <- setNames(as.character(cell_annotation[as.character(cells$cell_id)]),
    as.character(cells$cell_id))
  totals_before <- Matrix::colSums(counts_before)
  profiles <- bench_type_profiles(counts_before, cell_types)
  top_type <- colnames(profiles)[max.col(profiles, ties.method = "first")]

  pair_defs <- list()
  for (S in types) for (T_type in setdiff(types, S)) {
    T_cells <- intersect(names(cell_types)[!is.na(cell_types) & cell_types == T_type],
      intersect(rownames(expo$counts), colnames(counts_before)))
    if (length(T_cells) < 200) next
    e <- expo$counts[T_cells, S]
    if (sum(e == 0) < 100) next
    baseline_T0 <- bench_pseudobulk(counts_before, T_cells[e == 0])
    pool <- bench_marker_pool(profiles, S, T_type, baseline_T0, n_pool = n_pool)
    if (length(pool) < 6) next
    strict <- bench_pool_strict(pool, profiles, S, baseline_T0)
    bins <- bench_exposure_bins(e)
    det <- bench_detect(bench_bin_rates(ag_mcount(counts_before, pool, T_cells),
      totals_before[T_cells], bins))
    if (det$p > 1e-4 || det$excess_molecules < 200) next
    # Alternating split along the contrast ranking (pool is rank-ordered).
    idx <- seq_along(pool)
    A <- pool[(idx + seed_offset) %% 2 == 1]
    B <- pool[(idx + seed_offset) %% 2 == 0]
    rates_of <- function(genes, m) bench_bin_rates(ag_mcount(m, genes, T_cells),
      totals_before[T_cells], bins)
    pair_defs[[paste(S, T_type, sep = " -> ")]] <- list(
      S = S, T_type = T_type, T_cells = T_cells, bins = bins,
      e = setNames(as.numeric(e), T_cells),
      pool = pool, A = A, B = B, strict = strict, strictB = intersect(strict, B),
      excess_pool = det$excess_molecules,
      excess_B = bench_detect(rates_of(B, counts_before))$excess_molecules,
      rates_before = list(A = rates_of(A, counts_before),
        B = rates_of(B, counts_before),
        strictB = if (length(intersect(strict, B)) >= 2)
          rates_of(intersect(strict, B), counts_before) else NULL))
  }

  native <- lapply(types, function(t) {
    own <- rownames(profiles)[top_type == t]
    own[order(profiles[own, t], decreasing = TRUE)][seq_len(min(30, length(own)))]
  })
  names(native) <- types

  list(pairs = pair_defs, native = native, cell_types = cell_types,
    totals_before = totals_before, counts_before = counts_before,
    profiles = profiles)
}

# Evaluate a corrected count matrix against the split definitions.
# Returns list(pairs = per-pair data frame, types = per-type false removal).
eval_split <- function(defs, counts_after, tag = "correction") {
  rows <- lapply(names(defs$pairs), function(nm) {
    d <- defs$pairs[[nm]]
    rates_of <- function(genes) bench_bin_rates(
      ag_mcount(counts_after, genes, d$T_cells),
      defs$totals_before[d$T_cells], d$bins)
    data.frame(pair = nm, source = d$S, target = d$T_type,
      excess_pool = round(d$excess_pool), excess_B = round(d$excess_B),
      power_A = bench_power(d$rates_before$A, rates_of(d$A)),
      power_B = bench_power(d$rates_before$B, rates_of(d$B)),
      power_strictB = if (!is.null(d$rates_before$strictB))
        bench_power(d$rates_before$strictB, rates_of(d$strictB)) else NA_real_,
      # Baseline erosion: how much of the exposure-0 B-marker signal was
      # removed. ~0 = calibrated removal; ~1 = the pair's whole gene content
      # was stripped, baseline included (over-removal unless the baseline is
      # pure ambient, i.e. the strict tier).
      depth_B = bench_removal_depth(d$rates_before$B, rates_of(d$B)),
      stringsAsFactors = FALSE)
  })
  pairs <- do.call(rbind, rows)
  fr <- do.call(rbind, lapply(names(defs$native), function(t) {
    tc <- names(defs$cell_types)[!is.na(defs$cell_types) & defs$cell_types == t]
    genes <- defs$native[[t]]
    before <- sum(ag_mcount(defs$counts_before, genes, tc))
    after <- sum(ag_mcount(counts_after, genes, tc))
    data.frame(cell_type = t, own_marker_molecules = round(before),
      false_removal = 1 - after / max(before, 1), stringsAsFactors = FALSE)
  }))
  wb <- function(p, ex) { ok <- !is.na(p); sum(p[ok] * ex[ok]) / sum(ex[ok]) }
  summary <- data.frame(tag = tag,
    n_pairs = nrow(pairs),
    power_B_weighted = wb(pairs$power_B, pairs$excess_B),
    power_B_median = stats::median(pairs$power_B, na.rm = TRUE),
    false_removal_pooled = sum(fr$false_removal * fr$own_marker_molecules) /
      sum(fr$own_marker_molecules),
    worst_false_removal = max(fr$false_removal))
  list(pairs = pairs, types = fr, summary = summary)
}

# Guide-pool <-> native-marker overlap check (should be empty per type-as-target).
check_disjoint <- function(defs) {
  bad <- 0
  for (nm in names(defs$pairs)) {
    d <- defs$pairs[[nm]]
    ov <- intersect(d$pool, defs$native[[d$T_type]])
    if (length(ov)) { bad <- bad + 1; message(nm, ": pool overlaps target natives: ",
      paste(ov, collapse = ",")) }
  }
  if (!bad) message("all pools disjoint from their target's native markers")
}
