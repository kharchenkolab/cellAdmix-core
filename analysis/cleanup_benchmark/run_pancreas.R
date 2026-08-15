# Cleanup benchmark, pancreas Xenium dataset: all detectable source->target
# pairs, both NMF variants, membrane and bridge scoring.
#
# Usage: Rscript analysis/cleanup_benchmark/run_pancreas.R
# Writes results/pancreas_pairs.csv and prints a compact scorecard.

.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/tmp/celladmix_r_lib"), .libPaths()))
suppressMessages(library(cellAdmixCore))
suppressMessages(library(Matrix))
bench_dir <- Sys.getenv("BENCH_DIR", "analysis/cleanup_benchmark")
source(file.path(bench_dir, "metrics.R"))

data_dir <- "examples/xenium_pancreas_membrane_377_full"
out_dir <- file.path(bench_dir, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

annotation_df <- read.csv(file.path(data_dir, "annotations/annotation.csv.gz"),
  stringsAsFactors = FALSE)
cell_annotation <- setNames(annotation_df$merged_annotation, annotation_df$cell_id)
setwd(data_dir)
ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)

message("== shared setup")
fit0 <- ds$fit(verbose = FALSE)
counts_before <- fit0$counts()
cells <- fit0$cell_factors()
annotation <- cell_annotation[as.character(cells$cell_id)]
exp15 <- cellAdmixCore:::.celladmix_source_exposure_counts(cells, cell_annotation, 15L)
exp30 <- cellAdmixCore:::.celladmix_source_exposure_counts(cells, cell_annotation, 30L)
rownames(exp15$counts) <- rownames(exp30$counts) <- as.character(cells$cell_id)
types <- exp15$types
cell_types <- setNames(as.character(annotation), as.character(cells$cell_id))
totals_before <- Matrix::colSums(counts_before)
profiles <- bench_type_profiles(counts_before, cell_types)

# Isolated-cell reference profiles for the NNLS index (zero non-self exposure
# at k=15); fall back to all cells of the type when too few are isolated.
profiles_iso <- sapply(types, function(t) {
  members <- names(cell_types)[!is.na(cell_types) & cell_types == t]
  members <- intersect(members, rownames(exp15$counts))
  nonself <- rowSums(exp15$counts[members, setdiff(types, t), drop = FALSE])
  iso <- members[nonself == 0]
  if (length(iso) < 50) iso <- members
  bench_pseudobulk(counts_before, iso)
})
rownames(profiles_iso) <- rownames(counts_before)

# Native marker sets per type (no baseline filter; top specificity).
native_markers <- lapply(types, function(t) {
  others <- setdiff(types, t)
  max_other <- apply(profiles[, others, drop = FALSE], 1, max)
  spec <- profiles[, t] / pmax(profiles[, t] + max_other, 1e-9)
  ok <- profiles[, t] >= 200
  head(rownames(profiles)[ok][order(spec[ok], decreasing = TRUE)], 10)
})
names(native_markers) <- types

mcount <- function(m, genes, cs) {
  genes <- intersect(genes, rownames(m))
  cs <- intersect(cs, colnames(m))
  Matrix::colSums(m[genes, cs, drop = FALSE])
}

message("== enumerating pairs")
pair_defs <- list()
for (S in types) {
  for (T_type in setdiff(types, S)) {
    T_cells <- names(cell_types)[!is.na(cell_types) & cell_types == T_type]
    T_cells <- intersect(T_cells, intersect(rownames(exp15$counts), colnames(counts_before)))
    if (length(T_cells) < 200) next
    expo <- exp15$counts[T_cells, S]
    e0 <- T_cells[expo == 0]
    if (length(e0) < 100) next
    baseline_T0 <- bench_pseudobulk(counts_before, e0)
    pool <- bench_marker_pool(profiles, S, T_type, baseline_T0)
    if (length(pool) < 3) next
    bins <- bench_exposure_bins(expo)
    rates <- bench_bin_rates(mcount(counts_before, pool, T_cells), totals_before[T_cells], bins)
    det <- bench_detect(rates)
    pair_defs[[paste(S, T_type, sep = " -> ")]] <- list(S = S, T_type = T_type,
      T_cells = T_cells, expo = expo, bins = bins, pool = pool,
      p_detect = det$p, excess = det$excess_molecules, rates_before = rates)
  }
}
pvals <- vapply(pair_defs, function(d) d$p_detect, 0)
qvals <- p.adjust(pvals, "BH")
keep <- names(pair_defs)[qvals < 0.01 &
  vapply(pair_defs, function(d) d$excess, 0) >= 200]
message(sprintf("pairs considered: %d, detected: %d", length(pair_defs), length(keep)))
pair_defs <- pair_defs[keep]

arms <- expand.grid(variant = c("ls_nmf", "invsqrt_kl"),
  method = c("membrane", "bridge"), stringsAsFactors = FALSE)

results <- list()
for (ai in seq_len(nrow(arms))) {
  variant <- arms$variant[[ai]]
  method <- arms$method[[ai]]
  arm <- paste(variant, method, sep = "/")
  message("== arm: ", arm)
  fit <- ds$fit(nmf_variant = variant, verbose = FALSE)
  score <- if (method == "membrane") fit$score_membrane() else fit$score_bridge()
  score_annotation <- score$annotation(p_thresh = 0.1)
  rules <- score$rules(p_thresh = 0.1)
  correction <- score$correct(rules = rules, name = paste0("cmp_", method, "_", variant))
  counts_after <- correction$counts()

  # molecule-level labels for attribution (sampled for memory)
  mols <- tryCatch(fit$molecules(sample_n = 3000000L), error = function(e) NULL)
  source_calls <- score_annotation$source_calls
  factors_of_source <- function(S) {
    hits <- names(source_calls)[vapply(source_calls, function(v)
      identical(as.character(v), S), TRUE)]
    as.integer(sub("^f_", "", hits))
  }

  for (pname in names(pair_defs)) {
    d <- pair_defs[[pname]]
    exposed <- d$T_cells[d$expo > 0]
    e0 <- d$T_cells[d$expo == 0]
    rb <- d$rates_before
    ra <- bench_bin_rates(mcount(counts_after, d$pool, d$T_cells),
      totals_before[d$T_cells], d$bins)
    ctrl <- bench_control_genes(profiles, d$S, d$pool)
    crb <- bench_bin_rates(mcount(counts_before, ctrl, d$T_cells),
      totals_before[d$T_cells], d$bins)
    cra <- bench_bin_rates(mcount(counts_after, ctrl, d$T_cells),
      totals_before[d$T_cells], d$bins)
    dose <- exp30$counts[d$T_cells, d$S]
    rule_rows <- rules[rules$source_cell_type == d$S &
      rules$target_cell_type == d$T_type, , drop = FALSE]
    label_frac <- NA_real_
    if (!is.null(mols) && all(c("gene", "cell", "factor_label") %in% names(mols))) {
      sub <- mols[mols$cell %in% exposed & mols$gene %in% d$pool, , drop = FALSE]
      if (nrow(sub) >= 50) {
        label_frac <- mean(sub$factor_label %in% factors_of_source(d$S))
      }
    }
    mk_before <- sum(mcount(counts_before, d$pool, exposed))
    mk_after <- sum(mcount(counts_after, d$pool, exposed))
    results[[length(results) + 1]] <- data.frame(
      pair = pname, source = d$S, target = d$T_type, variant = variant,
      method = method,
      excess_molecules = round(d$excess),
      n_exposed = length(exposed), n_e0 = length(e0),
      power = round(bench_power(rb, ra), 3),
      power_dose = round(bench_dose_power(
        mcount(counts_before, d$pool, d$T_cells),
        mcount(counts_after, d$pool, d$T_cells),
        totals_before[d$T_cells], dose), 3),
      ctrl_power = round(bench_power(crb, cra), 3),
      safety_baseline = round(bench_safety_baseline(rb, ra), 3),
      safety_native = round(bench_safety_native(
        mcount(counts_before, native_markers[[d$T_type]], d$T_cells),
        mcount(counts_after, native_markers[[d$T_type]], d$T_cells),
        totals_before[d$T_cells], d$bins), 3),
      profile_cor_e0 = round(bench_profile_integrity(counts_before, counts_after, e0), 4),
      nnls_before = round(bench_nnls_contamination(counts_before, profiles_iso,
        exposed, d$S, d$T_type), 4),
      nnls_after = round(bench_nnls_contamination(counts_after, profiles_iso,
        exposed, d$S, d$T_type), 4),
      shift_before = round(bench_doublet_shift(counts_before, d$T_cells, d$expo,
        profiles[, d$S], profiles_iso[, d$T_type]), 3),
      shift_after = round(bench_doublet_shift(counts_after, d$T_cells, d$expo,
        profiles[, d$S], profiles_iso[, d$T_type]), 3),
      marker_removed_frac = round(1 - mk_after / max(mk_before, 1), 3),
      rule_called = nrow(rule_rows) > 0,
      rule_kept = any(rule_rows$keep),
      label_frac = round(label_frac, 3),
      stringsAsFactors = FALSE)
  }
}

res <- do.call(rbind, results)
out_csv <- file.path("../..", bench_dir, "results", "pancreas_pairs.csv")
write.csv(res, out_csv, row.names = FALSE)
message("wrote ", out_csv)

message("\n== scorecard (power / safety_baseline) ==")
for (m in unique(res$method)) {
  cat("\n--- ", m, "---\n")
  wide <- reshape(res[res$method == m,
    c("pair", "variant", "power", "safety_baseline", "label_frac")],
    idvar = "pair", timevar = "variant", direction = "wide")
  print(wide, row.names = FALSE)
}
message("BENCH PANCREAS DONE")
