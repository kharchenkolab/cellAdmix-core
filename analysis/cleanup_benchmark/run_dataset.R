# Shared cleanup-benchmark driver. A per-dataset config script builds the
# dataset object and calls bench_run(); see run_pancreas.R for an example.
#
# Arms evaluated:
#  - variants x methods x seeds (bare methods)
#  - consensus: per variant x method, the union of native-check-passing
#    (source, target) decisions across seeds, imported into each seed's fit
#    subject to that fit's own native check
#  - extra_fits: pre-built runs (e.g. anchor-based fixed-H fits), by label

suppressMessages(library(Matrix))

`%||%` <- function(a, b) if (is.null(a)) b else a

bench_run <- function(ds, cell_annotation, label,
                      methods = c("membrane", "bridge"),
                      variants = c("ls_nmf", "invsqrt_kl"),
                      seeds = NULL,
                      consensus = FALSE,
                      extra_fits = NULL,
                      p_thresh = 0.1,
                      out_dir = NULL) {
  out_dir <- out_dir %||% file.path(Sys.getenv("BENCH_DIR", "analysis/cleanup_benchmark"), "results")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  message("== ", label, ": shared setup")
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

  profiles_iso <- sapply(types, function(t) {
    members <- intersect(names(cell_types)[!is.na(cell_types) & cell_types == t],
      rownames(exp15$counts))
    nonself <- rowSums(exp15$counts[members, setdiff(types, t), drop = FALSE])
    iso <- members[nonself == 0]
    if (length(iso) < 50) iso <- members
    bench_pseudobulk(counts_before, iso)
  })
  rownames(profiles_iso) <- rownames(counts_before)

  top_type <- colnames(profiles)[max.col(profiles, ties.method = "first")]
  native_markers <- lapply(types, function(t) {
    cand <- rownames(profiles)[top_type == t & profiles[, t] >= 200]
    max_other <- apply(profiles[cand, setdiff(types, t), drop = FALSE], 1, max)
    spec <- profiles[cand, t] / pmax(profiles[cand, t] + max_other, 1e-9)
    head(cand[order(spec, decreasing = TRUE)], 10)
  })
  names(native_markers) <- types

  mcount <- function(m, genes, cs) {
    genes <- intersect(genes, rownames(m))
    cs <- intersect(cs, colnames(m))
    Matrix::colSums(m[genes, cs, drop = FALSE])
  }

  message("== ", label, ": enumerating pairs")
  pair_defs <- list()
  for (S in types) {
    for (T_type in setdiff(types, S)) {
      T_cells <- intersect(
        names(cell_types)[!is.na(cell_types) & cell_types == T_type],
        intersect(rownames(exp15$counts), colnames(counts_before)))
      if (length(T_cells) < 200) next
      expo <- exp15$counts[T_cells, S]
      e0 <- T_cells[expo == 0]
      if (length(e0) < 100) next
      baseline_T0 <- bench_pseudobulk(counts_before, e0)
      pool <- bench_marker_pool(profiles, S, T_type, baseline_T0)
      if (length(pool) < 3) next
      strict <- bench_pool_strict(pool, profiles, S, baseline_T0)
      bins <- bench_exposure_bins(expo)
      rates <- bench_bin_rates(mcount(counts_before, pool, T_cells),
        totals_before[T_cells], bins)
      det <- bench_detect(rates)
      rates_s <- bench_bin_rates(mcount(counts_before, strict, T_cells),
        totals_before[T_cells], bins)
      det_s <- bench_detect(rates_s)
      pair_defs[[paste(S, T_type, sep = " -> ")]] <- list(S = S, T_type = T_type,
        T_cells = T_cells, expo = expo, bins = bins, pool = pool, strict = strict,
        p_detect = det$p, excess = det$excess_molecules, rates_before = rates,
        excess_strict = det_s$excess_molecules, rates_before_strict = rates_s)
    }
  }
  qvals <- p.adjust(vapply(pair_defs, function(d) d$p_detect, 0), "BH")
  keep <- names(pair_defs)[qvals < 0.01 &
    vapply(pair_defs, function(d) d$excess, 0) >= 200]
  message(sprintf("pairs considered: %d, detected: %d", length(pair_defs), length(keep)))
  pair_defs <- pair_defs[keep]

  results <- list()
  gene_results <- list()

  eval_arm <- function(fit, score, rules, counts_after, variant, method, seed, arm) {
    score_annotation <- score$annotation(p_thresh = p_thresh)
    source_calls <- score_annotation$source_calls
    factors_of_source <- function(S) {
      hits <- names(source_calls)[vapply(source_calls, function(v)
        identical(as.character(v), S), TRUE)]
      as.integer(sub("^f_", "", hits))
    }
    mols <- tryCatch(fit$molecules(sample_n = 3000000L), error = function(e) NULL)
    for (pname in names(pair_defs)) {
      d <- pair_defs[[pname]]
      exposed <- d$T_cells[d$expo > 0]
      e0 <- d$T_cells[d$expo == 0]
      rb <- d$rates_before
      ra <- bench_bin_rates(mcount(counts_after, d$pool, d$T_cells),
        totals_before[d$T_cells], d$bins)
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
      results[[length(results) + 1]] <<- data.frame(
        dataset = label, pair = pname, source = d$S, target = d$T_type,
        variant = variant, method = method, seed = seed, arm = arm,
        excess_molecules = round(d$excess),
        n_exposed = length(exposed), n_e0 = length(e0),
        power = round(bench_power(rb, ra), 3),
        power_strict = round(if (length(d$strict) >= 2) bench_power(d$rates_before_strict,
          bench_bin_rates(mcount(counts_after, d$strict, d$T_cells),
            totals_before[d$T_cells], d$bins)) else NA_real_, 3),
        excess_strict = round(d$excess_strict),
        n_strict_genes = length(d$strict),
        removal_depth = round(bench_removal_depth(rb, ra), 3),
        safety_native = round(bench_safety_native(
          mcount(counts_before, native_markers[[d$T_type]], d$T_cells),
          mcount(counts_after, native_markers[[d$T_type]], d$T_cells),
          totals_before[d$T_cells], d$bins), 3),
        integrity_native = round(bench_profile_integrity(counts_before, counts_after,
          e0, profiles = profiles, T_type = d$T_type), 4),
        nnls_before = round(bench_nnls_contamination(counts_before, profiles_iso,
          exposed, d$S, d$T_type), 4),
        nnls_after = round(bench_nnls_contamination(counts_after, profiles_iso,
          exposed, d$S, d$T_type), 4),
        marker_removed_frac = round(1 - mk_after / max(mk_before, 1), 3),
        rule_called = nrow(rule_rows) > 0,
        rule_kept = any(rule_rows$keep),
        label_frac = round(label_frac, 3),
        stringsAsFactors = FALSE)
      pg <- bench_power_per_gene(counts_before, counts_after, d$pool,
        d$T_cells, totals_before[d$T_cells], d$bins)
      pg$dataset <- label; pg$pair <- pname
      pg$variant <- variant; pg$method <- method; pg$seed <- seed; pg$arm <- arm
      gene_results[[length(gene_results) + 1]] <<- pg
    }
  }

  # ---- bare arms ----
  arm_ctx <- list()
  arms <- expand.grid(variant = variants, method = methods,
    seed = seeds %||% NA_integer_, stringsAsFactors = FALSE)
  for (ai in seq_len(nrow(arms))) {
    variant <- arms$variant[[ai]]; method <- arms$method[[ai]]; seed <- arms$seed[[ai]]
    arm <- paste0(variant, "/", method, if (!is.na(seed)) paste0("/s", seed))
    message("== ", label, " arm: ", arm)
    fit <- if (is.na(seed)) ds$fit(nmf_variant = variant, verbose = FALSE) else
      ds$fit(nmf_variant = variant, seed = seed, verbose = FALSE,
        run_id = sprintf("bench_seed%d_%s", seed, variant))
    score <- if (method == "membrane") fit$score_membrane() else fit$score_bridge()
    rules <- score$rules(p_thresh = p_thresh)
    corr_name <- paste0("cmp_", method, "_", variant, if (!is.na(seed)) paste0("_s", seed))
    correction <- score$correct(rules = rules, name = corr_name)
    eval_arm(fit, score, rules, correction$counts(), variant, method, seed, arm)
    arm_ctx[[arm]] <- list(fit = fit, score = score, rules = rules,
      variant = variant, method = method, seed = seed)
  }

  # ---- consensus arms ----
  if (isTRUE(consensus) && !is.null(seeds) && length(seeds) >= 2) {
    for (variant in variants) for (method in methods) {
      ctxs <- Filter(function(cx) cx$variant == variant && cx$method == method, arm_ctx)
      if (length(ctxs) < 2) next
      kept <- unique(unlist(lapply(ctxs, function(cx)
        paste(cx$rules$source_cell_type[cx$rules$keep],
          cx$rules$target_cell_type[cx$rules$keep], sep = "||"))))
      for (cx in ctxs) {
        arm <- paste0(variant, "/", method, "/cons_s", cx$seed)
        message("== ", label, " arm: ", arm)
        rules <- cx$rules
        ann <- cx$score$annotation(p_thresh = p_thresh)
        src <- ann$source_calls
        have <- unique(paste(rules$source_cell_type[rules$keep],
          rules$target_cell_type[rules$keep], sep = "||"))
        add <- list()
        for (cp in setdiff(kept, have)) {
          parts <- strsplit(cp, "\\|\\|")[[1]]
          fs <- as.integer(sub("^f_", "", names(src)[vapply(src, function(v)
            identical(as.character(v), parts[[1]]), TRUE)]))
          for (f in fs) {
            row <- rules[1, , drop = FALSE]
            row$factor <- f
            row$source_cell_type <- parts[[1]]
            row$target_cell_type <- parts[[2]]
            add[[length(add) + 1]] <- row
          }
        }
        if (length(add)) {
          imported <- cellAdmixCore:::.celladmix_apply_native_check(
            do.call(rbind, add), cx$fit)
          rules <- rbind(rules, imported[imported$keep, , drop = FALSE])
        }
        correction <- cx$score$correct(rules = rules,
          name = paste0("cmp_cons_", method, "_", variant, "_s", cx$seed))
        eval_arm(cx$fit, cx$score, rules, correction$counts(),
          variant, method, cx$seed, arm)
      }
    }
  }

  # ---- extra fits (e.g. anchor-based) ----
  for (nm in names(extra_fits %||% list())) {
    fit <- ds$read_fit(file.path(ds$prep$paths$runs_dir, extra_fits[[nm]]))
    for (method in methods) {
      arm <- paste0(nm, "/", method)
      message("== ", label, " arm: ", arm)
      score <- if (method == "membrane") fit$score_membrane() else fit$score_bridge()
      rules <- score$rules(p_thresh = p_thresh)
      correction <- score$correct(rules = rules, name = paste0("cmp_", nm, "_", method))
      eval_arm(fit, score, rules, correction$counts(), nm, method, NA_integer_, arm)
    }
  }

  res <- do.call(rbind, results)
  genes <- do.call(rbind, gene_results)
  write.csv(res, file.path(out_dir, paste0(label, "_pairs.csv")), row.names = FALSE)
  write.csv(genes, file.path(out_dir, paste0(label, "_genes.csv")), row.names = FALSE)

  message("\n== ", label, " scorecard ==")
  card <- do.call(rbind, lapply(split(res, res$arm), function(g) {
    data.frame(arm = g$arm[[1]],
      leakage_removed_overall = round(sum(g$power * g$excess_molecules) /
        sum(g$excess_molecules), 3),
      leakage_removed_strict = round(sum(g$power_strict * g$excess_strict, na.rm = TRUE) /
        sum(g$excess_strict[!is.na(g$power_strict)]), 3),
      pairs_ge_09 = sprintf("%d/%d", sum(g$power >= 0.9, na.rm = TRUE), nrow(g)),
      rules_kept = round(mean(g$rule_kept), 2),
      identity_retention_med = round(stats::median(g$safety_native, na.rm = TRUE), 3),
      identity_retention_min = round(min(g$safety_native, na.rm = TRUE), 3))
  }))
  print(card, row.names = FALSE)
  invisible(list(pairs = res, genes = genes, scorecard = card))
}
