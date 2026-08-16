# Admixture audit: exposure-gradient estimates of leaked molecules per
# cell-type pair, and verification of corrections against them.
#
# The measure exploits the spatial structure of segmentation-driven
# admixture: contamination of a target cell by a source cell type requires
# physical adjacency, so source-marker content in target cells rises with
# the number of source-type neighbor cells, while unexposed target cells
# form an internal negative control. See docs/benchmarks.md for the
# methodology and its validation.

.celladmix_audit_pseudobulk <- function(counts, cells, scale = 1e6) {
  cells <- intersect(cells, colnames(counts))
  if (!length(cells)) {
    return(stats::setNames(rep(NA_real_, nrow(counts)), rownames(counts)))
  }
  sums <- Matrix::rowSums(counts[, cells, drop = FALSE])
  sums / max(sum(sums), 1) * scale
}

.celladmix_audit_profiles <- function(counts, cell_types) {
  types <- sort(unique(cell_types[!is.na(cell_types)]))
  out <- sapply(types, function(t) {
    .celladmix_audit_pseudobulk(counts,
      names(cell_types)[!is.na(cell_types) & cell_types == t])
  })
  rownames(out) <- rownames(counts)
  out
}

# Source-marker panel for an ordered pair: among genes whose top expresser is
# the source, the top n by contrast against the unexposed-target baseline.
# Rank-based: an absolute baseline cutoff would be skewed by contamination.
.celladmix_audit_marker_pool <- function(profiles, source, baseline_cpm,
                                         n_pool = 20L, min_source_cpm = 50,
                                         baseline_floor = 20) {
  top_type <- colnames(profiles)[max.col(profiles, ties.method = "first")]
  cand <- rownames(profiles)[top_type == source & profiles[, source] >= min_source_cpm]
  contrast <- profiles[cand, source] / pmax(baseline_cpm[cand], baseline_floor)
  cand[order(contrast, profiles[cand, source], decreasing = TRUE)][
    seq_len(min(n_pool, length(cand)))]
}

.celladmix_audit_pool_strict <- function(pool, profiles, source, baseline_cpm,
                                         baseline_frac = 0.05) {
  pool[baseline_cpm[pool] < baseline_frac * profiles[pool, source]]
}

.celladmix_audit_bins <- function(exposure) {
  cut(exposure, breaks = c(-0.5, 0.5, 1.5, 2.5, Inf), labels = c("0", "1", "2", "3+"))
}

# Pooled panel counts and rates per exposure bin, with exact Poisson CIs.
.celladmix_audit_bin_rates <- function(marker_counts, totals, bins) {
  m <- tapply(marker_counts, bins, sum, default = 0)
  M <- tapply(totals, bins, sum, default = 0)
  n <- tapply(rep(1, length(bins)), bins, sum, default = 0)
  M_safe <- pmax(as.numeric(M), 1)
  data.frame(bin = names(m), markers = as.numeric(m), totals = as.numeric(M),
    n_cells = as.numeric(n),
    rate = as.numeric(m) / M_safe,
    rate_lo = stats::qgamma(0.025, pmax(as.numeric(m), 1e-9)) / M_safe,
    rate_hi = stats::qgamma(0.975, as.numeric(m) + 1) / M_safe,
    stringsAsFactors = FALSE)
}

# Estimated leakage: exceedance of exposed-bin rates over the reference rate.
.celladmix_audit_excess <- function(rates) {
  r0 <- rates$rate[rates$bin == "0"]
  exposed <- rates[rates$bin != "0", , drop = FALSE]
  if (!length(r0) || !nrow(exposed)) {
    return(list(excess = 0, p = 1))
  }
  excess <- sum(pmax(exposed$rate - r0, 0) * exposed$totals)
  expected <- r0 * sum(exposed$totals)
  m <- sum(exposed$markers)
  p <- if (expected <= 0 || m <= expected) 1 else
    stats::ppois(m - 1, expected, lower.tail = FALSE)
  list(excess = excess, p = p)
}

.celladmix_audit_power <- function(rates_before, rates_after) {
  r0b <- rates_before$rate[rates_before$bin == "0"]
  r0a <- rates_after$rate[rates_after$bin == "0"]
  eb <- pmax(rates_before$rate[rates_before$bin != "0"] - r0b, 0) *
    rates_before$totals[rates_before$bin != "0"]
  ea <- pmax(rates_after$rate[rates_after$bin != "0"] - r0a, 0) *
    rates_before$totals[rates_before$bin != "0"]
  if (sum(eb) <= 0) NA_real_ else 1 - sum(ea) / sum(eb)
}

.celladmix_audit_mcount <- function(m, genes, cs) {
  genes <- intersect(genes, rownames(m))
  cs <- intersect(cs, colnames(m))
  Matrix::colSums(m[genes, cs, drop = FALSE])
}

#' Audit Admixture by Spatial Exposure
#'
#' Estimates, for every ordered cell-type pair, the number of molecules
#' leaked from the source type into cells of the target type, using the
#' rise of source-marker content with source-neighbor exposure and the
#' unexposed target cells as an internal negative control. The estimates
#' are conservative lower bounds and require only the annotation, cell
#' positions, and counts (no factorization).
#'
#' @param fit A `CellAdmixFit`.
#' @param annotation Optional annotation name (defaults to the fit's).
#' @param neighbor_k Number of nearest cells defining exposure.
#' @param n_pool Maximum source-marker panel size per pair.
#' @param q_thresh BH-adjusted detection threshold for calling a pair.
#' @param min_excess Minimum estimated leaked molecules for a detected pair.
#' @param min_target_cells,min_reference_cells Minimum cells required in the
#'   target type overall and in its unexposed reference bin.
#' @return A [CellAdmixAudit] object.
#' @export
celladmix_audit_admixture <- function(fit, annotation = NULL, neighbor_k = 15L,
    n_pool = 20L, q_thresh = 0.01, min_excess = 200,
    min_target_cells = 200L, min_reference_cells = 100L) {
  if (!inherits(fit, "CellAdmixFit")) {
    stop("fit must be a CellAdmixFit")
  }
  ann <- tryCatch(
    fit$dataset$annotation(annotation %||% fit$annotation_name, as_vector = TRUE),
    error = function(e) NULL)
  if (is.null(ann)) {
    stop("audit_admixture requires a cell-type annotation")
  }
  cells <- fit$cell_factors()
  counts <- fit$counts()
  exposure <- .celladmix_source_exposure_counts(cells, ann, neighbor_k)
  rownames(exposure$counts) <- as.character(cells$cell_id)
  cell_types <- stats::setNames(as.character(ann[as.character(cells$cell_id)]),
    as.character(cells$cell_id))
  totals <- Matrix::colSums(counts)
  profiles <- .celladmix_audit_profiles(counts, cell_types)
  top_type <- colnames(profiles)[max.col(profiles, ties.method = "first")]
  native_markers <- lapply(exposure$types, function(t) {
    cand <- rownames(profiles)[top_type == t & profiles[, t] >= 200]
    others <- setdiff(exposure$types, t)
    max_other <- apply(profiles[cand, others, drop = FALSE], 1, max)
    spec <- profiles[cand, t] / pmax(profiles[cand, t] + max_other, 1e-9)
    utils::head(cand[order(spec, decreasing = TRUE)], 10)
  })
  names(native_markers) <- exposure$types

  pair_defs <- list()
  for (source in exposure$types) {
    for (target in setdiff(exposure$types, source)) {
      T_cells <- intersect(
        names(cell_types)[!is.na(cell_types) & cell_types == target],
        intersect(rownames(exposure$counts), colnames(counts)))
      if (length(T_cells) < min_target_cells) next
      expo <- exposure$counts[T_cells, source]
      if (sum(expo == 0) < min_reference_cells) next
      baseline <- .celladmix_audit_pseudobulk(counts, T_cells[expo == 0])
      pool <- .celladmix_audit_marker_pool(profiles, source, baseline, n_pool = n_pool)
      if (length(pool) < 3) next
      strict <- .celladmix_audit_pool_strict(pool, profiles, source, baseline)
      bins <- .celladmix_audit_bins(expo)
      rates <- .celladmix_audit_bin_rates(
        .celladmix_audit_mcount(counts, pool, T_cells), totals[T_cells], bins)
      det <- .celladmix_audit_excess(rates)
      rates_strict <- if (length(strict) >= 2) .celladmix_audit_bin_rates(
        .celladmix_audit_mcount(counts, strict, T_cells), totals[T_cells], bins) else NULL
      det_strict <- if (!is.null(rates_strict)) .celladmix_audit_excess(rates_strict) else
        list(excess = NA_real_, p = NA_real_)
      pair_defs[[paste(source, target, sep = " -> ")]] <- list(
        source = source, target = target, T_cells = T_cells, bins = bins,
        pool = pool, strict = strict, rates = rates, rates_strict = rates_strict,
        excess = det$excess, p = det$p, excess_strict = det_strict$excess)
    }
  }
  qvals <- stats::p.adjust(vapply(pair_defs, function(d) d$p, 0), "BH")
  for (i in seq_along(pair_defs)) {
    pair_defs[[i]]$q <- qvals[[i]]
    pair_defs[[i]]$detected <- qvals[[i]] < q_thresh && pair_defs[[i]]$excess >= min_excess
  }
  CellAdmixAudit$new(fit = fit, pair_defs = pair_defs, counts = counts,
    totals = totals, cell_types = cell_types, native_markers = native_markers,
    params = list(neighbor_k = neighbor_k, n_pool = n_pool, q_thresh = q_thresh,
      min_excess = min_excess))
}

#' Admixture Audit Result
#'
#' Created by [celladmix_audit_admixture()] (or `fit$audit_admixture()`).
#' Holds per-cell-type-pair leakage estimates, plotting methods, and the
#' `evaluate()` method that verifies a correction against the same
#' measurements.
#' @export
CellAdmixAudit <- R6::R6Class(
  "CellAdmixAudit",
  public = list(
    fit = NULL,
    params = NULL,

    initialize = function(fit, pair_defs, counts, totals, cell_types,
                          native_markers, params = list()) {
      self$fit <- fit
      self$params <- params
      private$.pair_defs <- pair_defs
      private$.counts <- counts
      private$.totals <- totals
      private$.cell_types <- cell_types
      private$.native_markers <- native_markers
    },

    pairs = function(detected_only = FALSE) {
      out <- do.call(rbind, lapply(private$.pair_defs, function(d) {
        exposed <- d$T_cells[d$bins != "0"]
        data.frame(source = d$source, target = d$target,
          excess = round(d$excess), excess_strict = round(d$excess_strict),
          pct_of_exposed = 100 * d$excess /
            max(sum(private$.totals[exposed]), 1),
          q_value = d$q, detected = d$detected,
          n_exposed = length(exposed), n_reference = sum(d$bins == "0"),
          n_markers = length(d$pool), n_strict = length(d$strict),
          stringsAsFactors = FALSE)
      }))
      rownames(out) <- NULL
      if (detected_only) out[out$detected, , drop = FALSE] else out
    },

    markers = function(source, target) {
      d <- private$.pair(source, target)
      list(pool = d$pool, strict = d$strict)
    },

    plot_map = function(value = c("percent", "molecules"), detected_only = TRUE) {
      .celladmix_require_ggplot2()
      value <- match.arg(value)
      df <- self$pairs()
      if (detected_only) df <- df[df$detected, , drop = FALSE]
      if (!nrow(df)) {
        return(.celladmix_empty_plot("No detected admixture pairs"))
      }
      df$fill <- if (value == "percent") df$pct_of_exposed else df$excess
      df$label <- ifelse(df$excess >= 1000,
        sprintf("%.0fk", df$excess / 1000), sprintf("%d", df$excess))
      ggplot2::ggplot(df, ggplot2::aes(x = target, y = source, fill = fill)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.4) +
        ggplot2::geom_text(ggplot2::aes(label = label), size = 2.8) +
        ggplot2::scale_fill_gradient(low = "#fff5eb", high = "#d94801",
          name = if (value == "percent") "% of exposed\ntarget molecules"
                 else "leaked molecules") +
        ggplot2::labs(x = "target cell type", y = "source cell type",
          title = "Estimated admixture by cell-type pair",
          subtitle = "cell text: estimated leaked molecules (conservative)") +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
          panel.grid = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(fill = NA, color = "grey35",
            linewidth = 0.35))
    },

    plot_exposure = function(source = NULL, target = NULL, correction = NULL,
                             strict = TRUE) {
      .celladmix_require_ggplot2()
      counts_after <- if (!is.null(correction)) correction$counts() else NULL
      if (is.null(source) != is.null(target)) {
        stop("Provide both source and target, or neither for the cumulative view")
      }
      build <- function(counts_mat, state) {
        if (is.null(source)) {
          m_tot <- NULL
          for (d in private$.pair_defs) {
            if (!d$detected) next
            genes <- if (strict && length(d$strict) >= 2) d$strict else d$pool
            mk <- .celladmix_audit_mcount(counts_mat, genes, d$T_cells)
            part <- data.frame(bin = d$bins, mk = mk,
              tot = private$.totals[d$T_cells])
            m_tot <- rbind(m_tot, part)
          }
          if (is.null(m_tot)) {
            stop("No detected pairs to pool")
          }
          rates <- .celladmix_audit_bin_rates(m_tot$mk, m_tot$tot, m_tot$bin)
        } else {
          d <- private$.pair(source, target)
          genes <- if (strict && length(d$strict) >= 2) d$strict else d$pool
          rates <- .celladmix_audit_bin_rates(
            .celladmix_audit_mcount(counts_mat, genes, d$T_cells),
            private$.totals[d$T_cells], d$bins)
        }
        rates$state <- state
        rates
      }
      df <- build(private$.counts, "before")
      if (!is.null(counts_after)) {
        df <- rbind(df, build(counts_after, "after cleanup"))
        df$state <- factor(df$state, levels = c("before", "after cleanup"))
      }
      df$bin <- factor(df$bin, levels = c("0", "1", "2", "3+"))
      r0 <- df$rate[df$bin == "0" & df$state == "before"]
      title <- if (is.null(source)) {
        "Admixture exposure profile (pooled over detected pairs)"
      } else {
        sprintf("%s → %s", source, target)
      }
      cols <- c("before" = "#c0392b", "after cleanup" = "#2980b9")
      ggplot2::ggplot(df, ggplot2::aes(x = bin, group = state)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = rate_lo * 1e3, ymax = rate_hi * 1e3,
          fill = state), alpha = 0.25) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = rate_lo * 1e3,
          ymax = rate_hi * 1e3, color = state), width = 0.12, linewidth = 0.4) +
        ggplot2::geom_line(ggplot2::aes(y = rate * 1e3, color = state)) +
        ggplot2::geom_point(ggplot2::aes(y = rate * 1e3, color = state), size = 1.8) +
        ggplot2::geom_hline(yintercept = r0 * 1e3, linetype = "dotted",
          color = "grey40") +
        ggplot2::scale_color_manual(values = cols, name = NULL) +
        ggplot2::scale_fill_manual(values = cols, name = NULL) +
        ggplot2::labs(x = "source-type cells among nearest neighbors",
          y = "source-marker rate (per 1,000 molecules)",
          title = title,
          subtitle = "error bars: 95% Poisson intervals (often narrower than the symbols); dotted line: unexposed reference") +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(legend.position = if (is.null(counts_after)) "none" else "bottom")
    },

    plot_remaining = function(corrections = list()) {
      .celladmix_require_ggplot2()
      if (length(corrections) && is.null(names(corrections))) {
        stop("corrections must be a named list, e.g. list(membrane = corr)")
      }
      state_excess <- function(counts_mat) {
        tot <- 0; var_tot <- 0
        for (d in private$.pair_defs) {
          if (!d$detected) next
          rates <- .celladmix_audit_bin_rates(
            .celladmix_audit_mcount(counts_mat, d$pool, d$T_cells),
            private$.totals[d$T_cells], d$bins)
          r0 <- rates$rate[rates$bin == "0"]
          exposed <- rates[rates$bin != "0", , drop = FALSE]
          m_ref <- rates$markers[rates$bin == "0"]
          M_ref <- max(rates$totals[rates$bin == "0"], 1)
          M_exp <- sum(exposed$totals)
          tot <- tot + sum(pmax(exposed$rate - r0, 0) * exposed$totals)
          var_tot <- var_tot + sum(exposed$markers) + (M_exp / M_ref)^2 * m_ref
        }
        c(excess = tot, sd = sqrt(var_tot))
      }
      states <- c(list(uncorrected = private$.counts),
        lapply(corrections, function(x) x$counts()))
      est <- t(vapply(states, state_excess, c(excess = 0, sd = 0)))
      # fixed denominator: all molecules before correction, so bars read as
      # the estimated admixed fraction of the dataset's molecules
      denom <- max(sum(private$.totals), 1)
      df <- data.frame(state = factor(rownames(est), levels = rownames(est)),
        y = est[, "excess"] / denom,
        lo = pmax(est[, "excess"] - 1.96 * est[, "sd"], 0) / denom,
        hi = (est[, "excess"] + 1.96 * est[, "sd"]) / denom)
      ggplot2::ggplot(df, ggplot2::aes(x = state, y = 100 * y)) +
        ggplot2::geom_col(fill = "#34495e", width = 0.6, alpha = 0.9) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = 100 * lo, ymax = 100 * hi),
          width = 0.15, linewidth = 0.4) +
        ggplot2::labs(x = NULL,
          y = "estimated admixture (% of all molecules)",
          title = "Remaining admixture by correction",
          subtitle = "conservative estimate over detected cell-type pairs; error bars: 95% intervals") +
        ggplot2::theme_classic(base_size = 10)
    },

    evaluate = function(correction, warn_uncovered = TRUE) {
      counts_after <- correction$counts()
      rows <- list()
      for (d in private$.pair_defs) {
        if (!d$detected) next
        rates_after <- .celladmix_audit_bin_rates(
          .celladmix_audit_mcount(counts_after, d$pool, d$T_cells),
          private$.totals[d$T_cells], d$bins)
        sens_strict <- if (!is.null(d$rates_strict)) {
          .celladmix_audit_power(d$rates_strict, .celladmix_audit_bin_rates(
            .celladmix_audit_mcount(counts_after, d$strict, d$T_cells),
            private$.totals[d$T_cells], d$bins))
        } else NA_real_
        rows[[length(rows) + 1]] <- data.frame(
          source = d$source, target = d$target,
          excess = round(d$excess), excess_strict = round(d$excess_strict),
          sensitivity = .celladmix_audit_power(d$rates, rates_after),
          sensitivity_strict = sens_strict,
          stringsAsFactors = FALSE)
      }
      pairs <- do.call(rbind, rows)
      rownames(pairs) <- NULL
      fr <- do.call(rbind, lapply(names(private$.native_markers), function(t) {
        tc <- names(private$.cell_types)[!is.na(private$.cell_types) &
          private$.cell_types == t]
        genes <- private$.native_markers[[t]]
        before <- sum(.celladmix_audit_mcount(private$.counts, genes, tc))
        after <- sum(.celladmix_audit_mcount(counts_after, genes, tc))
        data.frame(cell_type = t, own_marker_molecules = before,
          false_removal = 1 - after / max(before, 1), stringsAsFactors = FALSE)
      }))
      if (warn_uncovered && !is.null(correction$rules)) {
        covered <- paste(correction$rules$source_cell_type,
          correction$rules$target_cell_type)
        missed <- pairs[!(paste(pairs$source, pairs$target) %in% covered) &
          pairs$excess >= self$params$min_excess * 5, , drop = FALSE]
        for (i in seq_len(nrow(missed))) {
          warning(sprintf(paste0(
            "Detected ~%s leaked molecules from %s into %s, but no removal ",
            "rule covers this pair"),
            format(missed$excess[[i]], big.mark = ","),
            missed$source[[i]], missed$target[[i]]), call. = FALSE)
        }
      }
      CellAdmixCleanupReport$new(audit = self, correction = correction,
        pairs = pairs, false_removal = fr)
    }
  ),
  private = list(
    .pair_defs = NULL, .counts = NULL, .totals = NULL,
    .cell_types = NULL, .native_markers = NULL,
    .pair = function(source, target) {
      d <- private$.pair_defs[[paste(source, target, sep = " -> ")]]
      if (is.null(d)) {
        stop("No audited pair ", source, " -> ", target,
          " (see audit$pairs() for available pairs)")
      }
      d
    }
  )
)

#' Cleanup Verification Report
#'
#' Created by `audit$evaluate(correction)`: per-pair cleanup sensitivity
#' against the audit's leakage estimates, plus own-marker false-removal
#' rates per cell type.
#' @export
CellAdmixCleanupReport <- R6::R6Class(
  "CellAdmixCleanupReport",
  public = list(
    audit = NULL,
    correction = NULL,

    initialize = function(audit, correction, pairs, false_removal) {
      self$audit <- audit
      self$correction <- correction
      private$.pairs <- pairs
      private$.false_removal <- false_removal
    },

    pairs = function() private$.pairs,
    false_removal = function() private$.false_removal,

    summary = function() {
      p <- private$.pairs
      fr <- private$.false_removal
      total <- sum(p$excess)
      removed <- sum(p$excess * pmax(p$sensitivity, 0), na.rm = TRUE)
      list(
        detected_pairs = nrow(p),
        estimated_leaked_molecules = total,
        leakage_removed_overall = removed / max(total, 1),
        median_pair_sensitivity = stats::median(p$sensitivity, na.rm = TRUE),
        own_marker_false_removal = sum(fr$false_removal * fr$own_marker_molecules) /
          max(sum(fr$own_marker_molecules), 1),
        worst_false_removal = fr$cell_type[which.max(fr$false_removal)]
      )
    },

    plot_cleanup = function() {
      .celladmix_require_ggplot2()
      p <- private$.pairs
      p$pair <- stats::reorder(paste(p$source, "→", p$target), p$excess)
      s <- self$summary()
      ggplot2::ggplot(p, ggplot2::aes(y = pair)) +
        ggplot2::geom_col(ggplot2::aes(x = pmax(sensitivity, 0)),
          fill = "#34495e", alpha = 0.9) +
        ggplot2::geom_point(ggplot2::aes(x = 1.06, size = excess),
          color = "#e67e22", alpha = 0.85) +
        ggplot2::scale_size_area(max_size = 5, name = "leaked\nmolecules",
          labels = scales_comma) +
        ggplot2::coord_cartesian(xlim = c(0, 1.12)) +
        ggplot2::labs(x = "estimated cleanup sensitivity", y = NULL,
          title = "Cleanup verification by cell-type pair",
          subtitle = sprintf(
            "leakage removed overall: %.0f%%; own-marker false removal: %.2f%%",
            100 * s$leakage_removed_overall, 100 * s$own_marker_false_removal)) +
        ggplot2::theme_classic(base_size = 9)
    }
  ),
  private = list(.pairs = NULL, .false_removal = NULL)
)

scales_comma <- function(x) format(x, big.mark = ",", scientific = FALSE)
