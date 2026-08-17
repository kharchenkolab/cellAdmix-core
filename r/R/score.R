#' High-Level cellAdmix Score
#' @export
CellAdmixScore <- R6::R6Class(
  "CellAdmixScore",
  public = list(
    name = NULL,
    method = NULL,
    fit = NULL,
    result = NULL,
    annotation_name = NULL,
    annotation_hash = NULL,
    paths = NULL,
    params = NULL,
    member_rules_cache = NULL,

    initialize = function(name, method, fit, result, annotation, params = list()) {
      self$name <- .celladmix_clean_name(name, "score")
      self$method <- method
      self$fit <- fit
      self$result <- result
      self$annotation_name <- annotation$name
      self$annotation_hash <- annotation$hash
      self$paths <- result$paths %||% list()
      self$params <- params
      self$member_rules_cache <- list()
      private$.annotation <- annotation
    },

    save_metadata = function() {
      score_dir <- self$paths$score_dir %||% .celladmix_score_dir(self$fit$run, self$name)
      dir.create(score_dir, recursive = TRUE, showWarnings = FALSE)
      .celladmix_write_json(self$metadata(), file.path(score_dir, "score_metadata.json"))
      invisible(self)
    },

    metadata = function() {
      list(
        name = self$name,
        method = self$method,
        run_id = self$fit$run_id,
        annotation_name = self$annotation_name,
        annotation_hash = self$annotation_hash,
        paths = self$paths,
        params = self$params
      )
    },

    summary = function() self$result$summary,
    pairs = function() self$result$scores,
    annotation_vector = function() private$.annotation$labels,

    annotation = function(p_thresh = 0.1, adjust_p = FALSE, ...) {
      if (identical(self$method, "membrane")) {
        celladmix_membrane_annotation(self$result, p_thresh = p_thresh, adjust_p = adjust_p)
      } else if (identical(self$method, "coherence")) {
        celladmix_coherence_annotation(self$result, p_thresh = p_thresh, adjust_p = adjust_p, ...)
      } else {
        celladmix_bridge_annotation(self$result, p_thresh = p_thresh, adjust_p = adjust_p)
      }
    },

    rules = function(p_thresh = 0.1, adjust_p = FALSE, targets = NULL,
                     native_check = TRUE, native_median_thresh = 0.1,
                     native_expr_thresh = 0.05, native_outlier_min_frac = 0.1,
                     neighbor_k = 15L, ...) {
      rules <- if (identical(self$method, "membrane")) {
        celladmix_membrane_rules(self$result, p_thresh = p_thresh,
          adjust_p = adjust_p, target_cell_types = targets)
      } else if (identical(self$method, "coherence")) {
        celladmix_coherence_rules(self$result, p_thresh = p_thresh,
          adjust_p = adjust_p, target_cell_types = targets, ...)
      } else {
        celladmix_bridge_rules(self$result, p_thresh = p_thresh,
          adjust_p = adjust_p, target_cell_types = targets)
      }
      if (isTRUE(native_check)) {
        rules <- .celladmix_apply_native_check(rules, fit = self$fit,
          median_thresh = native_median_thresh,
          expr_thresh = native_expr_thresh,
          outlier_min_frac = native_outlier_min_frac,
          neighbor_k = neighbor_k)
      }
      rules
    },

    plot_heatmap = function(...) {
      if (identical(self$method, "membrane")) {
        celladmix_plot_membrane_heatmap(self$result, ...)
      } else if (identical(self$method, "coherence")) {
        celladmix_plot_coherence_heatmap(self$result, ...)
      } else {
        celladmix_plot_bridge_heatmap(self$result, ...)
      }
    },

    plot_pairs = function(...) celladmix_plot_score_pairs(self, ...),
    plot_score_pairs = function(...) self$plot_pairs(...),
    examples = function(...) celladmix_select_example_cells(self, ...),
    plot_example = function(example, p_thresh = 0.1, adjust_p = FALSE,
                            score_annotation = NULL, ...) {
      score_annotation <- score_annotation %||% self$annotation(p_thresh = p_thresh,
        adjust_p = adjust_p)
      if (is.character(example)) {
        example <- celladmix_select_example_cells(self,
          score_annotation = score_annotation, p_thresh = p_thresh,
          adjust_p = adjust_p, cells = example)
      }
      celladmix_plot_cell_example(example, fit = self$fit,
        score_annotation = score_annotation, ...)
    },
    plot_examples = function(examples = NULL, targets = NULL, n_per_target = 4,
                             p_thresh = 0.1, adjust_p = FALSE, rules = NULL,
                             score_annotation = NULL, cell_data = NULL,
                             cells = NULL, ncol = 2, ...) {
      score_annotation <- score_annotation %||% self$annotation(p_thresh = p_thresh,
        adjust_p = adjust_p)
      if (is.character(examples) && is.null(cells)) {
        cells <- examples
        examples <- NULL
      }
      if (is.null(examples)) {
        examples <- self$examples(rules = rules, score_annotation = score_annotation,
          cell_data = cell_data, targets = targets, n_per_target = n_per_target,
          cells = cells,
          p_thresh = p_thresh, adjust_p = adjust_p)
      }
      if (is.null(examples) || !nrow(examples)) {
        return(.celladmix_empty_plot("No example cells selected"))
      }
      plots <- lapply(seq_len(nrow(examples)), function(i) {
        self$plot_example(examples[i, , drop = FALSE],
          p_thresh = p_thresh, adjust_p = adjust_p,
          score_annotation = score_annotation, cell_data = cell_data, ...)
      })
      if (length(plots) == 1L) {
        return(plots[[1]])
      }
      if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("The cowplot package is required for plot_examples()")
      }
      cowplot::plot_grid(plotlist = plots, ncol = max(1L, as.integer(ncol)),
        align = "hv", axis = "tblr")
    },

    plot = function(...) self$plot_heatmap(...),

    correct = function(name = NULL, p_thresh = 0.1, adjust_p = FALSE,
                       targets = NULL, rules = NULL, ensemble = NULL,
                       vote = 0.3, ...) {
      self$fit$correct(
        self,
        name = name,
        p_thresh = p_thresh,
        adjust_p = adjust_p,
        targets = targets,
        rules = rules,
        ensemble = ensemble,
        vote = vote,
        ...
      )
    }
  ),
  private = list(
    .annotation = NULL
  )
)

#' Count annotated cell types among each cell's nearest neighbor cells
#' @keywords internal
#' @noRd
.celladmix_source_exposure_counts <- function(cells, annotation, neighbor_k = 15L) {
  types <- sort(unique(as.character(annotation[!is.na(annotation)])))
  codes <- match(as.character(annotation[as.character(cells$cell_id)]), types)
  codes[is.na(codes)] <- 0L
  counts <- .celladmix_cell_neighbor_type_counts(
    cells$x, cells$y, codes - 1L, length(types), neighbor_k)
  colnames(counts) <- types
  list(counts = counts, types = types)
}

#' Flag likely-native factor/target rules using source-distant cells
#'
#' For each rule (factor f, target T, source S), target cells with zero
#' S-type cells among their nearest neighbor cells form the source-distant
#' group. A rule is flagged `keep = FALSE` when the factor persists in that
#' group (median fraction above `median_thresh`), when its source-distant
#' expression rate is a cross-type outlier, when there is no positive
#' exposure gradient, or when the check cannot be evaluated safely.
#'
#' @keywords internal
#' @noRd
.celladmix_apply_native_check <- function(
    rules, fit,
    median_thresh = 0.1,
    expr_thresh = 0.05,
    outlier_min_frac = 0.1,
    neighbor_k = 15L,
    cell_factors = NULL
  ) {
  empty_cols <- function(df) {
    df$keep <- logical(0)
    df$native_check <- character(0)
    df$native_distant_n <- numeric(0)
    df$native_exposed_n <- numeric(0)
    df$native_distant_median <- numeric(0)
    df$native_exposure_gradient <- numeric(0)
    df$native_distant_expr_frac <- numeric(0)
    df
  }
  if (is.null(rules) || !nrow(rules)) {
    return(if (is.null(rules)) rules else empty_cols(rules))
  }
  annotation <- tryCatch(
    fit$dataset$annotation(fit$annotation_name, as_vector = TRUE),
    error = function(e) NULL)
  cells <- cell_factors %||% fit$cell_factors()
  if (is.null(annotation) || !nrow(cells)) {
    stop("native_check requires an annotation and cell factors")
  }
  exposure <- .celladmix_source_exposure_counts(cells, annotation, neighbor_k)
  counts <- exposure$counts
  types <- exposure$types
  cell_types <- as.character(annotation[as.character(cells$cell_id)])

  n <- nrow(rules)
  keep <- rep(TRUE, n)
  reason <- rep("pass", n)
  distant_n <- rep(NA_real_, n)
  exposed_n <- rep(NA_real_, n)
  distant_median <- rep(NA_real_, n)
  gradient <- rep(NA_real_, n)
  expr_frac <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    target <- as.character(rules$target_cell_type[[i]])
    source <- rules$source_cell_type[[i]]
    if (is.na(source) || !(as.character(source) %in% types)) {
      keep[i] <- FALSE; reason[i] <- "source_not_in_annotation"; next
    }
    if (!(target %in% types)) {
      keep[i] <- FALSE; reason[i] <- "target_not_in_annotation"; next
    }
    source <- as.character(source)
    fraction_col <- paste0("factor_", as.integer(rules$factor[[i]]), "_fraction")
    if (!(fraction_col %in% names(cells))) {
      keep[i] <- FALSE; reason[i] <- "missing_factor_fractions"; next
    }
    fractions <- as.numeric(cells[[fraction_col]])
    distant_mask <- counts[, source] == 0
    in_target <- !is.na(cell_types) & cell_types == target
    distant <- in_target & distant_mask
    exposed <- in_target & !distant_mask
    distant_n[i] <- sum(distant)
    exposed_n[i] <- sum(exposed)
    if (sum(distant) <= 1) {
      keep[i] <- FALSE; reason[i] <- "no_distant_cells"; next
    }
    if (sum(exposed) <= 1) {
      keep[i] <- FALSE; reason[i] <- "no_exposed_cells"; next
    }
    fr_d <- fractions[distant]
    fr_e <- fractions[exposed]
    distant_median[i] <- stats::median(fr_d)
    gradient[i] <- mean(fr_e) - mean(fr_d)
    expr_frac[i] <- mean(fr_d > expr_thresh)
    if (all(fr_d == 0) && all(fr_e == 0)) {
      keep[i] <- FALSE; reason[i] <- "no_expression"; next
    }
    rates <- vapply(setdiff(types, source), function(cell_type) {
      members <- !is.na(cell_types) & cell_types == cell_type & distant_mask
      if (sum(members) > 1) mean(fractions[members] > expr_thresh) else 0
    }, numeric(1))
    upper <- stats::quantile(rates, 0.75) + 1.5 * stats::IQR(rates)
    is_outlier <- expr_frac[i] > upper && expr_frac[i] > outlier_min_frac
    if (distant_median[i] > median_thresh) {
      keep[i] <- FALSE; reason[i] <- "native_median"
    } else if (is_outlier) {
      keep[i] <- FALSE; reason[i] <- "native_outlier"
    } else if (gradient[i] < 0) {
      keep[i] <- FALSE; reason[i] <- "no_gradient"
    }
  }

  rules$keep <- keep
  rules$native_check <- reason
  rules$native_distant_n <- distant_n
  rules$native_exposed_n <- exposed_n
  rules$native_distant_median <- distant_median
  rules$native_exposure_gradient <- gradient
  rules$native_distant_expr_frac <- expr_frac
  rules
}

#' Score one ensemble member with the same method and parameters
#' @noRd
.celladmix_score_member <- function(fit, score_obj, member) {
  fn <- switch(score_obj$method,
    membrane = celladmix_score_membrane,
    coherence = celladmix_score_coherence,
    celladmix_score_bridge)
  do.call(fn, c(
    list(fit$run,
      annotation = score_obj$annotation_vector(),
      ensemble_member = as.integer(member)),
    score_obj$params))
}

#' Cell factor fractions under one ensemble member's labeling
#' @noRd
.celladmix_member_cell_factors <- function(fit, member) {
  cells <- fit$cell_factors()
  member_fr <- .celladmix_ensemble_member_fractions(fit$run$path, member)
  idx <- match(as.character(cells$cell_id), as.character(member_fr$cell_id))
  for (k in seq_len(ncol(member_fr$fractions))) {
    cells[[paste0("factor_", k, "_fraction")]] <-
      as.numeric(member_fr$fractions[idx, k])
  }
  cells
}

#' Build per-member kept-rule tables for the molecule-vote ensemble
#'
#' Each member's rules are derived exactly like the primary member's: score
#' with the member's labeling, threshold, and vet with the native-factor
#' check computed on the member's own cell factor fractions. Returns the
#' stacked rule rows, their member indices, and the kept source->target
#' pair set per member (for the support diagnostic).
#' @noRd
.celladmix_ensemble_member_rules <- function(
    fit, score_obj, members, selected, primary_rules,
    p_thresh = 0.1, adjust_p = FALSE, targets = NULL,
    restrict_pairs = NULL) {
  pair_key <- function(df) paste(df$source_cell_type, df$target_cell_type, sep = "\r")
  cache_key <- function(r) paste(r, p_thresh, adjust_p,
    paste(targets %||% "", collapse = ","), sep = "|")
  # Member scoring dominates ensemble cost, so kept member rules are also
  # persisted in the run's scores directory, keyed by the member index and
  # every input that could change them; re-corrections and re-renders then
  # reuse them instead of re-scoring.
  disk_key <- function(r) paste(
    score_obj$method, r, p_thresh, adjust_p,
    paste(targets %||% "", collapse = ","),
    score_obj$annotation_hash,
    paste(deparse(score_obj$params), collapse = ""),
    sep = "|")
  disk_path <- function(r) file.path(
    fit$run$paths$scores_dir,
    sprintf("ensemble_rules_%s_m%d.rds", score_obj$method, r))
  frames <- list()
  frame_members <- integer(0)
  pair_sets <- list()
  for (r in members) {
    rr <- if (r == selected) {
      primary_rules
    } else {
      cached <- score_obj$member_rules_cache[[cache_key(r)]]
      if (is.null(cached) && file.exists(disk_path(r))) {
        stored <- tryCatch(readRDS(disk_path(r)), error = function(e) NULL)
        if (!is.null(stored) && identical(stored$key, disk_key(r))) {
          cached <- stored$rules
        }
      }
      if (!is.null(cached)) {
        score_obj$member_rules_cache[[cache_key(r)]] <- cached
        cached
      } else {
        result_r <- .celladmix_score_member(fit, score_obj, r)
        rules_r <- switch(score_obj$method,
          membrane = celladmix_membrane_rules(result_r, p_thresh = p_thresh,
            adjust_p = adjust_p, target_cell_types = targets),
          coherence = celladmix_coherence_rules(result_r, p_thresh = p_thresh,
            adjust_p = adjust_p, target_cell_types = targets),
          celladmix_bridge_rules(result_r, p_thresh = p_thresh,
            adjust_p = adjust_p, target_cell_types = targets))
        rules_r <- .celladmix_apply_native_check(rules_r, fit = fit,
          cell_factors = .celladmix_member_cell_factors(fit, r))
        rules_r <- rules_r[is.na(rules_r$keep) | rules_r$keep, , drop = FALSE]
        score_obj$member_rules_cache[[cache_key(r)]] <- rules_r
        saveRDS(list(key = disk_key(r), rules = rules_r), disk_path(r))
        rules_r
      }
    }
    if (!is.null(restrict_pairs) && nrow(rr)) {
      rr <- rr[pair_key(rr) %in% pair_key(restrict_pairs), , drop = FALSE]
    }
    pair_sets[[as.character(r)]] <- unique(pair_key(rr))
    if (nrow(rr)) {
      frames[[length(frames) + 1L]] <- data.frame(
        factor = as.integer(rr$factor),
        target_cell_type = as.character(rr$target_cell_type),
        stringsAsFactors = FALSE)
      frame_members <- c(frame_members, rep(as.integer(r), nrow(rr)))
    }
  }
  list(
    rules = if (length(frames)) do.call(rbind, frames) else
      data.frame(factor = integer(0), target_cell_type = character(0),
        stringsAsFactors = FALSE),
    member = frame_members,
    pair_sets = pair_sets)
}

#' Fraction of ensemble members keeping a rule for each pair
#' @noRd
.celladmix_rule_support <- function(rules, pair_sets) {
  if (is.null(rules) || !nrow(rules)) {
    return(numeric(0))
  }
  key <- paste(rules$source_cell_type, rules$target_cell_type, sep = "\r")
  n <- max(length(pair_sets), 1L)
  vapply(key, function(k) {
    sum(vapply(pair_sets, function(s) k %in% s, logical(1))) / n
  }, numeric(1), USE.NAMES = FALSE)
}
