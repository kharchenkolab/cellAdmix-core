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

    initialize = function(name, method, fit, result, annotation, params = list()) {
      self$name <- .celladmix_clean_name(name, "score")
      self$method <- method
      self$fit <- fit
      self$result <- result
      self$annotation_name <- annotation$name
      self$annotation_hash <- annotation$hash
      self$paths <- result$paths %||% list()
      self$params <- params
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
                       targets = NULL, rules = NULL, ...) {
      self$fit$correct(
        self,
        name = name,
        p_thresh = p_thresh,
        adjust_p = adjust_p,
        targets = targets,
        rules = rules,
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
    neighbor_k = 15L
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
  cells <- fit$cell_factors()
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
