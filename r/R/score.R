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

    rules = function(p_thresh = 0.1, adjust_p = FALSE, targets = NULL, ...) {
      if (identical(self$method, "membrane")) {
        celladmix_membrane_rules(self$result, p_thresh = p_thresh,
          adjust_p = adjust_p, target_cell_types = targets)
      } else if (identical(self$method, "coherence")) {
        celladmix_coherence_rules(self$result, p_thresh = p_thresh,
          adjust_p = adjust_p, target_cell_types = targets, ...)
      } else {
        celladmix_bridge_rules(self$result, p_thresh = p_thresh,
          adjust_p = adjust_p, target_cell_types = targets)
      }
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

#' High-Level cellAdmix Correction
