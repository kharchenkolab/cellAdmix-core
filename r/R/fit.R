#' High-Level cellAdmix Fit
#' @export
CellAdmixFit <- R6::R6Class(
  "CellAdmixFit",
  public = list(
    dataset = NULL,
    run = NULL,
    run_id = NULL,
    run_dir = NULL,
    annotation_name = NULL,
    annotation_hash = NULL,
    rank = NULL,
    scores_registry = NULL,
    corrections_registry = NULL,
    params = NULL,

    initialize = function(dataset, run, annotation, rank, params = list()) {
      self$dataset <- dataset
      self$run <- run
      self$run_dir <- run$path
      self$run_id <- basename(run$path)
      self$annotation_name <- annotation$name
      self$annotation_hash <- annotation$hash
      self$rank <- rank
      self$params <- params
      self$scores_registry <- list()
      self$corrections_registry <- list()
    },

    print = function(...) {
      cat("cellAdmix fit\n")
      cat("  run:", self$run_id, "\n")
      cat("  annotation:", self$annotation_name, "\n")
      cat("  rank:", self$rank, "\n")
      opts <- self$run$pipeline_options %||% list()
      if (!is.null(opts$nmf_variant)) {
        cat("  NMF:", opts$nmf_variant, "\n")
      }
      if (!is.null(opts$molecule_scoring)) {
        cat("  molecule scoring:", opts$molecule_scoring, "\n")
      }
      cat("  cells:", self$run$n_cells %||% NA_integer_, "\n")
      cat("  transcripts:", self$run$n_transcripts %||% NA_integer_, "\n")
      invisible(self)
    },

    save_metadata = function() {
      .celladmix_write_json(self$metadata(), file.path(self$run_dir, "fit_metadata.json"))
      invisible(self)
    },

    metadata = function() {
      list(
        run_id = self$run_id,
        run_dir = self$run_dir,
        annotation_name = self$annotation_name,
        annotation_hash = self$annotation_hash,
        rank = self$rank,
        params = self$params,
        scores = names(self$scores_registry),
        corrections = names(self$corrections_registry)
      )
    },

    summary = function() {
      opts <- self$run$pipeline_options %||% list()
      data.frame(
        run_id = self$run_id,
        annotation = self$annotation_name,
        rank = self$rank,
        nmf_variant = opts$nmf_variant %||% NA_character_,
        molecule_scoring = opts$molecule_scoring %||% NA_character_,
        nmf_n_runs = opts$nmf_n_runs %||% NA_integer_,
        training_rows = self$run$n_training_rows %||% NA_integer_,
        cells = self$run$n_cells %||% NA_integer_,
        transcripts = self$run$n_transcripts %||% NA_integer_,
        stringsAsFactors = FALSE
      )
    },
    loadings = function() self$run$h,
    cell_factors = function() celladmix_collect_cells(self$run),
    counts = function(...) celladmix_collect_counts_sparse(self$run, ...),
    molecules = function(...) celladmix_collect_transcripts(self$run, ...),
    training_molecules = function() celladmix_collect_training_molecules(self$run),
    region = function(bbox, ...) celladmix_collect_region(self$run, bbox = bbox, ...),
    report_data = function(...) celladmix_report_data(self$run, ...),
    stain = function(stain = "membrane", ...) {
      celladmix_discover_stain_image(self$run, stain = stain, ...)
    },
    stain_image = function(...) self$stain(...),
    stains = function(names = c("dapi", "membrane")) {
      .celladmix_discover_stain_images(self$run, names)
    },
    stain_crop = function(image, bbox, ...) celladmix_read_stain_crop(image, bbox = bbox, ...),
    cell_boundaries = function(boundary_path = NULL, cells = NULL, bbox = NULL) {
      path <- celladmix_discover_cell_boundaries(self$run, boundary_path = boundary_path)
      celladmix_read_cell_boundaries(path, cells = cells, bbox = bbox)
    },
    prepare_cell_example = function(example, ...) {
      celladmix_prepare_cell_example(self, example, ...)
    },
    plot_cell_example = function(example, ...) {
      celladmix_plot_cell_example(example, fit = self, ...)
    },
    stability = function(...) celladmix_nmf_stability_data(self$run, ...),
    score_factor_sources = function(name = "factor_sources", annotation = NULL,
                                    counts = NULL, ...) {
      ann <- private$score_annotation(annotation)
      result <- celladmix_score_factor_sources(
        self$run,
        annotation = ann,
        counts = counts,
        out_dir = .celladmix_score_dir(self$run, name),
        ...
      )
      CellAdmixFactorSourceScore$new(name, self, result, ann, list(...))
    },

    score = function(method = c("membrane", "bridge", "coherence"),
                     name = NULL, annotation = NULL, ...) {
      method <- match.arg(method)
      name <- name %||% method
      if (identical(method, "membrane")) {
        self$score_membrane(name = name, annotation = annotation, ...)
      } else if (identical(method, "bridge")) {
        self$score_bridge(name = name, annotation = annotation, ...)
      } else {
        self$score_coherence(name = name, annotation = annotation, ...)
      }
    },

    score_membrane = function(name = "membrane", annotation = NULL, ...) {
      ann <- private$score_annotation(annotation)
      dots <- .celladmix_dots_with_default_threads(list(...), self$dataset$num_threads)
      result <- do.call(celladmix_score_membrane, c(list(self$run, annotation = ann$labels), dots))
      result <- .celladmix_copy_score_outputs(
        result,
        .celladmix_score_dir(self$run, name),
        method = "membrane"
      )
      score <- CellAdmixScore$new(name, "membrane", self, result, ann, dots)
      self$scores_registry[[score$name]] <- score
      score$save_metadata()
      score
    },

    score_bridge = function(name = "bridge", annotation = NULL, ...) {
      ann <- private$score_annotation(annotation)
      dots <- .celladmix_dots_with_default_threads(list(...), self$dataset$num_threads)
      result <- do.call(celladmix_score_bridge, c(list(self$run, annotation = ann$labels), dots))
      result <- .celladmix_copy_score_outputs(
        result,
        .celladmix_score_dir(self$run, name),
        method = "bridge"
      )
      score <- CellAdmixScore$new(name, "bridge", self, result, ann, dots)
      self$scores_registry[[score$name]] <- score
      score$save_metadata()
      score
    },

    score_coherence = function(name = "coherence", annotation = NULL, ...) {
      ann <- private$score_annotation(annotation)
      dots <- .celladmix_dots_with_default_threads(list(...), self$dataset$num_threads)
      result <- do.call(celladmix_score_coherence, c(list(self$run, annotation = ann$labels), dots))
      result <- .celladmix_copy_score_outputs(
        result,
        .celladmix_score_dir(self$run, name),
        method = "coherence"
      )
      score <- CellAdmixScore$new(name, "coherence", self, result, ann, dots)
      self$scores_registry[[score$name]] <- score
      score$save_metadata()
      score
    },

    scores = function() {
      data.frame(
        name = names(self$scores_registry),
        method = vapply(self$scores_registry, function(x) x$method, character(1)),
        annotation = vapply(self$scores_registry, function(x) x$annotation_name, character(1)),
        stringsAsFactors = FALSE
      )
    },

    get_score = function(name) {
      self$scores_registry[[name]] %||% stop("Unknown score: ", name)
    },

    correct = function(score, name = NULL, p_thresh = 0.1, adjust_p = FALSE,
                       targets = NULL, rules = NULL, ...) {
      score_obj <- if (inherits(score, "CellAdmixScore")) {
        score
      } else if (is.character(score) && length(score) == 1L) {
        self$get_score(score)
      } else {
        stop("score must be a CellAdmixScore or registered score name")
      }
      if (is.null(rules)) {
        rules <- score_obj$rules(p_thresh = p_thresh, adjust_p = adjust_p, targets = targets)
      }
      if (!is.null(rules) && "keep" %in% names(rules)) {
        # Rules flagged by the native-factor check are excluded from correction.
        rules <- rules[is.na(rules$keep) | rules$keep, , drop = FALSE]
      }
      name <- .celladmix_clean_name(name %||% paste0(score_obj$name, "_clean"), "correction")
      out_dir <- file.path(self$run$paths$corrected_dir, name)
      corrected <- celladmix_correct(
        self$run,
        rules = rules,
        out_dir = out_dir,
        annotation = score_obj$annotation_vector()
      )
      correction <- CellAdmixCorrection$new(name, self, score_obj, corrected,
        params = c(list(p_thresh = p_thresh, adjust_p = adjust_p, targets = targets), list(...)))
      self$corrections_registry[[name]] <- correction
      correction$save_metadata()
      correction
    },

    corrections = function() {
      data.frame(
        name = names(self$corrections_registry),
        score = vapply(self$corrections_registry, function(x) x$score_name, character(1)),
        stringsAsFactors = FALSE
      )
    },

    get_correction = function(name) {
      self$corrections_registry[[name]] %||% stop("Unknown correction: ", name)
    },

    plot_loadings = function(n_genes = 10L, ncol = 5L, ...) {
      h <- self$run$h
      if (is.null(h) || !length(h)) {
        stop("Run does not contain factor loadings")
      }
      n_genes <- max(1L, as.integer(n_genes))
      if (!is.null(self$run$genes)) {
        colnames(h) <- self$run$genes
      }
      if (is.null(colnames(h))) {
        colnames(h) <- paste0("gene_", seq_len(ncol(h)))
      }
      if (requireNamespace("ggplot2", quietly = TRUE) &&
          requireNamespace("cowplot", quietly = TRUE)) {
        plots <- lapply(seq_len(nrow(h)), function(factor) {
          vals <- h[factor, ]
          denom <- sum(vals)
          frac <- if (is.finite(denom) && denom > 0) vals / denom else vals
          idx <- head(order(frac, decreasing = TRUE), n_genes)
          df <- data.frame(gene = colnames(h)[idx], loading = frac[idx],
            stringsAsFactors = FALSE)
          df$gene <- factor(df$gene, levels = rev(df$gene))
          breaks <- round(seq(0, max(df$loading), length.out = 3), 2)
          ggplot2::ggplot(df, ggplot2::aes(x = loading, y = gene)) +
            ggplot2::geom_col(fill = "grey25", width = 0.75) +
            ggplot2::scale_x_continuous(breaks = breaks, position = "top",
              expand = ggplot2::expansion(mult = c(0, 0.04))) +
            ggplot2::labs(x = "Loading fraction", y = paste0("Factor ", factor)) +
            ggplot2::theme_bw(base_size = 9) +
            ggplot2::theme(
              axis.title.x = ggplot2::element_text(size = 7),
              axis.title.y = ggplot2::element_text(size = 9, face = "bold"),
              axis.text.x = ggplot2::element_text(size = 6.5),
              axis.text.y = ggplot2::element_text(size = 8),
              panel.grid.major.y = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              plot.margin = ggplot2::margin(2, 5, 2, 2)
            )
        })
        return(cowplot::plot_grid(plotlist = plots,
          nrow = ceiling(length(plots) / max(1L, as.integer(ncol))), align = "hv"))
      }
      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)
      n_factors <- nrow(h)
      ncol_plot <- min(4L, n_factors)
      graphics::par(mfrow = c(ceiling(n_factors / ncol_plot), ncol_plot),
        mar = c(3.5, 5, 2, 0.5), mgp = c(2, 0.65, 0), cex = 2/3)
      for (factor in seq_len(n_factors)) {
        vals <- h[factor, ]
        idx <- head(order(vals, decreasing = TRUE), n_genes)
        graphics::barplot(rev(vals[idx]), horiz = TRUE, las = 1,
          names.arg = rev(colnames(h)[idx]), main = paste0("F", factor), ...)
      }
      invisible(self)
    },

    plot_stability = function(...) celladmix_plot_nmf_stability(self$run, ...)
  ),
  private = list(
    score_annotation = function(annotation = NULL) {
      ann <- self$dataset$.__enclos_env__$private$resolve_annotation(annotation, register = TRUE)
      if (!identical(ann$hash, self$annotation_hash)) {
        warning(
          "Scoring annotation differs from fit annotation; source/target labels use the scoring annotation.",
          call. = FALSE
        )
      }
      ann
    }
  )
)

#' High-Level cellAdmix Score
