#' High-Level cellAdmix Correction
#' @export
.celladmix_compact_correction_summary <- function(summary) {
  if (!is.data.frame(summary) || !nrow(summary)) {
    return(data.frame())
  }
  required <- c("cell_type", "n_molecules_before", "n_molecules_after",
    "n_removed", "fraction_removed")
  missing <- setdiff(required, names(summary))
  if (length(missing)) {
    stop("correction summary is missing columns: ", paste(missing, collapse = ", "))
  }
  df <- summary
  df$cell_type <- ifelse(is.na(df$cell_type) | !nzchar(df$cell_type),
    "unknown", as.character(df$cell_type))
  df$modified <- df$n_removed > 0

  summarize_group <- function(x, label) {
    modified <- x[x$modified, , drop = FALSE]
    data.frame(
      cell_type = label,
      n_cells = nrow(x),
      n_modified_cells = sum(x$modified, na.rm = TRUE),
      fraction_cells_modified = sum(x$modified, na.rm = TRUE) / max(1L, nrow(x)),
      molecules_before = sum(x$n_molecules_before, na.rm = TRUE),
      molecules_after = sum(x$n_molecules_after, na.rm = TRUE),
      molecules_removed = sum(x$n_removed, na.rm = TRUE),
      fraction_molecules_removed = sum(x$n_removed, na.rm = TRUE) /
        max(1, sum(x$n_molecules_before, na.rm = TRUE)),
      median_removed_per_modified_cell = if (nrow(modified)) {
        stats::median(modified$n_removed, na.rm = TRUE)
      } else 0,
      median_fraction_removed_per_modified_cell = if (nrow(modified)) {
        stats::median(modified$fraction_removed, na.rm = TRUE)
      } else 0,
      stringsAsFactors = FALSE
    )
  }

  pieces <- c(
    list(summarize_group(df, "all")),
    lapply(split(df, df$cell_type), function(x) summarize_group(x, unique(x$cell_type)[[1]]))
  )
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

CellAdmixCorrection <- R6::R6Class(
  "CellAdmixCorrection",
  public = list(
    name = NULL,
    fit = NULL,
    score_name = NULL,
    run = NULL,
    params = NULL,
    rules = NULL,

    initialize = function(name, fit, score, run, params = list(), rules = NULL) {
      self$name <- .celladmix_clean_name(name, "correction")
      self$fit <- fit
      self$score_name <- score$name
      self$run <- run
      self$params <- params
      self$rules <- rules
    },

    save_metadata = function() {
      .celladmix_write_json(self$metadata(), file.path(self$run$path, "correction_metadata.json"))
      invisible(self)
    },

    metadata = function() {
      list(name = self$name, score_name = self$score_name,
        run_dir = self$run$path, params = self$params)
    },

    summary = function() .celladmix_compact_correction_summary(self$cell_summary()),
    cell_summary = function() celladmix_collect_correction_summary(self$run),
    ensemble = function() {
      histogram <- self$run$vote_histogram
      if (length(histogram)) {
        names(histogram) <- paste0("votes_", seq_along(histogram))
      }
      list(
        members = self$params$ensemble %||% 1L,
        vote = self$params$vote %||% NA_real_,
        min_votes = self$run$min_votes %||% 1L,
        n_removed = self$run$n_removed %||% NA_real_,
        vote_histogram = histogram
      )
    },
    counts = function(...) celladmix_collect_counts_sparse(self$run, ...),
    cells = function() celladmix_collect_cells(self$run),
    molecules = function(...) celladmix_collect_transcripts(self$run, ...),
    cell_state_umap = function(annotation = NULL, cells_max = 5000L,
                               min_molecules = 10L, min_genes = 5L,
                               n_variable_genes = 1000L, pca_dims = 30L,
                               graph_k = 15L, cluster_resolution = 1,
                               compute_umap = TRUE, umap_neighbors = 15L,
                               umap_epochs = 200L,
                               umap_parallel_optimization = TRUE,
                               normalization_scale = 5000, seed = 1L,
                               verbose = FALSE, ...) {
      if (is.null(annotation)) {
        labels <- self$fit$dataset$annotation(as_vector = TRUE)
      } else if (is.character(annotation) && length(annotation) == 1L &&
          annotation %in% self$fit$dataset$annotations()$name) {
        labels <- self$fit$dataset$annotation(annotation, as_vector = TRUE)
      } else {
        labels <- annotation
      }
      dots <- .celladmix_dots_with_default_threads(list(...), self$fit$dataset$num_threads)
      result <- .celladmix_cluster_run_counts(
        self$run$path,
        min_molecules = min_molecules,
        min_genes = min_genes,
        cells_max = if (is.null(cells_max)) NA_integer_ else as.integer(cells_max),
        n_variable_genes = n_variable_genes,
        pca_dims = pca_dims,
        graph_k = graph_k,
        cluster_resolution = cluster_resolution,
        compute_umap = compute_umap,
        umap_neighbors = umap_neighbors,
        umap_epochs = umap_epochs,
        num_threads = dots$num_threads,
        umap_parallel_optimization = umap_parallel_optimization,
        normalization_scale = normalization_scale,
        seed = seed,
        verbose = verbose
      )
      .celladmix_clustering_frame(result, annotation = labels)
    },
    plot_removed_molecules = function(...) celladmix_plot_molecule_removal(self$run, ...)
  )
)
