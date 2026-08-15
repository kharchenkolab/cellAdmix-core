#' High-Level cellAdmix Dataset
#'
#' R6 object storing source metadata, active annotations, fitted runs, scores,
#' and corrections. Construct with [cellAdmix()].
#' @export
CellAdmixDataset <- R6::R6Class(
  "CellAdmixDataset",
  public = list(
    prep = NULL,
    format = NULL,
    output_dir = NULL,
    source_info = NULL,
    active_annotation = NULL,
    active_domain = NULL,
    fits = NULL,
    num_threads = NULL,

    initialize = function(prep, format, source_info = NULL, annotation = NULL,
                          annotation_name = "manual", annotation_col = NULL,
                          cell_id_col = "cell_id", num_threads = NULL) {
      self$prep <- prep
      self$format <- format
      self$output_dir <- prep$project_dir
      self$source_info <- source_info %||% prep$source
      self$num_threads <- as.integer(num_threads %||% .celladmix_default_threads())
      self$fits <- list()
      private$.annotations <- list()
      private$.domains <- list()
      if (!is.null(annotation)) {
        self$set_annotation(
          annotation,
          name = annotation_name,
          annotation_col = annotation_col,
          cell_id_col = cell_id_col,
          active = TRUE
        )
      }
    },

    print = function(...) {
      file_label <- function(path) {
        if (is.null(path) || !length(path) || is.na(path[[1]]) || !nzchar(path[[1]])) {
          return(NULL)
        }
        basename(path[[1]])
      }
      cat("cellAdmix dataset\n")
      cat("  format:", self$format, "\n")
      cat("  output cache: configured\n")
      cat("  threads:", self$num_threads, "\n")
      probe <- self$prep$probe %||% list()
      if (!is.null(probe$transcripts$preferred_path)) {
        cat("  molecules:", file_label(probe$transcripts$preferred_path), "\n")
      } else if (!is.null(probe$molecules$path)) {
        cat("  molecules:", file_label(probe$molecules$path), "\n")
      }
      if (!is.null(probe$cells$preferred_path) && nzchar(probe$cells$preferred_path)) {
        cat("  cells:", file_label(probe$cells$preferred_path), "\n")
      }
      anns <- self$annotations()
      if (nrow(anns)) {
        active <- self$active_annotation %||% ""
        active_row <- match(active, anns$name)
        cat("  annotation:", active, " (", anns$n_cells[active_row],
            " cells, ", anns$n_labels[active_row], " labels)\n", sep = "")
      } else {
        cat("  annotation: none\n")
        cat("  note: provide annotation=, call $set_annotation(), or run $cluster()\n")
      }
      invisible(self)
    },

    summary = function() {
      list(
        format = self$format,
        output_dir = self$output_dir,
        source = self$source_info,
        annotations = self$annotations(),
        active_annotation = self$active_annotation,
        num_threads = self$num_threads,
        domains = self$domains(),
        active_domain = self$active_domain,
        fits = names(self$fits)
      )
    },

    annotations = function() {
      if (!length(private$.annotations)) {
        return(data.frame(name = character(), active = logical(), n_cells = integer(),
          n_labels = integer(), hash = character(), stringsAsFactors = FALSE))
      }
      out <- do.call(rbind, lapply(names(private$.annotations), function(name) {
        ann <- private$.annotations[[name]]
        data.frame(
          name = name,
          active = identical(name, self$active_annotation),
          n_cells = ann$n_cells,
          n_labels = ann$n_labels,
          hash = ann$hash,
          path = ann$path,
          stringsAsFactors = FALSE
        )
      }))
      rownames(out) <- NULL
      out
    },

    annotation = function(name = NULL, as_vector = FALSE) {
      name <- name %||% self$active_annotation
      if (is.null(name) || !(name %in% names(private$.annotations))) {
        stop("Unknown annotation: ", name %||% "<none>")
      }
      ann <- private$.annotations[[name]]
      if (isTRUE(as_vector)) ann$labels else ann$table
    },

    set_annotation = function(annotation, name = NULL, active = TRUE,
                              annotation_col = NULL, cell_id_col = "cell_id") {
      name <- .celladmix_clean_name(name %||% if (!length(private$.annotations)) "manual" else "annotation",
        "annotation")
      ann <- .celladmix_make_annotation(
        annotation,
        annotations_dir = self$prep$paths$annotations_dir,
        name = name,
        cell_id_col = cell_id_col,
        label_col = annotation_col
      )
      private$.annotations[[ann$name]] <- ann
      if (isTRUE(active) || is.null(self$active_annotation)) {
        self$active_annotation <- ann$name
      }
      invisible(ann)
    },

    use_annotation = function(name) {
      if (!(name %in% names(private$.annotations))) {
        stop("Unknown annotation: ", name)
      }
      self$active_annotation <- name
      invisible(self)
    },

    drop_annotation = function(name) {
      private$.annotations[[name]] <- NULL
      if (identical(self$active_annotation, name)) {
        self$active_annotation <- names(private$.annotations)[[1]] %||% NULL
      }
      invisible(self)
    },

    ensure_store = function(required = c("molecules", "counts"), ...) {
      required <- match.arg(required)
      dots <- .celladmix_dots_with_default_threads(list(...), self$num_threads)
      do.call(.celladmix_as_store, c(list(self$prep, required = required), dots))
    },

    store = function(required = c("molecules", "counts"), ...) {
      self$ensure_store(required = required, ...)
    },

    counts = function(...) {
      dots <- .celladmix_dots_with_default_threads(list(...), self$num_threads)
      do.call(celladmix_collect_store_counts, c(list(self$prep), dots))
    },

    cell_state_umap = function(annotation = NULL, cells_max = 5000L,
                               min_molecules = 10L, min_genes = 5L,
                               n_variable_genes = 1000L, pca_dims = 30L,
                               graph_k = 15L, cluster_resolution = 1,
                               compute_umap = TRUE, umap_neighbors = 15L,
                               umap_epochs = 200L,
                               umap_parallel_optimization = TRUE,
                               normalization_scale = 5000, seed = 1L,
                               verbose = FALSE, ...) {
      ann <- private$resolve_annotation(annotation, register = TRUE)
      dots <- .celladmix_dots_with_default_threads(list(...), self$num_threads)
      store <- self$ensure_store(required = "counts",
        store_mode = "counts", num_threads = dots$num_threads,
        verbose = verbose)
      result <- .celladmix_cluster_store(
        store$path,
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
        clusters_out = tempfile("celladmix_state_clusters_", fileext = ".parquet"),
        embedding_out = tempfile("celladmix_state_embedding_", fileext = ".parquet"),
        verbose = verbose
      )
      .celladmix_clustering_frame(result, annotation = ann$labels)
    },

    cluster = function(name = "louvain", active = NULL, ...) {
      dots <- .celladmix_dots_with_default_threads(list(...), self$num_threads)
      clust <- do.call(celladmix_cluster_cells, c(list(self$prep, cluster_id = name), dots))
      labels <- stats::setNames(paste0("cluster_", clust$clusters$cluster), clust$clusters$cell_id)
      should_activate <- if (is.null(active)) is.null(self$active_annotation) else isTRUE(active)
      self$set_annotation(labels, name = name, active = should_activate)
      clust
    },

    identify_domains = function(name = NULL, annotation = NULL, active = TRUE, ...) {
      ann <- private$resolve_annotation(annotation, register = TRUE)
      domain_id <- .celladmix_clean_name(name %||% paste0(ann$name, "_domains"), "domains")
      dots <- .celladmix_dots_with_default_threads(list(...), self$num_threads)
      domains <- do.call(celladmix_identify_domains, c(list(
        prep = self$prep,
        annotations = ann$table,
        annotation_col = "label",
        annotation_cell_id_col = "cell_id",
        domain_id = domain_id
      ), dots))
      private$.domains[[domain_id]] <- domains
      if (isTRUE(active) || is.null(self$active_domain)) {
        self$active_domain <- domain_id
      }
      domains
    },

    load_domains = function(name) {
      domain_id <- .celladmix_clean_name(name, "domains")
      domains <- celladmix_load_domains(self$prep, domain_id)
      private$.domains[[domain_id]] <- domains
      if (is.null(self$active_domain)) {
        self$active_domain <- domain_id
      }
      domains
    },

    domains = function(name = NULL) {
      if (is.null(name)) {
        if (!length(private$.domains)) {
          return(data.frame(name = character(), active = logical(),
            n_cells = integer(), n_domains = integer(), stringsAsFactors = FALSE))
        }
        out <- do.call(rbind, lapply(names(private$.domains), function(domain_id) {
          domains <- private$.domains[[domain_id]]
          data.frame(
            name = domain_id,
            active = identical(domain_id, self$active_domain),
            n_cells = nrow(domains$domains %||% data.frame()),
            n_domains = length(unique((domains$domains %||% data.frame())$domain_label)),
            stringsAsFactors = FALSE
          )
        }))
        rownames(out) <- NULL
        return(out)
      }
      domain_id <- .celladmix_clean_name(name, "domains")
      private$.domains[[domain_id]] %||% stop("Unknown domains: ", domain_id)
    },

    list_fits = function() {
      data.frame(
        run_id = names(self$fits),
        run_dir = vapply(self$fits, function(x) x$run_dir, character(1)),
        annotation = vapply(self$fits, function(x) x$annotation_name, character(1)),
        rank = vapply(self$fits, function(x) as.integer(x$rank), integer(1)),
        stringsAsFactors = FALSE
      )
    },

    get_fit = function(run_id) {
      self$fits[[run_id]] %||% stop("Unknown fit: ", run_id)
    },

    read_fit = function(run_dir, annotation = NULL, name = NULL) {
      ann <- private$resolve_annotation(annotation, register = TRUE)
      run <- celladmix_read_run(run_dir)
      rank <- as.integer(run$n_factors %||% nrow(run$h) %||% NA_integer_)
      fit <- CellAdmixFit$new(
        dataset = self,
        run = run,
        annotation = ann,
        rank = rank,
        params = list(read_existing = TRUE)
      )
      key <- .celladmix_clean_name(name %||% basename(run$path), "fit")
      self$fits[[key]] <- fit
      fit
    },

    fit = function(annotation = NULL, rank = "auto", rank_multiplier = 1.2,
                   rank_cap = 30L, run_id = NULL, overwrite = FALSE, ...) {
      ann <- private$resolve_annotation(annotation, register = TRUE)
      resolved_rank <- if (is.null(rank) || identical(rank, "auto")) {
        .celladmix_annotation_rank(ann, multiplier = rank_multiplier, cap = rank_cap)
      } else {
        as.integer(rank[[1]])
      }
      dots <- list(...)
      dots <- .celladmix_dots_with_default_threads(dots, self$num_threads)
      nmf_variant <- dots$nmf_variant %||% "invsqrt_kl"
      run_id <- run_id %||% .celladmix_run_id("fit", ann, resolved_rank, nmf_variant)
      run_dir <- file.path(self$prep$paths$runs_dir, run_id)
      if (dir.exists(run_dir) && file.exists(file.path(run_dir, "run.json")) && !isTRUE(overwrite)) {
        fit <- self$read_fit(run_dir, annotation = ann, name = run_id)
        fit$params <- c(fit$params, list(reused_existing = TRUE, requested_rank = resolved_rank))
        return(fit)
      }
      run <- do.call(celladmix_fit, c(
        list(
          prep = self$prep,
          training_labels = ann,
          rank = resolved_rank,
          run_id = run_id
        ),
        dots
      ))
      explicit_ncv_k <- dots$ncv_k
      if (!is.null(explicit_ncv_k) && is.finite(explicit_ncv_k[[1]]) &&
          explicit_ncv_k[[1]] > 0) {
        n_genes <- length(run$genes)
        if (n_genes > 100 * as.integer(explicit_ncv_k[[1]])) {
          warning(sprintf(paste0(
            "ncv_k=%d is small for a %d-gene panel: neighborhoods carry ",
            "almost no gene co-occurrence signal and NMF factors become ",
            "seed-dependent. Consider the automatic default (omit ncv_k; ",
            "~%d for this panel, subject to cell size) or gene subsetting."),
            as.integer(explicit_ncv_k[[1]]), n_genes,
            as.integer(round(20 * sqrt(n_genes / 400)))), call. = FALSE)
        }
      }
      fit <- CellAdmixFit$new(
        dataset = self,
        run = run,
        annotation = ann,
        rank = resolved_rank,
        params = c(list(rank_policy = "auto_annotation", rank_multiplier = rank_multiplier,
          rank_cap = rank_cap, overwrite = isTRUE(overwrite)), dots)
      )
      self$fits[[run_id]] <- fit
      fit$save_metadata()
      fit
    }
  ),
  private = list(
    .annotations = NULL,
    .domains = NULL,

    resolve_annotation = function(annotation = NULL, register = FALSE) {
      if (is.null(annotation)) {
        if (is.null(self$active_annotation)) {
          stop("No active annotation is available. Provide annotation=, call $set_annotation(), or run $cluster().")
        }
        return(private$.annotations[[self$active_annotation]])
      }
      if (inherits(annotation, "celladmix_annotation")) {
        return(annotation)
      }
      if (is.character(annotation) && length(annotation) == 1L &&
          annotation %in% names(private$.annotations)) {
        return(private$.annotations[[annotation]])
      }
      if (!isTRUE(register)) {
        return(.celladmix_make_annotation(
          annotation,
          annotations_dir = self$prep$paths$annotations_dir,
          name = "override"
        ))
      }
      name <- if (is.character(annotation) && length(annotation) == 1L && !file.exists(annotation)) {
        annotation
      } else {
        paste0("annotation_", length(private$.annotations) + 1L)
      }
      self$set_annotation(annotation, name = name, active = FALSE)
      private$.annotations[[.celladmix_clean_name(name, "annotation")]]
    }
  )
)
