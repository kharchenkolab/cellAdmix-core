#' @keywords internal
#' @noRd
.celladmix_domain_paths <- function(prep, domain_id) {
  if (!inherits(prep, "celladmix_prep")) {
    stop("Domain paths require a celladmix_prep object")
  }
  if (!is.character(domain_id) || length(domain_id) != 1L || !nzchar(domain_id)) {
    stop("domain_id must be a non-empty string")
  }
  if (grepl("[/\\\\]", domain_id)) {
    stop("domain_id must not contain path separators")
  }
  base_dir <- file.path(prep$paths$annotations_dir, domain_id)
  list(
    domain_id = domain_id,
    dir = base_dir,
    domains_path = file.path(base_dir, "cell_domains.csv"),
    domain_annotation_path = file.path(base_dir, "domain_annotation.csv"),
    composition_path = file.path(base_dir, "domain_composition.csv"),
    summary_path = file.path(base_dir, "domain_summary.json"),
    params_path = file.path(base_dir, "domain_params.json")
  )
}

#' @keywords internal
#' @noRd
.celladmix_normalize_domain_annotations <- function(
    annotations,
    annotation_col,
    cell_id_col = "cell_id"
  ) {
  if (is.null(annotations)) {
    return(NULL)
  }
  if (inherits(annotations, "celladmix_clusters")) {
    annotations <- annotations$clusters
  } else if (is.character(annotations) && length(annotations) == 1L && file.exists(annotations)) {
    annotations <- utils::read.csv(annotations, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(annotations)) {
    stop("annotations must be NULL, a data frame, a CSV path, or a celladmix_clusters object")
  }
  if (!(annotation_col %in% colnames(annotations))) {
    stop("annotations is missing annotation_col: ", annotation_col)
  }
  if (!(cell_id_col %in% colnames(annotations))) {
    if ("cell_id" %in% colnames(annotations)) {
      cell_id_col <- "cell_id"
    } else if ("cell" %in% colnames(annotations)) {
      cell_id_col <- "cell"
    } else if (!is.null(rownames(annotations)) &&
               !identical(rownames(annotations), as.character(seq_len(nrow(annotations))))) {
      annotations$cell_id <- rownames(annotations)
      cell_id_col <- "cell_id"
    } else {
      stop("annotations is missing cell id column: ", cell_id_col)
    }
  }
  cells <- as.character(annotations[[cell_id_col]])
  labels <- as.character(annotations[[annotation_col]])
  keep <- !is.na(cells) & nzchar(cells) & !is.na(labels) & nzchar(labels)
  cells <- cells[keep]
  labels <- labels[keep]
  if (!length(cells)) {
    stop("No usable cell annotations were found")
  }
  if (anyDuplicated(cells)) {
    keep_first <- !duplicated(cells)
    warning("Duplicate cell ids in annotations; using the first label for each cell")
    cells <- cells[keep_first]
    labels <- labels[keep_first]
  }
  stats::setNames(labels, cells)
}

#' Identify Spatial Domains from Cell Annotations
#'
#' Detects lightweight tissue domains from a prepared cellAdmix input source and
#' a cell-level annotation. Domains are inferred from multiscale local
#' annotation composition, followed by spatial Potts smoothing and connected
#' region extraction.
#'
#' @param prep A `celladmix_prep` object.
#' @param annotations Optional cell annotation data frame, CSV path, or
#'   `celladmix_clusters` object. If `NULL`, cell types already present in the
#'   input store are used.
#' @param annotation_col Column in `annotations` containing the biological
#'   annotation used to define local composition.
#' @param annotation_cell_id_col Cell id column in `annotations`.
#' @param domain_id Identifier used for persisted outputs under
#'   `output_dir/annotations/<domain_id>`.
#' @param clust Optional `celladmix_clusters` object retained in the returned
#'   object for UMAP plotting.
#' @param scales Integer kNN neighborhood sizes used for multiscale local
#'   annotation composition.
#' @param n_domains Number of hard domains to initialize with k-means.
#' @param transform Feature transform for each composition block.
#' @param include_self Whether to include the central cell's annotation in each
#'   local composition vector.
#' @param self_weight Weight assigned to the central cell when `include_self` is
#'   `TRUE`.
#' @param include_center_type Whether to append a down-weighted one-hot encoding
#'   of the central cell annotation.
#' @param center_type_weight Weight for the appended central annotation block.
#' @param smooth_k Number of spatial neighbors used for Potts smoothing and
#'   connected-region extraction.
#' @param smooth_lambda Smoothness penalty for neighbor label disagreement.
#' @param max_edge_distance Optional maximum physical edge length. `NULL` keeps
#'   all kNN edges.
#' @param min_region_size Connected components smaller than this are reassigned
#'   to the most common neighboring domain when possible.
#' @param force Logical; if `TRUE`, rebuild the input store if needed.
#' @param seed Integer RNG seed for k-means initialization.
#' @param num_threads Number of worker threads for native graph construction and
#'   k-means assignment.
#' @param verbose Logical; emit high-level timing messages.
#'
#' @return A `celladmix_domains` object.
#' @keywords internal
#' @noRd
celladmix_identify_domains <- function(
    prep,
    annotations = NULL,
    annotation_col = "merged_annotation",
    annotation_cell_id_col = "cell_id",
    domain_id = paste0(annotation_col, "_domains"),
    clust = NULL,
    scales = c(15L, 50L, 150L),
    n_domains = 8L,
    transform = c("clr", "hellinger", "none"),
    include_self = TRUE,
    self_weight = 1,
    include_center_type = FALSE,
    center_type_weight = 0.25,
    smooth_k = 20L,
    smooth_lambda = 0.5,
    max_edge_distance = NULL,
    min_region_size = 25L,
    force = FALSE,
    seed = 1L,
    num_threads = 1L,
    verbose = FALSE
  ) {
  if (!inherits(prep, "celladmix_prep")) {
    stop("celladmix_identify_domains() expects a celladmix_prep object")
  }
  transform <- match.arg(transform)
  if (!is.numeric(scales) || !length(scales) || any(scales <= 0 | is.na(scales))) {
    stop("scales must be positive integers")
  }
  if (!is.numeric(n_domains) || length(n_domains) != 1L || is.na(n_domains) || n_domains < 1L) {
    stop("n_domains must be a positive integer")
  }
  if (!is.numeric(smooth_k) || length(smooth_k) != 1L || is.na(smooth_k) || smooth_k < 1L) {
    stop("smooth_k must be a positive integer")
  }
  if (!is.null(max_edge_distance) &&
      (!is.numeric(max_edge_distance) || length(max_edge_distance) != 1L || is.na(max_edge_distance) || max_edge_distance <= 0)) {
    stop("max_edge_distance must be NULL or a positive number")
  }

  info <- .celladmix_info_logger(verbose)
  t_start <- Sys.time()
  info("Starting spatial domain detection")

  paths <- .celladmix_domain_paths(prep, domain_id)
  dir.create(paths$dir, recursive = TRUE, showWarnings = FALSE)

  cell_types <- .celladmix_normalize_domain_annotations(
    annotations,
    annotation_col = annotation_col,
    cell_id_col = annotation_cell_id_col
  )
  store <- .celladmix_as_store(
    prep,
    required = "counts",
    store_mode = "counts",
    force = force,
    num_threads = num_threads,
    verbose = verbose
  )

  raw <- .celladmix_identify_domains_store(
    store$path,
    cell_types = cell_types,
    scales = as.integer(scales),
    transform = transform,
    include_self = include_self,
    self_weight = self_weight,
    include_center_type = include_center_type,
    center_type_weight = center_type_weight,
    n_domains = as.integer(n_domains),
    smooth_k = as.integer(smooth_k),
    smooth_lambda = smooth_lambda,
    max_edge_distance = if (is.null(max_edge_distance)) NA_real_ else max_edge_distance,
    min_region_size = as.integer(min_region_size),
    seed = as.integer(seed),
    num_threads = as.integer(num_threads)
  )

  domains <- data.frame(
    cell_id = as.character(raw$cell_id),
    x = as.numeric(raw$x),
    y = as.numeric(raw$y),
    z = as.numeric(raw$z),
    sample_id = as.character(raw$sample_id),
    fov_id = as.character(raw$fov_id),
    annotation = as.character(raw$annotation),
    domain_index = as.integer(raw$domain),
    domain = paste0("D", as.integer(raw$domain)),
    region_index = as.integer(raw$region),
    region_id = paste0("R", as.integer(raw$region)),
    stringsAsFactors = FALSE
  )
  domains$domain_label <- domains$domain

  composition_matrix <- raw$domain_composition
  composition <- as.data.frame(as.table(composition_matrix), stringsAsFactors = FALSE)
  colnames(composition) <- c("domain", "annotation", "fraction")
  composition$n_cells <- as.integer(raw$domain_sizes[match(composition$domain, names(raw$domain_sizes))])

  params <- list(
    annotation_col = annotation_col,
    domain_id = domain_id,
    scales = as.integer(scales),
    n_domains = as.integer(n_domains),
    transform = transform,
    include_self = isTRUE(include_self),
    self_weight = as.numeric(self_weight),
    include_center_type = isTRUE(include_center_type),
    center_type_weight = as.numeric(center_type_weight),
    smooth_k = as.integer(smooth_k),
    smooth_lambda = as.numeric(smooth_lambda),
    max_edge_distance = max_edge_distance,
    min_region_size = as.integer(min_region_size),
    seed = as.integer(seed),
    num_threads = as.integer(num_threads)
  )
  diagnostics <- list(
    n_cells = nrow(domains),
    n_store_cells = as.integer(raw$n_store_cells %||% nrow(domains)),
    n_dropped_unannotated = as.integer(raw$n_dropped_unannotated %||% 0L),
    n_annotation_types = length(raw$cell_type_levels),
    n_domains = length(raw$domain_sizes),
    n_regions = length(raw$region_sizes),
    spatial_coherence = raw$spatial_coherence,
    boundary_fraction = raw$boundary_fraction,
    domain_sizes = raw$domain_sizes,
    region_sizes = raw$region_sizes,
    kmeans_iterations = raw$kmeans_iterations,
    smooth_iterations = raw$smooth_iterations,
    elapsed_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  )

  utils::write.csv(domains, paths$domains_path, row.names = FALSE)
  utils::write.csv(composition, paths$composition_path, row.names = FALSE)
  .celladmix_write_json(diagnostics, paths$summary_path)
  .celladmix_write_json(params, paths$params_path)

  info(sprintf(
    "Finished spatial domain detection: %d cells, %d domains, %d regions",
    diagnostics$n_cells, diagnostics$n_domains, diagnostics$n_regions
  ), t_start)

  structure(list(
    output_dir = prep$project_dir,
    domain_id = domain_id,
    domains = domains,
    composition = composition,
    diagnostics = diagnostics,
    params = params,
    clust = clust,
    paths = paths
  ), class = "celladmix_domains")
}

#' Load Persisted Spatial Domains
#'
#' @param prep A `celladmix_prep` object.
#' @param domain_id Domain identifier used in [celladmix_identify_domains()].
#'
#' @return A `celladmix_domains` object.
#' @keywords internal
#' @noRd
celladmix_load_domains <- function(prep, domain_id) {
  paths <- .celladmix_domain_paths(prep, domain_id)
  domains_path <- if (file.exists(paths$domain_annotation_path)) {
    paths$domain_annotation_path
  } else {
    paths$domains_path
  }
  if (!file.exists(domains_path)) {
    stop("No domain table found for domain_id: ", domain_id)
  }
  composition <- if (file.exists(paths$composition_path)) {
    utils::read.csv(paths$composition_path, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  diagnostics <- if (file.exists(paths$summary_path)) .celladmix_read_json(paths$summary_path) else list()
  params <- if (file.exists(paths$params_path)) .celladmix_read_json(paths$params_path) else list()
  structure(list(
    output_dir = prep$project_dir,
    domain_id = domain_id,
    domains = utils::read.csv(domains_path, stringsAsFactors = FALSE),
    composition = composition,
    diagnostics = diagnostics,
    params = params,
    paths = paths
  ), class = "celladmix_domains")
}

#' @keywords internal
#' @noRd
celladmix_domain_composition <- function(domains) {
  if (!inherits(domains, "celladmix_domains")) {
    stop("celladmix_domain_composition() expects a celladmix_domains object")
  }
  domains$composition
}

#' @keywords internal
#' @noRd
celladmix_domain_diagnostics <- function(domains) {
  if (!inherits(domains, "celladmix_domains")) {
    stop("celladmix_domain_diagnostics() expects a celladmix_domains object")
  }
  domains$diagnostics
}

#' @keywords internal
#' @noRd
.celladmix_domains_df <- function(domains) {
  if (inherits(domains, "celladmix_domains")) {
    return(domains$domains)
  }
  if (is.data.frame(domains)) {
    return(domains)
  }
  stop("Expected a celladmix_domains object or a domain data frame")
}

#' Plot Spatial Domain Assignments
#'
#' @param domains A `celladmix_domains` object or its `$domains` data frame.
#' @param color Column used for point color.
#' @param point_size Point size.
#' @param alpha Point alpha.
#'
#' @return A `ggplot2` object.
#' @keywords internal
#' @noRd
celladmix_plot_domains_spatial <- function(
    domains,
    color = "domain_label",
    point_size = 0.15,
    alpha = 0.7
  ) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for domain plotting")
  }
  df <- .celladmix_domains_df(domains)
  if (!(color %in% colnames(df))) {
    stop("Unknown color column: ", color)
  }
  ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, color = .data[[color]])) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::coord_equal() +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
    ggplot2::labs(title = "Spatial domains", x = "x", y = "y", color = color) +
    ggplot2::theme_classic(base_size = 10)
}

#' Plot Domain Composition
#'
#' @param domains A `celladmix_domains` object.
#'
#' @return A `ggplot2` object.
#' @keywords internal
#' @noRd
celladmix_plot_domain_composition <- function(domains) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for domain plotting")
  }
  comp <- celladmix_domain_composition(domains)
  if (!nrow(comp)) {
    stop("Domain composition table is empty")
  }
  ggplot2::ggplot(comp, ggplot2::aes(.data$domain, .data$annotation, fill = .data$fraction)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "grey95", high = "#B2182B", limits = c(0, 1)) +
    ggplot2::labs(title = "Domain annotation composition", x = "Domain", y = NULL, fill = "fraction") +
    ggplot2::theme_classic(base_size = 9)
}

#' Plot Domains on a Cell UMAP
#'
#' @param domains A `celladmix_domains` object.
#' @param clust A `celladmix_clusters` object. Defaults to `domains$clust` when
#'   available.
#' @param color Domain table column used for point color.
#'
#' @return A `ggplot2` object.
#' @keywords internal
#' @noRd
celladmix_plot_domain_umap <- function(domains, clust = NULL, color = "domain_label") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for domain plotting")
  }
  if (!inherits(domains, "celladmix_domains")) {
    stop("celladmix_plot_domain_umap() expects a celladmix_domains object")
  }
  clust <- clust %||% domains$clust
  if (!inherits(clust, "celladmix_clusters")) {
    stop("A celladmix_clusters object is required for UMAP plotting")
  }
  if (!(color %in% colnames(domains$domains))) {
    stop("Unknown color column: ", color)
  }
  emb <- merge(
    clust$embedding,
    domains$domains[, c("cell_id", color), drop = FALSE],
    by = "cell_id",
    all.x = FALSE,
    sort = FALSE
  )
  ggplot2::ggplot(emb, ggplot2::aes(.data$umap_1, .data$umap_2, color = .data[[color]])) +
    ggplot2::geom_point(size = 0.15, alpha = 0.7) +
    ggplot2::coord_equal() +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
    ggplot2::labs(title = "Cell UMAP colored by spatial domain", x = "Cell UMAP 1", y = "Cell UMAP 2", color = color) +
    ggplot2::theme_classic(base_size = 10)
}
