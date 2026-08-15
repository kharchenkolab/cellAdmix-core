# Public batch/project API.

#' Return the Installed cellAdmix Core Version
#'
#' Returns the version string reported by the compiled `cellAdmix-core`
#' backend.
#'
#' @return A length-one character vector.
#' @export
celladmix_core_version <- function() {
  .Call("_cellAdmixCore_celladmix_core_version", PACKAGE = "cellAdmixCore")
}

#' @keywords internal
#' @noRd
.celladmix_info_logger <- function(verbose) {
  start <- Sys.time()
  force(verbose)
  function(message, stage_start = NULL) {
    if (!isTRUE(verbose)) {
      return(invisible(NULL))
    }
    now <- Sys.time()
    elapsed <- as.numeric(difftime(now, start, units = "secs"))
    suffix <- ""
    if (!is.null(stage_start)) {
      stage_sec <- as.numeric(difftime(now, stage_start, units = "secs"))
      suffix <- sprintf(" (%.3fs)", stage_sec)
    }
    cat(sprintf("[INFO %s +%.3fs] %s%s\n", format(now, "%H:%M:%S"), elapsed, message, suffix))
    flush.console()
    invisible(NULL)
  }
}

#' Read a Xenium Experiment Manifest
#'
#' Reads a Xenium `experiment.xenium` manifest and returns the parsed metadata
#' without loading transcript or cell tables.
#'
#' @param path Path to an `experiment.xenium` file.
#'
#' @return A named list describing the Xenium experiment manifest.
#' @keywords internal
#' @noRd
celladmix_read_xenium_manifest <- function(path) {
  .Call("_cellAdmixCore_celladmix_read_xenium_manifest", path, PACKAGE = "cellAdmixCore")
}

#' Simulate a Small NSCLC-Like Test Dataset
#'
#' Internal synthetic-data helper used by tests and quick local debugging.
#' This is intentionally not part of the public workflow API.
#'
#' @param transcripts_per_cell Approximate number of transcripts to simulate
#'   per cell.
#' @param admixture_per_target_cell Approximate number of admixture transcripts
#'   to add to each target cell.
#' @param seed Integer RNG seed.
#'
#' @return A list with simulated transcript and cell tables plus marker-gene
#'   metadata.
#' @keywords internal
#' @noRd
celladmix_simulate_nsclc <- function(
    transcripts_per_cell = 80L,
    admixture_per_target_cell = 24L,
    seed = 1L
  ) {
  .Call(
    "_cellAdmixCore_celladmix_simulate_nsclc",
    as.integer(transcripts_per_cell),
    as.integer(admixture_per_target_cell),
    as.integer(seed),
    PACKAGE = "cellAdmixCore"
  )
}

#' @keywords internal
#' @noRd
.celladmix_prepare_warning_messages <- function(probe) {
  warnings <- character()
  if (!nzchar(probe$transcripts$preferred_path)) {
    warnings <- c(warnings, "No Xenium transcript file was found.")
  }
  if (!nzchar(probe$cells$preferred_path)) {
    warnings <- c(warnings, "No cell metadata file was found; clustering and training strata will rely on transcript-derived cells.")
  } else if (!("cell_type" %in% probe$cells$columns)) {
    warnings <- c(warnings, "The cell metadata file does not expose a cell_type column; clustering will likely be needed before fitting.")
  }
  warnings
}

#' Prepare a Xenium Output Directory for cellAdmix
#'
#' Performs a cheap Xenium preflight step: resolves the manifest, checks the
#' expected transcript and cell files, normalizes crop definitions, creates a
#' output directory, and returns a small prep object that can be passed to
#' later workflow steps.
#'
#' This step does not load the full transcript table into R memory.
#'
#' @param bundle_dir Path to a Xenium bundle directory or directly to an
#'   `experiment.xenium` manifest.
#' @param output_dir Optional output/cache directory for cellAdmix artifacts. If
#'   `NULL`, a temporary output directory is created.
#' @param analysis_bbox Optional convenience bounding box for a single crop,
#'   given as `c(xmin, xmax, ymin, ymax)`.
#' @param analysis_crops Optional crop specification. May be a data frame, a
#'   list, or a path to a CSV file understood by
#'   [utils::.celladmix_normalize_crops()].
#' @param cell_filter,gene_filter Optional cell and gene IDs to retain while
#'   building native Xenium input stores.
#' @param min_qv Optional minimum transcript QV threshold applied in later
#'   native load/fit steps.
#' @param keep_unassigned Logical; whether to keep unassigned transcripts in
#'   later native load/fit steps.
#' @param keep_non_gene Logical; whether Xenium control/codeword/non-gene
#'   features should be retained in the input store. The default `FALSE`
#'   keeps only biological gene features for downstream analysis.
#' @param read_cells Logical; whether later native steps should read Xenium cell
#'   metadata if present.
#' @param prefer_parquet Logical; whether later native steps should prefer
#'   Parquet files over CSV files when both are available.
#'
#' @return A `celladmix_prep` object describing the prepared project.
#' @keywords internal
#' @noRd
celladmix_prepare_xenium <- function(
    bundle_dir,
    output_dir = NULL,
    analysis_bbox = NULL,
    analysis_crops = NULL,
    cell_filter = NULL,
    gene_filter = NULL,
    min_qv = NULL,
    keep_unassigned = FALSE,
    keep_non_gene = FALSE,
    read_cells = TRUE,
    prefer_parquet = TRUE
  ) {
  manifest_path <- if (dir.exists(bundle_dir)) {
    file.path(bundle_dir, "experiment.xenium")
  } else {
    bundle_dir
  }
  if (!file.exists(manifest_path)) {
    stop("Could not find experiment.xenium at ", manifest_path)
  }

  if (!is.null(analysis_bbox) && is.null(analysis_crops)) {
    analysis_crops <- list(
      crop_id = "crop_1",
      xmin = analysis_bbox[[1]],
      xmax = analysis_bbox[[2]],
      ymin = analysis_bbox[[3]],
      ymax = analysis_bbox[[4]]
    )
  }
  crops <- .celladmix_normalize_crops(analysis_crops)

  output_dir <- output_dir %||% .celladmix_default_output_dir()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  annotations_dir <- file.path(output_dir, "annotations")
  runs_dir <- file.path(output_dir, "runs")
  dir.create(annotations_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(runs_dir, recursive = TRUE, showWarnings = FALSE)

  probe <- .celladmix_probe_xenium(manifest_path)
  warnings <- .celladmix_prepare_warning_messages(probe)
  prep <- .celladmix_make_prep_object(list(
    project_dir = output_dir,
    source = list(
      type = "xenium",
      bundle_dir = probe$bundle_dir,
      manifest_path = probe$manifest_path
    ),
    crops = crops,
    defaults = list(
      cell_filter = if (is.null(cell_filter)) character() else unique(as.character(cell_filter)),
      gene_filter = if (is.null(gene_filter)) character() else unique(as.character(gene_filter)),
      min_qv = if (is.null(min_qv)) NA_real_ else min_qv,
      keep_unassigned = keep_unassigned,
      keep_non_gene = keep_non_gene,
      read_cells = read_cells,
      prefer_parquet = prefer_parquet
    ),
    paths = list(
      annotations_dir = annotations_dir,
      runs_dir = runs_dir,
      project_json = file.path(output_dir, "prep.json")
    ),
    probe = probe,
    warnings = warnings
  ))
  .celladmix_write_json(prep, prep$paths$project_json)
  prep
}

#' Prepare a Generic Tabular Output Directory for cellAdmix
#'
#' Performs a cheap preflight step for a CSV, TSV, or Parquet molecules table
#' and returns a small prep object that can be passed to later workflow steps.
#'
#' Generic tabular input must provide molecule coordinates and genes, plus one
#' of:
#' - a molecule-level cell identifier column
#' - a labeled TIFF segmentation mask that can be sampled at molecule
#'   coordinates
#'
#' @param molecules_path Path to a CSV, TSV, gzip-compressed CSV/TSV, or Parquet
#'   molecules table.
#' @param output_dir Optional output/cache directory for cellAdmix artifacts. If
#'   `NULL`, a temporary output directory is created.
#' @param analysis_bbox Optional convenience bounding box for a single crop,
#'   given as `c(xmin, xmax, ymin, ymax)`.
#' @param analysis_crops Optional crop specification. May be a data frame, a
#'   list, or a path to a CSV file understood by
#'   [utils::.celladmix_normalize_crops()].
#' @param x_col,y_col,z_col Column names for spatial coordinates. `z_col` may be
#'   `NULL` for 2D datasets.
#' @param gene_col Column name containing the molecule gene identity.
#' @param qv_col Optional quality-value column name.
#' @param cell_id_col Optional explicit molecule-level cell identifier column.
#' @param segmentation_mask_path Optional labeled TIFF segmentation mask used to
#'   assign molecules to cells when `cell_id_col` is not supplied.
#' @param cell_type_col Optional molecule-level cell type column.
#' @param cell_metadata_path Optional small CSV/TSV/Parquet cell metadata
#'   sidecar. CSV/TSV sidecars may be gzip-compressed. Used to attach cell-type
#'   labels when labels are not repeated on molecule rows.
#' @param cell_metadata_cell_id_col Cell identifier column in
#'   `cell_metadata_path`.
#' @param cell_metadata_cell_type_col Optional cell-type label column in
#'   `cell_metadata_path`.
#' @param sample_id_col,fov_id_col Optional molecule-level sample/FOV columns.
#' @param sample_id,fov_id Optional constant sample/FOV values written onto all
#'   loaded molecules when corresponding columns are not supplied.
#' @param min_qv Optional minimum transcript QV threshold applied in later
#'   native load/fit steps.
#' @param keep_unassigned Logical; whether to keep molecules that lack a cell
#'   assignment from `cell_id_col` or the segmentation mask.
#'
#' @return A `celladmix_prep` object describing the prepared project.
#' @keywords internal
#' @noRd
celladmix_prepare_tabular <- function(
    molecules_path,
    output_dir = NULL,
    analysis_bbox = NULL,
    analysis_crops = NULL,
    x_col = "x",
    y_col = "y",
    z_col = NULL,
    gene_col = "gene",
    qv_col = NULL,
    cell_id_col = NULL,
    segmentation_mask_path = NULL,
    cell_type_col = NULL,
    cell_metadata_path = NULL,
    cell_metadata_cell_id_col = "cell",
    cell_metadata_cell_type_col = NULL,
    sample_id_col = NULL,
    fov_id_col = NULL,
    sample_id = NULL,
    fov_id = NULL,
    min_qv = NULL,
    keep_unassigned = FALSE
  ) {
  if (!file.exists(molecules_path)) {
    stop("Could not find molecules table at ", molecules_path)
  }
  if (!is.null(cell_id_col) && !is.null(segmentation_mask_path)) {
    stop("Specify either cell_id_col or segmentation_mask_path, not both")
  }
  if (is.null(cell_id_col) && is.null(segmentation_mask_path)) {
    stop("Generic tabular input requires either cell_id_col or segmentation_mask_path")
  }
  if (!is.null(segmentation_mask_path) && !file.exists(segmentation_mask_path)) {
    stop("Could not find segmentation mask at ", segmentation_mask_path)
  }
  if (!is.null(cell_metadata_path) && !file.exists(cell_metadata_path)) {
    stop("Could not find cell metadata sidecar at ", cell_metadata_path)
  }
  if (!is.null(cell_type_col) && !is.null(cell_metadata_cell_type_col)) {
    stop("Specify either molecule-level cell_type_col or cell_metadata_cell_type_col, not both")
  }
  if (!is.null(cell_metadata_cell_type_col) && is.null(cell_metadata_path)) {
    stop("cell_metadata_cell_type_col requires cell_metadata_path")
  }
  if (!is.null(sample_id_col) && !is.null(sample_id)) {
    stop("Specify either sample_id_col or sample_id, not both")
  }
  if (!is.null(fov_id_col) && !is.null(fov_id)) {
    stop("Specify either fov_id_col or fov_id, not both")
  }

  if (!is.null(analysis_bbox) && is.null(analysis_crops)) {
    analysis_crops <- list(
      crop_id = "crop_1",
      xmin = analysis_bbox[[1]],
      xmax = analysis_bbox[[2]],
      ymin = analysis_bbox[[3]],
      ymax = analysis_bbox[[4]]
    )
  }
  crops <- .celladmix_normalize_crops(analysis_crops)

  output_dir <- output_dir %||% .celladmix_default_output_dir()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  annotations_dir <- file.path(output_dir, "annotations")
  runs_dir <- file.path(output_dir, "runs")
  dir.create(annotations_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(runs_dir, recursive = TRUE, showWarnings = FALSE)

  probe <- .celladmix_probe_tabular(molecules_path, segmentation_mask_path = segmentation_mask_path)
  warnings <- character()
  if (!is.null(cell_type_col) && !(cell_type_col %in% probe$molecules$columns)) {
    warnings <- c(warnings, sprintf("Requested cell_type_col '%s' was not found; cell-type training will be unavailable.", cell_type_col))
  } else if (is.null(cell_type_col) && is.null(cell_metadata_cell_type_col)) {
    warnings <- c(warnings, "No cell_type_col was supplied; clustering fallback will likely be needed before fitting.")
  }

  prep <- .celladmix_make_prep_object(list(
    project_dir = output_dir,
    source = list(
      type = "tabular",
      molecules_path = normalizePath(molecules_path, winslash = "/", mustWork = FALSE),
      x_col = x_col,
      y_col = y_col,
      z_col = z_col,
      gene_col = gene_col,
      qv_col = qv_col,
      cell_id_col = cell_id_col,
      segmentation_mask_path = if (is.null(segmentation_mask_path)) NULL else normalizePath(segmentation_mask_path, winslash = "/", mustWork = FALSE),
      cell_type_col = cell_type_col,
      cell_metadata_path = if (is.null(cell_metadata_path)) NULL else normalizePath(cell_metadata_path, winslash = "/", mustWork = FALSE),
      cell_metadata_cell_id_col = cell_metadata_cell_id_col,
      cell_metadata_cell_type_col = cell_metadata_cell_type_col,
      sample_id_col = sample_id_col,
      fov_id_col = fov_id_col,
      sample_id = sample_id,
      fov_id = fov_id
    ),
    crops = crops,
    defaults = list(
      min_qv = if (is.null(min_qv)) NA_real_ else min_qv,
      keep_unassigned = keep_unassigned
    ),
    paths = list(
      annotations_dir = annotations_dir,
      runs_dir = runs_dir,
      project_json = file.path(output_dir, "prep.json")
    ),
    probe = probe,
    warnings = warnings
  ))
  .celladmix_write_json(prep, prep$paths$project_json)
  prep
}

#' Build or Reuse a Normalized Input Store
#'
#' Materializes a source-independent store under `output_dir/input_store`.
#' Clustering can use a counts-only store, while fitting requires molecule rows.
#'
#' @param prep A `celladmix_prep` object.
#' @param force Logical; if `TRUE`, rebuild the store even when a compatible
#'   store manifest already exists.
#' @param store_mode Store materialization mode. `"counts"` builds only
#'   cell-gene counts, `"full"` builds molecule rows and cell offsets, and
#'   `"auto"` resolves to the mode required by the caller.
#' @param materialize_molecules Logical; whether to store molecule-level rows.
#'   Deprecated in favor of `store_mode`; still accepted for compatibility.
#' @param num_threads Number of worker threads reserved for store construction.
#'   The current first implementation records the interface; source-specific
#'   streaming builders will use it for row-group/block parallelism.
#' @param parquet_row_group_size Target row-group size for generated Parquet.
#' @param verbose Logical; if `TRUE`, emit `[INFO]` timing messages.
#'
#' @return A lightweight `celladmix_store` handle.
#' @keywords internal
#' @noRd
celladmix_store <- function(
    prep,
    force = FALSE,
    materialize_molecules = NULL,
    store_mode = c("auto", "full", "counts"),
    num_threads = 1L,
    parquet_row_group_size = 65536L,
    verbose = FALSE
  ) {
  if (!inherits(prep, "celladmix_prep")) {
    stop("celladmix_store() expects a celladmix_prep object")
  }
  store_mode <- match.arg(store_mode)
  required <- if (identical(store_mode, "counts") || identical(materialize_molecules, FALSE)) "counts" else "molecules"
  resolved_store <- .celladmix_resolve_store_mode(
    store_mode = store_mode,
    materialize_molecules = materialize_molecules,
    required = required
  )
  store_dir <- file.path(prep$project_dir, "input_store")
  if (identical(prep$source$type, "xenium")) {
    store <- .celladmix_build_xenium_store(
      prep$source$manifest_path,
      crops = prep$crops,
      cell_filter = prep$defaults$cell_filter %||% NULL,
      gene_filter = prep$defaults$gene_filter %||% NULL,
      min_qv = prep$defaults$min_qv,
      keep_unassigned = prep$defaults$keep_unassigned,
      keep_non_gene = prep$defaults$keep_non_gene %||% FALSE,
      read_cells = prep$defaults$read_cells,
      prefer_parquet = prep$defaults$prefer_parquet,
      store_dir = store_dir,
      materialize_molecules = resolved_store$materialize_molecules,
      force = force,
      num_threads = num_threads,
      parquet_row_group_size = parquet_row_group_size,
      verbose = verbose
    )
  } else if (identical(prep$source$type, "tabular")) {
    store <- .celladmix_build_tabular_store(
      prep$source$molecules_path,
      crops = prep$crops,
      min_qv = prep$defaults$min_qv,
      keep_unassigned = prep$defaults$keep_unassigned,
      x_col = prep$source$x_col,
      y_col = prep$source$y_col,
      z_col = prep$source$z_col,
      gene_col = prep$source$gene_col,
      qv_col = prep$source$qv_col,
      cell_id_col = prep$source$cell_id_col,
      segmentation_mask_path = prep$source$segmentation_mask_path,
      cell_type_col = prep$source$cell_type_col,
      cell_metadata_path = prep$source$cell_metadata_path,
      cell_metadata_cell_id_col = prep$source$cell_metadata_cell_id_col,
      cell_metadata_cell_type_col = prep$source$cell_metadata_cell_type_col,
      sample_id_col = prep$source$sample_id_col,
      fov_id_col = prep$source$fov_id_col,
      sample_id = prep$source$sample_id,
      fov_id = prep$source$fov_id,
      store_dir = store_dir,
      materialize_molecules = resolved_store$materialize_molecules,
      force = force,
      num_threads = num_threads,
      parquet_row_group_size = parquet_row_group_size,
      verbose = verbose
    )
  } else {
    stop("Unsupported prep source type: ", prep$source$type)
  }
  store$store_mode <- resolved_store$store_mode
  store$num_threads <- as.integer(num_threads)
  structure(store, class = unique(c("celladmix_store", class(store))))
}

#' Collect Cell-by-Gene Counts from an Input Store
#'
#' Builds or reuses a counts-capable input store and returns the sparse
#' gene-by-cell count matrix. This is intended for notebook diagnostics and
#' pre-fit DE summaries.
#'
#' @param x A `celladmix_prep` or `celladmix_store` object.
#' @param force Logical; rebuild the counts store when `x` is a prep object.
#' @param num_threads Number of worker threads for store construction.
#' @param verbose Logical; if `TRUE`, emit store-construction timing messages.
#'
#' @return A `Matrix::dgCMatrix` with genes in rows and cells in columns. Cell
#'   types, when available, are attached as the `cell_type` attribute.
#' @keywords internal
#' @noRd
celladmix_collect_store_counts <- function(
    x,
    force = FALSE,
    num_threads = 1L,
    verbose = FALSE
  ) {
  store <- .celladmix_as_store(
    x,
    required = "counts",
    store_mode = "counts",
    force = force,
    num_threads = num_threads,
    verbose = verbose
  )
  raw <- .celladmix_collect_input_store_counts(store$path)
  .celladmix_sparse_counts_from_raw(raw)
}

.celladmix_sparse_counts_from_raw <- function(raw) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required for sparse count collection")
  }
  out <- methods::new(
    "dgCMatrix",
    p = as.integer(raw$p),
    i = as.integer(raw$i),
    x = as.numeric(raw$x),
    Dim = as.integer(c(length(raw$genes), length(raw$cells))),
    Dimnames = list(raw$genes, raw$cells)
  )
  attr(out, "cell_type") <- stats::setNames(raw$cell_type, raw$cells)
  out
}

#' @keywords internal
#' @noRd
.celladmix_resolve_store_mode <- function(
    store_mode = c("auto", "full", "counts"),
    materialize_molecules = NULL,
    required = c("molecules", "counts")
  ) {
  required <- match.arg(required)
  store_mode <- match.arg(store_mode[[1]], c("auto", "full", "counts"))

  if (!is.null(materialize_molecules)) {
    store_mode <- if (isTRUE(materialize_molecules)) "full" else "counts"
  } else if (identical(store_mode, "auto")) {
    store_mode <- if (identical(required, "molecules")) "full" else "counts"
  }

  if (identical(required, "molecules") && !identical(store_mode, "full")) {
    stop("This operation requires a full input store with molecule rows")
  }

  list(
    store_mode = store_mode,
    materialize_molecules = identical(store_mode, "full")
  )
}

#' @keywords internal
#' @noRd
.celladmix_check_store_capability <- function(store, required = c("counts", "molecules")) {
  required <- match.arg(required)
  if (!inherits(store, "celladmix_store")) {
    stop("Expected a celladmix_store object")
  }
  caps <- store$capabilities %||% list()
  if (identical(required, "counts") && !isTRUE(caps$has_cell_gene_counts)) {
    stop("Input store does not provide cell-gene counts")
  }
  if (identical(required, "molecules") && !isTRUE(caps$has_molecule_rows %||% caps$has_molecules)) {
    stop("Input store does not provide molecule rows")
  }
  invisible(store)
}

#' @keywords internal
#' @noRd
.celladmix_as_store <- function(
    x,
    required = c("counts", "molecules"),
    store_mode = c("auto", "full", "counts"),
    force = FALSE,
    num_threads = 1L,
    parquet_row_group_size = 65536L,
    verbose = FALSE
  ) {
  required <- match.arg(required)
  store_mode <- match.arg(store_mode[[1]], c("auto", "full", "counts"))
  if (inherits(x, "celladmix_store")) {
    .celladmix_check_store_capability(x, required = required)
    return(x)
  }
  if (!inherits(x, "celladmix_prep")) {
    stop("Expected a celladmix_prep or celladmix_store object")
  }
  resolved <- .celladmix_resolve_store_mode(store_mode = store_mode, required = required)
  celladmix_store(
    x,
    force = force,
    materialize_molecules = resolved$materialize_molecules,
    store_mode = resolved$store_mode,
    num_threads = num_threads,
    parquet_row_group_size = parquet_row_group_size,
    verbose = verbose
  )
}

#' Cluster Cells for Training-Strata Fallback
#'
#' Runs a cell-level clustering workflow on the prepared Xenium scope. The
#' result is written under the project directory and returned as a small object
#' suitable for plotting and for use as `training_labels` in [celladmix_fit()].
#'
#' The current pipeline performs sparse cell-by-gene counting, library
#' normalization, `log1p`, variable-gene selection, PCA, kNN graph
#' construction, plain Louvain clustering, and UMAP embedding.
#'
#' @param prep A `celladmix_prep` object created by
#'   [celladmix_prepare_xenium()] or [celladmix_prepare_tabular()].
#' @param cluster_id Character identifier used to name the persisted clustering
#'   outputs under the project directory.
#' @param min_molecules Minimum number of molecules required for a cell to be
#'   included in clustering.
#' @param min_genes Minimum number of detected non-control genes required for a
#'   cell to be included in clustering.
#' @param cells_max Optional maximum number of cells to include in clustering.
#'   If `NULL`, all eligible cells in scope are used.
#' @param n_variable_genes Number of variable genes to retain.
#' @param pca_dims Number of PCA dimensions to compute.
#' @param graph_k Number of nearest neighbors used for the cell graph.
#' @param resolution Louvain resolution parameter.
#' @param compute_umap Logical; whether to compute a cell UMAP embedding.
#' @param umap_neighbors Number of neighbors used by UMAP.
#' @param umap_epochs Number of UMAP optimization epochs.
#' @param num_threads Number of worker threads for parallelizable clustering
#'   steps.
#' @param store_mode Input-store mode for clustering. `"auto"` currently uses a
#'   counts store, `"counts"` forces count-only materialization, and `"full"`
#'   materializes molecule rows as well for downstream reuse.
#' @param umap_parallel_optimization Logical; if `TRUE`, allow parallel UMAP
#'   layout optimization in addition to parallel neighbor search. Defaults to
#'   `TRUE`.
#' @param normalization_scale Target transcript-count scale used for cell-level
#'   normalization before `log1p`.
#' @param seed Integer RNG seed.
#' @param verbose Logical; if `TRUE`, emit `[INFO]` timing messages around the
#'   clustering API call.
#'
#' @return A `celladmix_clusters` object with cluster assignments, UMAP
#'   embedding coordinates, and paths to persisted cluster outputs.
#' @keywords internal
#' @noRd
celladmix_cluster_cells <- function(
    prep,
    cluster_id = "default",
    min_molecules = 10L,
    min_genes = 5L,
    cells_max = NULL,
    n_variable_genes = 1000L,
    pca_dims = 30L,
    graph_k = 15L,
    resolution = 1,
    compute_umap = TRUE,
    umap_neighbors = 15L,
    umap_epochs = 200L,
    num_threads = 1L,
    store_mode = c("auto", "counts", "full"),
    umap_parallel_optimization = TRUE,
    normalization_scale = 5000,
    seed = 1L,
    verbose = FALSE
  ) {
  if (!inherits(prep, "celladmix_prep")) {
    stop("celladmix_cluster_cells() expects a celladmix_prep object")
  }

  info <- .celladmix_info_logger(verbose)
  info("Starting cell clustering")
  cluster_start <- Sys.time()
  store_mode <- match.arg(store_mode)

  paths <- .celladmix_cluster_paths(prep, cluster_id = cluster_id)
  dir.create(paths$dir, recursive = TRUE, showWarnings = FALSE)
  store <- .celladmix_as_store(
    prep,
    required = "counts",
    store_mode = store_mode,
    num_threads = num_threads,
    verbose = verbose
  )
  out <- .celladmix_cluster_store(
    store$path,
    min_molecules = min_molecules,
    min_genes = min_genes,
    cells_max = if (is.null(cells_max)) NA_integer_ else as.integer(cells_max),
    n_variable_genes = n_variable_genes,
    pca_dims = pca_dims,
    graph_k = graph_k,
    cluster_resolution = resolution,
    compute_umap = compute_umap,
    umap_neighbors = umap_neighbors,
    umap_epochs = umap_epochs,
    num_threads = num_threads,
    umap_parallel_optimization = umap_parallel_optimization,
    normalization_scale = normalization_scale,
    seed = seed,
    clusters_out = paths$clusters_path,
    embedding_out = paths$embedding_path,
    verbose = verbose
  )

  summary <- list(
    n_cells = out$n_cells,
    n_clusters = out$n_clusters,
    n_variable_genes = out$n_variable_genes,
    variable_genes = out$variable_genes,
    pca_variance_explained = out$pca_variance_explained,
    min_molecules = as.integer(min_molecules),
    min_genes = as.integer(min_genes),
    cells_max = if (is.null(cells_max)) NA_integer_ else as.integer(cells_max),
    compute_umap = isTRUE(compute_umap),
    num_threads = as.integer(num_threads),
    umap_parallel_optimization = isTRUE(umap_parallel_optimization),
    resolution = as.numeric(resolution)
  )
  .celladmix_write_json(summary, paths$summary_path)
  if (is.null(out$embedding_path)) {
    paths$embedding_path <- NULL
  }
  info(sprintf("Finished cell clustering: %d cells, %d clusters", out$n_cells, out$n_clusters), cluster_start)

  structure(list(
    output_dir = prep$project_dir,
    cluster_id = cluster_id,
    clusters = out$clusters,
    embedding = out$embedding,
    paths = paths,
    summary = summary
  ), class = "celladmix_clusters")
}

#' @keywords internal
#' @noRd
.celladmix_resolve_training_labels <- function(prep, training_labels) {
  if (inherits(training_labels, "celladmix_annotation")) {
    if (is.null(training_labels$path) || !file.exists(training_labels$path)) {
      stop("celladmix_annotation does not have a persisted label sidecar")
    }
    return(list(path = training_labels$path, use_cell_type_training = FALSE))
  }
  if (is.character(training_labels) && length(training_labels) == 1L && file.exists(training_labels)) {
    return(list(path = normalizePath(training_labels, winslash = "/", mustWork = TRUE), use_cell_type_training = FALSE))
  }
  if (inherits(training_labels, "celladmix_clusters")) {
    if (!identical(normalizePath(training_labels$output_dir, winslash = "/", mustWork = FALSE), prep$project_dir)) {
      stop("training_labels was created from a different output directory")
    }
    return(list(path = training_labels$paths$clusters_path, use_cell_type_training = FALSE))
  }

  if (is.null(training_labels) || identical(training_labels, "auto")) {
    default_paths <- .celladmix_cluster_paths(prep, cluster_id = "default")
    if (file.exists(default_paths$clusters_path)) {
      return(list(path = default_paths$clusters_path, use_cell_type_training = FALSE))
    }
    return(list(path = NULL, use_cell_type_training = TRUE))
  }
  if (identical(training_labels, "cluster")) {
    default_paths <- .celladmix_cluster_paths(prep, cluster_id = "default")
    if (!file.exists(default_paths$clusters_path)) {
      stop("No default cell clustering was found under ", default_paths$clusters_path)
    }
    return(list(path = default_paths$clusters_path, use_cell_type_training = FALSE))
  }
  if (identical(training_labels, "cell_type")) {
    return(list(path = NULL, use_cell_type_training = TRUE))
  }
  if (identical(training_labels, "none")) {
    return(list(path = NULL, use_cell_type_training = FALSE))
  }
  stop("Unsupported training_labels specification")
}

#' @keywords internal
#' @noRd
.celladmix_normalize_scope <- function(scope) {
  if (is.null(scope)) {
    return(NULL)
  }
  scope <- as.character(scope)
  scope <- unique(scope[!is.na(scope) & nzchar(scope)])
  if (!length(scope)) {
    return(NULL)
  }
  scope
}

#' @keywords internal
#' @noRd
.celladmix_resolve_nmf_n_runs <- function(nmf_n_runs, nmf_init, num_threads) {
  if (is.null(nmf_n_runs) || !length(nmf_n_runs) || is.na(nmf_n_runs[[1]])) {
    return(max(1L, as.integer(num_threads)))
  }
  nmf_n_runs <- as.integer(nmf_n_runs[[1]])
  if (is.na(nmf_n_runs) || nmf_n_runs < 1L) {
    stop("nmf_n_runs must be a positive integer or NA for the default policy")
  }
  nmf_n_runs
}

#' @keywords internal
#' @noRd
.celladmix_recommended_rank <- function(n_clusters, max_clusters = 20L, factor_multiplier = 1.2, factor_cap = 30L) {
  stopifnot(length(n_clusters) == 1L, is.finite(n_clusters), n_clusters > 0)
  capped_clusters <- min(as.integer(n_clusters), as.integer(max_clusters))
  min(as.integer(factor_cap), max(2L, ceiling(factor_multiplier * capped_clusters)))
}

#' @keywords internal
#' @noRd
.celladmix_resolve_rank <- function(prep, training_labels, rank, default_rank = 4L, max_clusters = 20L, factor_multiplier = 1.2, factor_cap = 30L) {
  if (!is.null(rank) && !is.na(rank)) {
    return(as.integer(rank))
  }

  n_clusters <- NA_integer_
  if (inherits(training_labels, "celladmix_annotation")) {
    n_clusters <- as.integer(training_labels$n_labels %||% NA_integer_)
  } else if (inherits(training_labels, "celladmix_clusters")) {
    n_clusters <- as.integer(training_labels$summary$n_clusters %||% NA_integer_)
  } else if (is.character(training_labels) && length(training_labels) == 1L && file.exists(training_labels)) {
    n_clusters <- NA_integer_
  } else if (is.null(training_labels) || identical(training_labels, "auto") || identical(training_labels, "cluster")) {
    cluster_paths <- .celladmix_cluster_paths(prep, cluster_id = "default")
    if (file.exists(cluster_paths$summary_path)) {
      cluster_summary <- .celladmix_read_json(cluster_paths$summary_path)
      n_clusters <- as.integer(cluster_summary$n_clusters %||% NA_integer_)
    }
  }

  if (is.na(n_clusters) || n_clusters <= 0L) {
    return(as.integer(default_rank))
  }
  resolved_rank <- .celladmix_recommended_rank(
    n_clusters,
    max_clusters = max_clusters,
    factor_multiplier = factor_multiplier,
    factor_cap = factor_cap
  )
  if (n_clusters > as.integer(max_clusters)) {
    warning(
      sprintf(
        "Cell clustering produced %d clusters; capping the automatically selected factor count at %d. Set rank explicitly to override.",
        n_clusters,
        resolved_rank
      ),
      call. = FALSE
    )
  }
  resolved_rank
}

#' @keywords internal
#' @noRd
.celladmix_fit_prep <- function(
    prep, training_labels = "auto", out_dir = NULL, run_id = NULL, analysis_crop = NULL,
    scope = NULL,
    ncv_k = 20L, rank = NULL, graph_k = 10L, same_label_ratio = 5, nmf_iterations = 150L,
    nmf_init = c("auto", "random", "cluster"),
    nmf_variant = c("invsqrt_kl", "kl", "sqrt_kl", "ls_nmf"),
    molecule_scoring = c("gene_loadings", "ncv_projection", "auto"),
    nmf_n_runs = NA_integer_, nmf_train_max_rows = 10000L,
    nmf_min_molecules = 10L, num_threads = 1L,
    tile_size = 100, parquet_row_group_size = 65536L, report_ncv_umap = FALSE,
    seed = 1L, verbose = FALSE) {
  if (!inherits(prep, "celladmix_prep")) {
    stop(".celladmix_fit_prep() expects a celladmix_prep object")
  }

  nmf_init <- match.arg(nmf_init)
  nmf_variant <- match.arg(nmf_variant)
  molecule_scoring <- match.arg(molecule_scoring)
  nmf_n_runs <- .celladmix_resolve_nmf_n_runs(nmf_n_runs, nmf_init, num_threads)
  scope <- .celladmix_normalize_scope(scope)
  training <- .celladmix_resolve_training_labels(prep, training_labels)
  rank <- .celladmix_resolve_rank(prep, training_labels, rank)
  out_dir <- out_dir %||% .celladmix_make_run_dir(prep, run_id = run_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  store <- .celladmix_as_store(prep,
    required = "molecules", store_mode = "auto", num_threads = num_threads,
    parquet_row_group_size = parquet_row_group_size, verbose = verbose)

  .celladmix_make_run_object(
    .celladmix_fit_store_run(
      store$path,
      analysis_crop = analysis_crop,
      ncv_k = ncv_k,
      rank = rank,
      graph_k = graph_k,
      same_label_ratio = same_label_ratio,
      nmf_iterations = nmf_iterations,
      nmf_init = nmf_init,
      nmf_variant = nmf_variant,
      molecule_scoring = molecule_scoring,
      nmf_n_runs = nmf_n_runs,
      nmf_train_max_rows = nmf_train_max_rows,
      nmf_min_molecules = nmf_min_molecules,
      num_threads = num_threads,
      training_labels_path = training$path,
      use_cell_type_training = training$use_cell_type_training,
      training_scope_cell_types = scope,
      seed = seed,
      out_dir = out_dir,
      tile_size = tile_size,
      parquet_row_group_size = parquet_row_group_size,
      report_ncv_umap = report_ncv_umap,
      verbose = verbose
    ),
    class = "celladmix_run"
  )
}

#' Fit cellAdmix on a Prepared Project
#'
#' Runs the core cellAdmix factorization and classification pipeline on a
#' prepared project scope and writes a persisted run directory.
#'
#' Training labels may be supplied explicitly from a
#' `celladmix_clusters` object, or resolved automatically from the project.
#'
#' @param prep A `celladmix_prep` object created by
#'   [celladmix_prepare_xenium()].
#' @param training_labels Training-strata specification. Supported values are a
#'   `celladmix_clusters` object, `"auto"`, `"cluster"`, `"cell_type"`, or
#'   `"none"`.
#' @param out_dir Optional run output directory. If `NULL`, a run directory is
#'   created under the prepared project.
#' @param run_id Optional run identifier used when creating a run directory
#'   under the prepared project.
#' @param analysis_crop Optional crop identifier selecting one crop from a
#'   multi-crop prepared scope.
#' @param scope Optional character vector of `cell_type` values used only for
#'   selecting NMF training rows. Factors are still projected to all cells in
#'   the analysis scope.
#' @param ncv_k Number of within-cell neighbors used for NCV construction.
#' @param rank Optional NMF rank. If `NULL`, the fit uses an automatically
#'   recommended value of about `1.2 * n_clusters` from the clustering result,
#'   capped as if there were at most `20` clusters, with a warning when clustering exceeds `20`
#'   clusters.
#' @param graph_k Number of neighbors used for within-cell smoothing.
#' @param same_label_ratio Potts-style smoothing weight.
#' @param nmf_iterations Number of NMF update iterations.
#' @param nmf_init NMF initialization mode. `"random"` uses the original random
#'   start. `"cluster"` seeds factor loadings from mean NCV profiles of the
#'   sampled training strata. `"auto"` uses cluster initialization when
#'   informative training strata are available and otherwise falls back to
#'   random initialization.
#' @param nmf_variant NMF matrix/objective variant. The default
#'   `"invsqrt_kl"` applies inverse-square-root gene-prevalence weighting,
#'   `"kl"` uses raw NCV counts, and
#'   `"sqrt_kl"` applies square-root NCV counts before KL-NMF. `"ls_nmf"` uses
#'   the legacy weighted least-squares NMF on untransformed reference-style NCV
#'   counts and gene-normalized loading potentials for molecule scoring.
#' @param molecule_scoring Molecule node-potential strategy before smoothing.
#'   The default `"gene_loadings"` assigns each molecule a factor potential
#'   from the observed gene's normalized factor loadings. `"ncv_projection"`
#'   projects each molecule's NCV row through fixed loadings. `"auto"` keeps the
#'   older mixed policy: `"gene_loadings"` for `"ls_nmf"` and
#'   `"ncv_projection"` for KL variants.
#' @param nmf_n_runs Number of independent NMF restarts. When `NA`, the default
#'   is `num_threads` restarts. When greater than `1`, the native fit runs multiple seeds,
#'   keeps the lowest-loss result, and returns loss/stability diagnostics for
#'   the candidates.
#' @param nmf_train_max_rows Optional cap on the number of NCV rows used for NMF
#'   training.
#' @param nmf_min_molecules Minimum number of molecules required for a cell to
#'   participate in NMF training.
#' @param num_threads Number of worker threads used for threaded parts of the
#'   native fit pipeline, including parallel NMF restarts.
#' @param tile_size Tile size stored in run metadata for downstream tiled access.
#' @param parquet_row_group_size Target Parquet row-group size for persisted run
#'   outputs.
#' @param report_ncv_umap Logical; if `TRUE`, build the developer-only NCV UMAP
#'   report sidecar during the fit. Leave `FALSE` for lean production runs and
#'   call [celladmix_report_data()] later when notebook/report data is needed.
#' @param seed Integer RNG seed.
#' @param verbose Logical; if `TRUE`, emit `[INFO]` fit-stage progress messages
#'   from the native backend.
#'
#' @return A persisted `celladmix_run` object.
#' @keywords internal
#' @noRd
celladmix_fit <- function(
    prep, training_labels = "auto", out_dir = NULL, run_id = NULL, analysis_crop = NULL,
    scope = NULL,
    ncv_k = 20L, rank = NULL, graph_k = 10L, same_label_ratio = 5,
    nmf_iterations = 150L, nmf_init = c("auto", "random", "cluster"),
    nmf_variant = c("invsqrt_kl", "kl", "sqrt_kl", "ls_nmf"),
    molecule_scoring = c("gene_loadings", "ncv_projection", "auto"), nmf_n_runs = NA_integer_,
    nmf_train_max_rows = 10000L, nmf_min_molecules = 10L, num_threads = 1L,
    tile_size = 100, parquet_row_group_size = 65536L, report_ncv_umap = FALSE, seed = 1L,
    verbose = FALSE) {
  nmf_init <- match.arg(nmf_init)
  nmf_variant <- match.arg(nmf_variant)
  molecule_scoring <- match.arg(molecule_scoring)
  .celladmix_fit_prep(prep, training_labels = training_labels, out_dir = out_dir,
    run_id = run_id, analysis_crop = analysis_crop, ncv_k = ncv_k, rank = rank,
    scope = scope, graph_k = graph_k, same_label_ratio = same_label_ratio, nmf_iterations = nmf_iterations,
    nmf_init = nmf_init,
    nmf_variant = nmf_variant,
    molecule_scoring = molecule_scoring,
    nmf_n_runs = nmf_n_runs, nmf_train_max_rows = nmf_train_max_rows, nmf_min_molecules = nmf_min_molecules,
    num_threads = num_threads, tile_size = tile_size,
    parquet_row_group_size = parquet_row_group_size, report_ncv_umap = report_ncv_umap,
    seed = seed, verbose = verbose)
}

#' Read a Persisted cellAdmix Run
#'
#' Reads a previously written run directory and returns a small run object that
#' can be used with the downstream collection, scoring, and correction helpers.
#'
#' @param path Path to a persisted cellAdmix run directory.
#'
#' @return A `celladmix_run` object.
#' @keywords internal
#' @noRd
celladmix_read_run <- function(path) {
  .celladmix_make_run_object(.celladmix_read_run_manifest(path), class = "celladmix_run")
}

#' Build NMF Restart Stability Diagnostics
#'
#' Returns one row per selected factor when a run was fitted with multiple NMF
#' restarts. Stability is the mean Hungarian-matched correlation of each
#' selected factor's gene-ownership profile (loadings normalized per gene
#' across factors) against the independent random restarts; unrelated factors
#' score near zero, and values above ~0.3 indicate a reproducibly re-found
#' factor.
#' Importance defaults to the fraction of assigned molecules represented by the
#' factor, estimated from per-cell factor fractions and molecule counts.
#'
#' @param run A `celladmix_run` object.
#' @param cell_factors Optional output from `fit$cell_factors()` or
#'   [celladmix_collect_cells()]. Supplying this avoids rereading the cell table.
#' @param importance Factor importance definition. `"molecule_fraction"` uses
#'   cell factor fractions weighted by per-cell molecule counts,
#'   `"cell_mean_fraction"` averages cell factor fractions equally, and
#'   `"loading_fraction"` uses factor loading row sums.
#'
#' @return A data frame with `factor`, `stability`, and `importance`.
#' @keywords internal
#' @noRd
celladmix_nmf_stability_data <- function(
    run,
    cell_factors = NULL,
    importance = c("molecule_fraction", "cell_mean_fraction", "loading_fraction")
  ) {
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_nmf_stability_data() expects a celladmix_run")
  }
  importance <- match.arg(importance)
  n_factors <- if (!is.null(run$n_factors)) {
    as.integer(run$n_factors[[1]])
  } else if (!is.null(run$h)) {
    nrow(run$h)
  } else {
    0L
  }
  if (is.na(n_factors) || n_factors <= 0L) {
    return(data.frame())
  }

  stability <- run$nmf_factor_stability
  if ((is.null(stability) || !length(stability)) && !is.null(run$nmf_diagnostics)) {
    stability <- run$nmf_diagnostics$selected_factor_stability
  }
  stability <- as.numeric(stability)
  if (!length(stability)) {
    stability <- rep(NA_real_, n_factors)
  } else if (length(stability) < n_factors) {
    stability <- c(stability, rep(NA_real_, n_factors - length(stability)))
  } else if (length(stability) > n_factors) {
    stability <- stability[seq_len(n_factors)]
  }

  factor_cols <- paste0("factor_", seq_len(n_factors), "_fraction")
  importance_values <- rep(NA_real_, n_factors)
  if (importance %in% c("molecule_fraction", "cell_mean_fraction")) {
    if (is.null(cell_factors)) {
      cell_factors <- tryCatch(celladmix_collect_cells(run), error = function(e) NULL)
    }
    if (!is.null(cell_factors) && all(factor_cols %in% names(cell_factors))) {
      fraction_matrix <- as.matrix(cell_factors[, factor_cols, drop = FALSE])
      fraction_matrix[!is.finite(fraction_matrix)] <- 0
      if (importance == "molecule_fraction" && "transcript_count" %in% names(cell_factors)) {
        weights <- pmax(as.numeric(cell_factors$transcript_count), 0)
        weighted <- colSums(fraction_matrix * weights, na.rm = TRUE)
        importance_values <- if (sum(weighted) > 0) weighted / sum(weighted) else weighted
      } else {
        means <- colMeans(fraction_matrix, na.rm = TRUE)
        importance_values <- if (sum(means) > 0) means / sum(means) else means
      }
    }
  }
  if (!all(is.finite(importance_values)) && !is.null(run$h)) {
    loading <- rowSums(run$h, na.rm = TRUE)
    importance_values <- if (sum(loading) > 0) loading / sum(loading) else loading
    importance <- "loading_fraction"
  }

  out <- data.frame(
    factor = seq_len(n_factors),
    factor_label = paste0("F", seq_len(n_factors)),
    stability = stability,
    importance = importance_values,
    importance_mode = importance,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

#' Plot NMF Factor Stability versus Factor Importance
#'
#' Draws a compact scatter plot for multiseed NMF runs. The x-axis is factor
#' importance and the y-axis is per-factor restart stability.
#'
#' @inheritParams celladmix_nmf_stability_data
#' @param label Logical; label points by factor id.
#' @param min_stability Reference stability threshold shown as a horizontal
#'   dashed line.
#'
#' @return A `ggplot2` object when ggplot2 is available; otherwise invisibly
#'   returns the plotted diagnostic data frame.
#' @keywords internal
#' @noRd
celladmix_plot_nmf_stability <- function(
    run,
    cell_factors = NULL,
    importance = c("molecule_fraction", "cell_mean_fraction", "loading_fraction"),
    label = TRUE,
    min_stability = 0.3
  ) {
  importance <- match.arg(importance)
  df <- celladmix_nmf_stability_data(run, cell_factors = cell_factors,
    importance = importance)
  n_runs <- if (!is.null(run$pipeline_options$nmf_n_runs)) {
    as.integer(run$pipeline_options$nmf_n_runs[[1]])
  } else {
    length(run$nmf_candidate_final_objectives)
  }
  if (!nrow(df) || all(!is.finite(df$stability)) || is.na(n_runs) || n_runs <= 1L) {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      return(.celladmix_empty_plot(
        "NMF restart stability is unavailable\n(run was not fitted with multiple seeds)"
      ))
    }
    graphics::plot.new()
    graphics::text(0.5, 0.5, "NMF restart stability is unavailable\n(run was not fitted with multiple seeds)")
    return(invisible(df))
  }

  x <- 100 * pmax(df$importance, 0)
  y <- pmax(-1, pmin(1, df$stability))
  max_importance <- max(df$importance, na.rm = TRUE)
  point_cex <- if (is.finite(max_importance) && max_importance > 0) {
    0.8 + 3.0 * sqrt(pmax(df$importance, 0) / max_importance)
  } else {
    rep(1.2, nrow(df))
  }
  point_col <- ifelse(y >= 0.6, "#1B9E77", ifelse(y >= min_stability, "#7570B3", "#D95F02"))
  x_label <- switch(unique(df$importance_mode)[[1]],
    molecule_fraction = "Factor importance (% of assigned molecules)",
    cell_mean_fraction = "Factor importance (% mean cell fraction)",
    loading_fraction = "Factor importance (% of loading mass)",
    "Factor importance (%)")
  subtitle <- sprintf("NMF restarts: %d; selected run: %s; selected seed: %s",
    n_runs,
    if (!is.null(run$nmf_selected_run)) as.character(run$nmf_selected_run[[1]]) else "NA",
    if (!is.null(run$nmf_selected_seed)) as.character(run$nmf_selected_seed[[1]]) else "NA")

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot_df <- df
    plot_df$importance_percent <- x
    plot_df$stability_plot <- y
    plot_df$stability_class <- factor(
      ifelse(y >= 0.6, "high", ifelse(y >= min_stability, "moderate", "low")),
      levels = c("high", "moderate", "low")
    )
    plot_df$point_size <- if (is.finite(max_importance) && max_importance > 0) {
      sqrt(pmax(plot_df$importance, 0) / max_importance)
    } else {
      1
    }
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(
      x = .data$importance_percent,
      y = .data$stability_plot,
      color = .data$stability_class,
      size = .data$point_size,
      label = .data$factor_label
    )) +
      ggplot2::geom_hline(yintercept = min_stability, linetype = "dashed",
        color = "grey55", linewidth = 0.35) +
      ggplot2::geom_hline(yintercept = 0.6, linetype = "dotted",
        color = "grey70", linewidth = 0.35) +
      ggplot2::geom_point(alpha = 0.9) +
      ggplot2::scale_color_manual(values = c(high = "#1B9E77",
        moderate = "#7570B3", low = "#D95F02")) +
      ggplot2::scale_size_continuous(range = c(1.8, 5.0), guide = "none") +
      ggplot2::labs(
        title = "NMF factor stability vs. importance",
        subtitle = subtitle,
        x = x_label,
        y = "Matched ownership correlation across restarts",
        color = "Stability"
      ) +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
        legend.position = "bottom")
    if (isTRUE(label)) {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(size = 3.0, show.legend = FALSE,
          max.overlaps = Inf)
      } else {
        p <- p + ggplot2::geom_text(vjust = -0.6, size = 2.8,
          show.legend = FALSE)
      }
    }
    return(p)
  }

  graphics::plot(x, y, pch = 19, cex = point_cex, col = point_col,
    xlab = x_label,
    ylab = "Matched ownership correlation across restarts",
    main = "NMF factor stability vs. importance",
    ylim = c(max(-0.1, min(y, na.rm = TRUE) - 0.05), 1.02))
  graphics::abline(h = min_stability, lty = 2, col = "grey55")
  graphics::abline(h = 0.6, lty = 3, col = "grey70")
  if (label) {
    graphics::text(x, y, labels = df$factor_label, pos = 3, cex = 0.72)
  }
  graphics::mtext(subtitle, side = 3, line = 0.2, cex = 0.75)
  invisible(df)
}

#' Collect Transcripts from a Persisted Run
#'
#' Reads transcript-level records from a persisted run, optionally restricting
#' to a crop or bounding box and optionally sampling a subset of rows.
#'
#' @param run A `celladmix_run` object.
#' @param sample_n Optional maximum number of transcript rows to return.
#' @param seed Integer RNG seed used when `sample_n` is specified.
#' @param analysis_crop Optional crop identifier.
#' @param analysis_bbox Optional spatial bounding box given as
#'   `c(xmin, xmax, ymin, ymax)`.
#'
#' @return A data frame of transcript-level records from the run.
#' @keywords internal
#' @noRd
celladmix_collect_transcripts <- function(
    run,
    sample_n = NULL,
    seed = 1L,
    analysis_crop = NULL,
    analysis_bbox = NULL
  ) {
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_collect_transcripts() expects a celladmix_run")
  }
  .celladmix_collect_run_transcripts(
    run$path,
    analysis_crop = analysis_crop,
    bbox = analysis_bbox,
    sample_n = if (is.null(sample_n)) NA_integer_ else sample_n,
    seed = seed
  )
}

#' Collect Cell Summaries from a Persisted Run
#'
#' @param run A `celladmix_run` object.
#'
#' @return A data frame of per-cell summaries stored in the run.
#' @keywords internal
#' @noRd
celladmix_collect_cells <- function(run) {
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_collect_cells() expects a celladmix_run")
  }
  .celladmix_collect_run_cells(run$path)
}

#' Collect the Sampled Training Molecules from a Persisted Run
#'
#' Reads the report-sidecar molecular UMAP built from the sampled NCV rows used
#' in NMF training. This includes per-molecule 2D UMAP coordinates, final CRF
#' labels, and normalized NMF component weights.
#'
#' If the sidecar has not been built yet, call [celladmix_report_data()] with
#' `what = "ncv_umap"` first, or fit with `report_ncv_umap = TRUE`.
#'
#' @param run A `celladmix_run` object.
#'
#' @return A data frame of sampled training molecules with molecular UMAP
#'   coordinates and component weights.
#' @keywords internal
#' @noRd
celladmix_collect_training_molecules <- function(run) {
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_collect_training_molecules() expects a celladmix_run")
  }
  .celladmix_collect_training_molecules(run$path)
}

#' Build and Collect Lightweight Report Data
#'
#' Builds developer/report-side data only when requested. This keeps
#' [celladmix_fit()] optimized for core persisted outputs while still allowing
#' notebooks to request cell-level and NCV-level UMAP views on demand.
#'
#' The cell UMAP view is assembled in R by merging a `celladmix_clusters`
#' embedding with per-cell run summaries. The NCV UMAP view is an opt-in
#' report-sidecar stored under `run/report/`.
#'
#' @param run A `celladmix_run` object.
#' @param clust Optional `celladmix_clusters` object. Required when
#'   `what` includes `"cell_umap"`.
#' @param what Character vector selecting which report views to build/return.
#'   Supported values are `"cell_umap"` and `"ncv_umap"`.
#' @param force Logical; if `TRUE`, rebuild requested report-sidecar data even
#'   if it already exists on disk.
#' @param ncv_umap_neighbors Number of neighbors used by the NCV-space UMAP.
#' @param ncv_umap_epochs Number of optimization epochs used by the NCV-space
#'   UMAP.
#' @param ncv_umap_normalization_scale Library-size normalization scale applied
#'   before `log1p` when embedding sampled NCV rows.
#' @param ncv_umap_pca_dims Maximum number of PCA dimensions used before UMAP.
#' @param num_threads Number of worker threads used while building report-only
#'   embeddings.
#' @param seed Integer RNG seed used for the NCV-space UMAP.
#' @param verbose Logical; if `TRUE`, emit `[INFO]` timing messages around each
#'   requested report-data block.
#'
#' @return A named list containing one or both of:
#'   `cell_umap`, `ncv_umap`.
#' @keywords internal
#' @noRd
celladmix_report_data <- function(
    run,
    clust = NULL,
    what = c("cell_umap", "ncv_umap"),
    force = FALSE,
    ncv_umap_neighbors = 15L,
    ncv_umap_epochs = 200L,
    ncv_umap_normalization_scale = 5000,
    ncv_umap_pca_dims = 30L,
    num_threads = 1L,
    seed = 1L,
    verbose = FALSE
  ) {
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_report_data() expects a celladmix_run")
  }
  what <- unique(as.character(what))
  invalid <- setdiff(what, c("cell_umap", "ncv_umap"))
  if (length(invalid) > 0) {
    stop("Unsupported report data request: ", paste(invalid, collapse = ", "))
  }

  info <- .celladmix_info_logger(verbose)
  info("Starting report data collection")
  out <- list()

  if ("cell_umap" %in% what) {
    stage <- Sys.time()
    if (!inherits(clust, "celladmix_clusters")) {
      stop("cell_umap report data requires a celladmix_clusters object")
    }
    cell_df <- celladmix_collect_cells(run)
    out$cell_umap <- merge(clust$embedding, cell_df, by = "cell_id", all.x = TRUE, sort = FALSE)
    info(sprintf("Built cell UMAP report data: %d rows", nrow(out$cell_umap)), stage)
  }

  if ("ncv_umap" %in% what) {
    stage <- Sys.time()
    ncv_path <- run$paths$training_molecule_umap_parquet %||%
      file.path(run$path, "report", "training_molecule_umap.parquet")
    if (!isTRUE(force) && file.exists(ncv_path)) {
      info("Reused existing NCV UMAP report sidecar", stage)
    } else {
      .celladmix_build_report_data(
        run$path,
        build_ncv_umap = TRUE,
        force = force,
        umap_neighbors = ncv_umap_neighbors,
        umap_epochs = ncv_umap_epochs,
        normalization_scale = ncv_umap_normalization_scale,
        pca_dims = ncv_umap_pca_dims,
        num_threads = num_threads,
        seed = seed,
        verbose = verbose
      )
      info("Built NCV UMAP report sidecar", stage)
    }
    ncv_df <- celladmix_collect_training_molecules(run)
    if (inherits(clust, "celladmix_clusters")) {
      ncv_df <- merge(
        ncv_df,
        clust$clusters[, c("cell_id", "cluster")],
        by = "cell_id",
        all.x = TRUE,
        sort = FALSE
      )
    }
    out$ncv_umap <- ncv_df
    info(sprintf("Collected NCV UMAP report data: %d rows", nrow(out$ncv_umap)), stage)
  }

  info("Finished report data collection")
  out
}

#' Collect Cell-by-Gene Counts from a Persisted Run
#'
#' @param run A `celladmix_run` object.
#' @param analysis_crop Optional crop identifier.
#' @param analysis_bbox Optional spatial bounding box given as
#'   `c(xmin, xmax, ymin, ymax)`.
#'
#' @return A numeric matrix with genes in rows and cells in columns.
#' @keywords internal
#' @noRd
celladmix_collect_counts <- function(run, analysis_crop = NULL, analysis_bbox = NULL) {
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_collect_counts() expects a celladmix_run")
  }
  .celladmix_collect_run_counts(run$path, analysis_crop = analysis_crop, bbox = analysis_bbox)
}

#' Collect Sparse Cell-by-Gene Counts from a Persisted Run
#'
#' Streams the run molecule parquet and returns a sparse gene-by-cell count
#' matrix. This is preferred for corrected runs and large panels where dense
#' materialization is unnecessary.
#'
#' @param run A `celladmix_run` object.
#' @param analysis_crop Optional crop identifier.
#' @param analysis_bbox Optional spatial bounding box given as
#'   `c(xmin, xmax, ymin, ymax)`.
#'
#' @return A `Matrix::dgCMatrix` with genes in rows and cells in columns. Cell
#'   types, when available, are attached as the `cell_type` attribute.
#' @keywords internal
#' @noRd
celladmix_collect_counts_sparse <- function(run, analysis_crop = NULL, analysis_bbox = NULL) {
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_collect_counts_sparse() expects a celladmix_run")
  }
  raw <- .celladmix_collect_run_counts_sparse(run$path, analysis_crop = analysis_crop, bbox = analysis_bbox)
  .celladmix_sparse_counts_from_raw(raw)
}

#' Collect a Spatial Region from a Persisted Run
#'
#' Convenience wrapper over [celladmix_collect_transcripts()] for a bounding-box
#' region query.
#'
#' @param run A `celladmix_run` object.
#' @param bbox Spatial bounding box given as `c(xmin, xmax, ymin, ymax)`.
#' @param analysis_crop Optional crop identifier.
#' @param sample_n Optional maximum number of transcript rows to return.
#' @param seed Integer RNG seed used when `sample_n` is specified.
#'
#' @return A data frame of transcript-level records from the selected region.
#' @keywords internal
#' @noRd
celladmix_collect_region <- function(
    run,
    bbox,
    analysis_crop = NULL,
    sample_n = NULL,
    seed = 1L
  ) {
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_collect_region() expects a celladmix_run")
  }
  celladmix_collect_transcripts(
    run,
    analysis_crop = analysis_crop,
    analysis_bbox = bbox,
    sample_n = sample_n,
    seed = seed
  )
}

#' Score Bridge Evidence for a Factor
#'
#' Computes bridge-style cross-cell evidence for fitted factor labels using the
#' persisted run outputs.
#'
#' @param run A `celladmix_run` object.
#' @param annotation Optional cell annotation as a named vector, data frame, or
#'   CSV path. If supplied, labels override run cell types for scoring.
#' @param annotation_col Candidate annotation columns to use when `annotation`
#'   is a data frame or CSV path.
#' @param cell_id_col Cell identifier column in `annotation`.
#' @param factor Optional one-based factor index to keep in returned tables.
#' @param candidate_mode Candidate discovery mode. `"molecule_global"` is the
#'   reference molecule-level path; `"cell_center"` first samples likely
#'   adjacent cell pairs from cell-centroid kNN and then validates them with
#'   exact two-cell molecular kNN.
#' @param candidate_k Number of molecule neighbors used to discover candidate
#'   ordered adjacent cell pairs.
#' @param cell_candidate_k Number of cell-centroid neighbors considered per
#'   cell when `candidate_mode = "cell_center"`.
#' @param candidate_pairs_per_type_pair Maximum cell-center candidate pairs kept
#'   for each ordered cell-type pair before molecular validation.
#' @param crossing_k Number of neighbors used for two-cell bridge-crossing
#'   estimation.
#' @param min_type_pair_contacts Minimum molecule-level contacts required for a
#'   cell-type pair to be tested.
#' @param min_factor_molecules Minimum target-cell molecules assigned to a
#'   factor before that cell contributes to the bridge test.
#' @param min_pairs Minimum cell pairs required for a cell-type/factor summary.
#' @param max_cells_per_type_pair Maximum candidate target cells sampled for one
#'   cell-type/factor summary.
#' @param null_iterations Number of matched null iterations per selected pair.
#' @param null_max_iterations Maximum source-cell movement iterations in each
#'   null placement.
#' @param cell_candidate_halo Optional cell-radius halo for filtering
#'   cell-center candidates. Negative values disable radius filtering.
#' @param null_step Fractional movement step for matched null source cells.
#' @param null_pool_fraction Fraction of most factor-similar same-type source
#'   cells used as the replacement-source sampling pool.
#' @param num_threads Number of worker threads for bridge scoring.
#' @param seed Random seed for deterministic pair/null sampling.
#' @param compute_null Whether to compute null scores and p-values.
#' @param verbose Logical; if `TRUE`, emit `[INFO]` bridge timing and profiling
#'   messages.
#' @param analysis_crop Optional crop identifier.
#' @param analysis_bbox Optional spatial bounding box given as
#'   `c(xmin, xmax, ymin, ymax)`.
#'
#' @return A `celladmix_bridge_result` list with per-pair `scores`, statistical
#'   `summary`, and output `paths`.
#' @keywords internal
#' @noRd
celladmix_score_bridge <- function(
    run,
    annotation = NULL,
    annotation_col = c("merged_annotation", "cell_type", "cluster_label", "cluster"),
    cell_id_col = "cell_id",
    factor = NULL,
    candidate_mode = c("cell_center", "molecule_global"),
    candidate_k = 10L,
    cell_candidate_k = 50L,
    candidate_pairs_per_type_pair = 1000L,
    crossing_k = 20L,
    min_type_pair_contacts = 5L,
    min_factor_molecules = 5L,
    min_pairs = 5L,
    max_cells_per_type_pair = 400L,
    null_iterations = 3L,
    null_max_iterations = 50L,
    cell_candidate_halo = -1,
    null_step = 0.02,
    null_pool_fraction = 0.1,
    num_threads = 1L,
    seed = 1L,
    compute_null = TRUE,
    verbose = FALSE,
    analysis_crop = NULL,
    analysis_bbox = NULL
  ) {
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_score_bridge() expects a celladmix_run")
  }
  candidate_mode <- match.arg(candidate_mode)
  cell_types <- .celladmix_normalize_cell_type_override(
    annotation,
    annotation_col = annotation_col,
    cell_id_col = cell_id_col
  )
  out <- .celladmix_bridge_scores_run(
    run$path,
    analysis_crop = analysis_crop,
    bbox = analysis_bbox,
    cell_types = cell_types,
    candidate_mode = candidate_mode,
    candidate_k = candidate_k,
    cell_candidate_k = cell_candidate_k,
    candidate_pairs_per_type_pair = candidate_pairs_per_type_pair,
    crossing_k = crossing_k,
    min_type_pair_contacts = min_type_pair_contacts,
    min_factor_molecules = min_factor_molecules,
    min_pairs = min_pairs,
    max_cells_per_type_pair = max_cells_per_type_pair,
    null_iterations = null_iterations,
    null_max_iterations = null_max_iterations,
    cell_candidate_halo = cell_candidate_halo,
    null_step = null_step,
    null_pool_fraction = null_pool_fraction,
    num_threads = num_threads,
    seed = seed,
    compute_null = compute_null,
    verbose = verbose
  )
  if (!is.null(factor)) {
    factor <- as.integer(factor)
    out$scores <- out$scores[out$scores$factor %in% factor, , drop = FALSE]
    out$summary <- out$summary[out$summary$factor %in% factor, , drop = FALSE]
  }
  out
}

#' Collapse Bridge Scores for Annotation Plots
#'
#' Converts a bridge-test summary table into a target-cell-type by factor matrix
#' analogous to the original `plot_annot_hmap()` input.
#'
#' @param bridge A `celladmix_bridge_result` from [celladmix_score_bridge()] or
#'   a bridge summary data frame.
#' @param p_thresh P-value threshold used to mark target cell types for removal.
#' @param adjust_p Whether to FDR-adjust p-values before plotting and thresholding.
#'
#' @return A list containing `score_matrix`, `source_calls`, `remove_calls`,
#'   `threshold`, and the normalized `summary` table.
#' @keywords internal
#' @noRd
celladmix_bridge_annotation <- function(bridge, p_thresh = 0.1, adjust_p = FALSE) {
  summary <- if (inherits(bridge, "celladmix_bridge_result")) bridge$summary else bridge
  if (!is.data.frame(summary)) {
    stop("bridge must be a celladmix_bridge_result or bridge summary data frame")
  }
  required <- c("target_cell_type", "source_cell_type", "factor", "p_value", "neg_log10_p")
  missing <- setdiff(required, colnames(summary))
  if (length(missing) > 0L) {
    stop("bridge summary is missing required columns: ", paste(missing, collapse = ", "))
  }
  summary <- summary[!is.na(summary$target_cell_type) & !is.na(summary$source_cell_type) &
    !is.na(summary$factor), , drop = FALSE]
  if (nrow(summary) == 0L) {
    empty <- matrix(numeric(0), nrow = 0L, ncol = 0L)
    return(list(score_matrix = empty, source_calls = list(), remove_calls = character(),
                threshold = -log10(p_thresh), summary = summary))
  }

  summary$factor <- as.integer(summary$factor)
  if (adjust_p) {
    summary$plot_p_value <- stats::p.adjust(summary$p_value, method = "fdr")
    summary$plot_score <- -log10(summary$plot_p_value)
  } else {
    summary$plot_p_value <- summary$p_value
    summary$plot_score <- summary$neg_log10_p
  }
  summary$plot_score[is.nan(summary$plot_score)] <- NA_real_

  cell_types <- sort(unique(c(summary$target_cell_type, summary$source_cell_type)))
  factors <- sort(unique(summary$factor))
  mat <- matrix(NA_real_, nrow = length(cell_types), ncol = length(factors),
                dimnames = list(cell_types, paste0("F", factors)))
  source_score <- matrix(NA_real_, nrow = length(cell_types), ncol = length(factors),
                         dimnames = list(cell_types, paste0("F", factors)))

  for (ct in cell_types) {
    for (factor in factors) {
      idx <- summary$target_cell_type == ct & summary$factor == factor
      vals <- summary$plot_score[idx]
      vals <- vals[is.finite(vals)]
      if (length(vals) > 0L) {
        mat[ct, paste0("F", factor)] <- max(vals)
      }

      sidx <- summary$source_cell_type == ct & summary$factor == factor
      svals <- summary$plot_score[sidx]
      svals <- svals[is.finite(svals)]
      if (length(svals) > 0L) {
        source_score[ct, paste0("F", factor)] <- mean(svals)
      }
    }
  }

  threshold <- -log10(p_thresh)
  source_calls <- list()
  remove_calls <- character()
  for (factor in factors) {
    f_col <- paste0("F", factor)
    src_vals <- source_score[, f_col]
    if (any(is.finite(src_vals))) {
      source_idx <- which.max(replace(src_vals, !is.finite(src_vals), -Inf))
      source_name <- names(src_vals)[source_idx]
      if (is.null(source_name) || is.na(source_name)) {
        source_name <- rownames(source_score)[source_idx]
      }
      source_calls[[paste0("f_", factor)]] <- source_name
    }
    target_vals <- mat[, f_col]
    targets <- names(target_vals)[is.finite(target_vals) & target_vals > threshold]
    source <- source_calls[[paste0("f_", factor)]]
    if (!is.null(source)) {
      targets <- setdiff(targets, source)
    }
    if (length(targets) > 0L) {
      remove_calls <- c(remove_calls, paste0(factor, "_", targets))
    }
  }

  list(score_matrix = mat, source_calls = source_calls, remove_calls = remove_calls,
       threshold = threshold, summary = summary)
}

#' Convert Bridge Annotations to Correction Rules
#'
#' Builds the `factor`/`target_cell_type` rule table consumed by
#' [celladmix_correct()]. Sources are retained for diagnostics, but molecules
#' are removed only from target cell types.
#'
#' @param bridge A `celladmix_bridge_result`, bridge summary data frame, or
#'   object returned by [celladmix_bridge_annotation()].
#' @param p_thresh P-value threshold for target removal calls.
#' @param adjust_p Whether to FDR-adjust p-values before thresholding.
#' @param target_cell_types Optional character vector limiting returned target
#'   cell types.
#'
#' @return Data frame of correction rules.
#' @keywords internal
#' @noRd
celladmix_bridge_rules <- function(
    bridge,
    p_thresh = 0.1,
    adjust_p = FALSE,
    target_cell_types = NULL
  ) {
  ann <- if (is.list(bridge) && all(c("score_matrix", "source_calls", "remove_calls") %in% names(bridge))) {
    bridge
  } else {
    celladmix_bridge_annotation(bridge, p_thresh = p_thresh, adjust_p = adjust_p)
  }
  empty <- data.frame(
    factor = integer(),
    source_cell_type = character(),
    target_cell_type = character(),
    p_value = numeric(),
    neg_log10_p = numeric(),
    rule_id = character(),
    stringsAsFactors = FALSE
  )
  if (!length(ann$remove_calls)) {
    return(empty)
  }
  summary <- ann$summary
  out <- lapply(ann$remove_calls, function(rule_id) {
    factor <- as.integer(sub("^([0-9]+)_.*$", "\\1", rule_id))
    target <- sub("^[0-9]+_", "", rule_id)
    if (!is.null(target_cell_types) && !(target %in% target_cell_types)) {
      return(NULL)
    }
    source <- ann$source_calls[[paste0("f_", factor)]] %||% NA_character_
    idx <- summary$factor == factor & summary$target_cell_type == target
    rows <- summary[idx, , drop = FALSE]
    if (nrow(rows) > 0L && "plot_p_value" %in% names(rows)) {
      best <- which.min(rows$plot_p_value)
    } else if (nrow(rows) > 0L) {
      best <- which.min(rows$p_value)
    } else {
      best <- integer()
    }
    p_value <- if (length(best)) rows$p_value[best] else NA_real_
    neg_log10_p <- if (length(best)) rows$neg_log10_p[best] else NA_real_
    data.frame(
      factor = factor,
      source_cell_type = source,
      target_cell_type = target,
      p_value = p_value,
      neg_log10_p = neg_log10_p,
      rule_id = rule_id,
      stringsAsFactors = FALSE
    )
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) {
    return(empty)
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' Plot Bridge Annotation Heatmap
#'
#' Visualizes bridge-test evidence by target cell type and factor. Cells are
#' colored by the maximum `-log10(p)` over source cell types. An `S` marks the
#' inferred source cell type for a factor and `*` marks target cell types passing
#' `p_thresh`.
#'
#' @param bridge A `celladmix_bridge_result`, bridge summary data frame, or
#'   object returned by [celladmix_bridge_annotation()].
#' @param p_thresh P-value threshold used for marking target cell types.
#' @param adjust_p Whether to FDR-adjust p-values before plotting.
#' @param main Plot title for the base-graphics fallback.
#' @param max_neg_log10 Upper color scale limit.
#' @param source_prior Optional factor source prior from
#'   `fit$score_factor_sources()`. When supplied, `G` marks gene-content source
#'   calls, `S` marks spatial source calls, and `SG` marks agreement.
#' @param use_complexheatmap If `TRUE` and `ComplexHeatmap`/`circlize` are
#'   installed, return a `ComplexHeatmap::Heatmap` object. Otherwise draw a
#'   base R heatmap.
#'
#' @return Invisibly returns the annotation list for base graphics, or a
#'   `ComplexHeatmap::Heatmap` object when `use_complexheatmap = TRUE` and
#'   optional dependencies are available.
#' @keywords internal
#' @noRd
celladmix_plot_bridge_heatmap <- function(
    bridge,
    p_thresh = 0.1,
    adjust_p = FALSE,
    main = NULL,
    max_neg_log10 = 6,
    source_prior = NULL,
    source_prior_min_score = 0.15,
    source_prior_min_margin = 0.03,
    use_complexheatmap = FALSE
  ) {
  main <- main %||% sprintf("Bridge factor annotation, p < %.2g", p_thresh)
  ann <- if (is.list(bridge) && all(c("score_matrix", "source_calls", "remove_calls") %in% names(bridge))) {
    bridge
  } else {
    celladmix_bridge_annotation(bridge, p_thresh = p_thresh, adjust_p = adjust_p)
  }
  gene_source_calls <- list()
  source_prior_ann <- .celladmix_factor_source_prior_annotation(source_prior,
    min_score = source_prior_min_score, min_margin = source_prior_min_margin)
  if (!is.null(source_prior_ann)) {
    source_prior_ann <- source_prior_ann[source_prior_ann$called &
      !is.na(source_prior_ann$source_cell_type) &
      nzchar(source_prior_ann$source_cell_type), , drop = FALSE]
    if (nrow(source_prior_ann)) {
      gene_source_calls <- stats::setNames(
        source_prior_ann$source_cell_type,
        paste0("f_", source_prior_ann$factor)
      )
    }
  }
  source_call <- function(calls, factor_id) {
    key <- paste0("f_", factor_id)
    if (key %in% names(calls)) calls[[key]] else NULL
  }
  mat <- ann$score_matrix
  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    graphics::plot.new()
    graphics::title(main)
    graphics::text(0.5, 0.5, "No bridge summaries to plot")
    return(invisible(ann))
  }

  if (isTRUE(use_complexheatmap) && requireNamespace("ComplexHeatmap", quietly = TRUE) &&
      requireNamespace("circlize", quietly = TRUE)) {
    threshold <- ann$threshold
    col_fun <- circlize::colorRamp2(c(0, threshold, max_neg_log10), c("white", "white", "red"))
    hmap <- ComplexHeatmap::Heatmap(
      pmin(mat, max_neg_log10), name = "-log10(p)", cluster_rows = FALSE,
      cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
      col = col_fun, border = TRUE, row_names_side = "left",
      column_title = "Factor", row_title = NULL, column_title_side = "bottom",
      column_names_side = "bottom", row_title_side = "left", na_col = "white",
      cell_fun = function(j, i, x, y, width, height, fill) {
        factor_id <- sub("^F", "", colnames(mat)[j])
        ct <- rownames(mat)[i]
        remove_id <- paste0(factor_id, "_", ct)
        source <- source_call(ann$source_calls, factor_id)
        gene_source <- source_call(gene_source_calls, factor_id)
        fontsize <- min(grid::convertWidth(width, "mm", valueOnly = TRUE),
                        grid::convertHeight(height, "mm", valueOnly = TRUE)) * 2
        label <- character()
        if (!is.null(source) && !is.null(gene_source) && identical(source, ct) &&
            identical(gene_source, ct)) {
          label <- "SG"
        } else if (!is.null(source) && identical(source, ct)) {
          label <- "S"
        } else if (!is.null(gene_source) && identical(gene_source, ct)) {
          label <- "G"
        }
        if (remove_id %in% ann$remove_calls) {
          label <- paste0(label, "*")
        }
        if (length(label) && nzchar(label)) {
          grid::grid.text(label, x, y, gp = grid::gpar(fontsize = fontsize,
            fontface = "bold"))
        }
      }
    )
    return(hmap)
  }

  values <- pmin(mat, max_neg_log10)
  threshold <- ann$threshold
  palette_fun <- grDevices::colorRampPalette(c("white", "red"))
  color_for <- function(value) {
    if (!is.finite(value) || value <= threshold) {
      return("white")
    }
    scaled <- (value - threshold) / max(max_neg_log10 - threshold, .Machine$double.eps)
    palette_fun(100L)[max(1L, min(100L, ceiling(scaled * 100)))]
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  nr <- nrow(values)
  nc <- ncol(values)
  left_margin <- max(8, min(18, 4 + 0.35 * max(nchar(rownames(values)))))
  graphics::par(mar = c(if (!is.null(source_prior_ann)) 5 else 4, left_margin, 3, 5), xpd = NA)
  graphics::plot(NA, xlim = c(0.5, nc + 1.45), ylim = c(0.5, nr + 0.5),
                 axes = FALSE, xlab = "Factor", ylab = "", main = main)
  graphics::axis(1, at = seq_len(nc), labels = colnames(values), las = 2)
  graphics::axis(2, at = rev(seq_len(nr)), labels = rownames(values), las = 2)
  for (i in seq_len(nr)) {
    y <- nr - i + 1L
    for (j in seq_len(nc)) {
      graphics::rect(j - 0.5, y - 0.5, j + 0.5, y + 0.5,
                     col = color_for(values[i, j]), border = "grey70")
      factor_id <- sub("^F", "", colnames(values)[j])
      ct <- rownames(values)[i]
      remove_id <- paste0(factor_id, "_", ct)
      source <- source_call(ann$source_calls, factor_id)
      gene_source <- source_call(gene_source_calls, factor_id)
      label <- character()
      if (!is.null(source) && !is.null(gene_source) && identical(source, ct) &&
          identical(gene_source, ct)) {
        label <- "SG"
      } else if (!is.null(source) && identical(source, ct)) {
        label <- "S"
      } else if (!is.null(gene_source) && identical(gene_source, ct)) {
        label <- "G"
      }
      if (remove_id %in% ann$remove_calls) {
        label <- paste0(label, "*")
      }
      if (length(label) && nzchar(label)) {
        graphics::text(j, y, label, font = 2,
          cex = if (nchar(label) > 1L) 0.9 else 1.15)
      }
    }
  }
  if (!is.null(source_prior_ann)) {
    graphics::mtext("S spatial source; G gene-content source; SG agreement; separate S/G in one column indicates conflict; * target",
      side = 1, line = 3.7, cex = 0.68)
  }

  key_n <- 80L
  key_values <- seq(threshold, max_neg_log10, length.out = key_n)
  key_y <- seq(0.8, nr + 0.2, length.out = key_n + 1L)
  key_x0 <- nc + 0.78
  key_x1 <- nc + 0.96
  for (k in seq_len(key_n)) {
    graphics::rect(key_x0, key_y[k], key_x1, key_y[k + 1L],
      col = color_for(key_values[k]), border = NA)
  }
  graphics::rect(key_x0, key_y[1L], key_x1, key_y[key_n + 1L],
    border = "grey40", col = NA)
  graphics::text(key_x1 + 0.08, key_y[1L], sprintf("%.1f", threshold),
    adj = c(0, 0.5), cex = 0.75)
  graphics::text(key_x1 + 0.08, key_y[key_n + 1L], sprintf("%.1f", max_neg_log10),
    adj = c(0, 0.5), cex = 0.75)
  graphics::text(key_x1 + 0.34, mean(range(key_y)), "-log10(p)",
    srt = 90, cex = 0.8)
  invisible(ann)
}

#' @keywords internal
#' @noRd
.celladmix_manifest_path_from_run <- function(run) {
  candidates <- c(
    run$source$manifest_path %||% NA_character_,
    run$source$path %||% NA_character_
  )
  files <- run$source$files
  if (is.data.frame(files) && all(c("role", "path") %in% colnames(files))) {
    candidates <- c(candidates, files$path[files$role %in% c("experiment_manifest", "manifest")])
  } else if (is.list(files)) {
    for (file in files) {
      if (is.list(file) && identical(file$role %||% NA_character_, "experiment_manifest")) {
        candidates <- c(candidates, file$path %||% NA_character_)
      }
    }
  }
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) {
    return(NULL)
  }
  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

#' @keywords internal
#' @noRd
.celladmix_focus_image_candidate <- function(focus_path, focus_index) {
  if (is.null(focus_path) || length(focus_path) < 1L || !nzchar(focus_path[[1]])) {
    return(NULL)
  }
  if (is.null(focus_index) || is.na(focus_index)) {
    return(focus_path)
  }
  focus_index <- as.integer(focus_index[[1]])
  if (is.na(focus_index) || focus_index < 0L) {
    stop("focus_index must be a non-negative integer or NULL")
  }
  replacement <- sprintf("_%04d\\1", focus_index)
  sub("_[0-9]{4}(\\.ome\\.tiff?$)", replacement, focus_path, ignore.case = TRUE)
}

#' @keywords internal
#' @noRd
.celladmix_normalize_stain_name <- function(stain) {
  stain <- tolower(as.character(stain %||% "membrane")[[1]])
  switch(stain,
    boundary = "membrane",
    membranes = "membrane",
    nuclei = "dapi",
    nucleus = "dapi",
    stain
  )
}

#' Discover a Xenium Stain Image
#'
#' Resolves a named Xenium morphology-focus image and pixel size. The discovery
#' logic is intentionally soft: explicit `image_path` always wins; otherwise
#' known stain names map to conventional focus-image indices for current Xenium
#' multimodal bundles (`dapi` = 0, `membrane` = 1); `morphology` selects the
#' full morphology image from the manifest.
#'
#' @param run A `celladmix_run` object.
#' @param stain Stain name: `"membrane"`, `"dapi"`, or `"morphology"`.
#' @param image_path Optional explicit image path. Relative paths are resolved
#'   relative to the Xenium bundle when possible.
#' @param focus_index Optional Xenium morphology-focus image index override.
#' @param pixel_size Optional microns-per-pixel override. If `NULL`, the Xenium
#'   manifest pixel size is used.
#' @param x_offset,y_offset Optional physical-coordinate offset in microns.
#'
#' @return A list with `path`, `pixel_size`, `x_offset`, `y_offset`, and
#'   `source`.
#' @export
celladmix_discover_stain_image <- function(
    run,
    stain = "membrane",
    image_path = NULL,
    focus_index = NULL,
    pixel_size = NULL,
    x_offset = 0,
    y_offset = 0
  ) {
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_discover_stain_image() expects a celladmix_run")
  }
  stain <- .celladmix_normalize_stain_name(stain)
  if (!(stain %in% c("membrane", "dapi", "morphology"))) {
    stop("Unknown stain: ", stain, ". Use 'membrane', 'dapi', or 'morphology'.")
  }
  manifest_path <- .celladmix_manifest_path_from_run(run)
  manifest <- NULL
  bundle_dir <- NULL
  if (!is.null(manifest_path)) {
    manifest <- celladmix_read_xenium_manifest(manifest_path)
    bundle_dir <- dirname(manifest_path)
  }

  resolve_path <- function(path) {
    if (is.null(path) || length(path) < 1L || !nzchar(path[[1]])) {
      return(NULL)
    }
    path <- path[[1]]
    candidates <- path
    if (!is.null(bundle_dir) && !grepl("^/", path)) {
      candidates <- c(file.path(bundle_dir, path), candidates)
    }
    candidates <- candidates[file.exists(candidates)]
    if (!length(candidates)) {
      return(NULL)
    }
    normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
  }

  source <- "explicit"
  resolved_path <- resolve_path(image_path)
  if (!is.null(image_path) && is.null(resolved_path)) {
    stop("Explicit image_path does not exist: ", image_path[[1]])
  }
  if (is.null(resolved_path) && !identical(stain, "morphology")) {
    source <- paste0("xenium_focus_", stain)
    if (is.null(focus_index) || length(focus_index) < 1L || is.na(focus_index[[1]])) {
      focus_index <- switch(stain, dapi = 0L, membrane = 1L)
    }
    focus_rel <- manifest$morphology_focus_filepath %||% NULL
    focus_path <- resolve_path(.celladmix_focus_image_candidate(focus_rel, focus_index))
    if (is.null(focus_path) && identical(stain, "dapi")) {
      focus_path <- resolve_path(focus_rel)
    }
    resolved_path <- focus_path
  }
  if (is.null(resolved_path)) {
    source <- "xenium_morphology"
    resolved_path <- resolve_path(manifest$morphology_filepath %||% NULL)
  }
  if (is.null(resolved_path)) {
    stop(
      "Could not discover a ", stain, " image. Provide image_path explicitly, ",
      "or check the Xenium manifest/focus-image paths for this bundle."
    )
  }

  if (is.null(pixel_size) || is.na(pixel_size[[1]])) {
    pixel_size <- manifest$pixel_size %||% NA_real_
  }
  pixel_size <- as.numeric(pixel_size[[1]])
  if (!is.finite(pixel_size) || pixel_size <= 0) {
    stop("A positive pixel_size is required for membrane scoring")
  }

  list(
    path = resolved_path,
    pixel_size = pixel_size,
    x_offset = as.numeric(x_offset[[1]]),
    y_offset = as.numeric(y_offset[[1]]),
    source = source,
    stain = stain,
    manifest_path = manifest_path
  )
}

#' Discover All Available Stain Images for a Run
#'
#' Quiet counterpart of [celladmix_discover_stain_image()]: resolves each of
#' the requested stains and silently skips ones the source bundle does not
#' provide. Returns an empty list for non-Xenium sources, which makes the
#' result directly usable as a default example-plot background specification.
#'
#' @keywords internal
#' @noRd
.celladmix_discover_stain_images <- function(run, names = c("dapi", "membrane")) {
  out <- list()
  for (name in names) {
    image <- tryCatch(
      celladmix_discover_stain_image(run, stain = name),
      error = function(e) NULL
    )
    if (!is.null(image)) {
      out[[name]] <- image
    }
  }
  out
}

#' Discover a Xenium Membrane Image for Scoring
#'
#' Backward-compatible wrapper around [celladmix_discover_stain_image()] used by
#' membrane and stain-aware coherence scoring.
#'
#' @keywords internal
#' @noRd
celladmix_discover_membrane_image <- function(
    run,
    image_path = NULL,
    focus_index = 1L,
    pixel_size = NULL,
    x_offset = 0,
    y_offset = 0
  ) {
  celladmix_discover_stain_image(
    run,
    stain = "membrane",
    image_path = image_path,
    focus_index = focus_index,
    pixel_size = pixel_size,
    x_offset = x_offset,
    y_offset = y_offset
  )
}

#' Read a Small Stain-Image Crop
#'
#' Reads a physical-coordinate window from a single-channel Xenium morphology
#' image. This is intended for diagnostics and report notebooks, not for bulk
#' image processing.
#'
#' @param image A list returned by [celladmix_discover_membrane_image()], or a
#'   character image path.
#' @param bbox Physical-coordinate bounding box `c(xmin, xmax, ymin, ymax)`.
#' @param pixel_size Microns per image pixel. Required when `image` is a path.
#' @param x_offset,y_offset Physical-coordinate offsets in microns.
#' @param max_pixels Maximum returned width or height in pixels; larger windows
#'   are downsampled during reading.
#'
#' @return A list containing a matrix of image values plus `xlim`, `ylim`,
#'   `pixel_size`, `width`, `height`, and `path`.
#' @keywords internal
#' @noRd
celladmix_read_stain_crop <- function(
    image,
    bbox,
    pixel_size = NULL,
    x_offset = 0,
    y_offset = 0,
    max_pixels = 512L
  ) {
  if (is.list(image)) {
    image_path <- image$path
    pixel_size <- pixel_size %||% image$pixel_size
    x_offset <- image$x_offset %||% x_offset
    y_offset <- image$y_offset %||% y_offset
  } else {
    image_path <- as.character(image[[1]])
  }
  if (is.null(image_path) || !file.exists(image_path)) {
    stop("Image path does not exist: ", image_path %||% "<NULL>")
  }
  if (is.null(pixel_size) || !is.finite(as.numeric(pixel_size[[1]])) || as.numeric(pixel_size[[1]]) <= 0) {
    stop("A positive pixel_size is required")
  }
  if (!is.numeric(bbox) || length(bbox) != 4L) {
    stop("bbox must be numeric c(xmin, xmax, ymin, ymax)")
  }
  .celladmix_read_stain_image_crop(
    image_path = normalizePath(image_path, winslash = "/", mustWork = TRUE),
    pixel_size = as.numeric(pixel_size[[1]]),
    x_offset = as.numeric(x_offset[[1]]),
    y_offset = as.numeric(y_offset[[1]]),
    bbox = as.numeric(bbox),
    max_pixels = as.integer(max_pixels[[1]])
  )
}

#' @keywords internal
#' @noRd
.celladmix_normalize_cell_type_override <- function(
    annotation,
    annotation_col = c("merged_annotation", "cell_type", "cluster_label", "cluster"),
    cell_id_col = "cell_id"
  ) {
  if (is.null(annotation)) {
    return(NULL)
  }
  if (is.character(annotation) && length(annotation) == 1L && file.exists(annotation)) {
    annotation <- utils::read.csv(annotation, stringsAsFactors = FALSE)
  }
  if (is.atomic(annotation) && !is.null(names(annotation))) {
    labels <- as.character(annotation)
    cells <- names(annotation)
  } else if (is.data.frame(annotation)) {
    if (!(cell_id_col %in% colnames(annotation))) {
      stop("annotation is missing cell_id column: ", cell_id_col)
    }
    requested_cols <- annotation_col
    annotation_col <- annotation_col[annotation_col %in% colnames(annotation)]
    if (!length(annotation_col)) {
      stop("annotation is missing one of annotation_col: ", paste(requested_cols, collapse = ", "))
    }
    annotation_col <- annotation_col[[1]]
    cells <- as.character(annotation[[cell_id_col]])
    labels <- as.character(annotation[[annotation_col]])
  } else {
    stop("annotation must be NULL, a named vector, a data frame, or a CSV path")
  }
  keep <- !is.na(cells) & nzchar(cells) & !is.na(labels) & nzchar(labels)
  cells <- cells[keep]
  labels <- labels[keep]
  if (!length(cells)) {
    stop("annotation did not contain any usable cell labels")
  }
  keep_first <- !duplicated(cells)
  stats::setNames(labels[keep_first], cells[keep_first])
}

#' Score Membrane Evidence for Factor Transfer
#'
#' Computes membrane-stain enrichment for ordered target/source cell-type pairs.
#' Candidate cell pairs are sampled from cell-center adjacency, then per-factor
#' target-cell molecules are compared with same-cell distance-matched controls.
#'
#' @param run A `celladmix_run` object.
#' @param annotation Optional cell annotation as a named vector, data frame, or
#'   CSV path. If supplied, labels override run cell types for scoring.
#' @param annotation_col Candidate annotation columns to use when `annotation`
#'   is a data frame or CSV path.
#' @param cell_id_col Cell identifier column in `annotation`.
#' @param image_path Optional explicit membrane image path.
#' @param focus_index Xenium focus-image index used when discovering the image.
#' @param pixel_size Optional microns-per-pixel override.
#' @param x_offset,y_offset Optional physical-coordinate offset in microns.
#' @param epsilon Positive pseudocount added before log signal ratios.
#' @param factor Optional one-based factor index to keep in returned tables.
#' @param cell_candidate_k Number of cell-center neighbors considered per cell.
#' @param candidate_pairs_per_type_pair Maximum candidate ordered cell pairs
#'   kept for each cell-type pair before factor scoring.
#' @param cell_candidate_halo Optional cell-radius halo for filtering
#'   cell-center candidates. Negative values disable radius filtering.
#' @param min_factor_molecules Minimum target-cell factor molecules required
#'   before a pair/factor is scored.
#' @param min_pairs Minimum scored cell pairs required for a summary row.
#' @param max_cells_per_type_pair Maximum scored cell pairs retained in a
#'   cell-type/factor summary.
#' @param control_distance_fraction Relative centroid-distance tolerance used
#'   when selecting same-cell control molecules.
#' @param line_samples Number of interpolation samples along each
#'   molecule-to-centroid line segment.
#' @param num_threads Number of worker threads for scoring.
#' @param seed Random seed for deterministic candidate/control tie-breaking.
#' @param verbose Logical; if `TRUE`, emit `[INFO]` timing messages.
#' @param analysis_crop Optional crop identifier.
#' @param analysis_bbox Optional spatial bounding box given as
#'   `c(xmin, xmax, ymin, ymax)`.
#'
#' @return A `celladmix_membrane_result` list with per-pair `scores`,
#'   statistical `summary`, image metadata, and output `paths`.
#' @keywords internal
#' @noRd
celladmix_score_membrane <- function(
    run,
    annotation = NULL,
    annotation_col = c("merged_annotation", "cell_type", "cluster_label", "cluster"),
    cell_id_col = "cell_id",
    image_path = NULL,
    focus_index = 1L,
    pixel_size = NULL,
    x_offset = 0,
    y_offset = 0,
    epsilon = 1,
    factor = NULL,
    cell_candidate_k = 50L,
    candidate_pairs_per_type_pair = 400L,
    cell_candidate_halo = -1,
    min_factor_molecules = 5L,
    min_pairs = 5L,
    max_cells_per_type_pair = 400L,
    control_distance_fraction = 0.025,
    line_samples = 16L,
    num_threads = 1L,
    seed = 1L,
    verbose = FALSE,
    analysis_crop = NULL,
    analysis_bbox = NULL
  ) {
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_score_membrane() expects a celladmix_run")
  }
  image <- celladmix_discover_membrane_image(
    run,
    image_path = image_path,
    focus_index = focus_index,
    pixel_size = pixel_size,
    x_offset = x_offset,
    y_offset = y_offset
  )
  cell_types <- .celladmix_normalize_cell_type_override(
    annotation,
    annotation_col = annotation_col,
    cell_id_col = cell_id_col
  )
  out <- .celladmix_membrane_scores_run(
    run$path,
    analysis_crop = analysis_crop,
    bbox = analysis_bbox,
    image_path = image$path,
    pixel_size = image$pixel_size,
    x_offset = image$x_offset,
    y_offset = image$y_offset,
    epsilon = epsilon,
    cell_types = cell_types,
    cell_candidate_k = cell_candidate_k,
    candidate_pairs_per_type_pair = candidate_pairs_per_type_pair,
    cell_candidate_halo = cell_candidate_halo,
    min_factor_molecules = min_factor_molecules,
    min_pairs = min_pairs,
    max_cells_per_type_pair = max_cells_per_type_pair,
    control_distance_fraction = control_distance_fraction,
    line_samples = line_samples,
    num_threads = num_threads,
    seed = seed,
    verbose = verbose
  )
  out$image$source <- image$source
  out$image$manifest_path <- image$manifest_path
  if (!is.null(factor)) {
    factor <- as.integer(factor)
    out$scores <- out$scores[out$scores$factor %in% factor, , drop = FALSE]
    out$summary <- out$summary[out$summary$factor %in% factor, , drop = FALSE]
  }
  out
}

#' Collapse Membrane Scores for Annotation Plots
#'
#' @param membrane A `celladmix_membrane_result` from
#'   [celladmix_score_membrane()] or a membrane summary data frame.
#' @param p_thresh P-value threshold used to mark target cell types.
#' @param adjust_p Whether to FDR-adjust p-values before plotting and
#'   thresholding.
#'
#' @return A list containing `score_matrix`, `source_calls`, `remove_calls`,
#'   `threshold`, and the normalized `summary` table.
#' @keywords internal
#' @noRd
celladmix_membrane_annotation <- function(membrane, p_thresh = 0.1, adjust_p = FALSE) {
  summary <- if (inherits(membrane, "celladmix_membrane_result")) membrane$summary else membrane
  celladmix_bridge_annotation(summary, p_thresh = p_thresh, adjust_p = adjust_p)
}

#' Convert Membrane Annotations to Correction Rules
#'
#' @inheritParams celladmix_membrane_annotation
#' @param target_cell_types Optional character vector limiting returned target
#'   cell types.
#'
#' @return Data frame of correction rules.
#' @keywords internal
#' @noRd
celladmix_membrane_rules <- function(
    membrane,
    p_thresh = 0.1,
    adjust_p = FALSE,
    target_cell_types = NULL
  ) {
  ann <- if (is.list(membrane) && all(c("score_matrix", "source_calls", "remove_calls") %in% names(membrane))) {
    membrane
  } else {
    celladmix_membrane_annotation(membrane, p_thresh = p_thresh, adjust_p = adjust_p)
  }
  celladmix_bridge_rules(
    ann,
    p_thresh = p_thresh,
    adjust_p = adjust_p,
    target_cell_types = target_cell_types
  )
}

#' Plot Membrane Annotation Heatmap
#'
#' Visualizes membrane-test evidence by target cell type and factor. Cells are
#' colored by maximum `-log10(p)` over source cell types; `S` marks inferred
#' source cell types and `*` marks target cell types passing `p_thresh`.
#'
#' @inheritParams celladmix_plot_bridge_heatmap
#' @param membrane A `celladmix_membrane_result`, membrane summary data frame,
#'   or object returned by [celladmix_membrane_annotation()].
#'
#' @return Invisibly returns the annotation list for base graphics, or a
#'   `ComplexHeatmap::Heatmap` object when optional dependencies are available.
#' @keywords internal
#' @noRd
celladmix_plot_membrane_heatmap <- function(
    membrane,
    p_thresh = 0.1,
    adjust_p = FALSE,
    main = NULL,
    max_neg_log10 = 6,
    source_prior = NULL,
    source_prior_min_score = 0.15,
    source_prior_min_margin = 0.03,
    use_complexheatmap = FALSE
  ) {
  main <- main %||% sprintf("Membrane factor annotation, p < %.2g", p_thresh)
  ann <- if (is.list(membrane) && all(c("score_matrix", "source_calls", "remove_calls") %in% names(membrane))) {
    membrane
  } else {
    celladmix_membrane_annotation(membrane, p_thresh = p_thresh, adjust_p = adjust_p)
  }
  celladmix_plot_bridge_heatmap(
    ann,
    p_thresh = p_thresh,
    adjust_p = adjust_p,
    main = main,
    max_neg_log10 = max_neg_log10,
    source_prior = source_prior,
    source_prior_min_score = source_prior_min_score,
    source_prior_min_margin = source_prior_min_margin,
    use_complexheatmap = use_complexheatmap
  )
}

#' Score Local Molecular Coherence for Factor Transfer
#'
#' Computes a target-cell-local coherence score for factor-assigned molecules.
#' The score asks whether source-enriched factor molecules form coherent patches
#' inside target cells. With `stain_mode = "membrane_barrier"`, molecule
#' neighbor support is downweighted across membrane-stain barriers.
#'
#' @param run A `celladmix_run` object.
#' @param annotation Optional cell annotation as a named vector, data frame, or
#'   CSV path. If supplied, labels override run cell types for scoring.
#' @param annotation_col Candidate annotation columns to use when `annotation`
#'   is a data frame or CSV path.
#' @param cell_id_col Cell identifier column in `annotation`.
#' @param stain_mode `"none"` for Euclidean molecular coherence, or
#'   `"membrane_barrier"` to downweight molecule-neighbor edges across membrane
#'   signal.
#' @param image_path Optional explicit membrane image path when
#'   `stain_mode = "membrane_barrier"`.
#' @param focus_index Xenium focus-image index used when discovering the image.
#' @param pixel_size Optional microns-per-pixel override.
#' @param x_offset,y_offset Optional physical-coordinate offset in microns.
#' @param stain_max_pixels Maximum width or height of the stain crop used for
#'   barrier sampling.
#' @param factor Optional one-based factor index to keep in returned tables.
#' @param k_neighbors Number of within-cell molecule neighbors.
#' @param min_factor_molecules Minimum molecules assigned to a factor in a cell
#'   before the cell contributes a coherence score.
#' @param min_cells Minimum scored cells required for a target/source/factor
#'   summary.
#' @param max_neighbor_distance Optional maximum molecule-neighbor distance.
#'   Negative values disable filtering.
#' @param distance_sigma Optional Gaussian distance scale for neighbor weights.
#'   Negative values use equal distance weights.
#' @param lambda_coherence Weight of local same-factor neighbor support.
#' @param beta_margin Weight of per-molecule factor margin in the raw score.
#' @param score_threshold Molecule-level threshold used to count active
#'   molecules in summaries.
#' @param source_pseudocount Pseudocount for factor-source enrichment estimates.
#' @param min_source_log_enrichment Minimum source log enrichment required for a
#'   source/factor pair to be scored.
#' @param include_self_source Whether to keep source types matching target cell
#'   type in the cell-level scores.
#' @param compute_null Whether to estimate a within-cell matched null score.
#' @param null_method Null model. `"label_permutation"` shuffles the candidate
#'   factor label within each target cell while preserving nucleus/density
#'   strata and the observed factor-margin distribution. `"matched_subset"`
#'   samples matched decoy molecule subsets and preserves the previous behavior.
#' @param null_iterations Number of random matched subsets per cell/factor/source
#'   score when `compute_null = TRUE`.
#' @param null_exclude_factor Whether null subsets should avoid the observed
#'   factor molecules when enough other molecules are available.
#' @param null_match_nucleus Whether null subsets should preserve available
#'   Xenium nucleus-overlap and nucleus-distance strata.
#' @param null_match_density Whether null subsets should preserve local molecule
#'   density strata within each cell.
#' @param null_nucleus_distance_bins,null_density_bins Number of quantile bins
#'   used for nucleus-distance and local-density matching.
#' @param seed Random seed for null subset sampling.
#' @param normalize_membrane Whether to robustly rescale the membrane crop to
#'   `[0, 1]` before applying the barrier.
#' @param membrane_low_quantile,membrane_high_quantile Quantiles used for
#'   robust membrane rescaling.
#' @param membrane_alpha Barrier strength for membrane-weighted edges.
#' @param line_samples Number of samples along molecule-neighbor segments.
#' @param patch_edge_weight_min Minimum molecule-neighbor edge weight used when
#'   connecting same-factor molecules into patch components.
#' @param num_threads Number of worker threads.
#' @param verbose Logical; if `TRUE`, emit `[INFO]` timing messages.
#' @param analysis_crop Optional crop identifier.
#' @param analysis_bbox Optional spatial bounding box given as
#'   `c(xmin, xmax, ymin, ymax)`.
#'
#' @return A `celladmix_coherence_result` list with per-cell `scores`,
#'   statistical `summary`, and output `paths`.
#' @keywords internal
#' @noRd
celladmix_score_coherence <- function(
    run,
    annotation = NULL,
    annotation_col = c("merged_annotation", "cell_type", "cluster_label", "cluster"),
    cell_id_col = "cell_id",
    stain_mode = c("membrane_barrier", "none"),
    image_path = NULL,
    focus_index = 1L,
    pixel_size = NULL,
    x_offset = 0,
    y_offset = 0,
    stain_max_pixels = 2048L,
    factor = NULL,
    k_neighbors = 10L,
    min_factor_molecules = 5L,
    min_cells = 5L,
    max_neighbor_distance = -1,
    distance_sigma = -1,
    lambda_coherence = 1,
    beta_margin = 1,
    score_threshold = 0,
    source_pseudocount = 0.5,
    min_source_log_enrichment = 0,
    include_self_source = FALSE,
    compute_null = TRUE,
    null_method = c("label_permutation", "matched_subset"),
    null_iterations = 20L,
    null_exclude_factor = TRUE,
    null_match_nucleus = TRUE,
    null_match_density = TRUE,
    null_nucleus_distance_bins = 3L,
    null_density_bins = 3L,
    seed = 1L,
    normalize_membrane = TRUE,
    membrane_low_quantile = 0.05,
    membrane_high_quantile = 0.995,
    membrane_alpha = 3,
    line_samples = 12L,
    patch_edge_weight_min = 0.1,
    num_threads = 1L,
    verbose = FALSE,
    analysis_crop = NULL,
    analysis_bbox = NULL
  ) {
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_score_coherence() expects a celladmix_run")
  }
  stain_mode <- match.arg(stain_mode)
  null_method <- match.arg(null_method)
  image <- NULL
  if (identical(stain_mode, "membrane_barrier")) {
    image <- celladmix_discover_membrane_image(
      run,
      image_path = image_path,
      focus_index = focus_index,
      pixel_size = pixel_size,
      x_offset = x_offset,
      y_offset = y_offset
    )
  }
  cell_types <- .celladmix_normalize_cell_type_override(
    annotation,
    annotation_col = annotation_col,
    cell_id_col = cell_id_col
  )
  out <- .celladmix_coherence_scores_run(
    run$path,
    analysis_crop = analysis_crop,
    bbox = analysis_bbox,
    cell_types = cell_types,
    stain_mode = stain_mode,
    image_path = image$path %||% "",
    pixel_size = image$pixel_size %||% 1,
    x_offset = image$x_offset %||% 0,
    y_offset = image$y_offset %||% 0,
    stain_max_pixels = stain_max_pixels,
    k_neighbors = k_neighbors,
    min_factor_molecules = min_factor_molecules,
    min_cells = min_cells,
    max_neighbor_distance = max_neighbor_distance,
    distance_sigma = distance_sigma,
    lambda_coherence = lambda_coherence,
    beta_margin = beta_margin,
    score_threshold = score_threshold,
    source_pseudocount = source_pseudocount,
    min_source_log_enrichment = min_source_log_enrichment,
    include_self_source = include_self_source,
    compute_null = compute_null,
    null_method = null_method,
    null_iterations = null_iterations,
    null_exclude_factor = null_exclude_factor,
    null_match_nucleus = null_match_nucleus,
    null_match_density = null_match_density,
    null_nucleus_distance_bins = null_nucleus_distance_bins,
    null_density_bins = null_density_bins,
    seed = seed,
    normalize_membrane = normalize_membrane,
    membrane_low_quantile = membrane_low_quantile,
    membrane_high_quantile = membrane_high_quantile,
    membrane_alpha = membrane_alpha,
    line_samples = line_samples,
    patch_edge_weight_min = patch_edge_weight_min,
    num_threads = num_threads,
    verbose = verbose
  )
  if (!is.null(image)) {
    out$image <- image
  }
  if (!is.null(factor)) {
    factor <- as.integer(factor)
    out$scores <- out$scores[out$scores$factor %in% factor, , drop = FALSE]
    out$summary <- out$summary[out$summary$factor %in% factor, , drop = FALSE]
  }
  out
}

#' Collapse Coherence Scores for Annotation Plots
#'
#' @param coherence A `celladmix_coherence_result` from
#'   [celladmix_score_coherence()] or a coherence summary data frame.
#' @param source_mode How to infer the source cell type for each factor.
#'   Coherence target evidence is still based on p-values, but the default
#'   source call uses the factor/source prior estimated from gene enrichment.
#' @param score_cap Finite cap used for plotting and rule aggregation when
#'   statistical p-values underflow to zero.
#' @param score_mode Which coherence evidence to use for target cleanup. `"patch"`
#'   uses patch-level p-values when available, `"mean"` uses the original
#'   mean-score p-values, and `"combined"` requires both.
#' @param min_delta_score Minimum observed-minus-null mean score required for a
#'   target cleanup call when `mean_delta_score` is available.
#' @param min_delta_patch_score Minimum observed-minus-null patch score required
#'   for a target cleanup call when patch columns are available.
#' @param min_active_fraction Minimum active molecule fraction required for a
#'   target cleanup call when `active_fraction` is available.
#' @param min_largest_patch_fraction Minimum largest-patch fraction required
#'   when patch columns are available.
#' @param min_cells,min_molecules Optional minimum summary support required for
#'   a target cleanup call when `n_cells`/`n_molecules` are available.
#' @inheritParams celladmix_bridge_annotation
#'
#' @return A list containing `score_matrix`, `source_calls`, `remove_calls`,
#'   `threshold`, and the normalized `summary` table.
#' @keywords internal
#' @noRd
celladmix_coherence_annotation <- function(
    coherence,
    p_thresh = 0.1,
    adjust_p = FALSE,
    source_mode = c("source_probability", "source_log_enrichment", "plot_score"),
    score_cap = 30,
    score_mode = c("patch", "mean", "combined"),
    min_delta_score = 0,
    min_delta_patch_score = 0,
    min_active_fraction = 0,
    min_largest_patch_fraction = 0,
    min_cells = NULL,
    min_molecules = NULL
  ) {
  source_mode <- match.arg(source_mode)
  score_mode <- match.arg(score_mode)
  summary <- if (inherits(coherence, "celladmix_coherence_result")) coherence$summary else coherence
  if (!is.data.frame(summary)) {
    stop("coherence must be a celladmix_coherence_result or coherence summary data frame")
  }
  required <- c("target_cell_type", "source_cell_type", "factor", "p_value", "neg_log10_p")
  missing <- setdiff(required, colnames(summary))
  if (length(missing) > 0L) {
    stop("coherence summary is missing required columns: ", paste(missing, collapse = ", "))
  }
  summary <- summary[!is.na(summary$target_cell_type) & !is.na(summary$source_cell_type) &
    !is.na(summary$factor), , drop = FALSE]
  if (nrow(summary) == 0L) {
    empty <- matrix(numeric(0), nrow = 0L, ncol = 0L)
    return(list(score_matrix = empty, source_calls = list(), remove_calls = character(),
                threshold = -log10(p_thresh), summary = summary, source_mode = source_mode))
  }
  summary$factor <- as.integer(summary$factor)
  has_patch_p <- all(c("patch_p_value", "patch_neg_log10_p") %in% colnames(summary))
  effective_score_mode <- if (has_patch_p) score_mode else "mean"
  mean_p <- summary$p_value
  patch_p <- if (has_patch_p) summary$patch_p_value else summary$p_value
  if (adjust_p) {
    mean_plot_p <- stats::p.adjust(mean_p, method = "fdr")
    patch_plot_p <- stats::p.adjust(patch_p, method = "fdr")
  } else {
    mean_plot_p <- mean_p
    patch_plot_p <- patch_p
  }
  mean_plot_score <- -log10(mean_plot_p)
  patch_plot_score <- -log10(patch_plot_p)
  if (identical(effective_score_mode, "patch")) {
    summary$plot_p_value <- patch_plot_p
    summary$plot_score <- patch_plot_score
  } else if (identical(effective_score_mode, "combined")) {
    summary$plot_p_value <- pmax(mean_plot_p, patch_plot_p, na.rm = TRUE)
    summary$plot_score <- pmin(mean_plot_score, patch_plot_score, na.rm = TRUE)
  } else {
    summary$plot_p_value <- mean_plot_p
    summary$plot_score <- mean_plot_score
  }
  summary$mean_plot_score <- mean_plot_score
  if (has_patch_p) summary$patch_plot_score <- patch_plot_score
  summary$plot_score[is.nan(summary$plot_score)] <- NA_real_
  summary$plot_score[is.infinite(summary$plot_score) & summary$plot_score < 0] <- NA_real_
  if (is.finite(score_cap) && score_cap > 0) {
    summary$plot_score[is.infinite(summary$plot_score) & summary$plot_score > 0] <- score_cap
    summary$plot_score <- pmin(summary$plot_score, score_cap)
  }

  source_field <- if (source_mode %in% colnames(summary)) source_mode else "plot_score"
  source_summary <- summary[!is.na(summary$source_cell_type) & !is.na(summary$factor), , drop = FALSE]
  source_summary$factor <- as.integer(source_summary$factor)
  source_summary$.source_score <- as.numeric(source_summary[[source_field]])
  source_summary <- source_summary[is.finite(source_summary$.source_score), , drop = FALSE]
  if (nrow(source_summary) == 0L) {
    empty <- matrix(numeric(0), nrow = 0L, ncol = 0L)
    return(list(score_matrix = empty, source_calls = list(), remove_calls = character(),
                threshold = -log10(p_thresh), summary = summary, source_mode = source_field))
  }

  source_score <- stats::aggregate(
    source_summary$.source_score,
    by = list(source_cell_type = source_summary$source_cell_type, factor = source_summary$factor),
    FUN = function(x) {
      if (identical(source_field, "plot_score")) mean(x, na.rm = TRUE) else max(x, na.rm = TRUE)
    }
  )
  names(source_score)[names(source_score) == "x"] <- "source_score"

  source_calls <- list()
  for (factor in sort(unique(source_score$factor))) {
    frows <- source_score[source_score$factor == factor, , drop = FALSE]
    if (nrow(frows) == 0L) next
    best <- which.max(frows$source_score)
    source_calls[[paste0("f_", factor)]] <- frows$source_cell_type[best]
  }

  cell_types <- sort(unique(c(summary$target_cell_type, summary$source_cell_type)))
  factors <- sort(unique(summary$factor))
  mat <- matrix(NA_real_, nrow = length(cell_types), ncol = length(factors),
                dimnames = list(cell_types, paste0("F", factors)))
  summary$coherence_rule_pass <- FALSE
  threshold <- -log10(p_thresh)
  remove_calls <- character()
  for (factor in factors) {
    source <- source_calls[[paste0("f_", factor)]]
    if (is.null(source)) {
      next
    }
    for (ct in cell_types) {
      idx <- summary$factor == factor & summary$target_cell_type == ct &
        summary$source_cell_type == source
      rows <- which(idx)
      if (!length(rows)) {
        next
      }
      vals <- summary$plot_score[rows]
      vals <- vals[is.finite(vals)]
      if (length(vals) > 0L) {
        mat[ct, paste0("F", factor)] <- max(vals)
      }
      pass <- idx & is.finite(summary$plot_score) & summary$plot_score > threshold
      if ("mean_delta_score" %in% names(summary) && !is.null(min_delta_score)) {
        pass <- pass & is.finite(summary$mean_delta_score) &
          summary$mean_delta_score > min_delta_score
      }
      if ("mean_delta_patch_score" %in% names(summary) && !is.null(min_delta_patch_score)) {
        pass <- pass & is.finite(summary$mean_delta_patch_score) &
          summary$mean_delta_patch_score > min_delta_patch_score
      }
      if ("active_fraction" %in% names(summary) && !is.null(min_active_fraction) &&
          min_active_fraction > 0) {
        pass <- pass & is.finite(summary$active_fraction) &
          summary$active_fraction >= min_active_fraction
      }
      if ("mean_largest_patch_fraction" %in% names(summary) &&
          !is.null(min_largest_patch_fraction) && min_largest_patch_fraction > 0) {
        pass <- pass & is.finite(summary$mean_largest_patch_fraction) &
          summary$mean_largest_patch_fraction >= min_largest_patch_fraction
      }
      if ("n_cells" %in% names(summary) && !is.null(min_cells)) {
        pass <- pass & summary$n_cells >= min_cells
      }
      if ("n_molecules" %in% names(summary) && !is.null(min_molecules)) {
        pass <- pass & summary$n_molecules >= min_molecules
      }
      pass <- pass & summary$target_cell_type != source
      if (any(pass)) {
        summary$coherence_rule_pass[pass] <- TRUE
        remove_calls <- c(remove_calls, paste0(factor, "_", ct))
      }
    }
  }

  list(score_matrix = mat, source_calls = source_calls,
       remove_calls = unique(remove_calls), threshold = threshold,
       summary = summary, source_mode = source_field, source_score = source_score,
       score_cap = score_cap, score_mode = effective_score_mode,
       min_delta_score = min_delta_score,
       min_delta_patch_score = min_delta_patch_score,
       min_active_fraction = min_active_fraction, min_cells = min_cells,
       min_largest_patch_fraction = min_largest_patch_fraction,
       min_molecules = min_molecules)
}

#' Convert Coherence Annotations to Correction Rules
#'
#' @inheritParams celladmix_coherence_annotation
#' @param target_cell_types Optional character vector limiting returned target
#'   cell types.
#'
#' @return Data frame of correction rules.
#' @keywords internal
#' @noRd
celladmix_coherence_rules <- function(
    coherence,
    p_thresh = 0.1,
    adjust_p = FALSE,
    target_cell_types = NULL,
    ...
  ) {
  ann <- if (is.list(coherence) && all(c("score_matrix", "source_calls", "remove_calls") %in% names(coherence))) {
    coherence
  } else {
    celladmix_coherence_annotation(coherence, p_thresh = p_thresh, adjust_p = adjust_p, ...)
  }
  empty <- data.frame(
    factor = integer(),
    source_cell_type = character(),
    target_cell_type = character(),
    p_value = numeric(),
    neg_log10_p = numeric(),
    patch_p_value = numeric(),
    patch_neg_log10_p = numeric(),
    mean_delta_score = numeric(),
    mean_delta_patch_score = numeric(),
    mean_largest_patch_fraction = numeric(),
    active_fraction = numeric(),
    n_cells = integer(),
    n_molecules = integer(),
    rule_id = character(),
    stringsAsFactors = FALSE
  )
  if (!length(ann$remove_calls)) {
    return(empty)
  }
  summary <- ann$summary
  out <- lapply(ann$remove_calls, function(rule_id) {
    factor <- as.integer(sub("^([0-9]+)_.*$", "\\1", rule_id))
    target <- sub("^[0-9]+_", "", rule_id)
    if (!is.null(target_cell_types) && !(target %in% target_cell_types)) {
      return(NULL)
    }
    source <- ann$source_calls[[paste0("f_", factor)]] %||% NA_character_
    idx <- summary$factor == factor & summary$target_cell_type == target &
      summary$source_cell_type == source
    rows <- summary[idx, , drop = FALSE]
    if (nrow(rows) > 0L && "coherence_rule_pass" %in% names(rows)) {
      rows <- rows[rows$coherence_rule_pass, , drop = FALSE]
    }
    if (nrow(rows) > 0L && "plot_score" %in% names(rows)) {
      best <- which.max(replace(rows$plot_score, !is.finite(rows$plot_score), -Inf))
    } else if (nrow(rows) > 0L) {
      best <- which.min(rows$p_value)
    } else {
      best <- integer()
    }
    data.frame(
      factor = factor,
      source_cell_type = source,
      target_cell_type = target,
      p_value = if (length(best)) rows$p_value[best] else NA_real_,
      neg_log10_p = if (length(best)) rows$neg_log10_p[best] else NA_real_,
      patch_p_value = if (length(best) && "patch_p_value" %in% names(rows)) rows$patch_p_value[best] else NA_real_,
      patch_neg_log10_p = if (length(best) && "patch_neg_log10_p" %in% names(rows)) rows$patch_neg_log10_p[best] else NA_real_,
      mean_delta_score = if (length(best) && "mean_delta_score" %in% names(rows)) rows$mean_delta_score[best] else NA_real_,
      mean_delta_patch_score = if (length(best) && "mean_delta_patch_score" %in% names(rows)) rows$mean_delta_patch_score[best] else NA_real_,
      mean_largest_patch_fraction = if (length(best) && "mean_largest_patch_fraction" %in% names(rows)) rows$mean_largest_patch_fraction[best] else NA_real_,
      active_fraction = if (length(best) && "active_fraction" %in% names(rows)) rows$active_fraction[best] else NA_real_,
      n_cells = if (length(best) && "n_cells" %in% names(rows)) rows$n_cells[best] else NA_integer_,
      n_molecules = if (length(best) && "n_molecules" %in% names(rows)) rows$n_molecules[best] else NA_integer_,
      rule_id = rule_id,
      stringsAsFactors = FALSE
    )
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) {
    return(empty)
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' Plot Coherence Annotation Heatmap
#'
#' Visualizes local-coherence evidence by target cell type and factor.
#'
#' @inheritParams celladmix_plot_bridge_heatmap
#' @param coherence A `celladmix_coherence_result`, coherence summary data
#'   frame, or object returned by [celladmix_coherence_annotation()].
#'
#' @return Invisibly returns the annotation list for base graphics, or a
#'   `ComplexHeatmap::Heatmap` object when optional dependencies are available.
#' @keywords internal
#' @noRd
celladmix_plot_coherence_heatmap <- function(
    coherence,
    p_thresh = 0.1,
    adjust_p = FALSE,
    main = NULL,
    max_neg_log10 = 6,
    source_prior = NULL,
    source_prior_min_score = 0.15,
    source_prior_min_margin = 0.03,
    use_complexheatmap = FALSE,
    ...
  ) {
  main <- main %||% sprintf("Coherence factor annotation, p < %.2g", p_thresh)
  ann <- if (is.list(coherence) && all(c("score_matrix", "source_calls", "remove_calls") %in% names(coherence))) {
    coherence
  } else {
    celladmix_coherence_annotation(coherence, p_thresh = p_thresh, adjust_p = adjust_p, ...)
  }
  celladmix_plot_bridge_heatmap(
    ann,
    p_thresh = p_thresh,
    adjust_p = adjust_p,
    main = main,
    max_neg_log10 = max_neg_log10,
    source_prior = source_prior,
    source_prior_min_score = source_prior_min_score,
    source_prior_min_margin = source_prior_min_margin,
    use_complexheatmap = use_complexheatmap
  )
}

#' Write a Corrected child Run
#'
#' Applies factor-removal rules to a persisted run and writes a corrected child
#' run directory.
#'
#' @param run A `celladmix_run` object.
#' @param rules Data frame of correction rules. It must contain `factor` and
#'   `target_cell_type` columns.
#' @param out_dir Optional output directory for the corrected child run. If
#'   `NULL`, a timestamped directory is created under the parent run's
#'   `corrected` directory.
#' @param annotation Optional cell annotation as a named vector, data frame, or
#'   CSV path. If supplied, labels override run cell types for rule targeting.
#' @param annotation_col Candidate annotation columns to use when `annotation`
#'   is a data frame or CSV path.
#' @param cell_id_col Cell identifier column in `annotation`.
#'
#' @return A `celladmix_correction_run` object.
#' @keywords internal
#' @noRd
celladmix_correct <- function(
    run,
    rules,
    out_dir = NULL,
    annotation = NULL,
    annotation_col = c("merged_annotation", "cell_type", "cluster_label", "cluster"),
    cell_id_col = "cell_id"
  ) {
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_correct() expects a celladmix_run")
  }
  stopifnot(is.data.frame(rules))
  cell_types <- .celladmix_normalize_cell_type_override(
    annotation,
    annotation_col = annotation_col,
    cell_id_col = cell_id_col
  )
  out_dir <- out_dir %||% file.path(run$paths$corrected_dir, paste0("corrected_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  .celladmix_make_run_object(
    .celladmix_correct_run(run$path, rules, out_dir = out_dir, cell_types = cell_types),
    class = "celladmix_correction_run"
  )
}

#' Collect Correction Summary
#'
#' Reads the per-cell molecule removal summary written by
#' [celladmix_correct()].
#'
#' @param run A corrected `celladmix_run` object.
#'
#' @return Data frame with molecule counts before and after correction.
#' @keywords internal
#' @noRd
celladmix_collect_correction_summary <- function(run) {
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_collect_correction_summary() expects a celladmix_run")
  }
  .celladmix_collect_correction_summary(run$path)
}

#' Plot Molecule Removals from a Corrected Run
#'
#' Draws compact boxplots for modified cells using the per-cell correction
#' summary. The default `"summary"` mode shows both before/after molecule
#' counts and the number of removed molecules.
#'
#' @param x A correction summary data frame or corrected `celladmix_run`.
#' @param mode Either `"summary"`, `"before_after"`, `"removed"`, or
#'   `"fraction_removed"`.
#' @param modified_only If `TRUE`, plot only cells with at least one removed
#'   molecule.
#' @param title Optional plot title.
#'
#' @return A `ggplot2` object when ggplot2 is available; otherwise the plotted
#'   summary data, invisibly.
#' @export
celladmix_plot_molecule_removal <- function(
    x,
    mode = c("summary", "before_after", "removed", "fraction_removed"),
    modified_only = TRUE,
    title = NULL
  ) {
  mode <- match.arg(mode)
  summary <- if (inherits(x, "celladmix_run")) celladmix_collect_correction_summary(x) else x
  if (!is.data.frame(summary)) {
    stop("x must be a correction summary data frame or celladmix_run")
  }
  required <- c("cell_type", "n_molecules_before", "n_molecules_after", "n_removed", "fraction_removed")
  missing <- setdiff(required, names(summary))
  if (length(missing) > 0L) {
    stop("correction summary is missing columns: ", paste(missing, collapse = ", "))
  }
  plot_df <- summary
  if (modified_only) {
    plot_df <- plot_df[plot_df$n_removed > 0, , drop = FALSE]
  }
  if (nrow(plot_df) == 0L) {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      return(.celladmix_empty_plot("No modified cells"))
    }
    graphics::plot.new()
    graphics::title(title %||% "No modified cells")
    return(invisible(plot_df))
  }
  plot_df$cell_type <- ifelse(is.na(plot_df$cell_type) | plot_df$cell_type == "",
    "unknown", plot_df$cell_type)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    if (identical(mode, "before_after") || identical(mode, "summary")) {
      counts <- rbind(
        data.frame(cell_type = plot_df$cell_type, state = "original",
          molecules = plot_df$n_molecules_before, stringsAsFactors = FALSE),
        data.frame(cell_type = plot_df$cell_type, state = "cleaned",
          molecules = plot_df$n_molecules_after, stringsAsFactors = FALSE)
      )
      counts$state <- factor(counts$state, levels = c("original", "cleaned"))
      p_counts <- ggplot2::ggplot(counts, ggplot2::aes(x = cell_type, y = molecules,
        fill = state)) +
        ggplot2::geom_boxplot(outlier.size = 0.25, linewidth = 0.25) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = c(original = "grey70", cleaned = "#2C7FB8")) +
        ggplot2::xlab(NULL) +
        ggplot2::ylab("Molecules per modified cell") +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(legend.title = ggplot2::element_blank())
      if (identical(mode, "before_after")) {
        return(p_counts + ggplot2::ggtitle(title %||% "Molecule counts in modified cells"))
      }
    }
    if (identical(mode, "removed") || identical(mode, "summary")) {
      p_removed <- ggplot2::ggplot(plot_df, ggplot2::aes(x = cell_type, y = n_removed)) +
        ggplot2::geom_boxplot(fill = "#F28E2B", outlier.size = 0.25, linewidth = 0.25) +
        ggplot2::coord_flip() +
        ggplot2::xlab(NULL) +
        ggplot2::ylab("Removed molecules per cell") +
        ggplot2::theme_classic(base_size = 10)
      if (identical(mode, "removed")) {
        return(p_removed + ggplot2::ggtitle(title %||% "Molecules removed from modified cells"))
      }
    }
    if (identical(mode, "fraction_removed")) {
      return(ggplot2::ggplot(plot_df, ggplot2::aes(x = cell_type, y = fraction_removed)) +
        ggplot2::geom_boxplot(fill = "#F28E2B", outlier.size = 0.25, linewidth = 0.25) +
        ggplot2::coord_flip() +
        ggplot2::xlab(NULL) +
        ggplot2::ylab("Fraction removed") +
        ggplot2::ggtitle(title %||% "Fraction of molecules removed") +
        ggplot2::theme_classic(base_size = 10))
    }
    if (requireNamespace("cowplot", quietly = TRUE)) {
      return(cowplot::plot_grid(p_counts, p_removed, ncol = 1, align = "v"))
    }
    return(p_counts)
  }

  if (identical(mode, "before_after") || identical(mode, "summary")) {
    values <- c(plot_df$n_molecules_before, plot_df$n_molecules_after)
    state <- rep(c("original", "cleaned"), each = nrow(plot_df))
    group <- interaction(rep(plot_df$cell_type, 2L), state, sep = ": ", lex.order = TRUE)
    graphics::boxplot(values ~ group, las = 2, ylab = "Molecules per modified cell",
      main = title %||% "Molecule counts in modified cells")
  } else if (identical(mode, "removed")) {
    graphics::boxplot(plot_df$n_removed ~ plot_df$cell_type, las = 2,
      ylab = "Removed molecules per cell",
      main = title %||% "Molecules removed from modified cells")
  } else {
    graphics::boxplot(plot_df$fraction_removed ~ plot_df$cell_type, las = 2,
      ylab = "Fraction removed", main = title %||% "Fraction of molecules removed")
  }
  invisible(plot_df)
}

#' Select a Factor by Marker Enrichment
#'
#' Chooses the factor with the largest summed loading across the provided marker
#' genes.
#'
#' @param run A `celladmix_run` object.
#' @param marker_genes Character vector of marker genes to evaluate.
#'
#' @return A one-based factor index.
#' @keywords internal
#' @noRd
celladmix_select_factor_by_markers <- function(run, marker_genes) {
  stopifnot(!is.null(run$h), !is.null(run$genes))
  marker_idx <- which(run$genes %in% marker_genes)
  if (length(marker_idx) == 0) {
    stop("No marker genes were found in run$genes")
  }
  scores <- rowSums(run$h[, marker_idx, drop = FALSE])
  which.max(scores)
}

#' Plot a Cell Clustering Result
#'
#' Quick base-R plot for a `celladmix_clusters` object.
#'
#' @param x A `celladmix_clusters` object.
#' @param color_by Either `"cluster"` or `"analysis_crop"`.
#' @param ... Additional arguments passed to [graphics::plot()].
#'
#' @return The input object, invisibly.
#' @keywords internal
#' @noRd
plot.celladmix_clusters <- function(x, color_by = c("cluster", "analysis_crop"), ...) {
  color_by <- match.arg(color_by)
  df <- x$embedding
  if (nrow(df) == 0) {
    plot.new()
    return(invisible(x))
  }
  value <- df[[if (identical(color_by, "cluster")) "cluster" else "analysis_crop"]]
  palette <- grDevices::rainbow(length(unique(value)))
  colors <- grDevices::adjustcolor(
    palette[match(value, unique(value))],
    alpha.f = 0.5
  )
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(3.5, 3.5, 2.0, 0.5), mgp = c(2, 0.65, 0), cex = 2/3)
  plot_args <- c(
    list(
      x = df$umap_1,
      y = df$umap_2,
      col = colors,
      pch = 16,
      cex = 0.4,
      asp = 1,
      xlab = "UMAP 1",
      ylab = "UMAP 2"
    ),
    list(...)
  )
  do.call(graphics::plot, plot_args)
  invisible(x)
}
