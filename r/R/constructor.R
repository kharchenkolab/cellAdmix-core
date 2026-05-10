#' Construct a cellAdmix Dataset Object
#'
#' Creates a high-level R6 dataset object from a Xenium bundle, tabular source,
#' or existing prepared project.
#'
#' @param source Xenium bundle path, `experiment.xenium` path, tabular molecule
#'   path when `format = "tabular"`, or an existing `celladmix_prep` object.
#' @param output_dir Output/cache directory.
#' @param annotation Optional initial cell annotation.
#' @param format Input format. `"auto"` detects Xenium bundles and prep objects.
#' @param schema Optional [celladmix_schema()] for tabular inputs.
#' @param annotation_name Name for the initial annotation.
#' @param annotation_col Optional one-time label column override for
#'   `annotation` data frames.
#' @param cell_id_col Cell ID column for annotation data frames.
#' @param num_threads Default number of worker threads used by dataset methods
#'   when a method call does not specify `num_threads`.
#' @param molecules Tabular molecule path when `source` is not used.
#' @param xenium_dir Explicit Xenium bundle path used when `source` is a Seurat
#'   object. The Seurat object supplies cells, features, and annotations; the
#'   bundle supplies molecule-complete transcript records.
#' @param assay,image Optional Seurat assay and FOV image names. When omitted,
#'   sensible Seurat defaults are used.
#' @param ... Passed to source-specific prepare functions.
#'
#' @return A `CellAdmixDataset` object.
#' @export
cellAdmix <- function(
    source = NULL,
    output_dir = NULL,
    annotation = NULL,
    format = c("auto", "xenium", "tabular", "prep", "store", "seurat"),
    schema = NULL,
    annotation_name = "manual",
    annotation_col = NULL,
    cell_id_col = "cell_id",
    num_threads = NULL,
    molecules = NULL,
    xenium_dir = NULL,
    assay = NULL,
    image = NULL,
    ...
  ) {
  format <- match.arg(format)
  detected <- .celladmix_infer_format(source, format = format, molecules = molecules)
  dots <- list(...)
  num_threads <- as.integer(num_threads %||% .celladmix_default_threads())
  if ("project_dir" %in% names(dots)) {
    stop("Use output_dir= instead of project_dir=")
  }

  if (identical(detected, "prep")) {
    prep <- source
  } else if (identical(detected, "xenium")) {
    prep <- do.call(celladmix_prepare_xenium, c(
      list(bundle_dir = source, output_dir = output_dir),
      dots
    ))
  } else if (identical(detected, "tabular")) {
    schema <- schema %||% celladmix_schema()
    molecules_path <- molecules %||% source
    if (is.null(molecules_path)) {
      stop("Tabular input requires a molecules path via source or molecules=")
    }
    prep <- do.call(celladmix_prepare_tabular, c(
      list(
        molecules_path = molecules_path,
        output_dir = output_dir,
        x_col = schema$x,
        y_col = schema$y,
        z_col = schema$z,
        gene_col = schema$gene,
        qv_col = schema$qv,
        cell_id_col = schema$cell,
        segmentation_mask_path = schema$segmentation_mask,
        cell_type_col = schema$cell_type,
        cell_metadata_path = schema$cell_metadata,
        cell_metadata_cell_id_col = schema$cell_metadata_cell,
        cell_metadata_cell_type_col = schema$cell_metadata_cell_type,
        sample_id_col = schema$sample_id,
        fov_id_col = schema$fov_id
      ),
      dots
    ))
    if (is.null(annotation) && !is.null(schema$cell_metadata) &&
        !is.null(schema$cell_metadata_cell_type)) {
      annotation <- schema$cell_metadata
      annotation_name <- schema$cell_metadata_cell_type
      annotation_col <- schema$cell_metadata_cell_type
      cell_id_col <- schema$cell_metadata_cell %||% cell_id_col
    }
  } else if (identical(detected, "seurat")) {
    prepared <- do.call(celladmix_prepare_seurat, c(
      list(
        object = source,
        xenium_dir = xenium_dir,
        output_dir = output_dir,
        annotation = annotation,
        assay = assay,
        image = image
      ),
      dots
    ))
    prep <- prepared$prep
    annotation <- prepared$annotation
    annotation_name <- prepared$annotation_name
    annotation_col <- NULL
    cell_id_col <- "cell_id"
  } else if (identical(detected, "store")) {
    stop("Direct R6 construction from an existing input store is planned but not implemented yet; pass the original source or a celladmix_prep object.")
  } else {
    stop("Unsupported input format: ", detected)
  }

  CellAdmixDataset$new(
    prep = prep,
    format = prep$source$type %||% detected,
    source_info = prep$source,
    annotation = annotation,
    annotation_name = annotation_name,
    annotation_col = annotation_col,
    cell_id_col = cell_id_col,
    num_threads = num_threads
  )
}
