# Optional Seurat integration helpers. These functions keep Seurat handling on
# the R side and pass a normal Xenium-backed prep/store to the native core.

.celladmix_require_seuratobject <- function() {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("The SeuratObject package is required for Seurat integration")
  }
}

.celladmix_seurat_assays <- function(object) {
  .celladmix_require_seuratobject()
  as.character(SeuratObject::Assays(object))
}

.celladmix_resolve_seurat_assay <- function(object, assay = NULL) {
  assays <- .celladmix_seurat_assays(object)
  if (!length(assays)) {
    stop("Seurat object does not contain any assays")
  }
  if (!is.null(assay)) {
    assay <- as.character(assay)[[1]]
    if (!(assay %in% assays)) {
      stop("Seurat assay '", assay, "' was not found. Available assays: ",
        paste(assays, collapse = ", "))
    }
    return(assay)
  }
  default <- tryCatch(SeuratObject::DefaultAssay(object), error = function(e) NULL)
  if (!is.null(default) && default %in% assays) {
    return(default)
  }
  if ("Xenium" %in% assays) {
    return("Xenium")
  }
  if (length(assays) == 1L) {
    return(assays[[1]])
  }
  stop("Could not infer Seurat assay. Pass assay=. Available assays: ",
    paste(assays, collapse = ", "))
}

.celladmix_seurat_images <- function(object, assay = NULL) {
  .celladmix_require_seuratobject()
  out <- tryCatch(SeuratObject::Images(object, assay = assay), error = function(e) NULL)
  if (is.null(out) || !length(out)) {
    out <- tryCatch(SeuratObject::Images(object), error = function(e) character())
  }
  as.character(out)
}

.celladmix_resolve_seurat_image <- function(object, assay = NULL, image = NULL) {
  images <- .celladmix_seurat_images(object, assay = assay)
  if (!length(images)) {
    return(NULL)
  }
  if (!is.null(image)) {
    image <- as.character(image)[[1]]
    if (!(image %in% images)) {
      stop("Seurat image/FOV '", image, "' was not found. Available images: ",
        paste(images, collapse = ", "))
    }
    return(image)
  }
  default_fov <- tryCatch(SeuratObject::DefaultFOV(object, assay = assay),
    error = function(e) NULL)
  if (!is.null(default_fov) && default_fov %in% images) {
    return(default_fov)
  }
  if ("fov" %in% images) {
    return("fov")
  }
  if (length(images) == 1L) {
    return(images[[1]])
  }
  NULL
}

.celladmix_seurat_cells <- function(object, assay) {
  cells <- colnames(object[[assay]])
  if (is.null(cells) || !length(cells)) {
    cells <- tryCatch(SeuratObject::Cells(object, assay = assay), error = function(e) NULL)
  }
  cells <- as.character(cells)
  cells[!is.na(cells) & nzchar(cells)]
}

.celladmix_seurat_features <- function(object, assay) {
  features <- rownames(object[[assay]])
  features <- as.character(features)
  features[!is.na(features) & nzchar(features)]
}

.celladmix_seurat_metadata <- function(object) {
  meta <- object[[]]
  if (!is.data.frame(meta)) {
    stop("Could not extract Seurat metadata")
  }
  meta
}

.celladmix_named_annotation_from_seurat <- function(object, annotation, cells) {
  meta <- .celladmix_seurat_metadata(object)
  make_labels <- function(labels, name) {
    labels <- stats::setNames(as.character(labels), rownames(meta))
    labels <- labels[cells]
    keep <- !is.na(labels) & nzchar(labels)
    list(labels = labels[keep], name = .celladmix_clean_name(name, "annotation"))
  }

  if (is.null(annotation)) {
    ids <- tryCatch(SeuratObject::Idents(object), error = function(e) NULL)
    if (!is.null(ids)) {
      labels <- stats::setNames(as.character(ids), names(ids))
      labels <- labels[cells]
      values <- unique(labels[!is.na(labels) & nzchar(labels)])
      if (length(values) > 1L) {
        return(list(labels = labels[!is.na(labels) & nzchar(labels)], name = "ident"))
      }
    }
    stop("No usable Seurat annotation was supplied. Pass annotation= as a metadata field name, ",
      "a named vector, or set informative Idents(object).")
  }

  if (is.character(annotation) && length(annotation) == 1L && !file.exists(annotation)) {
    field <- annotation[[1]]
    if (field %in% c("ident", "idents", "active.ident")) {
      ids <- SeuratObject::Idents(object)
      labels <- stats::setNames(as.character(ids), names(ids))
      labels <- labels[cells]
      keep <- !is.na(labels) & nzchar(labels)
      return(list(labels = labels[keep], name = "ident"))
    }
    if (!(field %in% colnames(meta))) {
      stop("Seurat metadata field '", field, "' was not found. Available fields: ",
        paste(colnames(meta), collapse = ", "))
    }
    return(make_labels(meta[[field]], field))
  }

  if (is.atomic(annotation) && !is.null(names(annotation))) {
    labels <- stats::setNames(as.character(annotation), names(annotation))
    labels <- labels[cells]
    keep <- !is.na(labels) & nzchar(labels)
    return(list(labels = labels[keep], name = "seurat_annotation"))
  }

  if (is.data.frame(annotation) || (is.character(annotation) && length(annotation) == 1L)) {
    ann <- .celladmix_normalize_annotation_table(annotation)
    labels <- ann$labels[cells]
    keep <- !is.na(labels) & nzchar(labels)
    return(list(labels = labels[keep], name = ann$name))
  }

  stop("Unsupported Seurat annotation specification")
}

#' Prepare a Seurat Object for cellAdmix
#'
#' Uses the Seurat object to define the active cells, active genes, and
#' annotation, while using an explicit Xenium bundle for molecule-complete input.
#'
#' @keywords internal
#' @noRd
celladmix_prepare_seurat <- function(
    object,
    xenium_dir = NULL,
    output_dir = NULL,
    annotation = NULL,
    assay = NULL,
    image = NULL,
    ...
  ) {
  .celladmix_require_seuratobject()
  if (is.null(xenium_dir)) {
    warning("cellAdmix Seurat integration works best with explicit xenium_dir=. ",
      "Seurat stores only a reduced molecule coordinate representation; full fitting ",
      "requires the original Xenium bundle in this implementation.", call. = FALSE)
    stop("Pass xenium_dir= to run cellAdmix from a Seurat object.")
  }

  assay <- .celladmix_resolve_seurat_assay(object, assay = assay)
  image <- .celladmix_resolve_seurat_image(object, assay = assay, image = image)
  cells <- .celladmix_seurat_cells(object, assay)
  features <- .celladmix_seurat_features(object, assay)
  if (!length(cells)) {
    stop("Resolved Seurat assay '", assay, "' does not contain any cells")
  }
  if (!length(features)) {
    stop("Resolved Seurat assay '", assay, "' does not contain any features")
  }
  ann <- .celladmix_named_annotation_from_seurat(object, annotation, cells = cells)
  if (!length(ann$labels)) {
    stop("Seurat annotation did not label any active cells")
  }

  prep <- celladmix_prepare_xenium(
    bundle_dir = xenium_dir,
    output_dir = output_dir,
    cell_filter = cells,
    gene_filter = features,
    ...
  )
  prep$source$seurat <- list(
    assay = assay,
    image = image %||% NA_character_,
    n_cells = length(cells),
    n_features = length(features),
    annotation = ann$name
  )
  prep$defaults$cell_filter_source <- "seurat"
  prep$defaults$gene_filter_source <- "seurat"
  .celladmix_write_json(prep, prep$paths$project_json)

  list(
    prep = prep,
    annotation = ann$labels,
    annotation_name = ann$name,
    assay = assay,
    image = image
  )
}

.celladmix_seurat_add_metadata <- function(object, metadata, prefix = "celladmix_") {
  .celladmix_require_seuratobject()
  if (!is.data.frame(metadata) || !nrow(metadata)) {
    return(object)
  }
  if (!("cell_id" %in% colnames(metadata))) {
    stop("metadata must contain a cell_id column")
  }
  rownames(metadata) <- as.character(metadata$cell_id)
  metadata$cell_id <- NULL
  common <- intersect(colnames(object), rownames(metadata))
  if (!length(common)) {
    warning("No Seurat cells matched cellAdmix metadata", call. = FALSE)
    return(object)
  }
  metadata <- metadata[common, , drop = FALSE]
  colnames(metadata) <- paste0(prefix, colnames(metadata))
  SeuratObject::AddMetaData(object, metadata = metadata)
}

#' Add cellAdmix Factor Metadata to a Seurat Object
#'
#' @param object A Seurat object.
#' @param fit A `CellAdmixFit` object.
#' @param prefix Metadata column prefix.
#'
#' @return The Seurat object with factor metadata columns added.
#' @export
celladmix_add_factors_to_seurat <- function(object, fit, prefix = "celladmix_") {
  if (!inherits(fit, "CellAdmixFit")) {
    stop("fit must be a CellAdmixFit object")
  }
  cells <- fit$cell_factors()
  keep <- c("cell_id", "dominant_factor", grep("^factor_[0-9]+_fraction$", names(cells), value = TRUE))
  .celladmix_seurat_add_metadata(object, cells[, keep, drop = FALSE], prefix = prefix)
}

.celladmix_reorder_counts_for_seurat <- function(counts, object) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required to add corrected counts to Seurat")
  }
  cells <- colnames(object)
  if (is.null(cells) || !length(cells)) {
    stop("Seurat object does not contain cells")
  }
  common <- intersect(cells, colnames(counts))
  if (!length(common)) {
    stop("Corrected counts do not contain any cells from the Seurat object")
  }
  if (!all(cells %in% colnames(counts))) {
    full <- Matrix::sparseMatrix(
      i = integer(),
      j = integer(),
      x = numeric(),
      dims = c(nrow(counts), length(cells)),
      dimnames = list(rownames(counts), cells)
    )
    full[, common] <- counts[, common, drop = FALSE]
    counts <- full
  } else {
    counts <- counts[, cells, drop = FALSE]
  }
  counts
}

#' Add Corrected cellAdmix Counts as a Seurat Assay
#'
#' @param object A Seurat object.
#' @param correction A `CellAdmixCorrection` object.
#' @param assay Name of the assay to create.
#'
#' @return The Seurat object with a corrected sparse counts assay.
#' @export
celladmix_add_corrected_assay <- function(object, correction, assay = "cellAdmix") {
  .celladmix_require_seuratobject()
  if (!inherits(correction, "CellAdmixCorrection")) {
    stop("correction must be a CellAdmixCorrection object")
  }
  counts <- .celladmix_reorder_counts_for_seurat(correction$counts(), object)
  create_assay <- if ("CreateAssay5Object" %in% getNamespaceExports("SeuratObject")) {
    SeuratObject::CreateAssay5Object
  } else {
    SeuratObject::CreateAssayObject
  }
  object[[assay]] <- create_assay(counts = counts)
  object
}

.celladmix_add_correction_metadata_to_seurat <- function(object, correction, prefix = "celladmix_") {
  summary <- correction$cell_summary()
  if (!is.data.frame(summary) || !nrow(summary) || !("cell_id" %in% names(summary))) {
    return(object)
  }
  keep <- intersect(c("cell_id", "n_removed", "fraction_removed",
    "n_molecules_before", "n_molecules_after"), names(summary))
  .celladmix_seurat_add_metadata(object, summary[, keep, drop = FALSE],
    prefix = paste0(prefix, correction$name, "_"))
}

.celladmix_store_seurat_provenance <- function(object, fit = NULL, score = NULL,
                                               correction = NULL, assay = NULL) {
  misc <- methods::slot(object, "misc")
  entry <- misc$cellAdmix %||% list()
  if (!is.null(fit)) {
    entry$last_fit <- fit$metadata()
  }
  if (!is.null(score)) {
    entry$last_score <- score$metadata()
  }
  if (!is.null(correction)) {
    entry$last_correction <- correction$metadata()
  }
  if (!is.null(assay)) {
    entry$corrected_assay <- assay
  }
  misc$cellAdmix <- entry
  methods::slot(object, "misc") <- misc
  object
}

#' Add cellAdmix Results Back to a Seurat Object
#'
#' Adds fitted cell factors, correction metadata, and optionally a corrected
#' counts assay to a Seurat object.
#'
#' @param object A Seurat object.
#' @param fit Optional `CellAdmixFit` object.
#' @param score Optional `CellAdmixScore` object.
#' @param correction Optional `CellAdmixCorrection` object.
#' @param corrected_assay Name of corrected assay to create, or `FALSE` to skip.
#' @param prefix Metadata column prefix.
#'
#' @return Updated Seurat object.
#' @export
celladmix_add_to_seurat <- function(
    object,
    fit = NULL,
    score = NULL,
    correction = NULL,
    corrected_assay = "cellAdmix",
    prefix = "celladmix_"
  ) {
  .celladmix_require_seuratobject()
  if (!is.null(fit)) {
    object <- celladmix_add_factors_to_seurat(object, fit = fit, prefix = prefix)
  }
  if (!is.null(correction)) {
    object <- .celladmix_add_correction_metadata_to_seurat(object, correction, prefix = prefix)
    if (!identical(corrected_assay, FALSE)) {
      object <- celladmix_add_corrected_assay(object, correction, assay = corrected_assay)
    }
  }
  .celladmix_store_seurat_provenance(object, fit = fit, score = score,
    correction = correction, assay = if (identical(corrected_assay, FALSE)) NULL else corrected_assay)
}
