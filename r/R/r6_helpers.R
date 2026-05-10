.celladmix_clean_name <- function(x, fallback = "item") {
  x <- as.character(x %||% fallback)[[1]]
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) fallback else x
}

.celladmix_first_existing_col <- function(x, candidates) {
  hits <- candidates[candidates %in% colnames(x)]
  if (length(hits)) hits[[1]] else NULL
}

.celladmix_annotation_hash <- function(labels) {
  labels <- stats::setNames(as.character(labels), names(labels))
  ord <- order(names(labels), labels, na.last = TRUE)
  tmp <- tempfile("celladmix_annotation_", fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(paste(names(labels)[ord], labels[ord], sep = "\t"), tmp, useBytes = TRUE)
  unname(tools::md5sum(tmp))
}

.celladmix_normalize_annotation_table <- function(
    annotation,
    name = NULL,
    cell_id_col = "cell_id",
    label_col = NULL
  ) {
  if (inherits(annotation, "celladmix_annotation")) {
    return(annotation)
  }
  if (is.character(annotation) && length(annotation) == 1L && file.exists(annotation)) {
    annotation <- utils::read.csv(annotation, stringsAsFactors = FALSE)
  }

  if (is.atomic(annotation) && !is.null(names(annotation))) {
    cells <- names(annotation)
    labels <- as.character(annotation)
  } else if (is.data.frame(annotation)) {
    if (!(cell_id_col %in% colnames(annotation))) {
      stop("annotation is missing cell ID column: ", cell_id_col)
    }
    if (is.null(label_col)) {
      label_col <- .celladmix_first_existing_col(
        annotation,
        c("annotation", "merged_annotation", "cell_type", "celltype", "label",
          "cluster_label", "cluster")
      )
      if (is.null(label_col)) {
        remaining <- setdiff(colnames(annotation), cell_id_col)
        if (length(remaining) == 1L) {
          label_col <- remaining[[1]]
        }
      }
    }
    if (is.null(label_col) || !(label_col %in% colnames(annotation))) {
      stop("Could not infer annotation label column; pass annotation_col once at construction or set_annotation().")
    }
    cells <- as.character(annotation[[cell_id_col]])
    labels <- as.character(annotation[[label_col]])
  } else {
    stop("annotation must be a named vector, data frame, CSV path, or celladmix_annotation")
  }

  keep <- !is.na(cells) & nzchar(cells) & !is.na(labels) & nzchar(labels)
  cells <- cells[keep]
  labels <- labels[keep]
  if (!length(cells)) {
    stop("annotation did not contain any usable cell labels")
  }
  keep_first <- !duplicated(cells)
  labels <- stats::setNames(labels[keep_first], cells[keep_first])
  levels <- sort(unique(unname(labels)))
  list(
    name = .celladmix_clean_name(name %||% "annotation", "annotation"),
    labels = labels,
    table = data.frame(cell_id = names(labels), label = unname(labels), stringsAsFactors = FALSE),
    n_cells = length(labels),
    n_labeled = length(labels),
    n_labels = length(levels),
    label_counts = sort(table(unname(labels)), decreasing = TRUE),
    levels = levels,
    hash = .celladmix_annotation_hash(labels),
    path = NULL
  )
}

.celladmix_write_annotation_sidecar <- function(annotation, annotations_dir, row_group_size = 65536L) {
  if (!inherits(annotation, "celladmix_annotation")) {
    annotation <- structure(annotation, class = "celladmix_annotation")
  }
  label_dir <- file.path(annotations_dir, "r6", .celladmix_clean_name(annotation$name, "annotation"))
  dir.create(label_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(label_dir, "cell_labels.parquet")
  .celladmix_write_cell_labels(path, annotation$labels, row_group_size = row_group_size)
  annotation$path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  metadata <- annotation[c("name", "n_cells", "n_labeled", "n_labels",
    "label_counts", "levels", "hash", "path")]
  metadata$labels_path <- annotation$path
  .celladmix_write_json(metadata, file.path(label_dir, "annotation.json"))
  structure(annotation, class = "celladmix_annotation")
}

.celladmix_make_annotation <- function(annotation, annotations_dir, name = NULL,
                                       cell_id_col = "cell_id", label_col = NULL) {
  out <- .celladmix_normalize_annotation_table(
    annotation,
    name = name,
    cell_id_col = cell_id_col,
    label_col = label_col
  )
  structure(.celladmix_write_annotation_sidecar(out, annotations_dir), class = "celladmix_annotation")
}

.celladmix_annotation_rank <- function(annotation, multiplier = 1.2, cap = 30L) {
  n_labels <- as.integer(annotation$n_labels %||% 0L)
  if (is.na(n_labels) || n_labels <= 0L) {
    stop("Cannot auto-select rank without annotation labels")
  }
  min(as.integer(cap), max(2L, ceiling(multiplier * n_labels)))
}

.celladmix_infer_format <- function(source, format = "auto", molecules = NULL) {
  if (!identical(format, "auto")) {
    return(format)
  }
  if (inherits(source, "celladmix_prep")) {
    return("prep")
  }
  if (inherits(source, "celladmix_store")) {
    return("store")
  }
  if (inherits(source, "Seurat")) {
    return("seurat")
  }
  if (!is.null(molecules)) {
    return("tabular")
  }
  if (is.character(source) && length(source) == 1L) {
    path <- source[[1]]
    if (dir.exists(path) && file.exists(file.path(path, "store.json"))) {
      return("store")
    }
    if (dir.exists(path) && file.exists(file.path(path, "experiment.xenium"))) {
      return("xenium")
    }
    if (file.exists(path) && identical(basename(path), "experiment.xenium")) {
      return("xenium")
    }
    if (file.exists(path)) {
      return("tabular")
    }
  }
  stop("Could not infer input format; pass format = 'xenium' or format = 'tabular'")
}

.celladmix_run_id <- function(prefix, annotation, rank, variant) {
  paste(
    .celladmix_clean_name(prefix, "fit"),
    .celladmix_clean_name(annotation$name, "annotation"),
    paste0("rank", rank),
    .celladmix_clean_name(variant, "nmf"),
    sep = "_"
  )
}

.celladmix_score_dir <- function(run, name) {
  base <- run$paths$scores_dir %||% file.path(run$path, "scores")
  file.path(base, .celladmix_clean_name(name, "score"))
}

.celladmix_copy_score_outputs <- function(result, score_dir, method) {
  dir.create(score_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- result$paths %||% list()
  if (!is.null(paths$scores) && file.exists(paths$scores)) {
    file.copy(paths$scores, file.path(score_dir, "pairs.parquet"), overwrite = TRUE)
    result$paths$scores <- normalizePath(file.path(score_dir, "pairs.parquet"), winslash = "/", mustWork = FALSE)
  }
  if (!is.null(paths$summary) && file.exists(paths$summary)) {
    file.copy(paths$summary, file.path(score_dir, "summary.parquet"), overwrite = TRUE)
    result$paths$summary <- normalizePath(file.path(score_dir, "summary.parquet"), winslash = "/", mustWork = FALSE)
  }
  result$paths$score_dir <- normalizePath(score_dir, winslash = "/", mustWork = FALSE)
  result$method <- method
  result
}

#' High-Level cellAdmix Dataset
#'
#' R6 object storing source metadata, active annotations, fitted runs, scores,
