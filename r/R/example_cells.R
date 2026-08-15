#' Example-Cell Diagnostics
#'
#' These helpers draw molecule-level views around selected target cells. They
#' are intended for notebook diagnostics where fitted factor labels, optional
#' stain crops, and cell boundaries should be shown together.

.celladmix_as_factor_id <- function(x) {
  if (is.null(x) || !length(x)) {
    return(integer())
  }
  if (is.character(x)) {
    x <- sub("^[Ff]_?", "", x)
  }
  out <- suppressWarnings(as.integer(x))
  out[is.finite(out)]
}

.celladmix_factor_label <- function(x) {
  paste0("F", .celladmix_as_factor_id(x))
}

.celladmix_factor_sources <- function(score_annotation) {
  calls <- score_annotation$source_calls
  if (is.null(calls) || !length(calls)) {
    return(data.frame(factor = integer(), source_cell_type = character(),
      stringsAsFactors = FALSE))
  }
  data.frame(factor = .celladmix_as_factor_id(names(calls)),
    source_cell_type = unname(unlist(calls, use.names = FALSE)),
    stringsAsFactors = FALSE)
}

#' Discover a Cell-Boundary File for a Fitted Run
#'
#' Currently this looks for the standard Xenium `cell_boundaries` files next to
#' the source bundle recorded in the run metadata.
#'
#' @param run A `celladmix_run` object.
#' @param boundary_path Optional explicit boundary path.
#'
#' @return Normalized path, or `NA_character_` when no candidate exists.
#' @export
celladmix_discover_cell_boundaries <- function(run, boundary_path = NULL) {
  if (!is.null(boundary_path) && length(boundary_path) && nzchar(boundary_path[[1]])) {
    if (!file.exists(boundary_path[[1]])) {
      stop("Cell boundary path does not exist: ", boundary_path[[1]])
    }
    return(normalizePath(boundary_path[[1]], winslash = "/", mustWork = TRUE))
  }
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_discover_cell_boundaries() expects a celladmix_run")
  }
  source_path <- run$source$path
  if (is.null(source_path) || !nzchar(source_path)) {
    return(NA_character_)
  }
  source_dir <- if (dir.exists(source_path)) source_path else dirname(source_path)
  candidates <- file.path(source_dir, c("cell_boundaries.csv.gz",
    "cell_boundaries.csv", "cell_boundaries.parquet"))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) {
    normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
  } else {
    NA_character_
  }
}

.celladmix_boundary_columns <- function(boundaries) {
  cell_col <- c("cell_id", "cell", "cellID")[c("cell_id", "cell", "cellID") %in% names(boundaries)][1]
  x_col <- c("vertex_x", "x", "X")[c("vertex_x", "x", "X") %in% names(boundaries)][1]
  y_col <- c("vertex_y", "y", "Y")[c("vertex_y", "y", "Y") %in% names(boundaries)][1]
  if (is.na(cell_col) || is.na(x_col) || is.na(y_col)) {
    stop("Cell boundary table must contain cell_id/cell, x/vertex_x, and y/vertex_y columns")
  }
  list(cell = cell_col, x = x_col, y = y_col)
}

#' Read Cell Boundaries
#'
#' Reads standard Xenium-style cell boundary files. If `cells` or `bbox` are
#' supplied, filtering is done by cell after reading so that complete polygon
#' paths are retained for plotting and clipped only by the ggplot coordinate
#' system.
#'
#' @param path Boundary CSV, CSV.GZ, or parquet file.
#' @param cells Optional cell IDs to keep.
#' @param bbox Optional bounding box `c(xmin, xmax, ymin, ymax)`. Cells with at
#'   least one boundary vertex inside the box are retained with all vertices.
#'
#' @return Data frame with `cell_id`, `x`, and `y`.
#' @export
celladmix_read_cell_boundaries <- function(path, cells = NULL, bbox = NULL) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  if (grepl("\\.parquet$", path, ignore.case = TRUE)) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("The arrow package is required to read parquet cell boundaries")
    }
    boundaries <- as.data.frame(arrow::read_parquet(path))
  } else {
    boundaries <- utils::read.csv(path, stringsAsFactors = FALSE)
  }
  cols <- .celladmix_boundary_columns(boundaries)
  out <- data.frame(cell_id = as.character(boundaries[[cols$cell]]),
    x = as.numeric(boundaries[[cols$x]]),
    y = as.numeric(boundaries[[cols$y]]),
    stringsAsFactors = FALSE)
  keep <- is.finite(out$x) & is.finite(out$y) & !is.na(out$cell_id) & nzchar(out$cell_id)
  out <- out[keep, , drop = FALSE]
  if (!is.null(cells)) {
    out <- out[out$cell_id %in% as.character(cells), , drop = FALSE]
  }
  if (!is.null(bbox)) {
    if (!is.numeric(bbox) || length(bbox) != 4L) {
      stop("bbox must be numeric c(xmin, xmax, ymin, ymax)")
    }
    hit <- out$x >= bbox[[1]] & out$x <= bbox[[2]] &
      out$y >= bbox[[3]] & out$y <= bbox[[4]]
    out <- out[out$cell_id %in% unique(out$cell_id[hit]), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

.celladmix_cell_hull_contour <- function(molecules, center, cell_id, role, min_radius = 4) {
  points <- molecules[molecules$cell == cell_id & is.finite(molecules$x) &
    is.finite(molecules$y), c("x", "y"), drop = FALSE]
  if (nrow(points) >= 3L) {
    idx <- grDevices::chull(points$x, points$y)
    idx <- c(idx, idx[[1]])
    return(data.frame(cell_id = cell_id, x = points$x[idx], y = points$y[idx],
      role = role, stringsAsFactors = FALSE))
  }
  theta <- seq(0, 2 * pi, length.out = 72L)
  data.frame(cell_id = cell_id, x = center$x + min_radius * cos(theta),
    y = center$y + min_radius * sin(theta), role = role, stringsAsFactors = FALSE)
}

.celladmix_square_bbox <- function(xrange, yrange, padding = 5, min_side = 24) {
  xmid <- mean(xrange)
  ymid <- mean(yrange)
  side <- max(diff(range(xrange)), diff(range(yrange)), min_side) + 2 * padding
  half <- side / 2
  c(xmid - half, xmid + half, ymid - half, ymid + half)
}

.celladmix_target_contours <- function(molecules, centers, target_cell, bbox,
                                       boundaries = NULL) {
  if (is.character(boundaries) && length(boundaries) == 1L) {
    boundaries <- celladmix_read_cell_boundaries(boundaries, bbox = bbox)
  }
  if (!is.null(boundaries) && nrow(boundaries)) {
    hit <- boundaries$x >= bbox[[1]] & boundaries$x <= bbox[[2]] &
      boundaries$y >= bbox[[3]] & boundaries$y <= bbox[[4]]
    keep_cells <- unique(c(target_cell, boundaries$cell_id[hit]))
    boundary_df <- boundaries[boundaries$cell_id %in% keep_cells, , drop = FALSE]
    if (nrow(boundary_df)) {
      boundary_df$role <- ifelse(boundary_df$cell_id == target_cell,
        "Target cell contour", "Nearby cell contours")
      return(boundary_df[, c("cell_id", "x", "y", "role"), drop = FALSE])
    }
  }
  center <- centers[centers$cell == target_cell, , drop = FALSE]
  if (!nrow(center)) {
    return(data.frame(cell_id = character(), x = numeric(), y = numeric(),
      role = character(), stringsAsFactors = FALSE))
  }
  .celladmix_cell_hull_contour(molecules, center, target_cell, "Target cell contour")
}

.celladmix_normalize_channel <- function(x, probs = c(0.01, 0.998), gamma = 0.85) {
  vals <- as.numeric(x)
  lim <- stats::quantile(vals[is.finite(vals)], probs = probs, na.rm = TRUE)
  if (!is.finite(lim[[1]]) || !is.finite(lim[[2]]) || lim[[2]] <= lim[[1]]) {
    lim <- range(vals[is.finite(vals)], na.rm = TRUE)
  }
  if (!is.finite(lim[[1]]) || !is.finite(lim[[2]]) || lim[[2]] <= lim[[1]]) {
    return(matrix(0, nrow = nrow(x), ncol = ncol(x)))
  }
  out <- pmax(0, pmin(1, ((x - lim[[1]]) / (lim[[2]] - lim[[1]])))) ^ gamma
  dim(out) <- dim(x)
  out
}

.celladmix_compose_stain_raster <- function(crops, colors = NULL) {
  crops <- crops[!vapply(crops, is.null, logical(1))]
  if (!length(crops)) {
    return(NULL)
  }
  dims <- dim(crops[[1]]$values)
  if (!all(vapply(crops, function(x) identical(dim(x$values), dims), logical(1)))) {
    stop("All stain crops must have matching dimensions")
  }
  crop_names <- names(crops)
  if (is.null(crop_names) || any(!nzchar(crop_names))) {
    crop_names <- paste0("stain_", seq_along(crops))
  }
  default_colors <- list(
    dapi = c(0.84, 0.73, 0.96),
    nucleus = c(0.84, 0.73, 0.96),
    membrane = c(0.75, 0.93, 0.74),
    polyA = c(0.95, 0.82, 0.56),
    poly_a = c(0.95, 0.82, 0.56)
  )
  colors <- colors %||% list()
  rgb_arr <- array(1, dim = c(dims[[1]], dims[[2]], 3))
  for (i in seq_along(crops)) {
    channel <- .celladmix_normalize_channel(crops[[i]]$values)
    key <- crop_names[[i]]
    color <- colors[[key]] %||% default_colors[[key]] %||% c(0.35, 0.35, 0.35)
    for (j in 1:3) {
      rgb_arr[, , j] <- rgb_arr[, , j] * (1 - channel * (1 - color[[j]]))
    }
  }
  grDevices::as.raster(rgb_arr)
}

.celladmix_resolve_example_stains <- function(fit, stains) {
  if (is.character(stains) && length(stains) == 1L && identical(stains[[1]], "auto")) {
    return(fit$stains())
  }
  if (is.null(stains) || !length(stains)) {
    return(list())
  }
  if (!is.list(stains)) {
    stains <- list(stain = stains)
  }
  if (!is.null(stains$path)) {
    # A single stain descriptor rather than a named list of descriptors.
    stains <- stats::setNames(list(stains), stains$stain %||% "stain")
  }
  stains
}

.celladmix_example_columns <- c(
  "target_cell", "target_cell_type", "top_admix_factor", "top_admix_source",
  "top_admix_score", "top_factor_molecules", "top_factor_fraction",
  "top_evidence", "n_candidate_factors", "total_factor_molecules",
  "transcript_count", "dominant_factor")

# Assemble the final example table for explicitly requested cells: keep the
# score-derived rows that exist, append NA-evidence rows for the remaining
# requested cells, and preserve the requested order.
.celladmix_requested_cell_examples <- function(score, cells, examples, cell_data) {
  have <- if (!is.null(examples) && nrow(examples)) {
    as.character(examples$target_cell)
  } else {
    character()
  }
  missing_cells <- setdiff(cells, have)
  if (length(missing_cells)) {
    type_map <- tryCatch(
      score$fit$dataset$annotation(score$fit$annotation_name, as_vector = TRUE),
      error = function(e) NULL)
    rows <- data.frame(target_cell = missing_cells, stringsAsFactors = FALSE)
    rows$target_cell_type <- if (!is.null(type_map)) {
      unname(type_map[missing_cells])
    } else {
      NA_character_
    }
    for (column in setdiff(.celladmix_example_columns,
        c("target_cell", "target_cell_type"))) {
      rows[[column]] <- NA_real_
    }
    if (is.data.frame(cell_data) && "cell_id" %in% names(cell_data)) {
      idx <- match(missing_cells, as.character(cell_data$cell_id))
      if ("transcript_count" %in% names(cell_data)) {
        rows$transcript_count <- cell_data$transcript_count[idx]
      }
      if ("dominant_factor" %in% names(cell_data)) {
        rows$dominant_factor <- cell_data$dominant_factor[idx]
      }
    }
    if (!is.null(examples) && nrow(examples)) {
      rows <- rows[, names(examples), drop = FALSE]
      examples <- rbind(examples, rows)
    } else {
      examples <- rows
    }
  }
  examples <- examples[match(cells, as.character(examples$target_cell)), , drop = FALSE]
  examples <- examples[!is.na(examples$target_cell), , drop = FALSE]
  rownames(examples) <- NULL
  examples
}

.celladmix_resolve_example_cell_types <- function(fit, cell_types) {
  if (is.character(cell_types) && length(cell_types) == 1L &&
      identical(cell_types[[1]], "auto")) {
    return(tryCatch(
      fit$dataset$annotation(fit$annotation_name, as_vector = TRUE),
      error = function(e) NULL))
  }
  if (is.null(cell_types) || !length(cell_types)) {
    return(NULL)
  }
  if (is.data.frame(cell_types)) {
    type_col <- intersect(c("cell_type", "merged_annotation", "cluster_label"),
      names(cell_types))
    if (!("cell_id" %in% names(cell_types)) || !length(type_col)) {
      stop("cell_types data frame must contain cell_id and cell_type columns")
    }
    return(stats::setNames(as.character(cell_types[[type_col[[1]]]]),
      as.character(cell_types$cell_id)))
  }
  cell_types
}

.celladmix_resolve_example_boundaries <- function(fit, boundaries) {
  if (is.character(boundaries) && length(boundaries) == 1L &&
      identical(boundaries[[1]], "auto")) {
    path <- tryCatch(celladmix_discover_cell_boundaries(fit$run),
      error = function(e) NA_character_)
    if (is.na(path)) {
      return(NULL)
    }
    return(path)
  }
  boundaries
}

.celladmix_collect_example_stain_crops <- function(fit, stains, bbox, max_pixels) {
  stains <- .celladmix_resolve_example_stains(fit, stains)
  if (!length(stains)) {
    return(list())
  }
  lapply(stains, function(image) {
    if (is.null(image)) {
      return(NULL)
    }
    fit$stain_crop(image, bbox = bbox, max_pixels = max_pixels)
  })
}

#' Select Example Cells from a Score
#'
#' Selects target cells with strong score-supported non-native factor evidence.
#' Native factors are excluded using score source calls and the cell's dominant
#' factor; the displayed top admixture factor is then chosen by molecule
#' abundance inside the target cell, with score used as a tie-breaker.
#' The returned table can be passed to [celladmix_plot_cell_example()] or
#' `score$plot_example()`.
#'
#' @param score A [CellAdmixScore] object.
#' @param rules Optional correction-rule table. If `NULL`, rules are generated
#'   from `score` unless `use_rules = FALSE`.
#' @param score_annotation Optional score annotation object. If `NULL`, it is
#'   generated from `score`.
#' @param cell_data Optional cell-level data frame with `cell_id`,
#'   `transcript_count`, and `dominant_factor` columns.
#' @param targets Optional target cell types to keep.
#' @param n_per_target Number of examples per target cell type.
#' @param p_thresh Score p-value threshold used when generating rules.
#' @param adjust_p Whether generated rules use adjusted p-values.
#' @param use_rules If `TRUE`, only factor/target combinations present in rules
#'   are considered.
#' @param min_molecules Minimum target-cell molecule count.
#' @param min_factor_molecules Minimum score-table `factor_count` for a
#'   target-cell/factor candidate.
#' @param cells Optional explicit target cell IDs. When supplied, selection is
#'   restricted to these cells, the evidence filters are relaxed so every
#'   requested cell is returned (with `NA` factor columns when the score table
#'   has no evidence for it), and the cells come back in the requested order.
#'
#' @return Data frame of selected examples.
#' @export
celladmix_select_example_cells <- function(
    score,
    rules = NULL,
    score_annotation = NULL,
    cell_data = NULL,
    targets = NULL,
    n_per_target = 4,
    p_thresh = 0.1,
    adjust_p = FALSE,
    use_rules = TRUE,
    min_molecules = 50,
    min_factor_molecules = 3,
    cells = NULL
  ) {
  if (!inherits(score, "CellAdmixScore")) {
    stop("score must be a CellAdmixScore")
  }
  if (!is.null(cells)) {
    cells <- as.character(cells)
    use_rules <- FALSE
    min_molecules <- 0
    min_factor_molecules <- 0
    targets <- NULL
  }
  pairs <- score$pairs()
  if (is.null(pairs) || !nrow(pairs)) {
    pairs <- data.frame()
  }
  if (!is.null(cells) && nrow(pairs)) {
    pairs <- pairs[as.character(pairs$target_cell) %in% cells, , drop = FALSE]
  }
  if (is.null(cell_data)) {
    cell_data <- score$fit$cell_factors()
  }
  if (!nrow(pairs)) {
    if (is.null(cells)) {
      return(data.frame())
    }
    return(.celladmix_requested_cell_examples(score, cells, NULL, cell_data))
  }
  score_annotation <- score_annotation %||% score$annotation(p_thresh = p_thresh,
    adjust_p = adjust_p)
  if (isTRUE(use_rules) && is.null(rules)) {
    rules <- score$rules(p_thresh = p_thresh, adjust_p = adjust_p, targets = targets)
  }

  required <- c("target_cell", "target_cell_type", "factor", "mean_score")
  missing <- setdiff(required, names(pairs))
  if (length(missing)) {
    stop("score pair table is missing required columns: ", paste(missing, collapse = ", "))
  }
  if (!is.null(targets)) {
    pairs <- pairs[pairs$target_cell_type %in% targets, , drop = FALSE]
  }
  if ("used_in_summary" %in% names(pairs) && is.null(cells)) {
    pairs <- pairs[pairs$used_in_summary %in% TRUE, , drop = FALSE]
  }
  pairs <- pairs[is.finite(pairs$mean_score), , drop = FALSE]
  if (isTRUE(use_rules) && !is.null(rules) && nrow(rules)) {
    rule_key <- paste(rules$target_cell_type, rules$factor, sep = "|")
    pairs <- pairs[paste(pairs$target_cell_type, pairs$factor, sep = "|") %in% rule_key,
      , drop = FALSE]
  }
  if ("factor_count" %in% names(pairs)) {
    pairs <- pairs[pairs$factor_count >= min_factor_molecules, , drop = FALSE]
  }
  if (!nrow(pairs)) {
    if (is.null(cells)) {
      return(data.frame())
    }
    return(.celladmix_requested_cell_examples(score, cells, NULL, cell_data))
  }

  factor_sources <- .celladmix_factor_sources(score_annotation)
  split_key <- paste(pairs$target_cell, pairs$target_cell_type, pairs$factor, sep = "|")
  per_factor <- do.call(rbind, lapply(split(pairs, split_key), function(x) {
    best <- x[which.max(x$mean_score), , drop = FALSE]
    data.frame(target_cell = best$target_cell[[1]],
      target_cell_type = best$target_cell_type[[1]],
      factor = best$factor[[1]],
      max_score = max(x$mean_score, na.rm = TRUE),
      factor_count = if ("factor_count" %in% names(x)) max(x$factor_count, na.rm = TRUE) else NA_real_,
      scored_molecules = if ("scored_molecules" %in% names(x)) max(x$scored_molecules, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE)
  }))
  per_factor <- merge(per_factor, factor_sources, by = "factor", all.x = TRUE,
    sort = FALSE)

  if (is.data.frame(cell_data) && "cell_id" %in% names(cell_data)) {
    optional <- intersect(c("cell_id", "transcript_count", "dominant_factor"), names(cell_data))
    per_factor <- merge(per_factor, cell_data[, optional, drop = FALSE],
      by.x = "target_cell", by.y = "cell_id", all.x = TRUE, sort = FALSE)
  }
  per_factor <- per_factor[is.na(per_factor$source_cell_type) |
    per_factor$source_cell_type != per_factor$target_cell_type, , drop = FALSE]
  if ("dominant_factor" %in% names(per_factor)) {
    per_factor <- per_factor[is.na(per_factor$dominant_factor) |
      per_factor$factor != per_factor$dominant_factor, , drop = FALSE]
  }
  if ("transcript_count" %in% names(per_factor)) {
    per_factor <- per_factor[is.na(per_factor$transcript_count) |
      per_factor$transcript_count >= min_molecules, , drop = FALSE]
  }
  if (is.null(cells)) {
    per_factor <- per_factor[is.finite(per_factor$max_score) & per_factor$max_score > 0,
      , drop = FALSE]
  }
  if (!nrow(per_factor)) {
    if (is.null(cells)) {
      return(data.frame())
    }
    return(.celladmix_requested_cell_examples(score, cells, NULL, cell_data))
  }

  examples <- do.call(rbind, lapply(split(per_factor, per_factor$target_cell), function(x) {
    if ("factor_count" %in% names(x) && any(is.finite(x$factor_count))) {
      x <- x[order(x$factor_count, x$max_score, decreasing = TRUE), , drop = FALSE]
    } else {
      x <- x[order(x$max_score, decreasing = TRUE), , drop = FALSE]
    }
    top <- x[1, , drop = FALSE]
    top_molecules <- if ("factor_count" %in% names(top)) top$factor_count[[1]] else NA_real_
    total_molecules <- sum(x$factor_count, na.rm = TRUE)
    transcript_count <- if ("transcript_count" %in% names(x)) x$transcript_count[[1]] else NA_real_
    top_evidence <- if (is.finite(top_molecules) && is.finite(top$max_score[[1]])) {
      top_molecules * max(0, top$max_score[[1]])
    } else {
      NA_real_
    }
    data.frame(target_cell = x$target_cell[[1]],
      target_cell_type = x$target_cell_type[[1]],
      top_admix_factor = top$factor[[1]],
      top_admix_source = top$source_cell_type[[1]],
      top_admix_score = top$max_score[[1]],
      top_factor_molecules = top_molecules,
      top_factor_fraction = if (is.finite(transcript_count) && transcript_count > 0) {
        top_molecules / transcript_count
      } else {
        NA_real_
      },
      top_evidence = top_evidence,
      n_candidate_factors = nrow(x),
      total_factor_molecules = total_molecules,
      transcript_count = transcript_count,
      dominant_factor = if ("dominant_factor" %in% names(x)) x$dominant_factor[[1]] else NA_integer_,
      stringsAsFactors = FALSE)
  }))
  if (!is.null(cells)) {
    return(.celladmix_requested_cell_examples(score, cells, examples, cell_data))
  }
  examples <- examples[order(examples$top_evidence, examples$top_factor_molecules,
    examples$total_factor_molecules, examples$top_admix_score,
    examples$n_candidate_factors, decreasing = TRUE), , drop = FALSE]
  target_order <- targets %||% unique(examples$target_cell_type)
  examples <- do.call(rbind, lapply(target_order, function(cell_type) {
    head(examples[examples$target_cell_type == cell_type, , drop = FALSE],
      max(1L, as.integer(n_per_target)))
  }))
  rownames(examples) <- NULL
  examples
}

#' Prepare an Example-Cell Plot
#'
#' @param fit A [CellAdmixFit] object.
#' @param example One row from [celladmix_select_example_cells()].
#' @param score_annotation Optional score annotation used to identify native
#'   factors from source calls.
#' @param cell_data Optional cell-level data frame with coordinates and factor
#'   summaries. Defaults to `fit$cell_factors()`.
#' @param boundaries Cell boundaries used for contours and cell-type coloring.
#'   The default `"auto"` discovers the standard Xenium boundary file next to
#'   the source bundle (`NULL` when unavailable). Pass a boundary data frame or
#'   path to override, or `NULL` to fall back to molecule hulls.
#' @param stains Stain images to compose as the plot background. The default
#'   `"auto"` discovers the available Xenium stain images (DAPI and membrane)
#'   from the fit's bundle and is empty for non-Xenium sources. Use `NULL` to
#'   disable backgrounds, or pass a named list of stain-image descriptors to
#'   override.
#' @param cell_types Cell-type labels used to shade cell polygons. The default
#'   `"auto"` uses the dataset's active annotation; pass a named vector
#'   (`cell_id` names, type values), a data frame with `cell_id`/`cell_type`,
#'   or `NULL` to disable shading.
#' @param markers Optional marker genes to size-emphasize inside the target cell.
#' @param padding,min_side Controls the square spatial window around the target
#'   cell.
#' @param max_pixels Maximum stain-crop width/height.
#'
#' @return A list consumed by [celladmix_plot_cell_example()].
#' @export
celladmix_prepare_cell_example <- function(
    fit,
    example,
    score_annotation = NULL,
    cell_data = NULL,
    boundaries = "auto",
    stains = "auto",
    cell_types = "auto",
    markers = NULL,
    padding = 5,
    min_side = 24,
    max_pixels = 512
  ) {
  if (!inherits(fit, "CellAdmixFit")) {
    stop("fit must be a CellAdmixFit")
  }
  if (!is.data.frame(example) || nrow(example) < 1L) {
    stop("example must be a non-empty data frame")
  }
  example <- example[1, , drop = FALSE]
  if (is.null(cell_data)) {
    cell_data <- fit$cell_factors()
  }
  if (!all(c("cell_id", "x", "y") %in% names(cell_data))) {
    stop("cell_data must contain cell_id, x, and y columns")
  }
  target_cell <- as.character(example$target_cell[[1]])
  target_type <- as.character(example$target_cell_type[[1]])
  target_center <- cell_data[cell_data$cell_id == target_cell, , drop = FALSE]
  if (!nrow(target_center)) {
    stop("Target cell is not present in cell_data: ", target_cell)
  }

  boundaries <- .celladmix_resolve_example_boundaries(fit, boundaries)
  if (is.character(boundaries) && length(boundaries) == 1L) {
    boundaries <- celladmix_read_cell_boundaries(boundaries)
  }
  boundary <- if (!is.null(boundaries) && nrow(boundaries)) {
    boundaries[boundaries$cell_id == target_cell, , drop = FALSE]
  } else {
    NULL
  }
  if (!is.null(boundary) && nrow(boundary)) {
    bbox <- .celladmix_square_bbox(range(boundary$x, na.rm = TRUE),
      range(boundary$y, na.rm = TRUE), padding = padding, min_side = min_side)
  } else {
    bbox <- .celladmix_square_bbox(target_center$x, target_center$y,
      padding = padding, min_side = min_side)
  }

  stain_crops <- .celladmix_collect_example_stain_crops(fit, stains, bbox, max_pixels)
  background <- .celladmix_compose_stain_raster(stain_crops)
  molecules <- fit$region(bbox = bbox)
  molecules$factor_label <- as.integer(molecules$factor_label)
  molecules$inside_target <- !is.na(molecules$cell) & molecules$cell == target_cell

  factor_sources <- if (!is.null(score_annotation)) .celladmix_factor_sources(score_annotation) else data.frame()
  native_factors <- integer()
  if (nrow(factor_sources) && !is.na(target_type)) {
    native_factors <- factor_sources$factor[factor_sources$source_cell_type == target_type]
  }
  if ("dominant_factor" %in% names(target_center) &&
      is.finite(target_center$dominant_factor[[1]])) {
    native_factors <- unique(c(native_factors, as.integer(target_center$dominant_factor[[1]])))
  }
  top_factor <- .celladmix_as_factor_id(example$top_admix_factor[[1]])
  top_factor <- if (length(top_factor)) top_factor[[1]] else NA_integer_
  native_factors <- setdiff(unique(native_factors), top_factor)

  native_label <- if (length(native_factors)) {
    paste0("native (", paste(.celladmix_factor_label(native_factors), collapse = ","), ")")
  } else {
    "native factors"
  }
  top_label <- if (is.finite(top_factor)) {
    paste0("top admixture (", .celladmix_factor_label(top_factor), ")")
  } else {
    "top admixture"
  }
  other_label <- "other admixture factors"
  role_levels <- c(native_label, top_label, other_label)
  molecules$role <- other_label
  molecules$role[molecules$factor_label %in% native_factors] <- native_label
  if (is.finite(top_factor)) {
    molecules$role[molecules$factor_label == top_factor] <- top_label
  }
  molecules$role <- factor(molecules$role, levels = role_levels)

  molecules$is_marker <- FALSE
  if (!is.null(markers) && length(markers)) {
    molecules$is_marker <- molecules$inside_target & molecules$gene %in% markers
  }
  centers <- data.frame(cell = target_center$cell_id, x = target_center$x,
    y = target_center$y, stringsAsFactors = FALSE)
  contours <- .celladmix_target_contours(molecules, centers, target_cell, bbox,
    boundaries = boundaries)

  type_map <- .celladmix_resolve_example_cell_types(fit, cell_types)
  contours$cell_type <- if (!is.null(type_map) && nrow(contours)) {
    unname(type_map[as.character(contours$cell_id)])
  } else {
    NA_character_
  }

  structure(list(example = example, target_cell = target_cell,
    target_cell_type = target_type, bbox = bbox, background = background,
    molecules = molecules,
    contours = contours, role_levels = role_levels),
    class = "celladmix_cell_example")
}

#' Plot a Molecule-Level Example Cell
#'
#' @param example Prepared example from [celladmix_prepare_cell_example()], or a
#'   one-row example table when `fit` is supplied.
#' @param fit Optional [CellAdmixFit] used to prepare an unprepared example.
#' @param color_cell_types Whether to color nearby cell contours by cell type
#'   when boundaries and cell-type labels are available. Only the contour is
#'   colored so that stain backgrounds stay visible.
#' @param cell_type_linewidth Line width of the cell-type contours.
#' @param score_annotation,cell_data,boundaries,stains,cell_types,markers,padding,min_side,max_pixels
#'   Passed to [celladmix_prepare_cell_example()] when `fit` is supplied.
#' @param outside_size,inside_size Molecule point sizes outside and inside the
#'   target cell.
#' @param marker_size_multiplier Marker-gene point-size multiplier inside the
#'   target cell.
#' @param non_marker_size_multiplier Non-marker point-size multiplier.
#' @param contour_color Target-cell contour color.
#'
#' @return A `ggplot2` object.
#' @export
celladmix_plot_cell_example <- function(
    example,
    fit = NULL,
    score_annotation = NULL,
    cell_data = NULL,
    boundaries = "auto",
    stains = "auto",
    cell_types = "auto",
    markers = NULL,
    padding = 5,
    min_side = 24,
    max_pixels = 512,
    outside_size = 0.54,
    inside_size = 1.18,
    marker_size_multiplier = 1.1,
    non_marker_size_multiplier = 0.9,
    contour_color = "#e85d04",
    color_cell_types = TRUE,
    cell_type_linewidth = 0.55,
    title = NULL,
    subtitle = NULL
  ) {
  .celladmix_require_ggplot2()
  if (!inherits(example, "celladmix_cell_example")) {
    if (is.null(fit)) {
      stop("fit must be supplied when example is not precomputed")
    }
    example <- celladmix_prepare_cell_example(fit, example,
      score_annotation = score_annotation, cell_data = cell_data,
      boundaries = boundaries, stains = stains, cell_types = cell_types,
      markers = markers,
      padding = padding, min_side = min_side, max_pixels = max_pixels)
  }
  bbox <- example$bbox
  molecules <- example$molecules
  if (!("is_marker" %in% names(molecules))) {
    molecules$is_marker <- FALSE
  }
  contours <- example$contours
  target_contours <- contours[contours$role == "Target cell contour", , drop = FALSE]
  nearby_contours <- contours[contours$role == "Nearby cell contours", , drop = FALSE]
  target_cell <- example$target_cell
  typed <- if (isTRUE(color_cell_types) && "cell_type" %in% names(contours)) {
    contours[!is.na(contours$cell_type) & contours$cell_id != target_cell, ,
      drop = FALSE]
  } else {
    contours[0, , drop = FALSE]
  }
  plain_nearby <- if (nrow(typed)) {
    nearby_contours[!(nearby_contours$cell_id %in% typed$cell_id), , drop = FALSE]
  } else {
    nearby_contours
  }
  role_cols <- stats::setNames(c("#2b6cb0", "#c92a2a", "#f08c00"),
    example$role_levels)
  title <- title %||% sprintf("%s: %s", example$target_cell_type, example$target_cell)
  subtitle <- subtitle %||% if (!is.null(example$background)) {
    "stain background; target contour in orange"
  } else {
    "target contour in orange"
  }

  p <- ggplot2::ggplot()
  if (!is.null(example$background)) {
    p <- p + ggplot2::annotation_raster(example$background,
      xmin = bbox[[1]], xmax = bbox[[2]], ymin = bbox[[3]], ymax = bbox[[4]])
  }
  if (nrow(typed)) {
    # Only contours are colored by cell type; filled polygons would occlude
    # the stain background.
    p <- p + ggplot2::geom_polygon(data = typed,
      ggplot2::aes(x, y, group = cell_id, fill = cell_type,
        colour = ggplot2::after_scale(fill)),
      alpha = 0, linewidth = cell_type_linewidth) +
      ggplot2::guides(fill = ggplot2::guide_legend(
        override.aes = list(alpha = 1, linewidth = 0))) +
      ggplot2::labs(fill = "Cell type")
  }
  p +
    ggplot2::geom_path(data = plain_nearby,
      ggplot2::aes(x, y, group = cell_id), color = "white",
      linewidth = 0.35, alpha = 0.35) +
    ggplot2::geom_path(data = target_contours,
      ggplot2::aes(x, y, group = cell_id), color = "white",
      linewidth = 1.2, alpha = 0.95) +
    ggplot2::geom_path(data = target_contours,
      ggplot2::aes(x, y, group = cell_id), color = contour_color,
      linewidth = 0.75, alpha = 0.98) +
    ggplot2::geom_point(data = molecules[!molecules$inside_target, , drop = FALSE],
      ggplot2::aes(x, y, color = role),
      size = outside_size * non_marker_size_multiplier, alpha = 0.38) +
    ggplot2::geom_point(data = molecules[molecules$inside_target & !molecules$is_marker, , drop = FALSE],
      ggplot2::aes(x, y, color = role),
      size = inside_size * non_marker_size_multiplier, alpha = 0.92) +
    ggplot2::geom_point(data = molecules[molecules$inside_target & molecules$is_marker, , drop = FALSE],
      ggplot2::aes(x, y, color = role),
      size = inside_size * marker_size_multiplier, alpha = 0.92) +
    ggplot2::scale_color_manual(values = role_cols, drop = FALSE) +
    ggplot2::scale_y_reverse(limits = rev(bbox[3:4]), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(limits = bbox[1:2], expand = c(0, 0)) +
    ggplot2::coord_equal() +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = 3.4, alpha = 1))) +
    ggplot2::labs(title = title, subtitle = subtitle, color = "Molecule class") +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.spacing.y = ggplot2::unit(1, "pt"),
      legend.text = ggplot2::element_text(size = 7),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8))
}
