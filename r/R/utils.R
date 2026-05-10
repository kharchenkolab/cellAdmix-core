`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

.celladmix_make_run_object <- function(x, class = "celladmix_run") {
  x$path <- normalizePath(x$path, winslash = "/", mustWork = FALSE)
  if (is.null(x$source)) {
    x$source <- list(type = NA_character_)
  }
  structure(x, class = c(class, "celladmix_run"))
}

.celladmix_default_output_dir <- function(prefix = "celladmix_output_") {
  path <- tempfile(pattern = prefix)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

.celladmix_default_threads <- function(max_threads = 10L) {
  env <- Sys.getenv("CELLADMIX_NUM_THREADS", unset = NA_character_)
  if (!is.na(env) && nzchar(env)) {
    parsed <- suppressWarnings(as.integer(env[[1]]))
    if (!is.na(parsed) && parsed > 0L) {
      return(parsed)
    }
  }
  cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (is.na(cores) || cores < 1L) {
    cores <- 1L
  }
  max(1L, min(as.integer(max_threads), as.integer(cores)))
}

.celladmix_dots_with_default_threads <- function(dots, num_threads) {
  if (!("num_threads" %in% names(dots)) || is.null(dots$num_threads) ||
      !length(dots$num_threads) || is.na(dots$num_threads[[1]])) {
    dots$num_threads <- as.integer(num_threads)
  }
  dots
}

.celladmix_clustering_frame <- function(result, annotation = NULL) {
  clusters <- result$clusters
  embedding <- result$embedding
  if (is.data.frame(embedding) && nrow(embedding)) {
    df <- merge(clusters, embedding[, c("cell_id", "umap_1", "umap_2"), drop = FALSE],
      by = "cell_id", all.x = TRUE, sort = FALSE)
  } else {
    df <- clusters
  }
  if (!is.null(annotation)) {
    labels <- annotation
    if (is.data.frame(labels)) {
      id_col <- if ("cell_id" %in% names(labels)) "cell_id" else names(labels)[[1]]
      label_col <- setdiff(names(labels), id_col)[[1]]
      labels <- stats::setNames(as.character(labels[[label_col]]), as.character(labels[[id_col]]))
    }
    if (is.atomic(labels) && !is.null(names(labels))) {
      df$cell_type <- unname(labels[as.character(df$cell_id)])
    }
  }
  df
}

.celladmix_json_safe <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.data.frame(x)) {
    out <- x
    for (nm in names(out)) {
      out[[nm]] <- .celladmix_json_safe(out[[nm]])
    }
    return(out)
  }
  if (inherits(x, "table")) {
    values <- as.list(as.integer(x))
    names(values) <- names(x)
    return(values)
  }
  if (is.factor(x)) {
    return(as.character(x))
  }
  if (is.numeric(x)) {
    x[!is.finite(x)] <- NA_real_
    return(x)
  }
  if (is.atomic(x)) {
    return(x)
  }
  if (is.list(x)) {
    return(lapply(x, .celladmix_json_safe))
  }
  x
}

.celladmix_write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    .celladmix_json_safe(x),
    path,
    auto_unbox = TRUE,
    dataframe = "rows",
    null = "null",
    na = "null",
    pretty = TRUE
  )
  invisible(normalizePath(path, winslash = "/", mustWork = FALSE))
}

.celladmix_read_json <- function(path) {
  jsonlite::read_json(path, simplifyVector = TRUE)
}

.celladmix_make_prep_object <- function(x) {
  x$project_dir <- normalizePath(x$project_dir, winslash = "/", mustWork = FALSE)
  x$paths$project_json <- normalizePath(x$paths$project_json, winslash = "/", mustWork = FALSE)
  structure(x, class = "celladmix_prep")
}

.celladmix_normalize_crops <- function(crops) {
  if (is.null(crops)) {
    return(NULL)
  }

  if (is.character(crops) && length(crops) == 1L && file.exists(crops)) {
    crops <- utils::read.csv(crops, stringsAsFactors = FALSE)
  } else if (is.numeric(crops)) {
    if (!(length(crops) %in% c(4, 6))) {
      stop("Numeric crop specification must have length 4 or 6")
    }
    crops <- as.list(crops)
    names(crops) <- c("xmin", "xmax", "ymin", "ymax", if (length(crops) == 6) c("zmin", "zmax"))
    crops <- as.data.frame(crops, stringsAsFactors = FALSE)
  } else if (is.list(crops) && !is.data.frame(crops)) {
    if (all(c("xmin", "xmax", "ymin", "ymax") %in% names(crops))) {
      crops <- as.data.frame(crops, stringsAsFactors = FALSE)
    } else {
      crop_rows <- lapply(seq_along(crops), function(i) {
        crop <- crops[[i]]
        if (!is.list(crop)) {
          stop("Each crop in a crop list must be a named list")
        }
        as.data.frame(crop, stringsAsFactors = FALSE)
      })
      crops <- do.call(rbind, crop_rows)
    }
  }

  if (!is.data.frame(crops)) {
    stop("analysis_crops must be a numeric vector, a named list, or a data.frame")
  }

  rename_map <- c(
    x_min = "xmin",
    x_max = "xmax",
    y_min = "ymin",
    y_max = "ymax",
    z_min = "zmin",
    z_max = "zmax",
    region_id = "crop_id"
  )
  for (nm in names(rename_map)) {
    if (nm %in% colnames(crops) && !(rename_map[[nm]] %in% colnames(crops))) {
      colnames(crops)[colnames(crops) == nm] <- rename_map[[nm]]
    }
  }
  if (!all(c("xmin", "xmax", "ymin", "ymax") %in% colnames(crops)) &&
      all(c("center_x", "center_y", "half_size") %in% colnames(crops))) {
    crops$xmin <- as.numeric(crops$center_x) - as.numeric(crops$half_size)
    crops$xmax <- as.numeric(crops$center_x) + as.numeric(crops$half_size)
    crops$ymin <- as.numeric(crops$center_y) - as.numeric(crops$half_size)
    crops$ymax <- as.numeric(crops$center_y) + as.numeric(crops$half_size)
  }

  required <- c("xmin", "xmax", "ymin", "ymax")
  missing <- setdiff(required, colnames(crops))
  if (length(missing) > 0) {
    stop("analysis_crops is missing required columns: ", paste(missing, collapse = ", "))
  }
  if (!("crop_id" %in% colnames(crops))) {
    crops$crop_id <- sprintf("crop_%d", seq_len(nrow(crops)))
  }
  if (!("zmin" %in% colnames(crops))) {
    crops$zmin <- -Inf
  }
  if (!("zmax" %in% colnames(crops))) {
    crops$zmax <- Inf
  }
  crops <- crops[, c("crop_id", "xmin", "xmax", "ymin", "ymax", "zmin", "zmax"), drop = FALSE]
  crops$crop_id <- as.character(crops$crop_id)
  for (nm in c("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")) {
    crops[[nm]] <- as.numeric(crops[[nm]])
  }
  crops
}

.celladmix_cluster_paths <- function(prep, cluster_id = "default") {
  base_dir <- file.path(prep$paths$annotations_dir, cluster_id)
  list(
    cluster_id = cluster_id,
    dir = base_dir,
    clusters_path = file.path(base_dir, "cell_clusters.parquet"),
    embedding_path = file.path(base_dir, "cell_embedding.parquet"),
    summary_path = file.path(base_dir, "summary.json")
  )
}

.celladmix_make_run_dir <- function(prep, run_id = NULL) {
  if (is.null(run_id)) {
    run_id <- sprintf("run_%s_%06d", format(Sys.time(), "%Y%m%d_%H%M%S"), sample.int(999999L, 1L))
  }
  file.path(prep$paths$runs_dir, run_id)
}
