#' Factor Source Gene-Content Scores
#'
#' These helpers annotate NMF factors by comparing factor gene loadings against
#' cell-type-specific gene signatures estimated from the same dataset. They are
#' diagnostic/source-prior utilities and do not modify spatial scoring or
#' correction rules by themselves.

#' High-Level cellAdmix Factor Source Score
#' @export
CellAdmixFactorSourceScore <- R6::R6Class(
  "CellAdmixFactorSourceScore",
  public = list(
    name = NULL,
    method = "factor_sources",
    fit = NULL,
    result = NULL,
    annotation_name = NULL,
    annotation_hash = NULL,
    paths = NULL,
    params = NULL,

    initialize = function(name, fit, result, annotation, params = list()) {
      self$name <- .celladmix_clean_name(name, "factor_sources")
      self$fit <- fit
      self$result <- result
      self$annotation_name <- annotation$name
      self$annotation_hash <- annotation$hash
      self$paths <- result$paths %||% list()
      self$params <- params
      private$.annotation <- annotation
    },

    print = function(...) {
      cat("cellAdmix factor source score\n")
      cat("  name:", self$name, "\n")
      cat("  run:", self$fit$run_id, "\n")
      cat("  annotation:", self$annotation_name, "\n")
      cat("  factors:", length(unique(self$result$scores$factor)), "\n")
      cat("  cell types:", length(unique(self$result$scores$cell_type)), "\n")
      invisible(self)
    },

    summary = function() self$result$scores,
    scores = function() self$result$scores,
    marker_weights = function() self$result$marker_weights,
    top_markers = function(n_genes = 10L) {
      top <- self$result$top_markers
      if (!nrow(top)) {
        return(top)
      }
      top[top$rank <= as.integer(n_genes), , drop = FALSE]
    },
    top_factor_genes = function(n_genes = 8L, sources = NULL) {
      celladmix_factor_source_top_genes(self$result, n_genes = n_genes,
        sources = sources)
    },

    annotation = function(min_score = 0.15, min_margin = 0.03) {
      celladmix_factor_source_annotation(self$result, min_score = min_score,
        min_margin = min_margin)
    },

    plot_heatmap = function(...) {
      celladmix_plot_factor_source_heatmap(self$result, ...)
    },
    plot_top_genes = function(...) {
      celladmix_plot_factor_source_top_genes(self$result, ...)
    },
    plot = function(...) self$plot_heatmap(...)
  ),
  private = list(
    .annotation = NULL
  )
)

#' Score NMF Factor Sources by Cell-Type Gene Content
#'
#' Builds pseudobulk cell-type profiles, estimates cell-type-specific gene
#' weights as one-versus-rest log-expression margins, and scores each NMF
#' factor against each cell type with a marker-weighted cosine similarity.
#'
#' @param run A `celladmix_run`.
#' @param annotation Named cell-type vector or `celladmix_annotation` object.
#' @param counts Optional sparse gene-by-cell count matrix. When omitted,
#'   counts are collected from `run`.
#' @param out_dir Optional directory for CSV/JSON sidecar output.
#' @param normalize Pseudobulk normalization used before computing cell-type
#'   marker weights.
#' @param scale_factor Library-size scale factor for CPM normalization.
#' @param specificity_power Exponent applied to positive one-vs-rest marker
#'   weights. Values above 1 emphasize highly specific genes.
#' @param factor_power Exponent applied to non-negative factor loadings before
#'   cosine scoring.
#' @param min_marker_logfc Minimum positive one-vs-rest log-expression margin
#'   for a gene to contribute to a cell-type marker vector.
#' @param top_markers_per_type Optional cap on marker genes retained per cell
#'   type.
#' @param top_genes Number of top marker genes stored per cell type.
#'
#' @return A list containing score, marker-weight, and top-gene tables.
#' @keywords internal
#' @noRd
celladmix_score_factor_sources <- function(
    run,
    annotation,
    counts = NULL,
    out_dir = NULL,
    normalize = c("log1p_cpm", "cpm", "none"),
    scale_factor = 10000,
    specificity_power = 1,
    factor_power = 1,
    min_marker_logfc = 0,
    top_markers_per_type = NULL,
    top_genes = 25L
  ) {
  if (!inherits(run, "celladmix_run")) {
    stop("celladmix_score_factor_sources() expects a celladmix_run")
  }
  if (is.null(run$h) || !length(run$h)) {
    stop("Run does not contain factor loadings")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required for factor source scoring")
  }
  normalize <- match.arg(normalize)
  specificity_power <- as.numeric(specificity_power)
  factor_power <- as.numeric(factor_power)
  if (!is.finite(specificity_power) || specificity_power <= 0) {
    stop("specificity_power must be positive")
  }
  if (!is.finite(factor_power) || factor_power <= 0) {
    stop("factor_power must be positive")
  }

  ann <- if (inherits(annotation, "celladmix_annotation")) {
    annotation$labels
  } else {
    annotation
  }
  if (is.null(ann) || is.null(names(ann))) {
    stop("annotation must be a named cell-type vector or celladmix_annotation")
  }
  ann <- stats::setNames(as.character(ann), names(ann))
  ann <- ann[!is.na(ann) & nzchar(ann)]
  if (!length(ann)) {
    stop("annotation does not contain usable labels")
  }

  counts <- counts %||% celladmix_collect_counts_sparse(run)
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("counts must have gene row names and cell column names")
  }

  h <- run$h
  genes <- run$genes %||% colnames(h) %||% paste0("gene_", seq_len(ncol(h)))
  colnames(h) <- as.character(genes)
  rownames(h) <- paste0("F", seq_len(nrow(h)))

  common_genes <- intersect(colnames(h), rownames(counts))
  if (!length(common_genes)) {
    stop("No overlap between factor-loading genes and count-matrix genes")
  }
  common_cells <- intersect(colnames(counts), names(ann))
  if (!length(common_cells)) {
    stop("No overlap between count-matrix cells and annotation cells")
  }
  h <- h[, common_genes, drop = FALSE]
  counts <- counts[common_genes, common_cells, drop = FALSE]
  labels <- ann[common_cells]
  keep_cells <- !is.na(labels) & nzchar(labels)
  counts <- counts[, keep_cells, drop = FALSE]
  labels <- labels[keep_cells]
  cell_types <- sort(unique(unname(labels)))
  if (length(cell_types) < 2L) {
    stop("Factor source scoring requires at least two annotated cell types")
  }

  type_index <- match(labels, cell_types)
  design <- Matrix::sparseMatrix(
    i = seq_along(type_index),
    j = type_index,
    x = 1,
    dims = c(length(type_index), length(cell_types)),
    dimnames = list(names(labels), cell_types)
  )
  bulk <- counts %*% design
  lib <- Matrix::colSums(bulk)
  if (!identical(normalize, "none")) {
    scale <- scale_factor / pmax(as.numeric(lib), 1)
    bulk <- bulk %*% Matrix::Diagonal(x = scale)
  }
  if (identical(normalize, "log1p_cpm")) {
    bulk@x <- log1p(bulk@x)
  }
  profile <- as.matrix(bulk)
  colnames(profile) <- cell_types
  rownames(profile) <- common_genes

  weights <- matrix(0, nrow = nrow(profile), ncol = ncol(profile),
    dimnames = dimnames(profile))
  for (j in seq_along(cell_types)) {
    other <- if (ncol(profile) == 2L) {
      profile[, -j]
    } else {
      apply(profile[, -j, drop = FALSE], 1L, max)
    }
    margin <- pmax(profile[, j] - other, 0)
    margin[margin < min_marker_logfc] <- 0
    if (!is.null(top_markers_per_type)) {
      n_keep <- max(1L, as.integer(top_markers_per_type))
      keep <- head(order(margin, decreasing = TRUE), n_keep)
      drop <- setdiff(seq_along(margin), keep)
      margin[drop] <- 0
    }
    weights[, j] <- margin
  }
  weights <- weights ^ specificity_power

  factor_loadings <- pmax(as.matrix(h), 0) ^ factor_power
  h_norm <- sqrt(rowSums(factor_loadings ^ 2))
  w_norm <- sqrt(colSums(weights ^ 2))
  numerator <- factor_loadings %*% weights
  denom <- outer(h_norm, w_norm, `*`)
  score_matrix <- numerator / denom
  score_matrix[!is.finite(score_matrix)] <- 0
  colnames(score_matrix) <- cell_types
  rownames(score_matrix) <- paste0("F", seq_len(nrow(score_matrix)))

  score_df <- as.data.frame(as.table(score_matrix), stringsAsFactors = FALSE)
  colnames(score_df) <- c("factor_label", "cell_type", "score")
  score_df$factor <- as.integer(sub("^F", "", score_df$factor_label))
  score_df <- score_df[, c("factor", "factor_label", "cell_type", "score")]
  score_df <- do.call(rbind, lapply(split(score_df, score_df$factor), function(x) {
    ord <- order(x$score, decreasing = TRUE)
    x$rank <- match(seq_len(nrow(x)), ord)
    sorted <- x$score[ord]
    x$best_score <- sorted[[1]]
    x$second_score <- sorted[[min(2L, length(sorted))]]
    x$margin <- x$score - x$second_score
    x$is_best <- x$rank == 1L
    x
  }))
  rownames(score_df) <- NULL
  score_df <- score_df[order(score_df$factor, score_df$rank), , drop = FALSE]

  top_markers <- .celladmix_factor_source_top_markers(profile, weights,
    n_genes = top_genes)
  marker_weights <- as.data.frame(as.table(weights), stringsAsFactors = FALSE)
  colnames(marker_weights) <- c("gene", "cell_type", "weight")
  marker_weights <- marker_weights[marker_weights$weight > 0, , drop = FALSE]

  result <- list(
    method = "factor_sources",
    scores = score_df,
    score_matrix = score_matrix,
    marker_weights = marker_weights,
    top_markers = top_markers,
    profile = profile,
    factor_loadings = factor_loadings,
    params = list(
      normalize = normalize,
      scale_factor = scale_factor,
      specificity_power = specificity_power,
      factor_power = factor_power,
      min_marker_logfc = min_marker_logfc,
      top_markers_per_type = top_markers_per_type,
      n_genes = length(common_genes),
      n_cells = length(labels),
      n_cell_types = length(cell_types)
    )
  )
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(score_df, file.path(out_dir, "summary.csv"), row.names = FALSE)
    utils::write.csv(top_markers, file.path(out_dir, "top_markers.csv"), row.names = FALSE)
    .celladmix_write_json(result$params, file.path(out_dir, "source_score_metadata.json"))
    result$paths <- list(score_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE))
  } else {
    result$paths <- list()
  }
  result
}

.celladmix_factor_source_top_markers <- function(profile, weights, n_genes = 25L) {
  out <- lapply(colnames(weights), function(cell_type) {
    vals <- weights[, cell_type]
    ord <- head(order(vals, decreasing = TRUE), as.integer(n_genes))
    ord <- ord[vals[ord] > 0]
    if (!length(ord)) {
      return(data.frame())
    }
    data.frame(
      cell_type = cell_type,
      rank = seq_along(ord),
      gene = rownames(weights)[ord],
      marker_weight = as.numeric(vals[ord]),
      expression = as.numeric(profile[ord, cell_type]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  if (is.null(out)) {
    out <- data.frame(cell_type = character(), rank = integer(), gene = character(),
      marker_weight = numeric(), expression = numeric(), stringsAsFactors = FALSE)
  }
  rownames(out) <- NULL
  out
}

celladmix_factor_source_annotation <- function(result, min_score = 0.15, min_margin = 0.03) {
  scores <- result$scores
  best <- scores[scores$is_best, , drop = FALSE]
  best$called <- best$score >= min_score & best$margin >= min_margin
  best$source_cell_type <- ifelse(best$called, best$cell_type, NA_character_)
  best[, c("factor", "factor_label", "source_cell_type", "score", "margin", "called")]
}

.celladmix_factor_source_prior_annotation <- function(
    source_prior = NULL,
    min_score = 0.15,
    min_margin = 0.03
  ) {
  if (is.null(source_prior)) {
    return(NULL)
  }
  if (inherits(source_prior, "CellAdmixFactorSourceScore")) {
    return(source_prior$annotation(min_score = min_score, min_margin = min_margin))
  }
  if (is.list(source_prior) && identical(source_prior$method %||% NULL, "factor_sources") &&
      is.data.frame(source_prior$scores)) {
    return(celladmix_factor_source_annotation(source_prior,
      min_score = min_score, min_margin = min_margin))
  }
  if (is.data.frame(source_prior)) {
    out <- source_prior
    if (!("factor" %in% names(out)) && "factor_label" %in% names(out)) {
      out$factor <- as.integer(sub("^F", "", out$factor_label))
    }
    if (!("factor_label" %in% names(out)) && "factor" %in% names(out)) {
      out$factor_label <- paste0("F", as.integer(out$factor))
    }
    if (!("source_cell_type" %in% names(out)) && "cell_type" %in% names(out)) {
      out$source_cell_type <- out$cell_type
    }
    missing <- setdiff(c("factor", "factor_label", "source_cell_type"), names(out))
    if (length(missing)) {
      stop("source_prior data frame is missing columns: ", paste(missing, collapse = ", "))
    }
    if (!("called" %in% names(out))) {
      out$called <- !is.na(out$source_cell_type) & nzchar(out$source_cell_type)
    }
    if (!("score" %in% names(out))) {
      out$score <- NA_real_
    }
    if (!("margin" %in% names(out))) {
      out$margin <- NA_real_
    }
    return(out[, c("factor", "factor_label", "source_cell_type", "score", "margin", "called")])
  }
  if (is.atomic(source_prior) && !is.null(names(source_prior))) {
    factor_label <- names(source_prior)
    factor <- suppressWarnings(as.integer(sub("^F", "", factor_label)))
    return(data.frame(
      factor = factor,
      factor_label = ifelse(is.na(factor), factor_label, paste0("F", factor)),
      source_cell_type = as.character(source_prior),
      score = NA_real_,
      margin = NA_real_,
      called = !is.na(source_prior) & nzchar(source_prior),
      stringsAsFactors = FALSE
    ))
  }
  stop("source_prior must be a CellAdmixFactorSourceScore, result list, data frame, or named vector")
}

celladmix_factor_source_top_genes <- function(result, n_genes = 8L, sources = NULL) {
  ann <- celladmix_factor_source_annotation(result)
  if (!is.null(sources)) {
    by_label <- sources[ann$factor_label]
    by_factor <- sources[as.character(ann$factor)]
    ann$source_cell_type <- ifelse(!is.na(by_label), by_label, by_factor)
  }
  factor_loadings <- result$factor_loadings
  marker_weights <- result$marker_weights
  out <- lapply(seq_len(nrow(ann)), function(i) {
    source <- ann$source_cell_type[[i]]
    if (is.na(source) || !nzchar(source)) {
      return(data.frame())
    }
    weights <- marker_weights[marker_weights$cell_type == source, , drop = FALSE]
    if (!nrow(weights)) {
      return(data.frame())
    }
    all_loadings <- as.numeric(factor_loadings[ann$factor[[i]], , drop = TRUE])
    names(all_loadings) <- colnames(factor_loadings)
    loading_sum <- sum(all_loadings, na.rm = TRUE)
    loading_fraction <- if (is.finite(loading_sum) && loading_sum > 0) {
      all_loadings / loading_sum
    } else {
      rep(0, length(all_loadings))
    }
    factor_max_loading_fraction <- max(loading_fraction, na.rm = TRUE)
    factor_max_loading_fraction <- if (is.finite(factor_max_loading_fraction)) {
      factor_max_loading_fraction
    } else {
      0
    }
    h <- all_loadings[weights$gene]
    h_fraction <- loading_fraction[weights$gene]
    contribution <- as.numeric(h) * weights$weight
    ord <- head(order(contribution, decreasing = TRUE), as.integer(n_genes))
    data.frame(
      factor = ann$factor[[i]],
      factor_label = ann$factor_label[[i]],
      source_cell_type = source,
      rank = seq_along(ord),
      gene = weights$gene[ord],
      contribution = contribution[ord],
      loading = as.numeric(h[ord]),
      loading_fraction = as.numeric(h_fraction[ord]),
      factor_max_loading_fraction = factor_max_loading_fraction,
      marker_weight = weights$weight[ord],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  if (is.null(out)) {
    out <- data.frame()
  }
  rownames(out) <- NULL
  out
}

celladmix_plot_factor_source_heatmap <- function(
    result,
    min_score = 0.15,
    min_margin = 0.03,
    main = NULL,
    max_score = NULL,
    show_source_label = TRUE
  ) {
  main <- main %||% "Factor source gene-content score"
  scores <- result$scores
  if (!nrow(scores)) {
    graphics::plot.new()
    graphics::title(main)
    graphics::text(0.5, 0.5, "No factor source scores to plot")
    return(invisible(result))
  }
  factors <- paste0("F", sort(unique(scores$factor)))
  cell_types <- sort(unique(scores$cell_type))
  mat <- matrix(0, nrow = length(cell_types), ncol = length(factors),
    dimnames = list(cell_types, factors))
  idx <- cbind(match(scores$cell_type, cell_types), match(scores$factor_label, factors))
  mat[idx] <- scores$score
  max_score <- max_score %||% max(mat, na.rm = TRUE)
  max_score <- if (is.finite(max_score) && max_score > 0) max_score else 1
  values <- pmin(mat, max_score)

  ann <- celladmix_factor_source_annotation(result,
    min_score = min_score, min_margin = min_margin)
  source_calls <- stats::setNames(ann$source_cell_type, ann$factor_label)

  palette_fun <- grDevices::colorRampPalette(c("white", "#2166AC"))
  palette <- palette_fun(100L)
  color_for <- function(value) {
    if (!is.finite(value) || value <= 0) {
      return("white")
    }
    scaled <- value / max(max_score, .Machine$double.eps)
    palette[max(1L, min(100L, ceiling(scaled * 100)))]
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  nr <- nrow(values)
  nc <- ncol(values)
  left_margin <- max(8, min(18, 4 + 0.35 * max(nchar(rownames(values)))))
  graphics::par(mar = c(4, left_margin, 3, 5), xpd = NA)
  graphics::plot(NA, xlim = c(0.5, nc + 1.45), ylim = c(0.5, nr + 0.5),
    axes = FALSE, xlab = "Factor", ylab = "", main = main)
  graphics::axis(1, at = seq_len(nc), labels = colnames(values), las = 2)
  graphics::axis(2, at = rev(seq_len(nr)), labels = rownames(values), las = 2)
  for (i in seq_len(nr)) {
    y <- nr - i + 1L
    for (j in seq_len(nc)) {
      graphics::rect(j - 0.5, y - 0.5, j + 0.5, y + 0.5,
        col = color_for(values[i, j]), border = "grey70")
      source <- source_calls[[colnames(values)[j]]]
      if (isTRUE(show_source_label) && !is.null(source) &&
          !is.na(source) && identical(source, rownames(values)[i])) {
        graphics::text(j, y, "S", font = 2, cex = 1.1)
      }
    }
  }

  key_n <- 80L
  key_values <- seq(0, max_score, length.out = key_n)
  key_y <- seq(0.8, nr + 0.2, length.out = key_n + 1L)
  key_x0 <- nc + 0.78
  key_x1 <- nc + 0.96
  for (k in seq_len(key_n)) {
    graphics::rect(key_x0, key_y[k], key_x1, key_y[k + 1L],
      col = color_for(key_values[k]), border = NA)
  }
  graphics::rect(key_x0, key_y[1L], key_x1, key_y[key_n + 1L],
    border = "grey40", col = NA)
  graphics::text(key_x1 + 0.08, key_y[1L], "0",
    adj = c(0, 0.5), cex = 0.75)
  graphics::text(key_x1 + 0.08, key_y[key_n + 1L], sprintf("%.2f", max_score),
    adj = c(0, 0.5), cex = 0.75)
  graphics::text(key_x1 + 0.34, mean(range(key_y)), "source score",
    srt = 90, cex = 0.8)
  invisible(result)
}

celladmix_plot_factor_source_top_genes <- function(
    result,
    n_genes = 6L,
    ncol = 4L,
    sources = NULL,
    title = "Marker genes supporting factor source calls"
  ) {
  .celladmix_require_ggplot2()
  df <- celladmix_factor_source_top_genes(result, n_genes = n_genes,
    sources = sources)
  if (!nrow(df)) {
    return(.celladmix_empty_plot("No marker-weighted factor genes available"))
  }
  df$panel <- paste0(df$factor_label, "\n", df$source_cell_type)
  df <- df[order(df$factor, df$rank), , drop = FALSE]
  df$gene_panel <- paste(df$gene, df$panel, sep = "\r")
  df$gene_panel <- factor(df$gene_panel, levels = rev(unique(df$gene_panel)))
  anchors <- df[!duplicated(df$panel), c("panel", "gene_panel", "factor_max_loading_fraction")]
  ggplot2::ggplot(df, ggplot2::aes(x = .data$loading_fraction, y = .data$gene_panel)) +
    ggplot2::geom_blank(data = anchors,
      ggplot2::aes(x = .data$factor_max_loading_fraction, y = .data$gene_panel),
      inherit.aes = FALSE) +
    ggplot2::geom_col(fill = "grey25", width = 0.75) +
    ggplot2::facet_wrap(~ panel, scales = "free", ncol = max(1L, as.integer(ncol))) +
    ggplot2::scale_x_continuous(position = "top",
      expand = ggplot2::expansion(mult = c(0, 0.04))) +
    ggplot2::scale_y_discrete(labels = function(x) sub("\r.*$", "", x)) +
    ggplot2::labs(title = title, x = "Factor loading fraction", y = NULL) +
    ggplot2::theme_bw(base_size = 8.5) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}
