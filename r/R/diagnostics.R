#' Lightweight Diagnostics for cellAdmix Fits
#'
#' These helpers provide reusable downstream summaries that are useful in
#' examples and notebooks, but are intentionally separate from the core C++
#' fitting, scoring, and correction path.
#' @noRd
.celladmix_require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for this plot")
  }
}

.celladmix_empty_plot <- function(label) {
  .celladmix_require_ggplot2()
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = label, size = 4) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_void()
}

.celladmix_matrix_col_sums <- function(x) {
  if (inherits(x, "Matrix")) Matrix::colSums(x) else colSums(x)
}

.celladmix_matrix_row_means <- function(x) {
  if (inherits(x, "Matrix")) Matrix::rowMeans(x) else rowMeans(x)
}

.celladmix_normalize_counts_matrix <- function(
    counts,
    normalize = c("cpm", "log1p_cpm", "none"),
    scale_factor = 10000
  ) {
  normalize <- match.arg(normalize)
  if (identical(normalize, "none")) {
    return(counts)
  }
  lib <- .celladmix_matrix_col_sums(counts)
  scale <- scale_factor / pmax(lib, 1)
  norm <- if (inherits(counts, "Matrix")) {
    counts %*% Matrix::Diagonal(x = scale)
  } else {
    sweep(counts, 2L, scale, `*`)
  }
  if (identical(normalize, "log1p_cpm")) {
    norm <- log1p(norm)
  }
  norm
}

.celladmix_named_groups <- function(groups, cell_ids) {
  if (is.null(groups)) {
    stop("groups must be provided")
  }
  groups <- as.character(groups)
  if (is.null(names(groups))) {
    if (length(groups) != length(cell_ids)) {
      stop("Unnamed groups must have one value per count-matrix column")
    }
    names(groups) <- cell_ids
  }
  groups
}

.celladmix_de_column <- function(de, candidates, required = TRUE) {
  hit <- candidates[candidates %in% colnames(de)][1]
  if (is.na(hit) && required) {
    stop("DE table is missing expected columns: ", paste(candidates, collapse = " or "))
  }
  hit
}

.celladmix_standardize_de <- function(de) {
  if (is.null(de) || !nrow(de)) {
    return(NULL)
  }
  gene_col <- .celladmix_de_column(de, c("gene", "Gene"))
  logfc_col <- .celladmix_de_column(de, c("logFC", "M", "logfc"))
  padj_col <- .celladmix_de_column(de, c("p_adj", "padj", "q_value"))
  p_col <- .celladmix_de_column(de, c("p_value", "pval", "p"), required = FALSE)
  out <- data.frame(
    gene = as.character(de[[gene_col]]),
    logFC = as.numeric(de[[logfc_col]]),
    p_adj = as.numeric(de[[padj_col]]),
    stringsAsFactors = FALSE
  )
  out$p_value <- if (!is.na(p_col)) as.numeric(de[[p_col]]) else NA_real_
  for (nm in setdiff(colnames(de), c(gene_col, logfc_col, padj_col, p_col))) {
    out[[nm]] <- de[[nm]]
  }
  out
}

.celladmix_marker_frame <- function(markers, default_set = "Marker") {
  if (is.null(markers)) {
    return(data.frame(set = character(), gene = character(), stringsAsFactors = FALSE))
  }
  if (is.character(markers) || is.factor(markers)) {
    return(data.frame(set = default_set, gene = as.character(markers),
      stringsAsFactors = FALSE))
  }
  if (is.data.frame(markers)) {
    set_col <- c("set", "marker_set", "label")[c("set", "marker_set", "label") %in% colnames(markers)][1]
    if (is.na(set_col)) {
      stop("markers data.frame must contain a set, marker_set, or label column")
    }
    if (!("gene" %in% colnames(markers))) {
      stop("markers data.frame must contain a gene column")
    }
    return(data.frame(set = as.character(markers[[set_col]]),
      gene = as.character(markers$gene), stringsAsFactors = FALSE))
  }
  if (!is.list(markers)) {
    stop("markers must be a character vector, named list, or data.frame")
  }
  if (is.null(names(markers)) || any(!nzchar(names(markers)))) {
    names(markers) <- paste0("set_", seq_along(markers))
  }
  do.call(rbind, lapply(names(markers), function(nm) {
    data.frame(set = nm, gene = as.character(markers[[nm]]), stringsAsFactors = FALSE)
  }))
}

.celladmix_marker_labels <- function(label, marker_genes) {
  if (isTRUE(label)) {
    return(marker_genes)
  }
  if (isFALSE(label) || is.null(label)) {
    return(character())
  }
  if (is.character(label)) {
    return(as.character(label))
  }
  stop("label must be TRUE, FALSE, NULL, or a character vector")
}

.celladmix_marker_palette <- function(levels) {
  marker_levels <- levels[levels != "Other gene"]
  palette <- c("#B2182B", "#2166AC", "#1B7837", "#762A83", "#E08214", "#4D4D4D")
  cols <- stats::setNames(rep(palette, length.out = max(1L, length(marker_levels))),
    marker_levels)
  c(cols, "Other gene" = "grey70")
}

.celladmix_title_label <- function(x) {
  x <- as.character(x %||% "Score")
  if (!length(x) || is.na(x[[1]]) || !nzchar(x[[1]])) {
    return("Score")
  }
  x <- x[[1]]
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

.celladmix_score_column <- function(df) {
  score_candidates <- c("mean_score", "score", "mean_delta_score",
    "mean_coherence_score", "mean_raw_score")
  score_col <- score_candidates[score_candidates %in% colnames(df)][1]
  if (is.na(score_col)) {
    stop("score table is missing a score column: ",
      paste(score_candidates, collapse = ", "))
  }
  score_col
}

.celladmix_factor_id <- function(factor) {
  ids <- .celladmix_factor_ids(factor)
  if (length(ids) != 1L) {
    stop("factor must be a single factor ID such as 1 or 'F1'")
  }
  ids[[1]]
}

.celladmix_factor_ids <- function(factor) {
  if (is.null(factor)) {
    return(NULL)
  }
  out <- suppressWarnings(as.integer(sub("^F", "", as.character(factor))))
  if (any(is.na(out) | out < 1L)) {
    stop("factor must contain positive factor IDs such as 1 or 'F1'")
  }
  unique(out)
}

.celladmix_score_summary_factor_ids <- function(summary) {
  ids <- suppressWarnings(as.integer(sub("^F", "", as.character(summary$factor))))
  sort(unique(ids[is.finite(ids) & ids > 0]))
}

.celladmix_triangle_vertices <- function(df) {
  if (is.null(df) || !nrow(df)) {
    return(data.frame())
  }
  pieces <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    x <- df$x[[i]]
    y <- df$y[[i]]
    if (identical(as.character(df$triangle[[i]]), "upper")) {
      xs <- c(x - 0.5, x + 0.5, x + 0.5)
      ys <- c(y + 0.5, y + 0.5, y - 0.5)
    } else {
      xs <- c(x - 0.5, x - 0.5, x + 0.5)
      ys <- c(y + 0.5, y - 0.5, y - 0.5)
    }
    pieces[[i]] <- cbind(df[rep(i, 3L), , drop = FALSE],
      x_vertex = xs, y_vertex = ys, polygon_id = i)
  }
  do.call(rbind, pieces)
}

.celladmix_extract_score_tables <- function(x, method_label = NULL) {
  if (inherits(x, "CellAdmixScore")) {
    return(list(pairs = x$pairs(), summary = x$summary(),
      method_label = method_label %||% x$name %||% x$method,
      method = x$method))
  }
  if (is.data.frame(x)) {
    return(list(pairs = x, summary = x, method_label = method_label,
      method = NULL))
  }
  if (is.list(x)) {
    return(list(pairs = x$scores %||% x$summary,
      summary = x$summary %||% x$scores,
      method_label = method_label %||% x$name %||% x$method,
      method = x$method %||% NULL))
  }
  stop("x must be a CellAdmixScore, score result list, or score-pairs data.frame")
}

.celladmix_significance_column <- function(df, method = NULL) {
  if (identical(method, "coherence")) {
    candidates <- c("patch_neg_log10_p", "neg_log10_p", "plot_score")
  } else {
    candidates <- c("neg_log10_p", "patch_neg_log10_p", "plot_score")
  }
  hit <- candidates[candidates %in% colnames(df)][1]
  if (!is.na(hit)) {
    return(hit)
  }
  if ("p_value" %in% colnames(df)) {
    return("p_value")
  }
  stop("score summary is missing a significance column: ",
    paste(c(candidates, "p_value"), collapse = ", "))
}

.celladmix_significance_values <- function(df, column) {
  values <- as.numeric(df[[column]])
  if (identical(column, "p_value")) {
    values <- -log10(values)
  }
  values[is.nan(values) | (is.infinite(values) & values < 0)] <- NA_real_
  values
}

.celladmix_source_score_summary <- function(summary, factor_id, method = NULL) {
  summary_factor <- suppressWarnings(as.integer(sub("^F", "", as.character(summary$factor))))
  df <- summary[summary_factor == factor_id, , drop = FALSE]
  if (!nrow(df)) {
    return(list(scores = data.frame(source_cell_type = character(),
      source_score = numeric(), stringsAsFactors = FALSE), label = "source score"))
  }
  if (identical(method, "coherence") && "source_probability" %in% names(df)) {
    source_col <- "source_probability"
    source_values <- as.numeric(df[[source_col]])
    fun <- function(x) max(x, na.rm = TRUE)
    label <- "source probability"
  } else if (identical(method, "coherence") && "source_log_enrichment" %in% names(df)) {
    source_col <- "source_log_enrichment"
    source_values <- as.numeric(df[[source_col]])
    fun <- function(x) max(x, na.rm = TRUE)
    label <- "source enrichment"
  } else {
    source_col <- .celladmix_significance_column(df, method = method)
    source_values <- .celladmix_significance_values(df, source_col)
    fun <- function(x) mean(x, na.rm = TRUE)
    label <- if (identical(source_col, "p_value")) "mean -log10(p)" else
      sprintf("mean %s", source_col)
  }
  source_df <- data.frame(source_cell_type = as.character(df$source_cell_type),
    source_score = source_values, stringsAsFactors = FALSE)
  source_df <- source_df[is.finite(source_df$source_score) &
    !is.na(source_df$source_cell_type), , drop = FALSE]
  if (!nrow(source_df)) {
    return(list(scores = data.frame(source_cell_type = character(),
      source_score = numeric(), stringsAsFactors = FALSE), label = label))
  }
  scores <- stats::aggregate(source_score ~ source_cell_type, source_df, fun)
  scores$source_score[!is.finite(scores$source_score)] <- NA_real_
  list(scores = scores, label = label)
}

.celladmix_source_score_plot_data <- function(source_scores, cell_levels, y) {
  base <- data.frame(source_cell_type = cell_levels, stringsAsFactors = FALSE)
  df <- merge(base, source_scores, by = "source_cell_type", all.x = TRUE, sort = FALSE)
  df <- df[match(cell_levels, df$source_cell_type), , drop = FALSE]
  df$x <- match(df$source_cell_type, cell_levels)
  finite <- df$source_score[is.finite(df$source_score)]
  best <- if (length(finite)) which.max(replace(df$source_score, !is.finite(df$source_score), -Inf)) else integer()
  df$source_score_label <- ifelse(is.finite(df$source_score),
    format(signif(df$source_score, 2), trim = TRUE), "")
  df$source_fill <- "#F2F2F2"
  if (length(finite)) {
    lo <- min(finite)
    hi <- max(finite)
    scaled <- if (hi > lo) (df$source_score - lo) / (hi - lo) else
      ifelse(is.finite(df$source_score), 1, NA_real_)
    pal <- grDevices::colorRampPalette(c("white", "#525252"))(101L)
    idx <- pmax(1L, pmin(101L, round(scaled * 100) + 1L))
    df$source_fill[is.finite(scaled)] <- pal[idx[is.finite(scaled)]]
  }
  df$source_best <- FALSE
  if (length(best) == 1L && is.finite(df$source_score[[best]])) {
    df$source_best[[best]] <- TRUE
  }
  df$y <- y
  df
}

.celladmix_plot_score_pair_matrix <- function(
    summary,
    factor,
    method_label,
    method = NULL,
    max_abs_significance = 6
  ) {
  .celladmix_require_ggplot2()
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("The cowplot package is required for factor-specific pair plots")
  }
  if (is.null(summary) || !is.data.frame(summary) || !nrow(summary)) {
    return(.celladmix_empty_plot(paste(method_label, "pair summaries unavailable")))
  }
  required <- c("target_cell_type", "source_cell_type", "factor")
  missing <- setdiff(required, colnames(summary))
  if (length(missing)) {
    stop("score summary is missing required columns: ", paste(missing, collapse = ", "))
  }
  factor_id <- .celladmix_factor_id(factor)
  summary_factor <- suppressWarnings(as.integer(sub("^F", "", as.character(summary$factor))))
  df <- summary[summary_factor == factor_id, , drop = FALSE]
  if (!nrow(df)) {
    return(.celladmix_empty_plot(paste0(method_label, " F", factor_id, " unavailable")))
  }
  score_col <- .celladmix_score_column(df)
  significance_col <- .celladmix_significance_column(df, method = method)
  df$score_value <- as.numeric(df[[score_col]])
  df$significance_value <- .celladmix_significance_values(df, significance_col)
  df$signed_significance <- sign(df$score_value) * df$significance_value
  df <- df[is.finite(df$signed_significance), , drop = FALSE]
  if (!nrow(df)) {
    return(.celladmix_empty_plot(paste0(method_label, " F", factor_id, " finite scores unavailable")))
  }
  mean_finite <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) mean(x) else NA_real_
  }
  agg <- stats::aggregate(signed_significance ~ target_cell_type + source_cell_type,
    df, mean_finite)
  cell_levels <- sort(unique(c(agg$target_cell_type, agg$source_cell_type)))
  base <- expand.grid(target_cell_type = cell_levels, source_cell_type = cell_levels,
    stringsAsFactors = FALSE)
  plot_df <- merge(base, agg, by = c("target_cell_type", "source_cell_type"),
    all.x = TRUE, sort = FALSE)
  factor_y_levels <- rev(cell_levels)
  plot_df$x <- match(plot_df$source_cell_type, cell_levels)
  plot_df$y <- match(plot_df$target_cell_type, factor_y_levels)
  source_y <- length(factor_y_levels) + 0.72
  sig_limit <- max_abs_significance
  if (is.null(sig_limit)) {
    finite <- plot_df$signed_significance[is.finite(plot_df$signed_significance)]
    sig_limit <- max(abs(finite), na.rm = TRUE)
    if (!is.finite(sig_limit) || sig_limit <= 0) sig_limit <- 1
  }
  sig_limit <- max(abs(sig_limit), .Machine$double.eps)
  plot_df$signed_significance_plot <- pmax(-sig_limit,
    pmin(sig_limit, plot_df$signed_significance))

  grid <- expand.grid(x = seq_along(cell_levels), y = seq_along(factor_y_levels))
  source_info <- .celladmix_source_score_summary(summary, factor_id, method = method)
  source_df <- .celladmix_source_score_plot_data(source_info$scores, cell_levels,
    y = source_y)
  heatmap <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = grid, ggplot2::aes(x, y), fill = "white",
      color = "grey86", linewidth = 0.22) +
    ggplot2::geom_tile(data = plot_df,
      ggplot2::aes(x, y, fill = signed_significance_plot),
      color = "white", linewidth = 0.18) +
    ggplot2::geom_tile(data = source_df, ggplot2::aes(x, y),
      fill = source_df$source_fill, color = "white", linewidth = 0.2,
      height = 0.34) +
    ggplot2::geom_text(data = source_df, ggplot2::aes(x, y, label = source_score_label),
      size = 2.3, color = "black") +
    ggplot2::geom_tile(data = source_df[source_df$source_best, , drop = FALSE],
      ggplot2::aes(x, y), fill = NA, color = "black", linewidth = 0.65,
      height = 0.34) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-sig_limit, sig_limit), na.value = "#F2F2F2",
      name = "Signed -log10(p)",
      guide = ggplot2::guide_colorbar(direction = "horizontal",
        title.position = "top", barwidth = grid::unit(2.6, "cm"),
        barheight = grid::unit(0.18, "cm"))) +
    ggplot2::scale_x_continuous(breaks = seq_along(cell_levels), labels = cell_levels,
      expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = c(seq_along(factor_y_levels), source_y),
      labels = c(factor_y_levels, "source evidence"),
      limits = c(0.5, source_y + 0.18), expand = c(0, 0)) +
    ggplot2::coord_fixed(expand = FALSE) +
    ggplot2::labs(title = sprintf("%s F%d cell-type pair scores", method_label, factor_id),
      x = "Source cell type", y = "Target cell type") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 7),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(5.5, 4, 5.5, 12))
  heatmap
}

.celladmix_plot_score_pair_matrices <- function(
    summary,
    factors = NULL,
    method_label,
    method = NULL,
    max_abs_significance = 6,
    ncol = 2
  ) {
  if (is.null(summary) || !is.data.frame(summary) || !nrow(summary)) {
    return(.celladmix_empty_plot(paste(method_label, "pair summaries unavailable")))
  }
  factor_ids <- if (is.null(factors)) .celladmix_score_summary_factor_ids(summary) else
    .celladmix_factor_ids(factors)
  if (!length(factor_ids)) {
    return(.celladmix_empty_plot(paste(method_label, "factor summaries unavailable")))
  }
  plots <- lapply(factor_ids, function(factor_id) {
    .celladmix_plot_score_pair_matrix(summary, factor = factor_id,
      method_label = method_label, method = method,
      max_abs_significance = max_abs_significance)
  })
  if (length(plots) == 1L) {
    return(plots[[1]])
  }
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("The cowplot package is required for multi-factor pair plots")
  }
  cowplot::plot_grid(plotlist = plots, ncol = ncol, align = "hv", axis = "tblr")
}

#' Plot score strength by factor and target-source cell-type pair.
#'
#' This diagnostic summarizes the pair-level score table returned by bridge,
#' membrane, or coherence scoring. It can be called directly on a
#' [CellAdmixScore] object, on a score result list with a `scores` table, or on
#' the score table itself.
#'
#' @param x A [CellAdmixScore], score result list, or score-pairs data.frame.
#' @param method_label Optional label used in the plot title.
#' @param top_pairs Number of target-source cell-type pairs to show.
#' @param used_only If `TRUE`, restrict to rows marked `used_in_summary` when
#' that column is available.
#' @param show_top_pairs If `TRUE`, draw the older compact summary of top
#'   target-source pairs across all factors. The default draws factor-specific
#'   directional target-by-source significance matrices.
#' @param factor Optional factor IDs such as `1`, `"F1"`, or `c("F1", "F2")`.
#'   When omitted, the default factor-specific view draws all factors.
#' @param max_abs_significance Symmetric color limit for factor-specific signed
#'   `-log10(p)` values.
#' @param ncol Number of columns for multi-factor factor-specific plots.
#' @export
celladmix_plot_score_pairs <- function(
    x,
    method_label = NULL,
    top_pairs = 12,
    used_only = TRUE,
    show_top_pairs = FALSE,
    factor = NULL,
    max_abs_significance = 6,
    ncol = 2
  ) {
  .celladmix_require_ggplot2()

  tables <- .celladmix_extract_score_tables(x, method_label = method_label)
  scores <- tables$pairs
  summary <- tables$summary
  method_label <- tables$method_label
  method <- tables$method
  method_label <- .celladmix_title_label(method_label %||% "Score")
  if (!isTRUE(show_top_pairs)) {
    return(.celladmix_plot_score_pair_matrices(summary, factors = factor,
      method_label = method_label, method = method,
      max_abs_significance = max_abs_significance, ncol = ncol))
  }

  if (is.null(scores) || !nrow(scores)) {
    return(.celladmix_empty_plot(paste(method_label, "scores unavailable")))
  }
  required <- c("target_cell_type", "source_cell_type", "factor")
  missing <- setdiff(required, colnames(scores))
  if (length(missing)) {
    stop("score table is missing required columns: ", paste(missing, collapse = ", "))
  }
  score_col <- .celladmix_score_column(scores)

  df <- scores
  if (isTRUE(used_only) && "used_in_summary" %in% names(df)) {
    df <- df[df$used_in_summary %in% TRUE, , drop = FALSE]
  }
  df$score_value <- as.numeric(df[[score_col]])
  df <- df[is.finite(df$score_value), , drop = FALSE]
  if (!nrow(df)) {
    return(.celladmix_empty_plot(paste(method_label, "finite scores unavailable")))
  }

  df$pair <- paste(df$target_cell_type, "<-", df$source_cell_type)
  pair_rank <- stats::aggregate(score_value ~ pair, df, max)
  pair_order <- pair_rank$pair[order(pair_rank$score_value, decreasing = TRUE)]
  keep_pairs <- head(pair_order, max(1, as.integer(top_pairs)))
  df <- df[df$pair %in% keep_pairs, , drop = FALSE]

  mean_df <- stats::aggregate(score_value ~ pair + factor, df, mean)
  n_df <- stats::aggregate(score_value ~ pair + factor, df, length)
  names(n_df)[names(n_df) == "score_value"] <- "n_pairs"
  plot_df <- merge(mean_df, n_df, by = c("pair", "factor"), all.x = TRUE)

  factor_chr <- as.character(plot_df$factor)
  factor_num <- suppressWarnings(as.numeric(factor_chr))
  factor_order <- unique(factor_chr[order(ifelse(is.na(factor_num), Inf, factor_num), factor_chr)])
  factor_labels <- ifelse(grepl("^F", factor_order), factor_order, paste0("F", factor_order))
  plot_df$factor_label <- ifelse(grepl("^F", factor_chr), factor_chr, paste0("F", factor_chr))
  plot_df$factor_label <- factor(plot_df$factor_label, levels = factor_labels)
  plot_df$pair <- factor(plot_df$pair, levels = rev(keep_pairs))

  ggplot2::ggplot(plot_df, ggplot2::aes(factor_label, pair, fill = score_value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = "Mean score") +
    ggplot2::labs(title = paste(method_label, "factor-pair score summary"),
      subtitle = sprintf("Top %d cell-type pairs by maximum factor score",
        length(keep_pairs)),
      x = "Factor", y = "Target <- source") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
      vjust = 0.5), plot.title = ggplot2::element_text(hjust = 0.5))
}

.celladmix_score_annotation <- function(x, p_thresh, adjust_p, ...) {
  if (inherits(x, "CellAdmixScore")) {
    return(x$annotation(p_thresh = p_thresh, adjust_p = adjust_p, ...))
  }
  if (is.list(x) && !is.null(x$score_matrix) && !is.null(x$source_calls)) {
    return(x)
  }
  stop("Each score must be a CellAdmixScore or an annotation list returned by score$annotation()")
}

.celladmix_score_factor_ids <- function(annotation) {
  ids <- integer()
  if (!is.null(annotation$score_matrix) && ncol(annotation$score_matrix)) {
    ids <- c(ids, as.integer(sub("^F", "", colnames(annotation$score_matrix))))
  }
  if (!is.null(annotation$source_calls) && length(annotation$source_calls)) {
    ids <- c(ids, as.integer(sub("^f_", "", names(annotation$source_calls))))
  }
  sort(unique(ids[is.finite(ids)]))
}

.celladmix_source_agreement <- function(annotations) {
  factor_ids <- sort(unique(unlist(lapply(annotations, .celladmix_score_factor_ids),
    use.names = FALSE)))
  source_for <- function(annotation, factor_ids) {
    out <- rep("unassigned", length(factor_ids))
    names(out) <- factor_ids
    calls <- annotation$source_calls
    if (!is.null(calls) && length(calls)) {
      call_ids <- as.integer(sub("^f_", "", names(calls)))
      for (i in seq_along(call_ids)) {
        value <- as.character(calls[[i]])
        if (is.finite(call_ids[[i]]) && length(value) && nzchar(value[[1]])) {
          out[as.character(call_ids[[i]])] <- value[[1]]
        }
      }
    }
    unname(out)
  }
  pieces <- lapply(names(annotations), function(method) {
    data.frame(factor = factor_ids, factor_label = paste0("F", factor_ids),
      method = method, source = source_for(annotations[[method]], factor_ids),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, pieces)
  out$method <- factor(out$method, levels = names(annotations))
  out
}

.celladmix_admixture_agreement <- function(annotations, p_thresh, max_neg_log10) {
  matrices <- lapply(annotations, `[[`, "score_matrix")
  factor_ids <- sort(unique(unlist(lapply(annotations, .celladmix_score_factor_ids),
    use.names = FALSE)))
  targets <- sort(unique(unlist(lapply(matrices, rownames), use.names = FALSE)))
  base <- expand.grid(target_cell_type = targets, factor = factor_ids,
    method = names(annotations), stringsAsFactors = FALSE)
  get_score <- function(mat, target, factor) {
    col <- paste0("F", factor)
    if (is.null(dim(mat)) || !(target %in% rownames(mat)) || !(col %in% colnames(mat))) {
      return(NA_real_)
    }
    as.numeric(mat[target, col])
  }
  base$logp <- mapply(function(method, target, factor) {
    get_score(matrices[[method]], target, factor)
  }, base$method, base$target_cell_type, base$factor)
  base$logp <- pmin(base$logp, max_neg_log10)
  threshold <- -log10(p_thresh)
  base$call <- is.finite(base$logp) & base$logp >= threshold
  base$factor_label <- paste0("F", base$factor)
  base$method <- factor(base$method, levels = names(annotations))
  base
}

.celladmix_agreement_cell_levels <- function(source_agreement, admixture_agreement) {
  vals <- c(source_agreement$source, admixture_agreement$target_cell_type)
  vals <- vals[!is.na(vals) & nzchar(vals) & !vals %in% c("unassigned", "__unannotated__")]
  sort(unique(vals))
}

.celladmix_agreement_factor_levels <- function(source_agreement, admixture_agreement) {
  vals <- c(source_agreement$factor, admixture_agreement$factor)
  paste0("F", sort(unique(vals[is.finite(vals)])))
}

.celladmix_source_pair <- function(source_agreement, method_a, method_b) {
  df <- source_agreement[source_agreement$method %in% c(method_a, method_b),
    c("factor", "factor_label", "method", "source"), drop = FALSE]
  if (!nrow(df)) {
    return(data.frame())
  }
  wide <- reshape(df, idvar = c("factor", "factor_label"), timevar = "method",
    direction = "wide")
  source_a <- paste0("source.", method_a)
  source_b <- paste0("source.", method_b)
  if (!(source_a %in% names(wide))) {
    wide[[source_a]] <- "unassigned"
  }
  if (!(source_b %in% names(wide))) {
    wide[[source_b]] <- "unassigned"
  }
  wide
}

.celladmix_plot_source_pair <- function(source_agreement, method_a, method_b,
                                        cell_levels, factor_levels) {
  pair_df <- .celladmix_source_pair(source_agreement, method_a, method_b)
  if (is.null(pair_df) || !nrow(pair_df)) {
    return(.celladmix_empty_plot("Source-call agreement unavailable"))
  }
  source_a <- paste0("source.", method_a)
  source_b <- paste0("source.", method_b)
  plot_df <- rbind(
    data.frame(factor = pair_df$factor, factor_label = pair_df$factor_label,
      method = method_a, triangle = "upper", cell_type = pair_df[[source_a]],
      stringsAsFactors = FALSE),
    data.frame(factor = pair_df$factor, factor_label = pair_df$factor_label,
      method = method_b, triangle = "lower", cell_type = pair_df[[source_b]],
      stringsAsFactors = FALSE)
  )
  plot_df <- plot_df[plot_df$cell_type %in% cell_levels, , drop = FALSE]
  factor_y_levels <- rev(factor_levels)
  plot_df$x <- match(plot_df$cell_type, cell_levels)
  plot_df$y <- match(plot_df$factor_label, factor_y_levels)
  triangles <- .celladmix_triangle_vertices(plot_df)
  grid <- expand.grid(x = seq_along(cell_levels), y = seq_along(factor_y_levels))
  method_cols <- stats::setNames(c("#D95F02", "#1F78B4"), c(method_a, method_b))

  ggplot2::ggplot() +
    ggplot2::geom_tile(data = grid, ggplot2::aes(x, y), fill = "white",
      color = "grey86", linewidth = 0.22) +
    ggplot2::geom_polygon(data = triangles,
      ggplot2::aes(x_vertex, y_vertex, group = polygon_id, fill = method),
      color = "white", linewidth = 0.18, alpha = 0.95) +
    ggplot2::scale_fill_manual(values = method_cols, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = seq_along(cell_levels), labels = cell_levels,
      expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_along(factor_y_levels),
      labels = factor_y_levels, expand = c(0, 0)) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::labs(title = "Source calls by factor",
      subtitle = sprintf("Upper-right = %s; lower-left = %s", method_a, method_b),
      x = "Cell type", y = "Factor", fill = "Method") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 8, hjust = 0.5),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))
}

.celladmix_admixture_pair <- function(admixture_agreement, method_a, method_b) {
  df <- admixture_agreement[admixture_agreement$method %in% c(method_a, method_b),
    c("target_cell_type", "factor", "factor_label", "method", "logp", "call"),
    drop = FALSE]
  if (!nrow(df)) {
    return(data.frame())
  }
  reshape(df, idvar = c("target_cell_type", "factor", "factor_label"),
    timevar = "method", direction = "wide")
}

.celladmix_plot_admixture_pair <- function(admixture_agreement, method_a, method_b,
                                           p_thresh, cell_levels, factor_levels,
                                           max_neg_log10) {
  pair_df <- .celladmix_admixture_pair(admixture_agreement, method_a, method_b)
  if (is.null(pair_df) || !nrow(pair_df)) {
    return(.celladmix_empty_plot("Admixture significance agreement unavailable"))
  }
  logp_a <- paste0("logp.", method_a)
  logp_b <- paste0("logp.", method_b)
  call_a <- paste0("call.", method_a)
  call_b <- paste0("call.", method_b)
  for (nm in c(logp_a, logp_b)) {
    if (!(nm %in% names(pair_df))) pair_df[[nm]] <- NA_real_
  }
  for (nm in c(call_a, call_b)) {
    if (!(nm %in% names(pair_df))) pair_df[[nm]] <- FALSE
  }
  plot_df <- rbind(
    data.frame(factor = pair_df$factor, factor_label = pair_df$factor_label,
      method = method_a, triangle = "upper", cell_type = pair_df$target_cell_type,
      logp = pair_df[[logp_a]], call = pair_df[[call_a]], stringsAsFactors = FALSE),
    data.frame(factor = pair_df$factor, factor_label = pair_df$factor_label,
      method = method_b, triangle = "lower", cell_type = pair_df$target_cell_type,
      logp = pair_df[[logp_b]], call = pair_df[[call_b]], stringsAsFactors = FALSE)
  )
  plot_df <- plot_df[plot_df$cell_type %in% cell_levels, , drop = FALSE]
  plot_df$logp_plot <- pmax(0, pmin(plot_df$logp, max_neg_log10), na.rm = FALSE)
  factor_y_levels <- rev(factor_levels)
  plot_df$x <- match(plot_df$cell_type, cell_levels)
  plot_df$y <- match(plot_df$factor_label, factor_y_levels)
  triangles <- .celladmix_triangle_vertices(plot_df)

  ggplot2::ggplot() +
    ggplot2::geom_tile(data = expand.grid(x = seq_along(cell_levels),
      y = seq_along(factor_y_levels)), ggplot2::aes(x, y), fill = "white",
      color = "grey86", linewidth = 0.22) +
    ggplot2::geom_polygon(data = triangles,
      ggplot2::aes(x_vertex, y_vertex, group = polygon_id, fill = logp_plot),
      color = "white", linewidth = 0.18) +
    ggplot2::scale_fill_gradient(low = "white", high = "#B2182B",
      limits = c(0, max_neg_log10), na.value = "white") +
    ggplot2::scale_x_continuous(breaks = seq_along(cell_levels), labels = cell_levels,
      expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_along(factor_y_levels),
      labels = factor_y_levels, expand = c(0, 0)) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::labs(title = "Admixture significance by factor",
      subtitle = sprintf("Upper-right = %s; lower-left = %s", method_a, method_b),
      x = "Cell type", y = "Factor", fill = "-log10(p)") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 8, hjust = 0.5),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))
}

.celladmix_plot_score_agreement_pair <- function(source_agreement, admixture_agreement,
                                                 method_a, method_b, p_thresh,
                                                 max_neg_log10) {
  pair_source <- source_agreement[source_agreement$method %in% c(method_a, method_b), ,
    drop = FALSE]
  pair_admixture <- admixture_agreement[admixture_agreement$method %in% c(method_a, method_b),
    , drop = FALSE]
  cell_levels <- .celladmix_agreement_cell_levels(pair_source, pair_admixture)
  factor_levels <- .celladmix_agreement_factor_levels(pair_source, pair_admixture)
  if (!length(cell_levels) || !length(factor_levels)) {
    return(.celladmix_empty_plot("Score agreement unavailable"))
  }
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("The cowplot package is required for score-agreement plots")
  }
  cowplot::plot_grid(
    .celladmix_plot_source_pair(source_agreement, method_a, method_b,
      cell_levels, factor_levels),
    .celladmix_plot_admixture_pair(admixture_agreement, method_a, method_b,
      p_thresh, cell_levels, factor_levels, max_neg_log10),
    ncol = 2,
    align = "hv",
    axis = "tblr"
  )
}

#' Compare source and admixture calls across scoring methods.
#'
#' This diagnostic compares score annotations across two or more scoring
#' methods. Each method can be supplied as a [CellAdmixScore] object or as the
#' annotation list returned by `score$annotation()`. For each pair of methods,
#' the left panel compares inferred factor source cell types and the right panel
#' compares target-cell admixture significance.
#'
#' @param scores Named list of [CellAdmixScore] objects or score annotation
#'   lists. Names are used as method labels.
#' @param p_thresh P-value threshold used to call significant admixture.
#' @param adjust_p Passed to `score$annotation()` when `scores` contains
#'   [CellAdmixScore] objects.
#' @param method_pairs Optional list of two-element method-name vectors to plot.
#'   Defaults to all pairwise comparisons.
#' @param max_neg_log10 Cap for plotted `-log10(p)` values.
#' @param ... Additional arguments passed to `score$annotation()` for
#'   [CellAdmixScore] inputs.
#' @export
celladmix_plot_score_agreement <- function(
    scores,
    p_thresh = 0.1,
    adjust_p = FALSE,
    method_pairs = NULL,
    max_neg_log10 = 30,
    ...
  ) {
  .celladmix_require_ggplot2()
  if (!is.list(scores) || length(scores) < 2L) {
    stop("scores must be a named list with at least two scoring methods")
  }
  if (is.null(names(scores)) || any(!nzchar(names(scores)))) {
    names(scores) <- paste0("method_", seq_along(scores))
  }
  annotations <- lapply(scores, .celladmix_score_annotation,
    p_thresh = p_thresh, adjust_p = adjust_p, ...)
  source_agreement <- .celladmix_source_agreement(annotations)
  admixture_agreement <- .celladmix_admixture_agreement(annotations,
    p_thresh = p_thresh, max_neg_log10 = max_neg_log10)

  methods <- names(annotations)
  if (is.null(method_pairs)) {
    method_pairs <- utils::combn(methods, 2, simplify = FALSE)
  }
  plots <- lapply(method_pairs, function(pair) {
    if (length(pair) != 2L || !all(pair %in% methods)) {
      stop("Each method pair must contain two names from scores")
    }
    .celladmix_plot_score_agreement_pair(source_agreement, admixture_agreement,
      pair[[1]], pair[[2]], p_thresh = p_thresh, max_neg_log10 = max_neg_log10)
  })
  if (length(plots) == 1L) {
    return(plots[[1]])
  }
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("The cowplot package is required for multi-pair score-agreement plots")
  }
  cowplot::plot_grid(plotlist = plots, ncol = 1, align = "v", axis = "lr")
}

#' Run a simple two-group differential-expression diagnostic.
#'
#' This helper is meant for small interpretability checks in examples, not as a
#' replacement for a full single-cell DE workflow. Counts are normalized per
#' cell, compared with a Wilcoxon test per gene, and returned with standardized
#' column names.
#'
#' @param counts Gene-by-cell count matrix.
#' @param groups Named vector of group labels, or an unnamed vector in column order.
#' @param contrast Two group labels; logFC is `contrast[1]` over `contrast[2]`.
#' @param method Currently only `"wilcox"` is implemented.
#' @param normalize Count normalization used before testing.
#' @param scale_factor Library-size normalization scale factor.
#' @param pseudocount Pseudocount used only for logFC calculation.
#' @param p_adjust_method Multiple-testing correction method passed to [p.adjust()].
#' @export
celladmix_de <- function(
    counts,
    groups,
    contrast = NULL,
    method = c("wilcox"),
    normalize = c("cpm", "log1p_cpm", "none"),
    scale_factor = 10000,
    pseudocount = 0.1,
    p_adjust_method = "BH"
  ) {
  method <- match.arg(method)
  normalize <- match.arg(normalize)
  if (!length(rownames(counts))) {
    rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
  }
  if (!length(colnames(counts))) {
    stop("counts must have cell IDs as column names")
  }

  groups <- .celladmix_named_groups(groups, colnames(counts))
  common <- intersect(colnames(counts), names(groups))
  counts <- counts[, common, drop = FALSE]
  groups <- groups[common]
  keep <- !is.na(groups) & nzchar(groups)
  counts <- counts[, keep, drop = FALSE]
  groups <- groups[keep]
  groups_seen <- unique(groups)
  if (is.null(contrast)) {
    if (length(groups_seen) != 2L) {
      stop("contrast must be provided unless groups contains exactly two labels")
    }
    contrast <- groups_seen
  }
  if (length(contrast) != 2L) {
    stop("contrast must contain exactly two group labels")
  }
  missing_groups <- setdiff(contrast, unique(groups))
  if (length(missing_groups)) {
    stop("contrast groups absent from groups: ", paste(missing_groups, collapse = ", "))
  }

  norm <- .celladmix_normalize_counts_matrix(counts, normalize = normalize,
    scale_factor = scale_factor)
  group_a <- groups == contrast[[1]]
  group_b <- groups == contrast[[2]]
  mean_a <- .celladmix_matrix_row_means(norm[, group_a, drop = FALSE])
  mean_b <- .celladmix_matrix_row_means(norm[, group_b, drop = FALSE])
  logfc <- log2((mean_a + pseudocount) / (mean_b + pseudocount))
  pval <- vapply(seq_len(nrow(norm)), function(i) {
    xa <- as.numeric(norm[i, group_a])
    xb <- as.numeric(norm[i, group_b])
    if (!length(xa) || !length(xb) || length(unique(c(xa, xb))) <= 1L) {
      return(1)
    }
    tryCatch(stats::wilcox.test(xa, xb, exact = FALSE)$p.value,
      error = function(e) NA_real_)
  }, numeric(1))

  data.frame(
    gene = rownames(norm),
    group_a = contrast[[1]],
    group_b = contrast[[2]],
    logFC = as.numeric(logfc),
    mean_a = as.numeric(mean_a),
    mean_b = as.numeric(mean_b),
    mean_expr = as.numeric((mean_a + mean_b) / 2),
    p_value = pval,
    p_adj = stats::p.adjust(pval, method = p_adjust_method),
    stringsAsFactors = FALSE
  )
}

.celladmix_metadata_cell_ids <- function(cell_data, cell_id_col) {
  if (!is.data.frame(cell_data)) {
    stop("cell_data must be a data frame")
  }
  if (cell_id_col %in% names(cell_data)) {
    return(as.character(cell_data[[cell_id_col]]))
  }
  rn <- rownames(cell_data)
  if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(cell_data))))) {
    return(as.character(rn))
  }
  stop("cell_data must contain cell_id_col, or meaningful row names")
}

.celladmix_metadata_vector <- function(value, cell_data, cells, what) {
  if (is.character(value) && length(value) == 1L && value %in% names(cell_data)) {
    out <- as.character(cell_data[[value]])
    names(out) <- cells
    return(out)
  }
  if (is.atomic(value) && !is.null(names(value))) {
    out <- as.character(value)
    return(out)
  }
  stop(what, " must be a cell_data column name or a named vector")
}

.celladmix_de_by_group_summary <- function(results) {
  summary <- attr(results, "summary")
  if (is.null(summary)) {
    data.frame()
  } else {
    summary
  }
}

.celladmix_counts_cells <- function(counts) {
  cells <- colnames(counts)
  if (is.null(cells) || !length(cells)) {
    stop("counts must have cell IDs as column names")
  }
  as.character(cells)
}

.celladmix_sparse_cell_gene_matrix <- function(counts) {
  if (inherits(counts, "Matrix")) {
    Matrix::t(counts)
  } else {
    Matrix::Matrix(t(as.matrix(counts)), sparse = TRUE)
  }
}

.celladmix_append_sccore_metrics <- function(
    de,
    counts,
    groups,
    contrast,
    low_expression_threshold = 0,
    append_auc = FALSE,
    strict = FALSE
  ) {
  if (!requireNamespace("sccore", quietly = TRUE)) {
    if (isTRUE(strict)) {
      stop("method = 'sccore' requires the sccore package")
    }
    return(list(de = de, used = FALSE, error = "sccore unavailable"))
  }
  if (isTRUE(append_auc) && !requireNamespace("pROC", quietly = TRUE)) {
    if (isTRUE(strict)) {
      stop("sccore_auc = TRUE requires the pROC package")
    }
    append_auc <- FALSE
  }
  common_cells <- intersect(colnames(counts), names(groups))
  counts <- counts[, common_cells, drop = FALSE]
  groups <- groups[common_cells]
  keep <- !is.na(groups) & nzchar(groups)
  counts <- counts[, keep, drop = FALSE]
  groups <- groups[keep]
  genes <- intersect(de$gene, rownames(counts))
  if (!length(genes)) {
    return(list(de = de, used = FALSE, error = "no common genes"))
  }
  p2_counts <- .celladmix_sparse_cell_gene_matrix(counts[genes, , drop = FALSE])
  tmp <- data.frame(Gene = genes, stringsAsFactors = FALSE)
  appended <- tryCatch(
    sccore::appendSpecificityMetricsToDE(
      tmp,
      clusters = groups,
      cluster.id = contrast[[1]],
      p2.counts = p2_counts,
      low.expression.threshold = low_expression_threshold,
      append.auc = append_auc
    ),
    error = function(e) e
  )
  if (inherits(appended, "error")) {
    if (isTRUE(strict)) {
      stop("sccore specificity metrics failed: ", conditionMessage(appended))
    }
    return(list(de = de, used = FALSE, error = conditionMessage(appended)))
  }
  idx <- match(de$gene, appended$Gene)
  if ("Specificity" %in% names(appended)) {
    de$sccore_specificity <- as.numeric(appended$Specificity[idx])
  }
  if ("Precision" %in% names(appended)) {
    de$sccore_precision <- as.numeric(appended$Precision[idx])
  }
  if ("ExpressionFraction" %in% names(appended)) {
    de$sccore_expression_fraction <- as.numeric(appended$ExpressionFraction[idx])
  }
  if ("AUC" %in% names(appended)) {
    de$sccore_auc <- as.numeric(appended$AUC[idx])
  }
  list(de = de, used = TRUE, error = NA_character_)
}

#' Run DE Diagnostics Within Cell Subsets
#'
#' `celladmix_de_by_group()` is a convenience layer around [celladmix_de()] for
#' common notebook diagnostics: compare two context groups within one or more
#' cell subsets. For example, in a tissue-domain analysis, one can compare
#' tumor-domain versus adjacent-normal-domain endothelial cells and immune cells
#' with a single call.
#'
#' This function is intentionally not a replacement for a full single-cell DE
#' workflow. The statistical test is the lightweight [celladmix_de()] Wilcoxon
#' path. When `method = "auto"` and `sccore` is installed, sccore specificity
#' metrics are appended to the returned DE tables; when `method = "sccore"`,
#' sccore metric annotation is required and failures are treated as errors.
#' The p-values and logFC values still come from [celladmix_de()] so results are
#' always returned in the same standardized format used by cellAdmix plotting
#' helpers.
#'
#' @param counts Gene-by-cell count matrix. Column names must be cell IDs.
#' @param cell_data Cell-level metadata data frame.
#' @param group_by Cell metadata column name, or named vector, defining the two
#'   groups to compare within each subset.
#' @param contrast Two group labels. `logFC` is `contrast[1] / contrast[2]`.
#'   If `NULL`, each subset must contain exactly two group labels.
#' @param subset_by Optional cell metadata column name, or named vector, used to
#'   split cells before DE. If `NULL`, all cells are analyzed as one subset.
#' @param subsets Optional subset labels to analyze. If `NULL`, all non-missing
#'   subset labels are used.
#' @param cell_id_col Cell ID column in `cell_data`. If absent, meaningful row
#'   names are used.
#' @param method `"wilcox"` runs only [celladmix_de()]. `"auto"` uses the same
#'   DE test and appends sccore specificity metrics when `sccore` is available.
#'   `"sccore"` requires successful sccore metric annotation.
#' @param min_cells_per_group Minimum number of cells required in each contrast
#'   group within a subset.
#' @param normalize,scale_factor,pseudocount,p_adjust_method Passed to
#'   [celladmix_de()].
#' @param sccore_low_expression_threshold Threshold used by
#'   `sccore::appendSpecificityMetricsToDE()` when sccore metrics are appended.
#' @param sccore_auc Whether to request sccore AUC values. This also requires
#'   the optional `pROC` package.
#'
#' @return A named list of DE tables, one per subset. The object has class
#'   `celladmix_de_grouped`; call `summary()` on it to see subset sizes,
#'   backend usage, and skipped subsets.
#' @export
celladmix_de_by_group <- function(
    counts,
    cell_data,
    group_by,
    contrast = NULL,
    subset_by = NULL,
    subsets = NULL,
    cell_id_col = "cell_id",
    method = c("auto", "wilcox", "sccore"),
    min_cells_per_group = 3,
    normalize = c("cpm", "log1p_cpm", "none"),
    scale_factor = 10000,
    pseudocount = 0.1,
    p_adjust_method = "BH",
    sccore_low_expression_threshold = 0,
    sccore_auc = FALSE
  ) {
  method <- match.arg(method)
  normalize <- match.arg(normalize)
  count_cells <- .celladmix_counts_cells(counts)
  cells <- .celladmix_metadata_cell_ids(cell_data, cell_id_col)
  group_values <- .celladmix_metadata_vector(group_by, cell_data, cells, "group_by")
  subset_values <- if (is.null(subset_by)) {
    stats::setNames(rep("all", length(cells)), cells)
  } else {
    .celladmix_metadata_vector(subset_by, cell_data, cells, "subset_by")
  }

  common <- Reduce(intersect, list(count_cells, names(group_values), names(subset_values)))
  if (!length(common)) {
    stop("No cells overlap between counts, group_by, and subset_by")
  }
  group_values <- group_values[common]
  subset_values <- subset_values[common]
  keep <- !is.na(group_values) & nzchar(group_values) &
    !is.na(subset_values) & nzchar(subset_values)
  group_values <- group_values[keep]
  subset_values <- subset_values[keep]
  common <- names(group_values)

  if (is.null(subsets)) {
    subsets <- unique(subset_values)
  }
  subsets <- as.character(subsets)
  subset_name <- if (is.null(subset_by)) "all" else if (is.character(subset_by) && length(subset_by) == 1L) subset_by else "subset"
  group_name <- if (is.character(group_by) && length(group_by) == 1L) group_by else "group"

  results <- list()
  summary_rows <- vector("list", length(subsets))
  for (i in seq_along(subsets)) {
    subset_label <- subsets[[i]]
    subset_cells <- names(subset_values)[subset_values == subset_label]
    subset_cells <- intersect(subset_cells, common)
    groups <- group_values[subset_cells]
    groups <- groups[!is.na(groups) & nzchar(groups)]
    subset_cells <- names(groups)
    subset_contrast <- contrast
    status <- "ok"
    sccore_used <- FALSE
    sccore_error <- NA_character_
    de <- NULL

    if (is.null(subset_contrast)) {
      present <- unique(groups)
      if (length(present) == 2L) {
        subset_contrast <- present
      } else {
        status <- "skipped_requires_two_groups"
      }
    }
    counts_by_group <- if (!is.null(subset_contrast)) {
      stats::setNames(vapply(subset_contrast, function(g) sum(groups == g),
        integer(1)), subset_contrast)
    } else {
      integer()
    }
    if (identical(status, "ok")) {
      if (!all(subset_contrast %in% groups)) {
        status <- "skipped_missing_contrast_group"
      } else if (any(counts_by_group < min_cells_per_group)) {
        status <- "skipped_too_few_cells"
      }
    }
    if (identical(status, "ok")) {
      analysis_cells <- names(groups)[groups %in% subset_contrast]
      analysis_groups <- groups[analysis_cells]
      de <- celladmix_de(
        counts[, analysis_cells, drop = FALSE],
        analysis_groups,
        contrast = subset_contrast,
        method = "wilcox",
        normalize = normalize,
        scale_factor = scale_factor,
        pseudocount = pseudocount,
        p_adjust_method = p_adjust_method
      )
      de$subset <- subset_label
      de$subset_by <- subset_name
      de$group_by <- group_name
      de$n_group_a <- unname(counts_by_group[[subset_contrast[[1]]]])
      de$n_group_b <- unname(counts_by_group[[subset_contrast[[2]]]])
      if (method %in% c("auto", "sccore")) {
        sccore <- .celladmix_append_sccore_metrics(
          de,
          counts[, analysis_cells, drop = FALSE],
          analysis_groups,
          subset_contrast,
          low_expression_threshold = sccore_low_expression_threshold,
          append_auc = sccore_auc,
          strict = identical(method, "sccore")
        )
        de <- sccore$de
        sccore_used <- isTRUE(sccore$used)
        sccore_error <- sccore$error
      }
      results[[subset_label]] <- de
    }
    summary_rows[[i]] <- data.frame(
      subset = subset_label,
      subset_by = subset_name,
      group_by = group_name,
      group_a = if (!is.null(subset_contrast)) subset_contrast[[1]] else NA_character_,
      group_b = if (!is.null(subset_contrast)) subset_contrast[[2]] else NA_character_,
      n_group_a = if (length(counts_by_group)) unname(counts_by_group[[subset_contrast[[1]]]]) else NA_integer_,
      n_group_b = if (length(counts_by_group)) unname(counts_by_group[[subset_contrast[[2]]]]) else NA_integer_,
      n_cells = length(subset_cells),
      n_genes = if (!is.null(de)) nrow(de) else 0L,
      method = "wilcox",
      sccore_metrics = sccore_used,
      sccore_error = sccore_error,
      status = status,
      stringsAsFactors = FALSE
    )
  }
  summary <- do.call(rbind, summary_rows)
  attr(results, "summary") <- summary
  class(results) <- c("celladmix_de_grouped", class(results))
  results
}

#' @export
summary.celladmix_de_grouped <- function(object, ...) {
  .celladmix_de_by_group_summary(object)
}

#' Plot a compact volcano diagnostic from a cellAdmix DE table.
#'
#' @param de Output from [celladmix_de()] or a compatible DE table.
#' @param markers Genes to highlight, supplied as a character vector, named list,
#'   or data.frame with `gene` and `set`/`marker_set`/`label` columns.
#' @param label Label policy for text annotations. `TRUE` labels all highlighted
#'   markers, `FALSE` labels none, and a character vector labels that subset.
#' @param marker_label Legend label used when `markers` is a character vector.
#' @param p_cap Cap for `-log10(adjusted p)`.
#' @param title Optional plot title.
#' @param subtitle Optional plot subtitle.
#' @export
celladmix_plot_volcano <- function(
    de,
    markers = NULL,
    label = TRUE,
    marker_label = "Marker",
    p_cap = 30,
    title = NULL,
    subtitle = NULL
  ) {
  .celladmix_require_ggplot2()
  de <- .celladmix_standardize_de(de)
  if (is.null(de) || !nrow(de)) {
    return(.celladmix_empty_plot("DE table unavailable"))
  }
  de$logp <- pmin(-log10(pmax(de$p_adj, .Machine$double.xmin)), p_cap)
  marker_df <- .celladmix_marker_frame(markers, default_set = marker_label)
  marker_df <- marker_df[!is.na(marker_df$gene) & nzchar(marker_df$gene), , drop = FALSE]
  marker_df <- marker_df[!duplicated(marker_df$gene), , drop = FALSE]
  marker_genes <- unique(marker_df$gene)
  if (nrow(marker_df)) {
    set_by_gene <- stats::setNames(marker_df$set, marker_df$gene)
    de$marker <- unname(set_by_gene[de$gene])
    de$marker[is.na(de$marker) | !nzchar(de$marker)] <- "Other gene"
    marker_levels <- c(unique(marker_df$set), "Other gene")
    de$marker <- factor(de$marker, levels = marker_levels)
    cols <- .celladmix_marker_palette(marker_levels)
  } else {
    de$marker <- factor("Other gene", levels = "Other gene")
    cols <- c("Other gene" = "grey70")
  }
  label_genes <- .celladmix_marker_labels(label, marker_genes)
  de$label <- ifelse(de$gene %in% label_genes, de$gene, NA_character_)
  label_df <- de[!is.na(de$label), , drop = FALSE]

  ymax <- max(de$logp, na.rm = TRUE)
  if (!is.finite(ymax)) {
    ymax <- 1
  }
  ymax <- ymax + 1

  p <- ggplot2::ggplot(de, ggplot2::aes(x = logFC, y = logp, color = marker, label = label)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, color = "grey45") +
    ggplot2::geom_vline(xintercept = c(-log2(1.5), log2(1.5)),
      linewidth = 0.3, linetype = "dashed", color = "red3") +
    ggplot2::geom_hline(yintercept = -log10(0.05),
      linewidth = 0.3, linetype = "dashed", color = "red3") +
    ggplot2::geom_point(data = de[de$marker == "Other gene", , drop = FALSE],
      alpha = 0.45, size = 1.1) +
    ggplot2::geom_point(data = de[de$marker != "Other gene", , drop = FALSE],
      alpha = 0.75, size = 1.5) +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE) +
    ggplot2::coord_cartesian(ylim = c(0, ymax)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "log2 fold change",
      y = "-log10 adjusted p", color = NULL) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "bottom")
  if (nrow(label_df) && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(data = label_df, size = 3.5,
      show.legend = FALSE, na.rm = TRUE, max.overlaps = Inf)
  } else if (nrow(label_df)) {
    p <- p + ggplot2::geom_text(data = label_df, size = 3,
      vjust = -0.4, show.legend = FALSE, na.rm = TRUE)
  }
  p
}

#' Plot mean marker expression before and optionally after correction.
#'
#' @param counts_before Gene-by-cell count matrix before correction.
#' @param counts_after Optional corrected count matrix.
#' @param cells Optional cell IDs to include.
#' @param markers Marker genes supplied as a character vector, named list, or
#'   data.frame with `gene` and `set`/`marker_set`/`label` columns.
#' @param normalize Count normalization used before averaging.
#' @param scale_factor Library-size normalization scale factor.
#' @param state_labels Labels for before/after states.
#' @export
celladmix_plot_marker_expression <- function(
    counts_before,
    counts_after = NULL,
    cells = NULL,
    markers,
    normalize = c("cpm", "log1p_cpm", "none"),
    scale_factor = 10000,
    state_labels = c("before", "after")
  ) {
  .celladmix_require_ggplot2()
  normalize <- match.arg(normalize)
  marker_df <- .celladmix_marker_frame(markers)
  if (!nrow(marker_df)) {
    return(.celladmix_empty_plot("No marker genes supplied"))
  }
  if (is.null(cells)) {
    cells <- colnames(counts_before)
  }
  cells <- intersect(cells, colnames(counts_before))
  if (!is.null(counts_after)) {
    cells <- intersect(cells, colnames(counts_after))
  }
  if (!length(cells)) {
    return(.celladmix_empty_plot("No cells available"))
  }
  marker_genes <- unique(marker_df$gene)
  marker_genes <- intersect(marker_genes, rownames(counts_before))
  if (!is.null(counts_after)) {
    marker_genes <- intersect(marker_genes, rownames(counts_after))
  }
  if (!length(marker_genes)) {
    return(.celladmix_empty_plot("No marker genes present"))
  }
  marker_df <- marker_df[marker_df$gene %in% marker_genes, , drop = FALSE]

  summarize_state <- function(counts, state) {
    norm <- .celladmix_normalize_counts_matrix(counts[marker_genes, cells, drop = FALSE],
      normalize = normalize, scale_factor = scale_factor)
    expr <- .celladmix_matrix_row_means(norm)
    data.frame(gene = names(expr), state = state, expression = as.numeric(expr),
      stringsAsFactors = FALSE)
  }
  states <- list(summarize_state(counts_before, state_labels[[1]]))
  if (!is.null(counts_after)) {
    states[[2]] <- summarize_state(counts_after, state_labels[[2]])
  }
  df <- merge(do.call(rbind, states), marker_df, by = "gene", all.x = TRUE)
  df$gene <- factor(df$gene, levels = unique(marker_df$gene))
  df$state <- factor(df$state, levels = state_labels[seq_along(unique(df$state))])
  fill_values <- stats::setNames(c("grey65", "#2C7FB8")[seq_along(levels(df$state))],
    levels(df$state))

  ggplot2::ggplot(df, ggplot2::aes(x = gene, y = expression, fill = state)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72), width = 0.64) +
    ggplot2::facet_wrap(~ set, scales = "free_x") +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Mean normalized expression") +
    ggplot2::scale_fill_manual(values = fill_values) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.title = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"))
}

#' Compare signed DE significance before and after correction.
#'
#' @param de_before DE table before correction.
#' @param de_after DE table after correction.
#' @param markers Genes to highlight, supplied as a character vector, named list,
#'   or data.frame with `gene` and `set`/`marker_set`/`label` columns.
#' @param label Label policy for text annotations. `TRUE` labels all highlighted
#'   markers, `FALSE` labels none, and a character vector labels that subset.
#' @param marker_label Legend label used when `markers` is a character vector.
#' @param signed If TRUE, multiply `-log10(p_adj)` by the sign of logFC.
#' @param p_cap Cap for plotted `-log10(p_adj)` values.
#' @param title Optional plot title.
#' @export
celladmix_plot_de_shift <- function(
    de_before,
    de_after,
    markers = NULL,
    label = TRUE,
    marker_label = "Marker",
    signed = TRUE,
    p_cap = 30,
    title = "DE significance before/after cleanup"
  ) {
  .celladmix_require_ggplot2()
  before <- .celladmix_standardize_de(de_before)
  after <- .celladmix_standardize_de(de_after)
  if (is.null(before) || is.null(after) || !nrow(before) || !nrow(after)) {
    return(.celladmix_empty_plot("DE p-value comparison unavailable"))
  }
  common <- intersect(before$gene, after$gene)
  if (!length(common)) {
    return(.celladmix_empty_plot("DE p-value comparison unavailable"))
  }
  before <- before[match(common, before$gene), , drop = FALSE]
  after <- after[match(common, after$gene), , drop = FALSE]
  transform_logp <- function(x) {
    y <- pmin(-log10(pmax(x$p_adj, .Machine$double.xmin)), p_cap)
    if (isTRUE(signed)) {
      y <- sign(x$logFC) * y
    }
    y
  }
  marker_df <- .celladmix_marker_frame(markers, default_set = marker_label)
  marker_df <- marker_df[!is.na(marker_df$gene) & nzchar(marker_df$gene), , drop = FALSE]
  marker_df <- marker_df[!duplicated(marker_df$gene), , drop = FALSE]
  marker_genes <- unique(marker_df$gene)
  if (nrow(marker_df)) {
    set_by_gene <- stats::setNames(marker_df$set, marker_df$gene)
    marker <- unname(set_by_gene[common])
    marker[is.na(marker) | !nzchar(marker)] <- "Other gene"
    marker_levels <- c(unique(marker_df$set), "Other gene")
  } else {
    marker <- rep("Other gene", length(common))
    marker_levels <- "Other gene"
  }
  label_genes <- .celladmix_marker_labels(label, marker_genes)
  df <- data.frame(
    gene = common,
    before = transform_logp(before),
    after = transform_logp(after),
    marker = marker,
    label = ifelse(common %in% label_genes, common, NA_character_),
    stringsAsFactors = FALSE
  )
  df$marker <- factor(df$marker, levels = marker_levels)
  label_df <- df[!is.na(df$label), , drop = FALSE]
  limit <- if (isTRUE(signed)) {
    max(abs(c(df$before, df$after)), na.rm = TRUE)
  } else {
    max(c(df$before, df$after), na.rm = TRUE)
  }
  if (!is.finite(limit) || limit <= 0) {
    limit <- 1
  }
  limit <- min(p_cap, ceiling(limit))
  lims <- if (isTRUE(signed)) c(-limit, limit) else c(0, limit)

  p <- ggplot2::ggplot(df, ggplot2::aes(before, after, color = marker, label = label)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
      color = "grey55", linewidth = 0.35) +
    ggplot2::geom_point(data = df[df$marker == "Other gene", , drop = FALSE],
      alpha = 0.35, size = 1.0) +
    ggplot2::geom_point(data = df[df$marker != "Other gene", , drop = FALSE],
      alpha = 0.85, size = 1.8) +
    ggplot2::scale_colour_manual(values = .celladmix_marker_palette(marker_levels),
      drop = FALSE) +
    ggplot2::coord_equal(xlim = lims, ylim = lims, expand = TRUE) +
    ggplot2::labs(title = title,
      x = if (signed) "Before cleanup: signed -log10 adjusted p" else "Before cleanup: -log10 adjusted p",
      y = if (signed) "After cleanup: signed -log10 adjusted p" else "After cleanup: -log10 adjusted p",
      color = NULL) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5),
      legend.position = "bottom")
  if (isTRUE(signed)) {
    p <- p +
      ggplot2::geom_hline(yintercept = 0, color = "grey45", linewidth = 0.25) +
      ggplot2::geom_vline(xintercept = 0, color = "grey45", linewidth = 0.25)
  }
  if (nrow(label_df) && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(data = label_df, size = 3.2,
      show.legend = FALSE, na.rm = TRUE, max.overlaps = Inf)
  } else if (nrow(label_df)) {
    p <- p + ggplot2::geom_text(data = label_df, size = 2.8, vjust = -0.4,
      show.legend = FALSE, na.rm = TRUE)
  }
  p
}

#' Collect corrected counts and recompute a DE diagnostic.
#'
#' @param correction A `CellAdmixCorrection` object.
#' @param counts_before Original gene-by-cell count matrix.
#' @param cells Cell IDs to include in the diagnostic.
#' @param groups Named group vector for the selected cells.
#' @param contrast Optional two-group contrast passed to [celladmix_de()].
#' @param ... Additional arguments passed to [celladmix_de()].
#' @export
celladmix_correction_impact <- function(
    correction,
    counts_before,
    cells,
    groups,
    contrast = NULL,
    ...
  ) {
  if (is.null(correction)) {
    return(NULL)
  }
  if (!inherits(correction, "CellAdmixCorrection")) {
    stop("correction must be a CellAdmixCorrection object")
  }
  counts_after <- correction$counts()
  cells <- intersect(intersect(cells, colnames(counts_before)), colnames(counts_after))
  groups <- .celladmix_named_groups(groups, cells)
  groups <- groups[cells]
  groups <- groups[!is.na(groups)]
  cells <- names(groups)
  de <- NULL
  if (length(unique(groups)) == 2L || (!is.null(contrast) && all(contrast %in% groups))) {
    de <- celladmix_de(counts_after[, cells, drop = FALSE], groups,
      contrast = contrast, ...)
  }
  list(counts = counts_after, cells = cells, groups = groups, de = de)
}
