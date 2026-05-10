# Cluster marker diagnostics used by annotation reports.

#' Score Cluster Marker Genes by Sparse Detection AUC
#'
#' Computes simple one-vs-rest marker statistics for cell clusters from a sparse
#' gene-by-cell count matrix. The AUC-like score is based on the difference in
#' detection fraction inside and outside each cluster, which is cheap to compute
#' on sparse count matrices and works well for quick annotation reports.
#'
#' @param counts Sparse gene-by-cell count matrix.
#' @param clusters Named vector of cluster labels, with names matching cell IDs.
#' @param normalization_scale Per-cell library-size normalization scale.
#' @param min_expression_fraction Minimum within-cluster detection fraction.
#' @param min_logfc Minimum normalized-expression log2 fold-change.
#'
#' @return A data.frame with one row per cluster-gene pair passing filters.
#' @export
celladmix_score_cluster_markers <- function(
    counts,
    clusters,
    normalization_scale = 10000,
    min_expression_fraction = 0.05,
    min_logfc = 0
  ) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required for cluster marker scoring")
  }
  if (!inherits(counts, "dgCMatrix")) {
    counts <- methods::as(counts, "dgCMatrix")
  }
  if (is.null(names(clusters))) {
    stop("clusters must be a named vector with cell IDs as names")
  }
  common_cells <- intersect(colnames(counts), names(clusters))
  if (!length(common_cells)) {
    stop("No common cells between counts and clusters")
  }
  counts <- counts[, common_cells, drop = FALSE]
  clusters <- factor(clusters[common_cells])

  lib_size <- Matrix::colSums(counts)
  scale <- normalization_scale / pmax(lib_size, 1)
  norm <- counts
  norm@x <- norm@x * rep.int(scale, diff(norm@p))

  detected <- counts
  detected@x <- rep(1, length(detected@x))

  total_detected <- Matrix::rowSums(detected)
  total_expression <- Matrix::rowSums(norm)
  n_cells <- ncol(counts)

  out <- lapply(levels(clusters), function(cluster) {
    in_cluster <- clusters == cluster
    n_in <- sum(in_cluster)
    n_out <- n_cells - n_in
    detected_in <- Matrix::rowSums(detected[, in_cluster, drop = FALSE])
    expression_in <- Matrix::rowSums(norm[, in_cluster, drop = FALSE])

    fraction_in <- detected_in / pmax(n_in, 1)
    fraction_out <- (total_detected - detected_in) / pmax(n_out, 1)
    auc <- 0.5 + 0.5 * (fraction_in - fraction_out)

    mean_in <- expression_in / pmax(n_in, 1)
    mean_out <- (total_expression - expression_in) / pmax(n_out, 1)
    logfc <- log2((mean_in + 1e-6) / (mean_out + 1e-6))
    precision <- detected_in / pmax(total_detected, 1)
    specificity <- (n_cells - total_detected) / pmax(n_cells - detected_in, 1)

    data.frame(cluster = cluster, gene = rownames(counts), auc = as.numeric(auc),
      logfc = as.numeric(logfc), expression_fraction = as.numeric(fraction_in),
      outside_expression_fraction = as.numeric(fraction_out),
      precision = as.numeric(precision), specificity = as.numeric(specificity),
      mean_norm_in = as.numeric(mean_in), mean_norm_out = as.numeric(mean_out),
      stringsAsFactors = FALSE)
  })

  markers <- do.call(rbind, out)
  markers <- markers[markers$expression_fraction >= min_expression_fraction &
    markers$logfc > min_logfc, , drop = FALSE]
  markers <- markers[order(markers$cluster, -markers$auc, -markers$logfc), , drop = FALSE]
  rownames(markers) <- NULL
  markers
}

#' Select Top Cluster Markers
#'
#' @param markers Output of [celladmix_score_cluster_markers()].
#' @param n Number of marker genes to keep per cluster.
#'
#' @return A data.frame with up to `n` marker rows per cluster.
#' @export
celladmix_top_cluster_markers <- function(markers, n = 5) {
  if (!nrow(markers)) {
    return(markers)
  }
  cluster_levels <- unique(as.character(markers$cluster))
  numeric_levels <- suppressWarnings(as.integer(cluster_levels))
  if (all(!is.na(numeric_levels))) {
    cluster_levels <- as.character(sort(numeric_levels))
  }
  pieces <- split(markers, factor(markers$cluster, levels = cluster_levels))
  out <- do.call(rbind, lapply(pieces, utils::head, n))
  rownames(out) <- NULL
  out
}

.celladmix_cluster_marker_genes <- function(markers, top_markers = NULL, n = 5) {
  if (!is.null(top_markers) && nrow(top_markers)) {
    return(unique(top_markers$gene))
  }
  unique(celladmix_top_cluster_markers(markers, n = n)$gene)
}

.celladmix_cluster_marker_levels <- function(markers) {
  levels <- unique(as.character(markers$cluster))
  numeric_levels <- suppressWarnings(as.integer(levels))
  if (all(!is.na(numeric_levels))) {
    levels <- as.character(sort(numeric_levels))
  }
  levels
}

#' Plot Cluster Marker AUC Heatmap
#'
#' @param markers Output of [celladmix_score_cluster_markers()].
#' @param top_markers Optional output of [celladmix_top_cluster_markers()].
#' @param n Top markers per cluster if `top_markers` is not supplied.
#'
#' @return A `ggplot2` object.
#' @export
celladmix_plot_cluster_marker_heatmap <- function(markers, top_markers = NULL, n = 5) {
  .celladmix_require_ggplot2()
  if (!nrow(markers)) {
    return(.celladmix_empty_plot("No cluster markers available"))
  }
  plot_markers <- .celladmix_cluster_marker_genes(markers, top_markers, n)
  heat_df <- markers[markers$gene %in% plot_markers, , drop = FALSE]
  heat_df$cluster <- factor(heat_df$cluster,
    levels = .celladmix_cluster_marker_levels(markers))
  heat_df$gene <- factor(heat_df$gene, levels = rev(plot_markers))

  ggplot2::ggplot(heat_df, ggplot2::aes(cluster, gene, fill = auc)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "grey95", high = "#B2182B",
      limits = c(0.5, 1)) +
    ggplot2::labs(title = "Detection AUC for top cluster markers",
      x = "Cluster", y = NULL, fill = "AUC") +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))
}

#' Plot Cluster Marker Dotplot
#'
#' @param markers Output of [celladmix_score_cluster_markers()].
#' @param top_markers Optional output of [celladmix_top_cluster_markers()].
#' @param n Top markers per cluster if `top_markers` is not supplied.
#'
#' @return A `ggplot2` object.
#' @export
celladmix_plot_cluster_marker_dotplot <- function(markers, top_markers = NULL, n = 5) {
  .celladmix_require_ggplot2()
  if (!nrow(markers)) {
    return(.celladmix_empty_plot("No cluster markers available"))
  }
  plot_markers <- .celladmix_cluster_marker_genes(markers, top_markers, n)
  dot_df <- markers[markers$gene %in% plot_markers, , drop = FALSE]
  dot_df$cluster <- factor(dot_df$cluster,
    levels = .celladmix_cluster_marker_levels(markers))
  dot_df$gene <- factor(dot_df$gene, levels = rev(plot_markers))

  ggplot2::ggplot(dot_df, ggplot2::aes(cluster, gene)) +
    ggplot2::geom_point(ggplot2::aes(size = expression_fraction, color = logfc),
      alpha = 0.85) +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "grey90",
      high = "#B2182B", midpoint = 0) +
    ggplot2::scale_size(range = c(0.2, 3.2), limits = c(0, 1)) +
    ggplot2::labs(title = "Top marker expression fraction and logFC",
      x = "Cluster", y = NULL, color = "logFC", size = "fraction") +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))
}

.celladmix_normalized_marker_expression <- function(counts, genes, normalization_scale) {
  if (!inherits(counts, "dgCMatrix")) {
    counts <- methods::as(counts, "dgCMatrix")
  }
  genes <- intersect(genes, rownames(counts))
  if (!length(genes)) {
    return(NULL)
  }
  lib_size <- Matrix::colSums(counts)
  scale <- normalization_scale / pmax(lib_size, 1)
  expr <- counts[genes, , drop = FALSE]
  expr@x <- log1p(expr@x * rep.int(scale, diff(expr@p)))
  as.matrix(expr)
}

#' Plot Top Marker Expression on a Cluster UMAP
#'
#' @param embedding Data.frame with `cell_id`, `umap_1`, and `umap_2` columns.
#' @param counts Sparse gene-by-cell count matrix.
#' @param top_markers Output of [celladmix_top_cluster_markers()].
#' @param normalization_scale Per-cell library-size normalization scale.
#' @param max_panels Maximum number of cluster marker panels to draw.
#' @param max_cells Maximum number of cells to draw. Large UMAP marker panels
#'   are subsampled deterministically because faceted full-tissue plots can
#'   otherwise create tens of millions of points.
#' @param seed Random seed used when `max_cells` triggers subsampling.
#' @param point_size UMAP point size.
#'
#' @return A `ggplot2` object.
#' @export
celladmix_plot_cluster_marker_umap <- function(
    embedding,
    counts,
    top_markers,
    normalization_scale = 10000,
    max_panels = 25,
    max_cells = 100000,
    seed = 1,
    point_size = 0.08
  ) {
  .celladmix_require_ggplot2()
  if (is.null(embedding) || !is.data.frame(embedding) || !nrow(embedding)) {
    return(.celladmix_empty_plot("No cluster UMAP available"))
  }
  required <- c("cell_id", "umap_1", "umap_2")
  if (!all(required %in% colnames(embedding))) {
    stop("embedding must contain: ", paste(required, collapse = ", "))
  }
  if (is.null(top_markers) || !nrow(top_markers)) {
    return(.celladmix_empty_plot("No top cluster markers available"))
  }
  top_gene <- top_markers[!duplicated(top_markers$cluster),
    c("cluster", "gene", "auc"), drop = FALSE]
  top_gene$cluster <- as.character(top_gene$cluster)
  numeric_clusters <- suppressWarnings(as.integer(top_gene$cluster))
  if (all(!is.na(numeric_clusters))) {
    top_gene <- top_gene[order(numeric_clusters), , drop = FALSE]
  }
  if (nrow(top_gene) > max_panels) {
    top_gene <- top_gene[seq_len(max_panels), , drop = FALSE]
  }

  common_cells <- intersect(colnames(counts), embedding$cell_id)
  if (!length(common_cells)) {
    return(.celladmix_empty_plot("No common cells between counts and UMAP"))
  }
  if (is.finite(max_cells) && max_cells > 0 && length(common_cells) > max_cells) {
    set.seed(seed)
    common_cells <- sample(common_cells, max_cells)
  }
  counts <- counts[, common_cells, drop = FALSE]
  embedding <- embedding[match(common_cells, embedding$cell_id), , drop = FALSE]

  expr <- .celladmix_normalized_marker_expression(counts, unique(top_gene$gene),
    normalization_scale)
  if (is.null(expr)) {
    return(.celladmix_empty_plot("Top marker genes are not present"))
  }

  frames <- lapply(seq_len(nrow(top_gene)), function(i) {
    gene <- top_gene$gene[[i]]
    values <- as.numeric(expr[gene, ])
    positive <- values[values > 0]
    cap <- if (length(positive)) {
      stats::quantile(positive, probs = 0.995, na.rm = TRUE, names = FALSE)
    } else {
      0
    }
    cap <- if (is.finite(cap) && cap > 0) cap else max(values, na.rm = TRUE)
    plot_value <- if (is.finite(cap) && cap > 0) pmin(values / cap, 1) else values
    data.frame(cell_id = common_cells, umap_1 = embedding$umap_1,
      umap_2 = embedding$umap_2, cluster = top_gene$cluster[[i]],
      gene = gene, auc = top_gene$auc[[i]], expression = values,
      expression_scaled = plot_value,
      panel = paste0("C", top_gene$cluster[[i]], ": ", gene),
      stringsAsFactors = FALSE)
  })
  plot_df <- do.call(rbind, frames)
  plot_df <- plot_df[order(plot_df$panel, plot_df$expression), , drop = FALSE]
  plot_df$panel <- factor(plot_df$panel, levels = unique(plot_df$panel))

  ggplot2::ggplot(plot_df, ggplot2::aes(umap_1, umap_2, color = expression_scaled)) +
    ggplot2::geom_point(size = point_size, alpha = 0.85) +
    ggplot2::facet_wrap(~ panel) +
    ggplot2::coord_equal() +
    ggplot2::scale_colour_gradientn(colors = c("#d9d9d9", "#b9d9ea",
      "#3c8dbc", "#08306b", "#67000d"), limits = c(0, 1),
      name = "scaled\nexpression") +
    ggplot2::labs(title = "Top marker expression by cluster",
      x = "Cell UMAP 1", y = "Cell UMAP 2") +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 7),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "right")
}
