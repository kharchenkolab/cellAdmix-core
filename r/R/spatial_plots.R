#' Generic Spatial Cell Plotting Helpers
#'
#' These helpers operate on plain cell-level data frames. They intentionally do
#' not depend on a specific input backend, so the same plotting call can be used
#' for tabular examples, Xenium stores, fitted cell-factor tables, and external
#' annotations merged by the user.

.celladmix_auto_point_size <- function(n) {
  n <- as.integer(n %||% 0L)
  if (n <= 5000L) {
    return(0.35)
  }
  if (n <= 50000L) {
    return(0.12)
  }
  if (n <= 250000L) {
    return(0.045)
  }
  0.02
}

.celladmix_is_integer_like <- function(x) {
  x <- x[!is.na(x)]
  length(x) > 0L && is.numeric(x) && all(is.finite(x)) && all(abs(x - round(x)) < 1e-8)
}

.celladmix_is_discrete_color <- function(x, color_by = NULL) {
  if (is.factor(x) || is.character(x) || is.logical(x)) {
    return(TRUE)
  }
  if (!is.null(color_by) &&
      grepl("(^|_)(factor|cluster|label|type|annotation|domain)($|_)", color_by)) {
    return(.celladmix_is_integer_like(x))
  }
  FALSE
}

.celladmix_geom_spatial_point <- function(data, mapping, raster, point_size, alpha) {
  if (isTRUE(raster) && requireNamespace("ggrastr", quietly = TRUE)) {
    return(ggrastr::geom_point_rast(data = data, mapping = mapping,
      size = point_size, alpha = alpha))
  }
  ggplot2::geom_point(data = data, mapping = mapping, size = point_size, alpha = alpha)
}

#' Plot cell-level values in tissue coordinates.
#'
#' `celladmix_plot_spatial()` is a data-frame based helper for visualizing cell
#' annotations, dominant factors, factor fractions, and other cell-level
#' summaries. It expects one row per cell and coordinate columns, but otherwise
#' makes no assumptions about the upstream data source.
#'
#' @param data Cell-level data frame.
#' @param color_by Column used for point color.
#' @param x,y Coordinate column names.
#' @param point_size Point size. When `NULL`, a size is chosen from `nrow(data)`.
#' @param alpha Point alpha.
#' @param colors Optional named color vector for discrete labels.
#' @param legend_title Legend title. Defaults to `color_by`.
#' @param title Optional plot title.
#' @param reverse_y If `TRUE`, reverse the y-axis. Useful for image-coordinate
#'   datasets such as Xenium when aligning to stain images.
#' @param raster If `TRUE`, use `ggrastr::geom_point_rast()` when available.
#' @param theme Plot theme style.
#' @param continuous_palette Two colors used for continuous gradients.
#'
#' @return A `ggplot2` object.
#' @export
celladmix_plot_spatial <- function(
    data,
    color_by,
    x = "x",
    y = "y",
    point_size = NULL,
    alpha = 0.8,
    colors = NULL,
    legend_title = NULL,
    title = NULL,
    reverse_y = FALSE,
    raster = TRUE,
    theme = c("void", "classic"),
    continuous_palette = c("grey90", "#2166AC")
  ) {
  .celladmix_require_ggplot2()
  theme <- match.arg(theme)
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  required <- c(x, y, color_by)
  missing <- setdiff(required, colnames(data))
  if (length(missing)) {
    stop("data is missing required columns: ", paste(missing, collapse = ", "))
  }
  plot_df <- data[is.finite(as.numeric(data[[x]])) & is.finite(as.numeric(data[[y]])), ,
    drop = FALSE]
  if (!nrow(plot_df)) {
    stop("No rows with finite spatial coordinates are available")
  }
  point_size <- point_size %||% .celladmix_auto_point_size(nrow(plot_df))
  legend_title <- legend_title %||% color_by
  color_values <- plot_df[[color_by]]
  is_discrete <- .celladmix_is_discrete_color(color_values, color_by = color_by)
  if (is_discrete) {
    plot_df[[color_by]] <- factor(color_values)
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[x]], y = .data[[y]],
    color = .data[[color_by]])) +
    .celladmix_geom_spatial_point(
      data = plot_df,
      mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]],
        color = .data[[color_by]]),
      raster = raster,
      point_size = point_size,
      alpha = alpha
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, x = x, y = y, color = legend_title) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = 3, alpha = 1)
    ))

  if (isTRUE(reverse_y)) {
    p <- p + ggplot2::scale_y_reverse()
  }
  if (is_discrete) {
    if (!is.null(colors)) {
      p <- p + ggplot2::scale_colour_manual(values = colors)
    }
  } else {
    p <- p + ggplot2::scale_colour_gradient(low = continuous_palette[[1]],
      high = continuous_palette[[2]], na.value = "grey85")
  }
  if (theme == "void") {
    p <- p + ggplot2::theme_void()
  } else {
    p <- p + ggplot2::theme_classic(base_size = 10)
  }
  p + ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.text = ggplot2::element_text(size = 8),
    legend.title = ggplot2::element_text(size = 10)
  )
}
