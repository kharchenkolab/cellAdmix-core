# R6 workflow API layered over the lower-level batch functions.

#' Define a Tabular Input Schema
#'
#' Convenience helper for [cellAdmix()] when `format = "tabular"`.
#'
#' @return A named list of tabular column mappings.
#' @export
celladmix_schema <- function(
    x = "x",
    y = "y",
    z = "z",
    gene = "gene",
    qv = NULL,
    cell = "cell",
    cell_type = NULL,
    sample_id = NULL,
    fov_id = NULL,
    segmentation_mask = NULL,
    cell_metadata = NULL,
    cell_metadata_cell = "cell",
    cell_metadata_cell_type = NULL
  ) {
  list(
    x = x,
    y = y,
    z = z,
    gene = gene,
    qv = qv,
    cell = cell,
    cell_type = cell_type,
    sample_id = sample_id,
    fov_id = fov_id,
    segmentation_mask = segmentation_mask,
    cell_metadata = cell_metadata,
    cell_metadata_cell = cell_metadata_cell,
    cell_metadata_cell_type = cell_metadata_cell_type
  )
}
