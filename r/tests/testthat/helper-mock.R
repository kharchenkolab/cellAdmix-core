write_mock_xenium_bundle <- function(td, sim, include_cell_type = TRUE) {
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      "{",
      '  "run_name": "Mock Xenium",',
      '  "pixel_size": 0.2125,',
      '  "z_step_size": 3.0',
      "}"
    ),
    file.path(td, "experiment.xenium")
  )

  tx <- data.frame(
    transcript_id = sim$transcripts$transcript_id,
    cell_id = sim$transcripts$cell,
    overlaps_nucleus = seq_len(nrow(sim$transcripts)) %% 3L == 0L,
    feature_name = sim$transcripts$gene,
    x_location = sim$transcripts$x,
    y_location = sim$transcripts$y,
    z_location = sim$transcripts$z,
    qv = rep(30, nrow(sim$transcripts))
  )
  tx$nucleus_distance <- as.numeric(seq_len(nrow(tx)) %% 11L)
  utils::write.csv(tx, gzfile(file.path(td, "transcripts.csv.gz")), row.names = FALSE)

  cells <- data.frame(
    cell_id = sim$cells$cell_id,
    x_centroid = sim$cells$x,
    y_centroid = sim$cells$y,
    z_centroid = sim$cells$z
  )
  if (include_cell_type) {
    cells$cell_type <- sim$cells$cell_type
  }
  utils::write.csv(cells, gzfile(file.path(td, "cells.csv.gz")), row.names = FALSE)
}

write_mock_tabular_bundle <- function(path, sim, include_cell_type = TRUE) {
  tx <- data.frame(
    x = sim$transcripts$x,
    y = sim$transcripts$y,
    z = sim$transcripts$z,
    gene = sim$transcripts$gene,
    cell_id = sim$transcripts$cell
  )
  if (include_cell_type) {
    tx$cell_type <- sim$transcripts$cell_type
  }
  utils::write.csv(tx, path, row.names = FALSE)
}
