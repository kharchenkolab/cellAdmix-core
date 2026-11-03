# Atera output bundles carry an `experiment.spatial` manifest and zipped Zarr
# stores. The store builder reads the Xenium layout, so an Atera bundle is read
# through Xenium-layout Parquet tables that the Python module celladmix.atera
# writes once under <output_dir>/atera_xenium. The conversion is not duplicated
# in R or C++: R reuses an existing conversion whose recorded source files are
# unchanged, and otherwise runs the Python command.

.celladmix_atera_stores <- c("transcripts.zarr.zip", "cells.zarr.zip", "cell_feature_matrix.zarr.zip")
.celladmix_atera_tables <- c(
  "transcripts.parquet", "cells.parquet",
  "cell_boundaries.parquet", "nucleus_boundaries.parquet"
)

# Bundle directory for an Atera output directory or its experiment.spatial;
# NULL for anything else, including directories that also have experiment.xenium.
.celladmix_atera_bundle_dir <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    return(NULL)
  }
  if (!dir.exists(path) && file.exists(path) && identical(basename(path), "experiment.spatial")) {
    path <- dirname(path)
  }
  if (dir.exists(path) && file.exists(file.path(path, "experiment.spatial")) &&
      !file.exists(file.path(path, "experiment.xenium"))) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  NULL
}

# Whether `cache_dir` holds a complete conversion of the current bundle files,
# judged by the sizes and modification times recorded in atera_source.json.
.celladmix_atera_tables_current <- function(bundle_dir, cache_dir) {
  record_path <- file.path(cache_dir, "atera_source.json")
  needed <- file.path(cache_dir, c(.celladmix_atera_tables, "experiment.xenium"))
  if (!file.exists(record_path) || !all(file.exists(needed))) {
    return(FALSE)
  }
  record <- tryCatch(jsonlite::read_json(record_path), error = function(e) NULL)
  if (is.null(record) || is.null(record$files)) {
    return(FALSE)
  }
  for (name in c("experiment.spatial", .celladmix_atera_stores)) {
    info <- file.info(file.path(bundle_dir, name))
    recorded <- record$files[[name]]
    if (is.na(info$size) || is.null(recorded) ||
        !isTRUE(as.numeric(recorded$size) == info$size) ||
        !isTRUE(as.numeric(recorded$mtime) == floor(as.numeric(info$mtime)))) {
      return(FALSE)
    }
  }
  TRUE
}

# Python interpreter used for the conversion: option cellAdmix.python, then the
# CELLADMIX_PYTHON environment variable, then python3 on the PATH.
.celladmix_atera_python <- function() {
  python <- getOption("cellAdmix.python", Sys.getenv("CELLADMIX_PYTHON", ""))
  if (!nzchar(python)) {
    python <- Sys.which("python3")
  }
  unname(python)
}

# Return the Xenium-layout directory for an Atera bundle, converting when the
# cached tables are missing or were written from different source files.
.celladmix_atera_xenium_dir <- function(bundle_dir, output_dir, force = FALSE) {
  cache_dir <- file.path(output_dir, "atera_xenium")
  if (!isTRUE(force) && .celladmix_atera_tables_current(bundle_dir, cache_dir)) {
    return(normalizePath(cache_dir, winslash = "/", mustWork = TRUE))
  }
  args <- c("-m", "celladmix.atera", if (isTRUE(force)) "--force", shQuote(bundle_dir), shQuote(cache_dir))
  python <- .celladmix_atera_python()
  command <- paste(c(if (nzchar(python)) shQuote(python) else "python3", args), collapse = " ")
  status <- if (nzchar(python) && nzchar(Sys.which(python))) {
    tryCatch(suppressWarnings(system2(python, args)), error = function(e) 127L)
  } else {
    127L
  }
  if (!identical(as.integer(status), 0L) || !.celladmix_atera_tables_current(bundle_dir, cache_dir)) {
    stop(
      "Reading the Atera bundle ", bundle_dir, " needs its Zarr stores written as ",
      "Xenium-layout tables, which the Python package celladmix (with zarr >= 3) does. ",
      "The conversion command failed or no Python interpreter was found:\n  ", command, "\n",
      "Install celladmix for that interpreter, or set options(cellAdmix.python = ...) ",
      "or CELLADMIX_PYTHON to one that has it, or run the command yourself and retry; ",
      "the tables in ", cache_dir, " are then reused.",
      call. = FALSE
    )
  }
  normalizePath(cache_dir, winslash = "/", mustWork = TRUE)
}
