#!/usr/bin/env Rscript

# Prepare tabular inputs for the CosMx NSCLC cellAdmix example.
#
# Download the raw example files from:
#   http://pklab.org/peterk/cellAdmix/examples/cosmx_nsclc/
#
# Place the files in a local input directory, normally `data/`, then run:
#   Rscript prepare_cosmx_tabular.R --input data --output prepared
#
# Expected input files:
#   data/SMI_Giotto_Object.RData
#   data/Lung5_Rep1_tx_file.csv.gz  (or .csv)
#   data/annotation_adj.csv.gz      (or .csv)
#
# Outputs:
#   prepared/molecules_all.csv.gz
#   prepared/cell_metadata_all.csv.gz
#
# The script does not require the Giotto package. It loads the saved Giotto
# object and reads only object attributes containing cell metadata and spatial
# cell centroids.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grepl(file_arg, args)], fixed = TRUE)
example_dir <- if (length(script_path)) {
  dirname(normalizePath(script_path[[1]], winslash = "/", mustWork = TRUE))
} else {
  getwd()
}

celladmix_prepare_usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript prepare_cosmx_tabular.R --input data --output prepared [options]\n\n",
    "Options:\n",
    "  --input DIR       Raw input directory. Default: data\n",
    "  --output DIR      Prepared output directory. Default: prepared\n",
    "  --force           Recreate prepared files even if they already exist\n",
    "  --fov IDS         Optional comma-separated FOV subset, e.g. 1,2,3\n",
    "  --chunk-size N    Transcript rows per streaming chunk. Default: 500000\n",
    "  --help            Show this message\n\n",
    "Download the raw files from:\n",
    "  http://pklab.org/peterk/cellAdmix/examples/cosmx_nsclc/\n\n",
    "Expected input files:\n",
    "  SMI_Giotto_Object.RData\n",
    "  Lung5_Rep1_tx_file.csv.gz  (or .csv)\n",
    "  annotation_adj.csv.gz      (or .csv)\n",
    sep = ""
  )
}

celladmix_parse_prepare_args <- function(argv) {
  opts <- list(
    input = file.path(example_dir, "data"),
    output = file.path(example_dir, "prepared"),
    force = Sys.getenv("CELLADMIX_FORCE_PREPARE", unset = "false") %in%
      c("1", "true", "TRUE", "yes", "YES"),
    fov = Sys.getenv("CELLADMIX_FOV_KEEP", unset = ""),
    chunk_size = Sys.getenv("CELLADMIX_PREPARE_CHUNK_SIZE", unset = "500000")
  )
  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (arg %in% c("--help", "-h")) {
      opts$help <- TRUE
    } else if (arg == "--force") {
      opts$force <- TRUE
    } else if (arg %in% c("--input", "-i")) {
      i <- i + 1L
      if (i > length(argv)) stop("--input requires a directory", call. = FALSE)
      opts$input <- argv[[i]]
    } else if (startsWith(arg, "--input=")) {
      opts$input <- sub("^--input=", "", arg)
    } else if (arg %in% c("--output", "-o")) {
      i <- i + 1L
      if (i > length(argv)) stop("--output requires a directory", call. = FALSE)
      opts$output <- argv[[i]]
    } else if (startsWith(arg, "--output=")) {
      opts$output <- sub("^--output=", "", arg)
    } else if (arg == "--fov") {
      i <- i + 1L
      if (i > length(argv)) stop("--fov requires comma-separated IDs", call. = FALSE)
      opts$fov <- argv[[i]]
    } else if (startsWith(arg, "--fov=")) {
      opts$fov <- sub("^--fov=", "", arg)
    } else if (arg == "--chunk-size") {
      i <- i + 1L
      if (i > length(argv)) stop("--chunk-size requires an integer", call. = FALSE)
      opts$chunk_size <- argv[[i]]
    } else if (startsWith(arg, "--chunk-size=")) {
      opts$chunk_size <- sub("^--chunk-size=", "", arg)
    } else {
      stop("Unknown argument: ", arg, "\nRun with --help for usage.", call. = FALSE)
    }
    i <- i + 1L
  }
  opts
}

celladmix_abs_path <- function(path) {
  if (grepl("^/", path)) {
    return(path)
  }
  normalizePath(file.path(getwd(), path), winslash = "/", mustWork = FALSE)
}

celladmix_gzip_file <- function(path, force = TRUE) {
  gz_path <- paste0(path, ".gz")
  if (file.exists(gz_path) && !force) {
    return(gz_path)
  }
  status <- suppressWarnings(system2("gzip", c(if (force) "-f" else NULL, path)))
  if (!identical(as.integer(status), 0L)) {
    input <- file(path, "rb")
    output <- gzfile(gz_path, "wb")
    on.exit({
      try(close(input), silent = TRUE)
      try(close(output), silent = TRUE)
    }, add = TRUE)
    repeat {
      chunk <- readBin(input, what = raw(), n = 1024 * 1024)
      if (!length(chunk)) {
        break
      }
      writeBin(chunk, output)
    }
    unlink(path)
  }
  gz_path
}

celladmix_compress_existing_csv <- function(gz_path, force = FALSE) {
  legacy_path <- sub("\\.gz$", "", gz_path)
  if (file.exists(gz_path) && !force) {
    return(gz_path)
  }
  if (file.exists(legacy_path) && !force) {
    return(celladmix_gzip_file(legacy_path, force = TRUE))
  }
  gz_path
}

celladmix_find_input_file <- function(input_dir, candidates, label) {
  paths <- file.path(input_dir, candidates)
  hit <- paths[file.exists(paths)]
  if (length(hit)) {
    return(hit[[1]])
  }
  stop(
    "Missing ", label, " in ",
    normalizePath(input_dir, winslash = "/", mustWork = FALSE), ".\n",
    "Expected one of: ", paste(candidates, collapse = ", "), call. = FALSE
  )
}

celladmix_load_nsclc_cell_metadata <- function(input_dir) {
  giotto_path <- celladmix_find_input_file(input_dir,
    "SMI_Giotto_Object.RData", "Giotto object")
  annotation_path <- celladmix_find_input_file(input_dir,
    c("annotation_adj.csv.gz", "annotation_adj.csv"), "cell annotation CSV")

  env <- new.env(parent = emptyenv())
  load(giotto_path, envir = env)
  if (!exists("gem", envir = env, inherits = FALSE)) {
    stop("Expected object named 'gem' in ", giotto_path, call. = FALSE)
  }
  gem_attr <- attributes(get("gem", envir = env))
  cell_meta <- as.data.frame(gem_attr$cell_metadata$rna)
  cell_locs <- as.data.frame(gem_attr$spatial_locs$raw)
  cell_meta <- cbind(cell_locs[, c("sdimx", "sdimy")], cell_meta)
  cell_meta <- cell_meta[cell_meta$Run_Tissue_name == "Lung5_Rep1", , drop = FALSE]
  cell_meta$cell <- cell_meta$cell_ID
  names(cell_meta)[names(cell_meta) == "sdimx"] <- "x"
  names(cell_meta)[names(cell_meta) == "sdimy"] <- "y"
  cell_meta$z <- 1

  annotation <- read.csv(annotation_path, stringsAsFactors = FALSE)
  required <- c("cell", "cell_type")
  missing <- setdiff(required, names(annotation))
  if (length(missing)) {
    stop("annotation_adj.csv is missing columns: ", paste(missing, collapse = ", "),
      call. = FALSE)
  }
  if (nrow(annotation) != nrow(cell_meta) ||
      !identical(annotation$cell, cell_meta$cell)) {
    idx <- match(cell_meta$cell, annotation$cell)
    if (anyNA(idx)) {
      stop("annotation_adj.csv does not cover all Lung5_Rep1 Giotto cells",
        call. = FALSE)
    }
    annotation <- annotation[idx, , drop = FALSE]
  }

  cell_meta$celltype <- annotation$cell_type
  cell_meta$cell_type <- annotation$cell_type
  if ("fov" %in% names(annotation)) {
    cell_meta$fov <- annotation$fov
  }
  if ("cell_id" %in% names(annotation)) {
    cell_meta$cell_id <- annotation$cell_id
  }

  immune <- c(
    "B-cell", "B cells", "NK", "T CD4 memory", "T CD4 naive",
    "T CD8 memory", "T CD8 naive", "Treg", "plasmablast", "mast",
    "mDC", "monocyte", "pDC", "neutrophil", "DC", "CD4+ T cells",
    "CD8+ T cells"
  )
  malignant <- startsWith(cell_meta$celltype, "tumor")
  cell_meta$celltype[malignant] <- "malignant"
  cell_meta$cell_type[startsWith(cell_meta$cell_type, "tumor")] <- "malignant"
  cell_meta$cell_type_coarse <- ifelse(cell_meta$celltype %in% immune,
    "immune other", cell_meta$celltype)
  cell_meta$regions_compare <- vapply(
    cell_meta$niche,
    function(x) {
      if (identical(x, "tumor interior")) return("tumor")
      if (identical(x, "stroma")) return("stroma")
      NA_character_
    },
    character(1)
  )
  rownames(cell_meta) <- cell_meta$cell
  cell_meta
}

celladmix_write_prepared_inputs <- function(input_dir, prepared_dir, force = FALSE,
                                            fov_keep = NULL, chunk_size = 500000) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("The readr package is required to stream the raw CosMx transcript CSV",
      call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The data.table package is required to write normalized transcript chunks",
      call. = FALSE)
  }

  tx_path <- celladmix_find_input_file(input_dir,
    c("Lung5_Rep1_tx_file.csv.gz", "Lung5_Rep1_tx_file.csv"),
    "CosMx transcript CSV")
  dir.create(prepared_dir, recursive = TRUE, showWarnings = FALSE)

  cell_meta <- celladmix_load_nsclc_cell_metadata(input_dir)
  if (!is.null(fov_keep)) {
    cell_meta <- cell_meta[cell_meta$fov %in% fov_keep, , drop = FALSE]
  }
  if (nrow(cell_meta) == 0L) {
    stop("No cells remain after applying the FOV filter", call. = FALSE)
  }
  cell_types <- stats::setNames(cell_meta$celltype, cell_meta$cell)

  scope_id <- if (is.null(fov_keep)) {
    "all"
  } else {
    paste0("fov_", paste(sort(unique(as.integer(fov_keep))), collapse = "_"))
  }

  cell_metadata_path <- file.path(prepared_dir,
    paste0("cell_metadata_", scope_id, ".csv.gz"))
  cell_metadata_cols <- intersect(
    c(
      "cell", "cell_type", "celltype", "cell_type_coarse", "regions_compare",
      "niche", "x", "y", "z", "fov", "cell_id"
    ),
    names(cell_meta)
  )
  cell_metadata_path <- celladmix_compress_existing_csv(cell_metadata_path,
    force = force)
  if (force || !file.exists(cell_metadata_path)) {
    con <- gzfile(cell_metadata_path, "wt")
    on.exit(try(if (isOpen(con)) close(con), silent = TRUE), add = TRUE)
    utils::write.csv(cell_meta[, cell_metadata_cols, drop = FALSE], con,
      row.names = FALSE)
    close(con)
  }

  molecules_path <- file.path(prepared_dir, paste0("molecules_", scope_id, ".csv.gz"))
  if (force && file.exists(molecules_path)) {
    unlink(molecules_path)
  }
  molecules_path <- celladmix_compress_existing_csv(molecules_path, force = force)
  molecules_tmp_path <- sub("\\.gz$", "", molecules_path)

  if (!file.exists(molecules_path)) {
    if (file.exists(molecules_tmp_path)) {
      unlink(molecules_tmp_path)
    }
    wrote_header <- FALSE
    callback <- readr::SideEffectChunkCallback$new(function(chunk, pos) {
      chunk <- as.data.frame(chunk)
      cell <- paste0("c_1_", chunk$fov, "_", chunk$cell_ID)
      keep <- chunk$cell_ID != 0 &
        !grepl("^NegPrb", chunk$target) &
        cell %in% names(cell_types)
      if (!is.null(fov_keep)) {
        keep <- keep & chunk$fov %in% fov_keep
      }
      if (!any(keep)) {
        return(invisible(NULL))
      }
      out <- data.frame(
        x = chunk$x_global_px[keep] * 0.18,
        y = chunk$y_global_px[keep] * 0.18,
        z = chunk$z[keep] * 0.8,
        gene = chunk$target[keep],
        cell = cell[keep],
        stringsAsFactors = FALSE
      )
      data.table::fwrite(out, molecules_tmp_path, append = wrote_header,
        col.names = !wrote_header)
      wrote_header <<- TRUE
      invisible(NULL)
    })
    readr::read_csv_chunked(
      tx_path,
      callback = callback,
      chunk_size = as.integer(chunk_size),
      col_types = readr::cols(
        fov = readr::col_integer(),
        cell_ID = readr::col_integer(),
        x_global_px = readr::col_double(),
        y_global_px = readr::col_double(),
        z = readr::col_double(),
        target = readr::col_character(),
        .default = readr::col_skip()
      ),
      progress = interactive()
    )
    if (file.exists(molecules_tmp_path)) {
      celladmix_gzip_file(molecules_tmp_path, force = TRUE)
    }
  }

  list(
    cell_meta = cell_meta,
    molecules_path = molecules_path,
    cell_metadata_path = cell_metadata_path,
    input_dir = input_dir,
    prepared_dir = prepared_dir
  )
}

argv <- commandArgs(trailingOnly = TRUE)
opts <- celladmix_parse_prepare_args(argv)
if (isTRUE(opts$help)) {
  celladmix_prepare_usage()
  quit(status = 0L)
}

input_dir <- celladmix_abs_path(opts$input)
prepared_dir <- celladmix_abs_path(opts$output)
fov_keep <- if (nzchar(opts$fov)) {
  as.integer(strsplit(opts$fov, ",", fixed = TRUE)[[1]])
} else {
  NULL
}
chunk_size <- as.integer(opts$chunk_size)
if (is.na(chunk_size) || !is.finite(chunk_size) || chunk_size <= 0L) {
  stop("--chunk-size must be a positive integer", call. = FALSE)
}

if (!dir.exists(input_dir)) {
  cat("Input directory does not exist: ", input_dir, "\n\n", sep = "")
  celladmix_prepare_usage()
  quit(status = 1L)
}

prepared <- tryCatch(
  celladmix_write_prepared_inputs(
    input_dir = input_dir,
    prepared_dir = prepared_dir,
    force = opts$force,
    fov_keep = fov_keep,
    chunk_size = chunk_size
  ),
  error = function(e) {
    cat(conditionMessage(e), "\n\n", sep = "")
    celladmix_prepare_usage()
    quit(status = 1L)
  }
)

cat("Prepared CosMx tabular files:\n")
cat("  input: ", prepared$input_dir, "\n", sep = "")
cat("  output: ", prepared$prepared_dir, "\n", sep = "")
cat("  molecules: ", prepared$molecules_path, "\n", sep = "")
cat("  cell metadata: ", prepared$cell_metadata_path, "\n", sep = "")
cat("  cells: ", nrow(prepared$cell_meta), "\n", sep = "")
cat("  metadata columns: ",
  paste(names(read.csv(prepared$cell_metadata_path, nrows = 1)), collapse = ", "),
  "\n", sep = "")
