#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(cellAdmixCore))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

usage <- function(status = 0L) {
  cat("
Usage:
  Rscript scripts/celladmix_batch.R --input PATH --output DIR [options]

Core options:
  --input PATH                         Xenium bundle, tabular molecule table, or prepared source
  --output DIR                         Output/cache directory
  --format auto|xenium|tabular         Input format [auto]
  --threads N                          Worker threads [parallel::detectCores()]
  --verbose                            Print progress messages
  --overwrite                          Recompute fit even if cached run exists
  --keep-non-gene                      Keep Xenium control/codeword/non-gene features

Annotation:
  --annotation PATH                    Annotation CSV/CSV.GZ
  --annotation-col COL                 Annotation label column
  --cell-id-col COL                    Annotation cell id column [cell_id]
  --annotation-name NAME               Annotation name [manual]
  --auto-annotate                      Run clustering if no annotation is supplied or inferred
  --annotation-report                  Render clustering report without main --report
  --cluster-name NAME                  Generated cluster annotation name [cluster]
  --cluster-min-molecules N            Minimum molecules per clustered cell [10]
  --cluster-min-genes N                Minimum detected non-control genes per clustered cell [5]
  --cluster-resolution X               Louvain resolution [1]
  --cluster-umap true|false            Compute UMAP for annotation report [true]
  --cluster-markers true|false         Add marker genes to annotation report [true]
  --cluster-top-markers N              Top marker genes per cluster [5]
  --cluster-marker-min-fraction X      Minimum within-cluster detection fraction [0.05]
  --cluster-marker-min-logfc X         Minimum marker log2 fold-change [0]
  --cluster-marker-max-panels N        Maximum marker UMAP panels [25]

Tabular input:
  --molecules PATH                     Molecule table path; defaults to --input for tabular
  --cell-metadata PATH                 Cell metadata CSV/CSV.GZ
  --x-col COL                          Molecule x column [x]
  --y-col COL                          Molecule y column [y]
  --z-col COL                          Molecule z column [z]
  --gene-col COL                       Molecule gene column [gene]
  --cell-col COL                       Molecule cell column [cell]
  --qv-col COL                         Molecule quality column
  --cell-type-col COL                  Molecule table cell-type column
  --cell-metadata-cell COL             Cell metadata cell id column [cell]
  --cell-metadata-cell-type COL        Cell metadata annotation column
  --segmentation-mask PATH             Optional segmentation mask for tabular input
  --sample-id-col COL                  Optional sample id column
  --fov-id-col COL                     Optional FOV id column

Fit:
  --rank auto|N                        Factor rank [auto]
  --rank-multiplier X                  Auto-rank multiplier [1.2]
  --rank-cap N                         Auto-rank cap [30]
  --nmf-variant invsqrt_kl|kl|sqrt_kl|ls_nmf
                                      NMF variant [invsqrt_kl]
  --nmf-init auto|random|cluster       NMF initialization [auto]
  --nmf-runs auto|N                    NMF random starts [auto = threads]
  --nmf-iterations N                   NMF iterations [package default]
  --seed N                             Random seed [1]
  --run-id NAME                        Optional fit run id

Scoring and correction:
  --score auto|membrane|bridge|coherence
                                      Score method [auto]
  --score-name NAME                    Score registry name [method]
  --p-thresh X                         Rule p-value threshold [0.1]
  --adjust-p                           Use adjusted p-values for rules
  --targets A,B,C                      Restrict correction rules to target cell types
  --no-correct                         Score only; skip correction
  --correction-name NAME               Correction run name [<score>_clean]
  --max-cells-per-type-pair N          Pair summary cap for membrane/bridge [400]
  --candidate-pairs-per-type-pair N    Membrane candidate-pair cap [400]
  --min-factor-molecules N             Minimum factor molecules per pair/cell [5]
  --min-pairs N                        Minimum scored pairs for summary [5]
  --no-null                            Skip bridge null computation
  --null-iterations N                  Bridge null iterations [3]

Reports:
  --report                             Write HTML report(s)
  --report-file PATH                   Report output path
  --top-genes N                        Top loading genes per factor [10]

Examples:
  Rscript scripts/celladmix_batch.R --input data --output out \\
    --annotation annotations/annotation.csv.gz --annotation-col merged_annotation \\
    --score auto --report --threads 10

  Rscript scripts/celladmix_batch.R --input data --output out \\
    --auto-annotate --score auto --report --threads 10

  Rscript scripts/celladmix_batch.R --format tabular --input molecules.csv.gz \\
    --output out --cell-metadata cells.csv.gz --cell-metadata-cell-type cell_type
")
  quit(save = "no", status = status)
}

parse_args <- function(args) {
  out <- list()
  i <- 1L
  flags <- c("help", "verbose", "overwrite", "adjust-p", "no-correct",
    "annotation-report", "auto-annotate", "report", "no-null", "keep-non-gene")
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!startsWith(arg, "--")) {
      stop("Unexpected positional argument: ", arg)
    }
    key_value <- sub("^--", "", arg)
    if (grepl("=", key_value, fixed = TRUE)) {
      parts <- strsplit(key_value, "=", fixed = TRUE)[[1]]
      key <- parts[[1]]
      value <- paste(parts[-1], collapse = "=")
    } else {
      key <- key_value
      if (key %in% flags) {
        value <- TRUE
      } else {
        if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
          stop("Missing value for --", key)
        }
        i <- i + 1L
        value <- args[[i]]
      }
    }
    out[[gsub("-", "_", key, fixed = TRUE)]] <- value
    i <- i + 1L
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  if (is.logical(x)) {
    return(isTRUE(x))
  }
  tolower(as.character(x[[1]])) %in% c("1", "true", "yes", "y", "on")
}

as_int <- function(x, default = NULL) {
  if (is.null(x) || !nzchar(as.character(x[[1]]))) {
    return(default)
  }
  value <- suppressWarnings(as.integer(x[[1]]))
  if (is.na(value)) {
    stop("Expected integer, got: ", x[[1]])
  }
  value
}

as_num <- function(x, default = NULL) {
  if (is.null(x) || !nzchar(as.character(x[[1]]))) {
    return(default)
  }
  value <- suppressWarnings(as.numeric(x[[1]]))
  if (is.na(value)) {
    stop("Expected numeric value, got: ", x[[1]])
  }
  value
}

as_chr <- function(x, default = NULL) {
  if (is.null(x) || length(x) == 0L || !nzchar(as.character(x[[1]]))) {
    return(default)
  }
  as.character(x[[1]])
}

as_csv <- function(x) {
  value <- as_chr(x, NULL)
  if (is.null(value)) {
    return(NULL)
  }
  value <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  value[nzchar(value)]
}

info <- function(message) {
  cat(sprintf("[batch %s] %s\n", format(Sys.time(), "%H:%M:%S"), message))
  flush.console()
}

require_file <- function(path, label) {
  if (is.null(path) || !file.exists(path)) {
    stop(label, " does not exist: ", path %||% "<missing>")
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

write_json_flat <- function(x, path) {
  quote_json <- function(value) {
    if (is.null(value) || length(value) == 0L || is.na(value[[1]])) {
      return("null")
    }
    if (is.logical(value)) {
      return(if (isTRUE(value[[1]])) "true" else "false")
    }
    if (is.numeric(value)) {
      return(as.character(value[[1]]))
    }
    paste0("\"", gsub("\"", "\\\\\"", as.character(value[[1]]), fixed = TRUE), "\"")
  }
  lines <- vapply(names(x), function(name) {
    sprintf("  \"%s\": %s", name, quote_json(x[[name]]))
  }, character(1))
  writeLines(c("{", paste(lines, collapse = ",\n"), "}"), path)
}

factor_top_genes <- function(fit, n = 10L) {
  h <- fit$loadings()
  genes <- fit$run$genes %||% colnames(h) %||% paste0("gene_", seq_len(ncol(h)))
  colnames(h) <- genes
  pieces <- lapply(seq_len(nrow(h)), function(factor) {
    vals <- as.numeric(h[factor, ])
    denom <- sum(vals, na.rm = TRUE)
    frac <- if (is.finite(denom) && denom > 0) vals / denom else vals
    idx <- head(order(frac, decreasing = TRUE), max(1L, as.integer(n)))
    data.frame(
      factor = paste0("F", factor),
      rank = seq_along(idx),
      gene = genes[idx],
      loading = vals[idx],
      loading_fraction = frac[idx],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

safe_write_csv <- function(x, path) {
  if (is.null(x)) {
    x <- data.frame()
  }
  if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }
  utils::write.csv(x, path, row.names = FALSE)
}

membrane_available <- function(fit) {
  ok <- tryCatch({
    img <- fit$stain("membrane")
    is.list(img) && file.exists(img$path)
  }, error = function(e) FALSE)
  isTRUE(ok)
}

resolve_score_method <- function(score_method, fit) {
  if (!identical(score_method, "auto")) {
    if (identical(score_method, "membrane") && !membrane_available(fit)) {
      stop("--score membrane was requested, but no membrane stain could be discovered")
    }
    return(score_method)
  }
  if (membrane_available(fit)) {
    "membrane"
  } else {
    "bridge"
  }
}

render_minimal_report <- function(path, env, top_genes_n, score_pairs_height) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    warning("The rmarkdown package is not available; skipping HTML report")
    return(FALSE)
  }
  rmd <- file.path(dirname(path), paste0(tools::file_path_sans_ext(basename(path)), ".Rmd"))
  writeLines(c(
    "---",
    "title: \"cellAdmix Batch Report\"",
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_depth: 2",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.align = 'center', out.width = '100%', dpi = 150)",
    "```",
    "",
    "## Summary",
    "",
    "```{r summary}",
    "fit$summary()",
    "```",
    "",
    "```{r score-summary}",
    "head(score$summary(), 20)",
    "```",
    "",
    "```{r correction-summary}",
    "if (!is.null(correction)) correction$summary() else data.frame(note = 'correction skipped')",
    "```",
    "",
    "## Factor Loadings",
    "",
    "```{r top-genes}",
    "head(top_genes, max(20, fit$rank * top_genes_n))",
    "```",
    "",
    "```{r loadings-plot, fig.width=9, fig.height=max(4.5, ceiling(fit$rank / 5) * 2.8)}",
    "fit$plot_loadings(n_genes = top_genes_n)",
    "```",
    "",
    "```{r stability, fig.width=6.2, fig.height=4.2}",
    "fit$plot_stability()",
    "```",
    "",
    "## Score",
    "",
    "```{r score-rules}",
    "head(rules[, intersect(c('factor', 'source_cell_type', 'target_cell_type', 'p_value'), names(rules)), drop = FALSE], 20)",
    "```",
    "",
    "```{r score-heatmap, fig.width=8.5, fig.height=5.4}",
    "score$plot_heatmap(p_thresh = p_thresh)",
    "```",
    "",
    "```{r score-pairs, fig.width=9.6, fig.height=score_pairs_height}",
    "score$plot_pairs()",
    "```",
    "",
    "## Correction",
    "",
    "```{r removal, fig.width=6.5, fig.height=5.2}",
    "if (!is.null(correction)) correction$plot_removed_molecules()",
    "```"
  ), rmd)
  env$top_genes_n <- top_genes_n
  env$score_pairs_height <- score_pairs_height
  rmarkdown::render(rmd, output_file = basename(path), output_dir = dirname(path),
    envir = env, quiet = TRUE)
  TRUE
}

render_annotation_report <- function(path, clust, markers = NULL, top_markers = NULL,
                                     counts = NULL, top_n = 5L, max_panels = 25L) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    warning("The rmarkdown package is not available; skipping annotation report")
    return(FALSE)
  }
  rmd <- file.path(dirname(path), paste0(tools::file_path_sans_ext(basename(path)), ".Rmd"))
  env <- new.env(parent = globalenv())
  env$clust <- clust
  env$markers <- markers
  env$top_markers <- top_markers
  env$counts <- counts
  env$top_n <- top_n
  env$max_panels <- max_panels
  writeLines(c(
    "---",
    "title: \"cellAdmix Cluster Annotation Report\"",
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_depth: 2",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.align = 'center', out.width = '100%', dpi = 150)",
    "```",
    "",
    "This report shows unsupervised cell clusters used as temporary labels. They are not curated biological annotations.",
    "",
    "## Cluster Summary",
    "",
    "```{r cluster-summary}",
    "clust$summary",
    "table(clust$clusters$cluster)",
    "```",
    "",
    "```{r cluster-umap, fig.width=6, fig.height=5}",
    "if (is.data.frame(clust$embedding) && nrow(clust$embedding)) {",
    "  df <- clust$embedding",
    "  value <- factor(df$cluster)",
    "  cols <- grDevices::adjustcolor(grDevices::rainbow(nlevels(value))[value], alpha.f = 0.55)",
    "  graphics::par(mar = c(3.5, 3.5, 2.0, 0.5), mgp = c(2, 0.65, 0), cex = 2/3)",
    "  graphics::plot(df$umap_1, df$umap_2, col = cols, pch = 16, cex = 0.35, asp = 1, xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'Cell UMAP by cluster')",
    "} else {",
    "  data.frame(note = 'UMAP was not computed for this clustering run')",
    "}",
    "```",
    "",
    "## Cluster Markers",
    "",
    "Marker scoring uses a sparse one-vs-rest detection AUC on cell-by-gene counts. These are provisional markers intended to help assign biological labels to the generated clusters.",
    "",
    "```{r marker-table}",
    "if (is.data.frame(top_markers) && nrow(top_markers)) {",
    "  display_cols <- head(names(top_markers), max(0, length(names(top_markers)) - 2))",
    "  top_markers[, display_cols, drop = FALSE]",
    "} else {",
    "  data.frame(note = 'cluster marker scoring was not run')",
    "}",
    "```",
    "",
    "```{r marker-expression-umap, fig.width=12, fig.height=12}",
    "if (is.data.frame(top_markers) && nrow(top_markers) && !is.null(counts)) {",
    "  celladmix_plot_cluster_marker_umap(clust$embedding, counts, top_markers, max_panels = max_panels)",
    "} else {",
    "  data.frame(note = 'top marker UMAP panels unavailable')",
    "}",
    "```",
    "",
    "```{r marker-auc-heatmap, fig.width=9, fig.height=max(5, length(unique(top_markers$gene)) * 0.16)}",
    "if (is.data.frame(markers) && nrow(markers)) {",
    "  celladmix_plot_cluster_marker_heatmap(markers, top_markers = top_markers, n = top_n)",
    "} else {",
    "  data.frame(note = 'marker AUC heatmap unavailable')",
    "}",
    "```",
    "",
    "```{r marker-dotplot, fig.width=9, fig.height=max(5, length(unique(top_markers$gene)) * 0.16)}",
    "if (is.data.frame(markers) && nrow(markers)) {",
    "  celladmix_plot_cluster_marker_dotplot(markers, top_markers = top_markers, n = top_n)",
    "} else {",
    "  data.frame(note = 'marker dotplot unavailable')",
    "}",
    "```"
  ), rmd)
  rmarkdown::render(rmd, output_file = basename(path), output_dir = dirname(path),
    envir = env, quiet = TRUE)
  TRUE
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
opt <- function(name) args[[name]]

if (as_bool(opt("help"), FALSE)) {
  usage(0L)
}

input <- as_chr(opt("input"))
output_dir <- as_chr(opt("output"))
if (is.null(input) || is.null(output_dir)) {
  usage(1L)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

format <- match.arg(as_chr(opt("format"), "auto"), c("auto", "xenium", "tabular"))
threads <- as_int(opt("threads"), max(1L, parallel::detectCores(logical = FALSE) %||% 1L))
verbose <- as_bool(opt("verbose"), FALSE)
overwrite <- as_bool(opt("overwrite"), FALSE)
seed <- as_int(opt("seed"), 1L)
p_thresh <- as_num(opt("p_thresh"), 0.1)
adjust_p <- as_bool(opt("adjust_p"), FALSE)
correct_enabled <- !as_bool(opt("no_correct"), FALSE)
targets <- as_csv(opt("targets"))
report_enabled <- as_bool(opt("report"), FALSE)

annotation <- as_chr(opt("annotation"), NULL)
if (!is.null(annotation)) {
  annotation <- require_file(annotation, "annotation")
}
annotation_col <- as_chr(opt("annotation_col"), NULL)
cell_id_col <- as_chr(opt("cell_id_col"), "cell_id")
annotation_name <- as_chr(opt("annotation_name"), "manual")

schema <- NULL
source <- input
if (identical(format, "tabular")) {
  molecules <- as_chr(opt("molecules"), input)
  schema <- celladmix_schema(
    x = as_chr(opt("x_col"), "x"),
    y = as_chr(opt("y_col"), "y"),
    z = as_chr(opt("z_col"), "z"),
    gene = as_chr(opt("gene_col"), "gene"),
    qv = as_chr(opt("qv_col"), NULL),
    cell = as_chr(opt("cell_col"), "cell"),
    cell_type = as_chr(opt("cell_type_col"), NULL),
    sample_id = as_chr(opt("sample_id_col"), NULL),
    fov_id = as_chr(opt("fov_id_col"), NULL),
    segmentation_mask = as_chr(opt("segmentation_mask"), NULL),
    cell_metadata = as_chr(opt("cell_metadata"), NULL),
    cell_metadata_cell = as_chr(opt("cell_metadata_cell"), "cell"),
    cell_metadata_cell_type = as_chr(opt("cell_metadata_cell_type"), NULL)
  )
  source <- molecules
}

info("Constructing dataset")
dataset_args <- list(
  source = source,
  output_dir = output_dir,
  format = format,
  schema = schema,
  annotation = annotation,
  annotation_name = annotation_name,
  annotation_col = annotation_col,
  cell_id_col = cell_id_col,
  num_threads = threads
)
if (!identical(format, "tabular")) {
  dataset_args$keep_non_gene <- as_bool(opt("keep_non_gene"), FALSE)
}
ds <- do.call(cellAdmix, dataset_args)
print(ds)

clust <- NULL
if (is.null(ds$active_annotation)) {
  if (!as_bool(opt("auto_annotate"), FALSE)) {
    stop("No annotation was supplied or inferred. Provide --annotation/--annotation-col, ",
      "configure tabular --cell-metadata-cell-type, or pass --auto-annotate ",
      "to generate provisional cluster labels.")
  }
  info("Running automatic cell clustering annotation")
  clust <- ds$cluster(
    name = as_chr(opt("cluster_name"), "cluster"),
    active = TRUE,
    min_molecules = as_int(opt("cluster_min_molecules"), 10L),
    min_genes = as_int(opt("cluster_min_genes"), 5L),
    resolution = as_num(opt("cluster_resolution"), 1),
    compute_umap = as_bool(opt("cluster_umap"), TRUE),
    seed = seed,
    verbose = verbose
  )
  generated_annotation <- ds$annotation(as_vector = FALSE)
  safe_write_csv(generated_annotation, file.path(output_dir, "generated_annotation.csv"))
  render_cluster_report <- isTRUE(report_enabled) || as_bool(opt("annotation_report"), FALSE)
  if (isTRUE(render_cluster_report)) {
    marker_enabled <- as_bool(opt("cluster_markers"), TRUE)
    markers <- NULL
    top_markers <- NULL
    marker_counts <- NULL
    top_n <- as_int(opt("cluster_top_markers"), 5L)
    if (isTRUE(marker_enabled)) {
      info("Collecting counts and scoring cluster markers")
      marker_counts <- ds$counts(force = FALSE, num_threads = threads, verbose = verbose)
      cluster_labels <- stats::setNames(as.factor(clust$clusters$cluster), clust$clusters$cell_id)
      common_cells <- intersect(colnames(marker_counts), names(cluster_labels))
      marker_counts <- marker_counts[, common_cells, drop = FALSE]
      cluster_labels <- cluster_labels[colnames(marker_counts)]
      markers <- celladmix_score_cluster_markers(
        marker_counts,
        cluster_labels,
        min_expression_fraction = as_num(opt("cluster_marker_min_fraction"), 0.05),
        min_logfc = as_num(opt("cluster_marker_min_logfc"), 0)
      )
      top_markers <- celladmix_top_cluster_markers(markers, n = top_n)
      safe_write_csv(markers, file.path(output_dir, "cluster_markers_auc.csv"))
      safe_write_csv(top_markers, file.path(output_dir, "cluster_markers_auc_top.csv"))
    }
    render_annotation_report(file.path(output_dir, "annotation_report.html"), clust,
      markers = markers, top_markers = top_markers, counts = marker_counts,
      top_n = top_n,
      max_panels = as_int(opt("cluster_marker_max_panels"), 25L))
  }
}

rank_arg <- as_chr(opt("rank"), "auto")
rank <- if (identical(rank_arg, "auto")) "auto" else as_int(rank_arg)
nmf_runs_arg <- as_chr(opt("nmf_runs"), "auto")
nmf_n_runs <- if (identical(nmf_runs_arg, "auto")) threads else as_int(nmf_runs_arg)
nmf_variant <- match.arg(as_chr(opt("nmf_variant"), "ls_nmf"),
  c("ls_nmf", "invsqrt_kl", "kl", "sqrt_kl"))
nmf_init <- match.arg(as_chr(opt("nmf_init"), "auto"), c("auto", "random", "cluster"))

fit_args <- list(
  rank = rank,
  rank_multiplier = as_num(opt("rank_multiplier"), 1.2),
  rank_cap = as_int(opt("rank_cap"), 30L),
  nmf_variant = nmf_variant,
  nmf_init = nmf_init,
  nmf_n_runs = nmf_n_runs,
  seed = seed,
  overwrite = overwrite,
  verbose = verbose
)
nmf_iterations <- as_int(opt("nmf_iterations"), NULL)
if (!is.null(nmf_iterations)) {
  fit_args$nmf_iterations <- nmf_iterations
}
run_id <- as_chr(opt("run_id"), NULL)
if (!is.null(run_id)) {
  fit_args$run_id <- run_id
}

info("Fitting factors")
fit <- do.call(ds$fit, fit_args)
print(fit)

score_method <- resolve_score_method(match.arg(as_chr(opt("score"), "auto"),
  c("auto", "membrane", "bridge", "coherence")), fit)
score_name <- as_chr(opt("score_name"), score_method)

score_args <- list(
  name = score_name,
  num_threads = threads,
  min_factor_molecules = as_int(opt("min_factor_molecules"), 5L),
  verbose = verbose
)
if (score_method %in% c("membrane", "bridge")) {
  score_args$max_cells_per_type_pair <- as_int(opt("max_cells_per_type_pair"), 400L)
  score_args$min_pairs <- as_int(opt("min_pairs"), 5L)
}
if (score_method %in% c("membrane", "bridge")) {
  score_args$candidate_pairs_per_type_pair <- as_int(opt("candidate_pairs_per_type_pair"), 400L)
}
if (identical(score_method, "bridge")) {
  score_args$compute_null <- !as_bool(opt("no_null"), FALSE)
  score_args$null_iterations <- as_int(opt("null_iterations"), 3L)
}

info(sprintf("Scoring with method: %s", score_method))
score <- do.call(fit$score, c(list(method = score_method), score_args))
rules <- score$rules(p_thresh = p_thresh, adjust_p = adjust_p, targets = targets)

correction <- NULL
if (isTRUE(correct_enabled)) {
  if (!nrow(rules)) {
    warning("No correction rules were called; skipping correction")
  } else {
    info("Applying correction rules")
    correction <- score$correct(
      rules = rules,
      name = as_chr(opt("correction_name"), paste0(score$name, "_clean")),
      p_thresh = p_thresh,
      adjust_p = adjust_p,
      targets = targets
    )
    print(correction$summary())
  }
}

top_genes <- factor_top_genes(fit, n = as_int(opt("top_genes"), 10L))
safe_write_csv(fit$summary(), file.path(output_dir, "fit_summary.csv"))
safe_write_csv(top_genes, file.path(output_dir, "factor_top_genes.csv"))
safe_write_csv(score$summary(), file.path(output_dir, "score_summary.csv"))
safe_write_csv(rules, file.path(output_dir, "score_rules.csv"))
if (!is.null(correction)) {
  safe_write_csv(correction$summary(), file.path(output_dir, "correction_summary.csv"))
}

summary <- list(
  input = normalizePath(input, winslash = "/", mustWork = FALSE),
  output_dir = output_dir,
  format = ds$format,
  annotation = ds$active_annotation %||% NA_character_,
  rank = fit$rank,
  nmf_variant = nmf_variant,
  nmf_n_runs = nmf_n_runs,
  score_method = score_method,
  score_name = score$name,
  n_rules = nrow(rules),
  corrected = !is.null(correction),
  correction_name = if (!is.null(correction)) correction$name else NA_character_,
  threads = threads
)
write_json_flat(summary, file.path(output_dir, "batch_summary.json"))
safe_write_csv(as.data.frame(summary, stringsAsFactors = FALSE), file.path(output_dir, "batch_summary.csv"))

if (isTRUE(report_enabled)) {
  report_file <- as_chr(opt("report_file"), file.path(output_dir, "report_minimal.html"))
  if (!grepl("\\.html?$", report_file, ignore.case = TRUE)) {
    stop("--report-file must end in .html")
  }
  if (!grepl("^/", report_file)) {
    report_file <- file.path(output_dir, report_file)
  }
  score_pairs_height <- max(6, ceiling(fit$rank / 2) * 4.9)
  env <- new.env(parent = globalenv())
  env$fit <- fit
  env$score <- score
  env$rules <- rules
  env$correction <- correction
  env$p_thresh <- p_thresh
  env$top_genes <- top_genes
  render_minimal_report(report_file, env,
    top_genes_n = as_int(opt("top_genes"), 10L),
    score_pairs_height = score_pairs_height)
}

info("Batch run complete")
cat("Output directory:", output_dir, "\n")
