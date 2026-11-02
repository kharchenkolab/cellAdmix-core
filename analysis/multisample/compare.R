# Compare the processing arms of a multi-sample dataset.
#
# Arms are discovered under <runs_root>/<dataset>/: one directory per sample
# (separate processing, run_separate.R), `consensus`, `joint`,
# `joint_profiles`, and `reference_<sample>` (run_collection.R / run_joint.R).
# For every arm and sample the script collects the audit's cleanup
# verification, the removal fractions, and the rule sets, and measures the
# divergence between the samples' corrected expression: per cell type, the
# correlation of pseudobulk profiles between samples, and the number of genes
# differentially expressed between the samples' spatial-block pseudobulks
# (Welch t-test, |log2 fold change| >= 0.585, BH FDR 0.05).
#
# Usage: Rscript compare.R <dataset> <runs_root> <results_dir> <sample>=<annotation.csv.gz> ...
#          [--groups=<name>:<sample>,<sample>,...;<name>:<sample>,...]
#
# The divergence is measured between two groups of samples (pooled counts,
# cells, and annotations): by default the first two samples; with --groups,
# the two named groups (for a cohort of tiles cut from two replicates, the
# tiles of each replicate).
.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/home/pkharchenko/cellAdmix/wt-dev_multisample/.Rlib"), .libPaths()))
suppressMessages(library(cellAdmixCore))
suppressMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("usage: compare.R <dataset> <runs_root> <results_dir> <sample>=<annotation> ...")
dataset <- args[[1]]
root <- file.path(args[[2]], dataset)
results_dir <- args[[3]]
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
rest <- args[-(1:3)]
groups_arg <- rest[startsWith(rest, "--groups=")]
rest <- rest[!startsWith(rest, "--groups=")]
specs <- strsplit(rest, "=", fixed = TRUE)
samples <- vapply(specs, `[[`, character(1), 1)
groups <- if (length(groups_arg)) {
  parts <- strsplit(strsplit(sub("^--groups=", "", groups_arg[[1]]), ";", fixed = TRUE)[[1]], ":", fixed = TRUE)
  stats::setNames(lapply(parts, function(x) strsplit(x[[2]], ",", fixed = TRUE)[[1]]), vapply(parts, `[[`, character(1), 1))
} else {
  stats::setNames(list(samples[[1]], samples[[2]]), samples[1:2])
}
if (length(groups) != 2L) stop("exactly two groups are compared")
annotation_paths <- stats::setNames(vapply(specs, `[[`, character(1), 2), samples)

read_annotation <- function(path) {
  ann <- read.csv(path, stringsAsFactors = FALSE)
  ann <- ann[!is.na(ann$cell_type) & ann$cell_type != "Unlabeled", ]
  stats::setNames(as.character(ann$cell_type), as.character(ann$cell_id))
}
annotations <- lapply(annotation_paths, read_annotation)
read_run_counts <- function(dir) {
  run <- cellAdmixCore:::celladmix_read_run(dir)
  cellAdmixCore:::celladmix_collect_counts_sparse(run)
}
read_run_cells <- function(dir) {
  cellAdmixCore:::celladmix_collect_cells(cellAdmixCore:::celladmix_read_run(dir))
}
strip_prefix <- function(ids, sample) sub(paste0("^", sample, "_"), "", ids)

# ---- locate arms --------------------------------------------------------------
# Each arm: per sample, the uncorrected run directory and the corrected run
# directory, and the arm's summary CSV directory.
arms <- list()
if (all(dir.exists(file.path(root, samples)))) {
  arms$separate <- lapply(samples, function(s) {
    p <- readRDS(file.path(root, s, "separate_pointers.rds"))
    list(run = p$run_dir, corrected = p$correction_dir, summary_dir = file.path(root, s), joint = FALSE)
  })
  names(arms$separate) <- samples
}
candidate_arms <- setdiff(basename(list.dirs(root, recursive = FALSE)), c(samples, "joint"))
for (arm in sort(candidate_arms)) {
  pointers_path <- file.path(root, arm, "pointers.csv")
  if (file.exists(pointers_path)) {
    pointers <- read.csv(pointers_path, stringsAsFactors = FALSE)
    arms[[arm]] <- lapply(samples, function(s) {
      row <- pointers[pointers$sample == s, , drop = FALSE]
      list(run = row$run_dir, corrected = row$correction_dir, summary_dir = file.path(root, arm), joint = FALSE)
    })
    names(arms[[arm]]) <- samples
  }
}
if (file.exists(file.path(root, "joint", "joint_pointers.rds"))) {
  p <- readRDS(file.path(root, "joint", "joint_pointers.rds"))
  arms$joint <- lapply(samples, function(s) {
    list(run = p$run_dir, corrected = p$correction_dir, summary_dir = file.path(root, "joint"), joint = TRUE)
  })
  names(arms$joint) <- samples
}
cat("arms found:", paste(names(arms), collapse = ", "), "\n")

# ---- counts per arm and sample --------------------------------------------------
# Joint runs hold all samples in one matrix; per-sample views are cut by the
# cell-id prefix and the prefix dropped so cells match the annotations.
load_counts <- function(entry, sample) {
  before <- read_run_counts(entry$run)
  after <- read_run_counts(entry$corrected)
  cells <- read_run_cells(entry$run)
  if (entry$joint) {
    keep <- startsWith(colnames(before), paste0(sample, "_"))
    before <- before[, keep, drop = FALSE]
    after <- after[, colnames(before), drop = FALSE]
    colnames(before) <- strip_prefix(colnames(before), sample)
    colnames(after) <- colnames(before)
    cells <- cells[startsWith(as.character(cells$cell_id), paste0(sample, "_")), , drop = FALSE]
    cells$cell_id <- strip_prefix(as.character(cells$cell_id), sample)
  } else {
    after <- after[, colnames(before), drop = FALSE]
  }
  list(before = before, after = after, cells = cells)
}

# ---- divergence between samples --------------------------------------------------
log_cpm <- function(v) log2(v / max(sum(v), 1) * 1e6 + 1)
pseudobulk_by_type <- function(counts, annotation) {
  ids <- intersect(colnames(counts), names(annotation))
  types <- annotation[ids]
  mat <- counts[, ids, drop = FALSE]
  groups <- split(seq_along(ids), types)
  sapply(groups, function(idx) rowSums(mat[, idx, drop = FALSE]))
}
block_pseudobulks <- function(counts, cells, annotation, type, block_um = 1000, min_cells = 30) {
  ids <- intersect(colnames(counts), names(annotation)[annotation == type])
  if (length(ids) < min_cells) return(NULL)
  pos <- cells[match(ids, as.character(cells$cell_id)), c("x", "y")]
  block <- paste(floor(pos$x / block_um), floor(pos$y / block_um))
  groups <- split(ids, block)
  groups <- groups[vapply(groups, length, integer(1)) >= min_cells]
  if (length(groups) < 3) return(NULL)
  sapply(groups, function(g) log_cpm(rowSums(counts[, g, drop = FALSE])))
}
# Pool the counts, cells, and annotations of a group of samples.
group_data <- function(members, counts_by_sample, which) {
  genes <- Reduce(intersect, lapply(members, function(s) rownames(counts_by_sample[[s]][[which]])))
  counts <- do.call(cbind, lapply(members, function(s) counts_by_sample[[s]][[which]][genes, , drop = FALSE]))
  cells <- do.call(rbind, lapply(members, function(s) {
    cells <- counts_by_sample[[s]]$cells
    cells[, intersect(c("cell_id", "x", "y"), names(cells)), drop = FALSE]
  }))
  annotation <- do.call(c, lapply(members, function(s) annotations[[s]]))
  list(counts = counts, cells = cells, annotation = annotation)
}
divergence_table <- function(arm, counts_by_sample, which = c("before", "after")) {
  which <- match.arg(which)
  ga <- group_data(groups[[1]], counts_by_sample, which)
  gb <- group_data(groups[[2]], counts_by_sample, which)
  a <- "a"; b <- "b"
  counts_by_sample <- list(a = list(cells = ga$cells), b = list(cells = gb$cells))
  counts_by_sample$a[[which]] <- ga$counts
  counts_by_sample$b[[which]] <- gb$counts
  annotations <- list(a = ga$annotation, b = gb$annotation)
  pb_a <- pseudobulk_by_type(counts_by_sample[[a]][[which]], annotations[[a]])
  pb_b <- pseudobulk_by_type(counts_by_sample[[b]][[which]], annotations[[b]])
  types <- intersect(colnames(pb_a), colnames(pb_b))
  genes <- intersect(rownames(pb_a), rownames(pb_b))
  do.call(rbind, lapply(types, function(t) {
    la <- log_cpm(pb_a[genes, t]); lb <- log_cpm(pb_b[genes, t])
    expressed <- (pb_a[genes, t] + pb_b[genes, t]) >= 20
    correlation <- if (sum(expressed) > 5) stats::cor(la[expressed], lb[expressed]) else NA_real_
    blocks_a <- block_pseudobulks(counts_by_sample[[a]][[which]], counts_by_sample[[a]]$cells, annotations[[a]], t)
    blocks_b <- block_pseudobulks(counts_by_sample[[b]][[which]], counts_by_sample[[b]]$cells, annotations[[b]], t)
    n_de <- NA_integer_; n_blocks_a <- 0L; n_blocks_b <- 0L
    if (!is.null(blocks_a) && !is.null(blocks_b)) {
      n_blocks_a <- ncol(blocks_a); n_blocks_b <- ncol(blocks_b)
      common <- intersect(rownames(blocks_a), rownames(blocks_b))
      blocks_a <- blocks_a[common, , drop = FALSE]; blocks_b <- blocks_b[common, , drop = FALSE]
      tested <- common[expressed[common]]
      if (length(tested) > 5) {
        p <- vapply(tested, function(g) {
          x <- blocks_a[g, ]; y <- blocks_b[g, ]
          if (stats::sd(x) == 0 && stats::sd(y) == 0) return(1)
          tryCatch(stats::t.test(x, y)$p.value, error = function(e) 1)
        }, numeric(1))
        lfc <- rowMeans(blocks_a[tested, , drop = FALSE]) - rowMeans(blocks_b[tested, , drop = FALSE])
        q <- stats::p.adjust(p, method = "BH")
        n_de <- sum(q < 0.05 & abs(lfc) >= log2(1.5))
      }
    }
    data.frame(arm = arm, stage = which, cell_type = t,
      cells_a = sum(annotations[[a]] == t), cells_b = sum(annotations[[b]] == t),
      correlation = correlation, n_de = n_de, n_blocks_a = n_blocks_a, n_blocks_b = n_blocks_b,
      stringsAsFactors = FALSE)
  }))
}

# ---- per-arm, per-sample summaries ---------------------------------------------
per_sample_rows <- list()
divergence_rows <- list()
rule_rows <- list()
for (arm in names(arms)) {
  cat("arm", arm, "\n")
  counts_by_sample <- lapply(samples, function(s) load_counts(arms[[arm]][[s]], s))
  names(counts_by_sample) <- samples
  for (s in samples) {
    entry <- arms[[arm]][[s]]
    evaluation <- read.csv(file.path(entry$summary_dir, "audit_evaluation.csv"), stringsAsFactors = FALSE)
    if ("sample" %in% names(evaluation)) evaluation <- evaluation[evaluation$sample == s, , drop = FALSE]
    correction <- read.csv(file.path(entry$summary_dir, "correction_summary.csv"), stringsAsFactors = FALSE)
    if ("sample" %in% names(correction)) correction <- correction[correction$sample == s, , drop = FALSE]
    all_row <- correction[correction$cell_type == "all", , drop = FALSE]
    rules <- read.csv(file.path(entry$summary_dir, "rules.csv"), stringsAsFactors = FALSE)
    kept <- rules[rules$keep %in% TRUE, , drop = FALSE]
    per_sample_rows[[length(per_sample_rows) + 1]] <- data.frame(
      arm = arm, sample = s,
      detected_pairs = evaluation$detected_pairs[[1]],
      leakage_removed_overall = evaluation$leakage_removed_overall[[1]],
      median_pair_sensitivity = evaluation$median_pair_sensitivity[[1]],
      own_marker_false_removal = evaluation$own_marker_false_removal[[1]],
      fraction_molecules_removed = all_row$fraction_molecules_removed[[1]],
      rules_kept = nrow(kept),
      rules_fully_supported = if ("sample_support" %in% names(kept)) sum(kept$sample_support >= 1) else NA_integer_,
      stringsAsFactors = FALSE)
  }
  # Rule-set agreement: kept source -> target pairs of the two samples (for
  # the separate arm these differ; shared arms agree by construction).
  pair_sets <- lapply(samples, function(s) {
    rules <- read.csv(file.path(arms[[arm]][[s]]$summary_dir, "rules.csv"), stringsAsFactors = FALSE)
    kept <- rules[rules$keep %in% TRUE, , drop = FALSE]
    unique(paste(kept$source_cell_type, kept$target_cell_type))
  })
  # Mean Jaccard over all sample pairs (one pair for two samples).
  pair_index <- utils::combn(length(pair_sets), 2)
  jaccard <- mean(apply(pair_index, 2, function(ij) {
    length(intersect(pair_sets[[ij[1]]], pair_sets[[ij[2]]])) /
      max(1, length(union(pair_sets[[ij[1]]], pair_sets[[ij[2]]])))
  }))
  rule_rows[[length(rule_rows) + 1]] <- data.frame(arm = arm, pairs_a = length(pair_sets[[1]]),
    pairs_b = length(pair_sets[[2]]), jaccard = jaccard, stringsAsFactors = FALSE)
  if (identical(arm, names(arms)[[1]])) {
    divergence_rows[[length(divergence_rows) + 1]] <- divergence_table("uncorrected", counts_by_sample, "before")
  }
  divergence_rows[[length(divergence_rows) + 1]] <- divergence_table(arm, counts_by_sample, "after")
}
per_sample <- do.call(rbind, per_sample_rows)
divergence <- do.call(rbind, divergence_rows)
rule_overlap <- do.call(rbind, rule_rows)
write.csv(per_sample, file.path(results_dir, paste0(dataset, "_arms_per_sample.csv")), row.names = FALSE)
write.csv(divergence, file.path(results_dir, paste0(dataset, "_divergence.csv")), row.names = FALSE)
write.csv(rule_overlap, file.path(results_dir, paste0(dataset, "_rule_overlap.csv")), row.names = FALSE)

# Divergence summary per arm: cell-weighted mean correlation and total DE genes
# over the types measured in every arm.
measured <- Reduce(intersect, lapply(split(divergence$cell_type, divergence$arm), unique))
summary <- do.call(rbind, lapply(split(divergence, divergence$arm), function(d) {
  d <- d[d$cell_type %in% measured, , drop = FALSE]
  w <- d$cells_a + d$cells_b
  data.frame(arm = d$arm[[1]], n_types = nrow(d),
    mean_correlation = sum(d$correlation * w, na.rm = TRUE) / sum(w[!is.na(d$correlation)]),
    total_de_genes = sum(d$n_de, na.rm = TRUE),
    types_with_de_test = sum(!is.na(d$n_de)), stringsAsFactors = FALSE)
}))
rownames(summary) <- NULL
write.csv(summary, file.path(results_dir, paste0(dataset, "_divergence_summary.csv")), row.names = FALSE)
print(per_sample, digits = 3)
print(rule_overlap, digits = 3)
print(summary, digits = 4)
