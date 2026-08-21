# Export the gene-split evaluation definitions and the "before" count matrix
# for the python-side generative model. Uses the audit harness verbatim
# (build_split_defs from analysis/audit_guided/01_split_eval.R) so that the
# counts, pair pools, A/B splits and exposure values are identical to those
# the regression corrector was evaluated on.
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(cellAdmixCore))
suppressMessages(library(Matrix))
ag_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/audit_guided"
gm_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
source(file.path(ag_dir, "01_split_eval.R"))

dataset <- Sys.getenv("GM_DATASET", "pancreas")
root <- "/home/pkharchenko/cellAdmix/cellAdmix-core"
if (dataset == "pancreas") {
  setwd(file.path(root, "examples/xenium_pancreas_membrane_377_full"))
  annotation <- read.csv("annotations/annotation.csv.gz", stringsAsFactors = FALSE)
  cell_annotation <- setNames(annotation$merged_annotation, annotation$cell_id)
  ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)
  fit <- ds$fit(nmf_variant = "invsqrt_kl", verbose = FALSE)
} else if (dataset == "nsclc") {
  setwd(file.path(root, "examples/cosmx_nsclc_giotto"))
  meta_ann <- read.csv("prepared/cell_metadata_all.csv.gz", stringsAsFactors = FALSE)
  cell_annotation <- setNames(meta_ann$cell_type_coarse, meta_ann$cell)
  ds <- cellAdmix("prepared/molecules_all.csv.gz", output_dir = "out",
    annotation = cell_annotation)
  fit <- ds$fit(verbose = FALSE)
} else if (dataset == "breast") {
  setwd(file.path(root, "examples/xenium_breast_membrane_5k_full"))
  annotation <- read.csv("annotations/annotation.csv.gz", stringsAsFactors = FALSE)
  cell_annotation <- setNames(annotation$merged_annotation, annotation$cell_id)
  cell_annotation <- cell_annotation[cell_annotation != "Ambiguous / low-quality"]
  ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation,
    num_threads = 10)
  fit <- ds$fit(verbose = FALSE)
} else stop("unknown dataset")

defs <- build_split_defs(fit, cell_annotation)
message("pairs: ", length(defs$pairs))
check_disjoint(defs)

data_dir <- file.path(gm_dir,
  if (dataset == "pancreas") "data" else paste0("data_", dataset))
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# Full defs for the R-side evaluation script (03) - avoids a rebuild.
saveRDS(defs, file.path(data_dir, "defs.rds"))

# Count matrix (gene x cell) as MatrixMarket plus dimnames.
counts <- defs$counts_before
Matrix::writeMM(counts, file.path(data_dir, "counts.mtx"))
writeLines(rownames(counts), file.path(data_dir, "genes.txt"))
writeLines(colnames(counts), file.path(data_dir, "cells.txt"))

# Pair definitions with A/B splits and strict tiers.
meta <- do.call(rbind, lapply(names(defs$pairs), function(nm) {
  d <- defs$pairs[[nm]]
  data.frame(pair = nm, source = d$S, target = d$T_type,
    pool = paste(d$pool, collapse = ";"),
    A = paste(d$A, collapse = ";"),
    B = paste(d$B, collapse = ";"),
    strict = paste(d$strict, collapse = ";"),
    excess_pool = round(d$excess_pool), excess_B = round(d$excess_B),
    stringsAsFactors = FALSE)
}))
write.csv(meta, file.path(data_dir, "pairs.csv"), row.names = FALSE)

# Per-pair per-target-cell exposure (number of source cells among the 15
# nearest neighbors).
con <- gzfile(file.path(data_dir, "exposure.csv.gz"), "w")
writeLines("pair,cell_id,exposure", con)
for (nm in names(defs$pairs)) {
  d <- defs$pairs[[nm]]
  writeLines(sprintf("%s,%s,%d", nm, names(d$e), as.integer(d$e)), con)
}
close(con)

# Cell types, native markers, and type cpm profiles.
write.csv(data.frame(cell_id = names(defs$cell_types),
  cell_type = as.character(defs$cell_types), stringsAsFactors = FALSE),
  file.path(data_dir, "cell_types.csv"), row.names = FALSE)
native <- do.call(rbind, lapply(names(defs$native), function(t) {
  data.frame(cell_type = t, gene = defs$native[[t]],
    rank = seq_along(defs$native[[t]]), stringsAsFactors = FALSE)
}))
write.csv(native, file.path(data_dir, "native_markers.csv"), row.names = FALSE)
prof <- as.data.frame(defs$profiles)
prof$gene <- rownames(defs$profiles)
write.csv(prof, file.path(data_dir, "type_profiles_cpm.csv"), row.names = FALSE)

# Sanity: counts should match the run's molecule table (cell_idx >= 0).
message("total molecules in counts: ", sum(counts))
message("cells: ", ncol(counts), "  genes: ", nrow(counts))
message("EXPORT DONE")
