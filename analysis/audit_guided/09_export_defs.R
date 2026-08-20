# Export pair definitions (pools, strict tiers, per-cell exposure) for the
# molecule-level diagnostics in python.
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(cellAdmixCore))
ag_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/audit_guided"
source(file.path(ag_dir, "01_split_eval.R"))

setwd("/home/pkharchenko/cellAdmix/cellAdmix-core/examples/xenium_pancreas_membrane_377_full")
annotation <- read.csv("annotations/annotation.csv.gz", stringsAsFactors = FALSE)
cell_annotation <- setNames(annotation$merged_annotation, annotation$cell_id)
ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)
fit <- ds$fit(nmf_variant = "invsqrt_kl", verbose = FALSE)
defs <- build_split_defs(fit, cell_annotation)

meta <- do.call(rbind, lapply(names(defs$pairs), function(nm) {
  d <- defs$pairs[[nm]]
  data.frame(pair = nm, source = d$S, target = d$T_type,
    pool = paste(d$pool, collapse = ";"),
    strict = paste(d$strict, collapse = ";"),
    stringsAsFactors = FALSE)
}))
write.csv(meta, file.path(ag_dir, "results", "pair_meta.csv"), row.names = FALSE)

con <- gzfile(file.path(ag_dir, "results", "pair_exposure.csv.gz"), "w")
writeLines("pair,cell_id,exposure", con)
for (nm in names(defs$pairs)) {
  d <- defs$pairs[[nm]]
  writeLines(sprintf("%s,%s,%d", nm, names(d$e), as.integer(d$e)), con)
}
close(con)

ann_out <- data.frame(cell_id = names(defs$cell_types),
  cell_type = as.character(defs$cell_types), stringsAsFactors = FALSE)
write.csv(ann_out, file.path(ag_dir, "results", "cell_types.csv"), row.names = FALSE)
message("EXPORT DONE: ", nrow(meta), " pairs")
