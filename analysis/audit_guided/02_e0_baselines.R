# E0: baseline corrections evaluated under the A/B gene split (pancreas).
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(cellAdmixCore))
ag_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/audit_guided"
source(file.path(ag_dir, "01_split_eval.R"))
dir.create(file.path(ag_dir, "results"), showWarnings = FALSE)

setwd("/home/pkharchenko/cellAdmix/cellAdmix-core/examples/xenium_pancreas_membrane_377_full")
annotation <- read.csv("annotations/annotation.csv.gz", stringsAsFactors = FALSE)
cell_annotation <- setNames(annotation$merged_annotation, annotation$cell_id)
ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)

fit_inv <- ds$fit(nmf_variant = "invsqrt_kl", verbose = FALSE)
fit_ls <- ds$fit(verbose = FALSE)

message("building split definitions ...")
defs <- build_split_defs(fit_inv, cell_annotation)
message("pairs: ", length(defs$pairs))
check_disjoint(defs)
saveRDS(defs[c("pairs", "native")], file.path(ag_dir, "results", "split_defs.rds"))

run_arm <- function(fit, method, arm_tag, name) {
  score <- switch(method,
    membrane = fit$score_membrane(verbose = FALSE),
    bridge = fit$score_bridge())
  corr <- suppressWarnings(fit$correct(score, name = name))
  res <- eval_split(defs, corr$counts(), arm_tag)
  write.csv(res$pairs, file.path(ag_dir, "results",
    sprintf("e0_%s_pairs.csv", arm_tag)), row.names = FALSE)
  write.csv(res$types, file.path(ag_dir, "results",
    sprintf("e0_%s_types.csv", arm_tag)), row.names = FALSE)
  print(res$summary, row.names = FALSE)
  res
}

message("\n== membrane ensemble on invsqrt ==")
m <- run_arm(fit_inv, "membrane", "membrane_inv", "ag_e0_membrane")
message("\n== bridge ensemble on ls ==")
b1 <- run_arm(fit_ls, "bridge", "bridge_ls", "ag_e0_bridge_ls")
message("\n== bridge ensemble on invsqrt ==")
b2 <- run_arm(fit_inv, "bridge", "bridge_inv", "ag_e0_bridge_inv")

message("\n-- per-pair power_B (membrane baseline), uncovered pairs highlighted --")
p <- m$pairs[order(-m$pairs$excess_B), ]
print(p[, c("pair", "excess_B", "power_A", "power_B", "power_strictB")], row.names = FALSE)
message("E0 DONE")
