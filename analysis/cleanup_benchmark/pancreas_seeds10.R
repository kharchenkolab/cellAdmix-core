# Extend the pancreas seed ensemble to 10 external seed fits and materialize
# per-seed corrections for the molecule-vote threshold sweep.
.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/tmp/celladmix_r_lib"), .libPaths()))
suppressMessages(library(cellAdmixCore))

annotation_df <- read.csv("examples/xenium_pancreas_membrane_377_full/annotations/annotation.csv.gz",
  stringsAsFactors = FALSE)
cell_annotation <- setNames(annotation_df$merged_annotation, annotation_df$cell_id)
setwd("examples/xenium_pancreas_membrane_377_full")
ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)

args <- commandArgs(trailingOnly = TRUE)
shard_variant <- args[[1]]
shard_seeds <- eval(parse(text = args[[2]]))
for (variant in shard_variant) {
  for (s in shard_seeds) {
    run_id <- sprintf("bench_seed%d_%s", s, variant)
    fit <- ds$fit(nmf_variant = variant, seed = s, verbose = FALSE, run_id = run_id)
    for (method in c("membrane", "bridge")) {
      corr_name <- sprintf("cmp_%s_%s_s%d", method, variant, s)
      corr_dir <- file.path("out", "runs", run_id, "corrected", corr_name)
      if (dir.exists(corr_dir) && file.exists(file.path(corr_dir, "molecules.parquet"))) {
        message("exists: ", corr_name)
        next
      }
      score <- if (method == "membrane") fit$score_membrane() else fit$score_bridge()
      rules <- score$rules(p_thresh = 0.1)
      invisible(score$correct(rules = rules, name = corr_name))
      message("built: ", corr_name)
    }
  }
}
message("SEEDS10 DONE")
