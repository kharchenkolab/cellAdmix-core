# Build per-seed corrections for the molecule-vote ensemble.
# Usage: Rscript analysis/cleanup_benchmark/seeds10.R <dataset> <variant> <seed_expr>
#   e.g. seeds10.R breast_crop ls_nmf 4:7
.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/tmp/celladmix_r_lib"), .libPaths()))
suppressMessages(library(cellAdmixCore))

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[[1]]
variant <- args[[2]]
seeds <- eval(parse(text = args[[3]]))

cfg <- switch(dataset,
  pancreas = list(dir = "examples/xenium_pancreas_membrane_377_full",
    ann = "annotations/annotation.csv.gz", ann_col = "merged_annotation",
    ann_id = "cell_id", out = "out", methods = c("membrane", "bridge"),
    build = function(ca) cellAdmix("data", output_dir = "out", annotation = ca)),
  breast_crop = list(dir = "examples/xenium_breast_membrane_5k_full",
    ann = "annotations/annotation.csv.gz", ann_col = "merged_annotation",
    ann_id = "cell_id", out = "out_medium_crop", methods = c("membrane", "bridge"),
    build = function(ca) cellAdmix("data", output_dir = "out_medium_crop",
      annotation = ca[ca != "Ambiguous / low-quality"], num_threads = 10,
      analysis_bbox = c(5000, 7000, 4500, 6500))),
  nsclc = list(dir = "examples/cosmx_nsclc_giotto",
    ann = "prepared/cell_metadata_all.csv.gz", ann_col = "cell_type_coarse",
    ann_id = "cell", out = "out", methods = "bridge",
    build = function(ca) cellAdmix(file.path("prepared", "molecules_all.csv.gz"),
      output_dir = "out", annotation = ca)))

ann_df <- read.csv(file.path(cfg$dir, cfg$ann), stringsAsFactors = FALSE)
cell_annotation <- setNames(ann_df[[cfg$ann_col]], ann_df[[cfg$ann_id]])
setwd(cfg$dir)
ds <- cfg$build(cell_annotation)

for (s in seeds) {
  run_id <- sprintf("bench_seed%d_%s", s, variant)
  fit <- ds$fit(nmf_variant = variant, seed = s, verbose = FALSE, run_id = run_id)
  for (method in cfg$methods) {
    corr_name <- sprintf("cmp_%s_%s_s%d", method, variant, s)
    corr_dir <- file.path(cfg$out, "runs", run_id, "corrected", corr_name)
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
message("SEEDS DONE ", dataset, " ", variant)
