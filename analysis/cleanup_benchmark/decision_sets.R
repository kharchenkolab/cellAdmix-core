# Kept (source, target) removal decisions per seed, for the decision-overlap
# panel of the stochasticity figure (pancreas).
.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/tmp/celladmix_r_lib"), .libPaths()))
suppressMessages(library(cellAdmixCore))
ann <- read.csv("examples/xenium_pancreas_membrane_377_full/annotations/annotation.csv.gz")
ca <- setNames(ann$merged_annotation, ann$cell_id)
out_csv <- file.path(normalizePath(Sys.getenv("BENCH_DIR", "analysis/cleanup_benchmark")),
  "results", "pancreas_decisions.csv")
setwd("examples/xenium_pancreas_membrane_377_full")
ds <- cellAdmix("data", output_dir = "out", annotation = ca)
rows <- list()
for (variant in c("ls_nmf", "invsqrt_kl")) {
  for (s in 1:10) {
    fit <- ds$read_fit(sprintf("out/runs/bench_seed%d_%s", s, variant))
    for (method in c("membrane", "bridge")) {
      score <- if (method == "membrane") fit$score_membrane() else fit$score_bridge()
      r <- score$rules(p_thresh = 0.1)
      r <- r[r$keep, , drop = FALSE]
      if (nrow(r)) {
        rows[[length(rows) + 1]] <- data.frame(variant = variant, seed = s,
          method = method, source = r$source_cell_type, target = r$target_cell_type)
      }
      message(variant, " s", s, " ", method, ": ", nrow(r), " kept decisions")
    }
  }
}
write.csv(do.call(rbind, rows), out_csv, row.names = FALSE)
message("DECISIONS DONE")
