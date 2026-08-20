# E2: full-replacement audit-guided anchors on the pancreas invsqrt/bridge
# failure regime, head-to-head with the historical co-occurrence (SPA)
# anchor runs (bench_anchor_s1..s3), all under the A/B split harness.
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(cellAdmixCore))
ag_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/audit_guided"
source(file.path(ag_dir, "01_split_eval.R"))

setwd("/home/pkharchenko/cellAdmix/cellAdmix-core/examples/xenium_pancreas_membrane_377_full")
annotation <- read.csv("annotations/annotation.csv.gz", stringsAsFactors = FALSE)
cell_annotation <- setNames(annotation$merged_annotation, annotation$cell_id)
ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)
fit_inv <- ds$fit(nmf_variant = "invsqrt_kl", verbose = FALSE)

defs <- build_split_defs(fit_inv, cell_annotation)
profiles <- defs$profiles
types <- colnames(profiles)

margin_row <- function(S) {
  others <- setdiff(types, S)
  m <- pmax(profiles[, S] - apply(profiles[, others, drop = FALSE], 1, max), 0)
  m + 0.05 * profiles[, S]
}

H0 <- fit_inv$loadings()
scale_to <- stats::median(rowSums(H0))
gene_order <- rownames(profiles)
H_types <- do.call(rbind, lapply(types, function(S) {
  v <- margin_row(S)[gene_order]; v / sum(v) * scale_to
}))
colnames(H_types) <- gene_order

fit_g <- ds$fit(nmf_fixed_h = H_types, nmf_variant = "kl", nmf_n_runs = 1L,
  run_id = "ag_anchor_alltypes_inv", overwrite = TRUE, verbose = FALSE)
fs <- fit_g$score_factor_sources()
best <- fs$scores()[fs$scores()$is_best, c("factor", "cell_type", "score", "margin")]
message("audit-guided full-anchor factor alignment:")
print(best, row.names = FALSE)

eval_arm <- function(fit_obj, tag) {
  score <- fit_obj$score_bridge()
  corr <- suppressWarnings(fit_obj$correct(score, ensemble = 1,
    name = paste0("ag_e2_", gsub("[^a-z0-9]", "", tolower(tag)))))
  res <- eval_split(defs, corr$counts(), tag)
  write.csv(res$pairs, file.path(ag_dir, "results",
    sprintf("e2_%s_pairs.csv", gsub("[^a-z0-9_]", "", tolower(tag)))),
    row.names = FALSE)
  print(res$summary, row.names = FALSE)
  invisible(res)
}

message("\n== audit-guided type anchors (bridge) ==")
eval_arm(fit_g, "guided_types")

for (s in 1:3) {
  message("\n== historical SPA anchors s", s, " (bridge) ==")
  fit_h <- ds$read_fit(sprintf("out/runs/bench_anchor_s%d", s))
  eval_arm(fit_h, sprintf("spa_s%d", s))
}
message("E2 DONE")
