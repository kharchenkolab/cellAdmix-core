# E1a: audit-guided Mural/pericyte anchor on the pancreas invsqrt fit,
# with a junk-anchor negative control and matched single-fit baselines.
#
# Anchor row = one-vs-rest positive margin of the source type's cpm profile
# (plus a small profile floor for shared genes): target-side blind, sharp,
# and covering B-half genes through the source's own expression. The A/B
# harness still validates on B halves. Fixed-H runs are single-fit, so the
# comparators here are single-fit corrections of the base fit.
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

margin_row <- function(S) {
  others <- setdiff(colnames(profiles), S)
  m <- pmax(profiles[, S] - apply(profiles[, others, drop = FALSE], 1, max), 0)
  m + 0.05 * profiles[, S]
}

# Junk control: an "anchor" built the same way from a fake profile - the
# abundance-matched unowned control genes get uniform weight.
junk_genes <- bench_control_genes(profiles, "Mural / pericyte",
  defs$pairs[["Mural / pericyte -> Endothelial"]]$A, n_per_marker = 2L)
junk_row <- setNames(numeric(nrow(profiles)), rownames(profiles))
junk_row[junk_genes] <- mean(margin_row("Mural / pericyte")[
  margin_row("Mural / pericyte") > 0])

H0 <- fit_inv$loadings()
gene_order <- rownames(profiles)
stopifnot(ncol(H0) == length(gene_order))
scale_to <- stats::median(rowSums(H0))
as_row <- function(v) { v <- v[gene_order]; v / sum(v) * scale_to }

anchored_fit <- function(row_vec, run_id) {
  H_aug <- rbind(H0, as_row(row_vec))
  colnames(H_aug) <- gene_order
  fit_a <- ds$fit(nmf_fixed_h = H_aug, nmf_variant = "kl", nmf_n_runs = 1L,
    run_id = run_id, overwrite = TRUE, verbose = FALSE)
  # The run reorders factors by mass - identify the anchor by loading match.
  Ha <- fit_a$loadings()
  sim <- apply(Ha, 1, function(r) stats::cor(r, as_row(row_vec)))
  anchor_f <- which.max(sim)
  message(run_id, ": anchor is factor ", anchor_f,
    sprintf(" (cor %.2f)", sim[anchor_f]))
  fs <- fit_a$score_factor_sources()
  b <- fs$scores()[fs$scores()$is_best & fs$scores()$factor == anchor_f, ]
  message("  anchor aligns to: ", b$cell_type, sprintf(" (score %.2f margin %.2f)",
    b$score, b$margin))
  list(fit = fit_a, anchor_f = anchor_f)
}

mural_watch <- c("Fibroblast / CAF -> Mural / pericyte",
  "Mural / pericyte -> Endothelial", "Mural / pericyte -> Fibroblast / CAF",
  "Mural / pericyte -> Immune")

eval_arm <- function(fit_obj, method, tag, single = TRUE, anchor_f = NULL) {
  score <- if (method == "membrane") fit_obj$score_membrane(verbose = FALSE)
    else fit_obj$score_bridge()
  rules <- score$rules(p_thresh = 0.1)
  if (!is.null(anchor_f)) {
    ar <- rules[rules$factor == anchor_f, , drop = FALSE]
    message("  [", tag, " ", method, "] anchor-factor rules: ", nrow(ar))
    if (nrow(ar)) print(ar[, c("factor", "source_cell_type", "target_cell_type",
      "p_value", "keep", "native_check")], row.names = FALSE)
  }
  corr <- suppressWarnings(if (single) fit_obj$correct(score, ensemble = 1,
    name = paste0("ag_", tag, "_", method)) else fit_obj$correct(score,
    name = paste0("ag_", tag, "_", method)))
  res <- eval_split(defs, corr$counts(), paste(tag, method))
  write.csv(res$pairs, file.path(ag_dir, "results",
    sprintf("e1b_%s_%s_pairs.csv", tag, method)), row.names = FALSE)
  print(res$summary, row.names = FALSE)
  print(res$pairs[res$pairs$pair %in% mural_watch,
    c("pair", "excess_B", "power_B", "depth_B")], row.names = FALSE)
  invisible(res)
}

message("\n==== base single-fit comparators ====")
for (m in c("membrane", "bridge")) eval_arm(fit_inv, m, "base_single")

message("\n==== mural anchor ====")
am <- anchored_fit(margin_row("Mural / pericyte"), "ag_anchor_mural_inv")
for (m in c("membrane", "bridge")) eval_arm(am$fit, m, "mural", anchor_f = am$anchor_f)

message("\n==== mural anchor + audit-attributed rules ====")
# The audit tells us the anchor factor's source identity and which targets
# receive mural admixture; bypass the scoring-level source call for the
# anchor factor only, keep auto-derived rules for everything else.
score_b <- am$fit$score_bridge()
auto <- score_b$rules(p_thresh = 0.1)
auto <- auto[auto$keep, , drop = FALSE]
guided <- do.call(rbind, lapply(c("Endothelial", "Fibroblast / CAF", "Immune"),
  function(tt) {
    row <- auto[1, , drop = FALSE]
    row$factor <- am$anchor_f
    row$source_cell_type <- "Mural / pericyte"
    row$target_cell_type <- tt
    row$keep <- TRUE
    row$native_check <- "audit_guided"
    row
  }))
rules_g <- rbind(auto, guided)
corr_g <- suppressWarnings(am$fit$correct(score_b, rules = rules_g,
  ensemble = 1, name = "ag_mural_guided"))
res_g <- eval_split(defs, corr_g$counts(), "mural guided bridge")
write.csv(res_g$pairs, file.path(ag_dir, "results",
  "e1b_mural_guided_bridge_pairs.csv"), row.names = FALSE)
print(res_g$summary, row.names = FALSE)
print(res_g$pairs[res_g$pairs$pair %in% mural_watch,
  c("pair", "excess_B", "power_B", "depth_B")], row.names = FALSE)

message("\n==== junk anchor ====")
aj <- anchored_fit(junk_row, "ag_anchor_junk_inv")
for (m in c("membrane", "bridge")) eval_arm(aj$fit, m, "junk", anchor_f = aj$anchor_f)

message("E1 MURAL DONE")
