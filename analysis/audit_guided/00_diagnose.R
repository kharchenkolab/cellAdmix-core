# Why are the 8 stromal pairs uncovered by the pancreas membrane correction?
# Failure modes per pair: (a) no factor aligned with the source type,
# (b) rule attempted but insignificant, (c) rule vetoed by the native check.
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(cellAdmixCore))
setwd("/home/pkharchenko/cellAdmix/cellAdmix-core/examples/xenium_pancreas_membrane_377_full")

annotation <- read.csv("annotations/annotation.csv.gz", stringsAsFactors = FALSE)
cell_annotation <- setNames(annotation$merged_annotation, annotation$cell_id)
ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)

uncovered <- rbind(
  c("Exocrine epithelial", "Immune"),
  c("Fibroblast / CAF", "Endothelial"),
  c("Fibroblast / CAF", "Immune"),
  c("Fibroblast / CAF", "Mural / pericyte"),
  c("Immune", "Fibroblast / CAF"),
  c("Mural / pericyte", "Endothelial"),
  c("Mural / pericyte", "Fibroblast / CAF"),
  c("Mural / pericyte", "Immune"))
colnames(uncovered) <- c("S", "T")

report_fit <- function(fit, tag) {
  cat("\n=====", tag, "=====\n")
  fs <- fit$score_factor_sources()
  best <- fs$scores()[fs$scores()$is_best, c("factor", "cell_type", "score", "margin")]
  cat("factor -> best type alignment:\n")
  print(best, row.names = FALSE)
  types_with_factor <- unique(best$cell_type[best$margin > 0.05])
  cat("types lacking an aligned factor:",
    paste(setdiff(sort(unique(cell_annotation)), types_with_factor), collapse = ", "), "\n")
}

report_rules <- function(score, tag) {
  cat("\n---- rules status for uncovered pairs [", tag, "] ----\n")
  rules <- score$rules(p_thresh = 0.1)
  for (i in seq_len(nrow(uncovered))) {
    S <- uncovered[i, 1]; T_type <- uncovered[i, 2]
    hit <- rules[rules$source_cell_type == S & rules$target_cell_type == T_type, , drop = FALSE]
    if (!nrow(hit)) {
      cat(sprintf("%-22s -> %-18s : NO RULE ROW\n", S, T_type))
    } else {
      for (j in seq_len(nrow(hit))) {
        cat(sprintf("%-22s -> %-18s : f%-2d p=%.3g keep=%s (%s)\n", S, T_type,
          hit$factor[j], hit$p_value[j], hit$keep[j],
          if ("native_check" %in% names(hit)) hit$native_check[j] else "-"))
      }
    }
  }
  invisible(rules)
}

fit_inv <- ds$fit(nmf_variant = "invsqrt_kl", verbose = FALSE)
fit_ls <- ds$fit(verbose = FALSE)

report_fit(fit_inv, "invsqrt_kl rank9 (membrane fit)")
report_fit(fit_ls, "ls_nmf rank9 (bridge fit)")

memb <- fit_inv$score_membrane(verbose = FALSE)
report_rules(memb, "membrane on invsqrt")

br_ls <- fit_ls$score_bridge()
report_rules(br_ls, "bridge on ls")

br_inv <- fit_inv$score_bridge()
report_rules(br_inv, "bridge on invsqrt")
cat("\nDIAGNOSE DONE\n")
