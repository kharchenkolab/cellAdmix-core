# Figure: how the audit's ambient reference is constructed and what it
# changes. Left: the reference ladder for one pancreas pair - the marker
# rate among target cells with zero source neighbors as the neighborhood
# grows. Right: the same pair's exposure profile with the zero-neighbor
# rate (dotted) and the chosen ambient reference (dashed).
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(cellAdmixCore)); suppressMessages(library(ggplot2))
setwd("/home/pkharchenko/cellAdmix/cellAdmix-core/examples/xenium_pancreas_membrane_377_full")
annotation <- read.csv("annotations/annotation.csv.gz", stringsAsFactors = FALSE)
ds <- cellAdmix("data", output_dir = "out",
  annotation = setNames(annotation$merged_annotation, annotation$cell_id))
fit <- ds$fit(nmf_variant = "invsqrt_kl", verbose = FALSE)
audit <- fit$audit_admixture()
p <- audit$pairs(detected_only = TRUE)
p <- p[p$reference_kind == "k240", ]
sel <- p[which.max(p$reference_inflation), ]
message("selected pair: ", sel$source, " -> ", sel$target,
  " (inflation ", round(sel$reference_inflation, 2), ")")
p1 <- audit$plot_reference(sel$source, sel$target) +
  labs(title = sprintf("(a) %s \u2192 %s: reference ladder", sel$source, sel$target),
    subtitle = NULL, x = "neighborhood size K") +
  theme(legend.position = "none")
p2 <- audit$plot_exposure(sel$source, sel$target, strict = FALSE) +
  labs(title = "(b) exposure profile (log scale)",
    subtitle = "dotted: zero-neighbor rate; dashed: ambient reference") +
  scale_y_log10()
fig <- cowplot::plot_grid(p1, p2, ncol = 2, align = "hv", axis = "tblr")
ggsave("/home/pkharchenko/cellAdmix/cellAdmix-core/docs/figures/benchmark_fig5.png",
  fig, width = 10.5, height = 3.8, dpi = 150)
message("FIG REFERENCE DONE")
