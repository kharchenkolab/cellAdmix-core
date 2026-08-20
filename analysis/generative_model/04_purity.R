# Downstream check: same-label kNN purity of corrected counts on a 5,000-cell
# subset, matching the diagnostic in analysis/audit_guided/08_phase_a.R
# (function copied verbatim so the numbers are comparable; the audit reported
# original 0.804, membrane ensemble 0.968, exposure-regression v3 0.967).
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(Matrix))
gm_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
defs <- readRDS(file.path(gm_dir, "data", "defs.rds"))
counts0 <- defs$counts_before
totals0 <- defs$totals_before

knn_purity <- function(m, n_cells = 5000L, n_genes = 1000L, k = 15L, seed = 7L) {
  set.seed(seed)
  ct <- defs$cell_types[colnames(m)]
  keep <- names(ct)[!is.na(ct) & totals0[names(ct)] >= 30]
  cells <- sample(keep, min(n_cells, length(keep)))
  sub <- m[, cells, drop = FALSE]
  norm <- log1p(t(t(as.matrix(sub)) / pmax(Matrix::colSums(sub), 1)) * 1e4)
  vars <- apply(norm, 1, stats::var)
  norm <- norm[order(-vars)[seq_len(min(n_genes, nrow(norm)))], ]
  pc <- stats::prcomp(t(norm), rank. = 30, center = TRUE, scale. = FALSE)$x
  dm <- as.matrix(stats::dist(pc))
  diag(dm) <- Inf
  lab <- as.character(ct[cells])
  same <- vapply(seq_along(cells), function(i) {
    nn <- order(dm[i, ])[seq_len(k)]
    mean(lab[nn] == lab[i])
  }, 0)
  base <- sum((table(lab) / length(lab))^2)
  (mean(same) - base) / (1 - base)
}

rem <- as(Matrix::readMM(file.path(gm_dir, "data", "gm_removed_production.mtx")),
  "CsparseMatrix")
dimnames(rem) <- dimnames(counts0)
cp <- counts0 - rem
cp@x <- pmax(cp@x, 0)

p0 <- knn_purity(counts0)
pg <- knn_purity(cp)
message(sprintf("original %.3f | generative production %.3f", p0, pg))
writeLines(c("arm,knn_purity",
  sprintf("original,%.4f", p0),
  sprintf("gm_production,%.4f", pg),
  "ensemble_audit_reported,0.9680",
  "exposure_v3_audit_reported,0.9670"),
  file.path(gm_dir, "results", "gm_purity.csv"))
message("PURITY DONE")
