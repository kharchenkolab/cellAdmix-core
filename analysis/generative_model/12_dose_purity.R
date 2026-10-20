# Downstream check for the dose-prior arms: same-label kNN purity of the
# corrected counts on a 5,000-cell subset (the diagnostic of 04_purity.R,
# same function and seed, so the numbers are comparable: original 0.804,
# generative production 0.963, membrane ensemble 0.968).
#
# Usage: Rscript 12_dose_purity.R [arm ...]; default: the production arms.
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(Matrix))
gm_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
dataset <- Sys.getenv("GM_DATASET", "pancreas")
data_dir <- file.path(gm_dir, if (dataset == "pancreas") "data" else paste0("data_", dataset))
prefix <- if (dataset == "pancreas") "" else paste0(dataset, "_")
defs <- readRDS(file.path(data_dir, "defs.rds"))
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

arms <- commandArgs(trailingOnly = TRUE)
if (!length(arms)) arms <- c("marginal_z0", "marginal", "joint", "ladder", "joint_ladder")
rows <- list(data.frame(arm = "original", knn_purity = knn_purity(counts0)))
for (arm in arms) {
  f <- file.path(data_dir, sprintf("gm_removed_dose_%s.mtx", arm))
  if (!file.exists(f)) { message("missing: ", f); next }
  rem <- as(Matrix::readMM(f), "CsparseMatrix")
  dimnames(rem) <- dimnames(counts0)
  cp <- counts0 - rem
  cp@x <- pmax(cp@x, 0)
  rows[[length(rows) + 1]] <- data.frame(arm = arm, knn_purity = knn_purity(cp))
  message(sprintf("%s: %.4f", arm, rows[[length(rows)]]$knn_purity))
}
out <- do.call(rbind, rows)
print(out, row.names = FALSE)
write.csv(out, file.path(gm_dir, "results", paste0(prefix, "gm_dose_purity.csv")),
  row.names = FALSE)
message("DOSE PURITY DONE")
