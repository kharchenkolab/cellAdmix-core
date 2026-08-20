# Ground-truth validation of the neighborhood-ladder reference.
#
# Synthetic tissue with a known amount of planted contamination in three
# components: contact admixture (proportional to in-plane source neighbors),
# hidden contamination (in target cells near the source region but with
# zero in-plane source neighbors - the analogue of out-of-section
# contamination), and uniform ambient background. The audit stores both the
# gradient-only excess (base zero-neighbor reference; blind to the hidden
# component) and the ladder-referenced excess; the test asks which recovers
# the true planted amount (contact + hidden; ambient excluded).
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(cellAdmixCore))
ag_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/audit_guided"

set.seed(21)
n_x <- 60; n_y <- 40; spacing <- 30
cells <- expand.grid(ix = seq_len(n_x), iy = seq_len(n_y))
cells$cell_id <- sprintf("c%04d", seq_len(nrow(cells)))
cells$cell_type <- ifelse(cells$ix <= 10, "A", "B")
cells$x <- cells$ix * spacing; cells$y <- cells$iy * spacing

genes_a <- paste0("a", 1:4); genes_b <- paste0("b", 1:4)
lam_ambient <- 0.6      # uniform background of source markers, per gene
hidden_scale <- 4       # hidden contamination at the source boundary
hidden_range <- 6       # columns over which hidden contamination decays
tx <- list()
truth_contact <- truth_hidden <- 0
for (i in seq_len(nrow(cells))) {
  ci <- cells[i, ]
  own <- if (ci$cell_type == "A") genes_a else genes_b
  g <- c(sample(own, 60L, replace = TRUE), rep("h1", 12L))
  if (ci$cell_type == "B") {
    dist_cols <- ci$ix - 10
    # contact admixture: only within ~1 grid step of the boundary
    n_contact <- if (dist_cols <= 1) 14L else 0L
    # hidden contamination: decays over hidden_range columns, present in
    # cells whose in-plane neighborhoods contain no A cells
    n_hidden <- rpois(1, hidden_scale * exp(-(dist_cols - 1) / hidden_range) *
      (dist_cols > 1))
    n_amb <- rpois(1, lam_ambient * length(genes_a))
    truth_contact <- truth_contact + n_contact
    truth_hidden <- truth_hidden + n_hidden
    add <- c(sample(genes_a, n_contact + n_hidden, replace = TRUE),
      sample(genes_a, n_amb, replace = TRUE))
    g <- c(g, add)
  }
  tx[[i]] <- data.frame(gene = g,
    x = ci$x + stats::runif(length(g), -10, 10),
    y = ci$y + stats::runif(length(g), -10, 10),
    z = 0, cell_id = ci$cell_id, cell_type = ci$cell_type,
    stringsAsFactors = FALSE)
}
tx <- do.call(rbind, tx)
csv_path <- tempfile("ladder_truth_", fileext = ".csv")
utils::write.csv(tx, csv_path, row.names = FALSE)

ds <- cellAdmix(csv_path, output_dir = tempfile("ladder_truth_out_"),
  schema = celladmix_schema(x = "x", y = "y", z = "z", gene = "gene",
    cell = "cell_id", cell_type = "cell_type"),
  annotation = stats::setNames(cells$cell_type, cells$cell_id),
  annotation_name = "manual")
fit <- ds$fit(rank = 2L, run_id = "truth", ncv_k = 6L, graph_k = 4L,
  nmf_iterations = 20L, nmf_n_runs = 1L, nmf_train_max_rows = 300L,
  num_threads = 4L, seed = 21L)
audit <- fit$audit_admixture(min_target_cells = 100L,
  min_reference_cells = 50L, min_excess = 50)
p <- audit$pairs()
ab <- p[p$source == "A" & p$target == "B", ]

# Gradient-only estimate under the base zero-neighbor reference, computed
# from the same marker panel with the harness helpers - what the audit
# would report without the ladder (blind to the hidden component).
pool <- audit$markers("A", "B")$pool
types <- stats::setNames(cells$cell_type, cells$cell_id)
cf <- fit$cell_factors()
expo <- cellAdmixCore:::.celladmix_source_exposure_counts(cf, types, 15L)
rownames(expo$counts) <- as.character(cf$cell_id)
counts <- fit$counts()
T_cells <- intersect(names(types)[types == "B"], colnames(counts))
totals <- Matrix::colSums(counts)
bins <- cellAdmixCore:::.celladmix_audit_bins(expo$counts[T_cells, "A"])
rates <- cellAdmixCore:::.celladmix_audit_bin_rates(
  cellAdmixCore:::.celladmix_audit_mcount(counts, pool, T_cells),
  totals[T_cells], bins)
grad <- cellAdmixCore:::.celladmix_audit_excess(rates)

truth_total <- truth_contact + truth_hidden
res <- data.frame(truth_contact = truth_contact, truth_hidden = truth_hidden,
  truth_total = truth_total,
  excess_gradient_base = round(grad$excess),
  excess_ladder = round(ab$excess),
  reference_kind = ab$reference_kind,
  reference_inflation = round(ab$reference_inflation, 2),
  recovery_gradient = round(grad$excess / truth_total, 3),
  recovery_ladder = round(ab$excess / truth_total, 3))
write.csv(res, file.path(ag_dir, "results", "ladder_truth.csv"),
  row.names = FALSE)
print(res, row.names = FALSE)
cat("LADDER TRUTH DONE\n")
