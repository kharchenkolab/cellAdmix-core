# Ground-truth tests of the dose prior's two estimation modes.
#
# Synthetic tissue with known planted contamination, run through the package
# (tabular input -> dataset audit -> generative fit -> correction), under the
# four combinations of dose_mode (marginal / joint) and dose_ladder (off /
# 30-60-120-240), plus the built-in covariate-permutation control.
#
# Scenario A - hidden contamination beyond the exposure window. Source A on
#   the left, target B elsewhere. B cells receive A-gene molecules from three
#   planted components: contact admixture (cells adjacent to the boundary),
#   hidden contamination decaying over several columns (the analogue of
#   out-of-section neighbors - present in cells with no A cell among their
#   15 nearest), and a uniform ambient background. Two of A's four genes are
#   also expressed natively by B at a low level, so they are pool genes but
#   not strict ones: the ambient tier (strict genes only) cannot reach their
#   contamination, and only the dose prior can. The question: how much of
#   the hidden component the correction removes with and without the ladder,
#   on the strict and on the shared genes, without removing native content.
#
# Scenario B - a co-located source sharing a gene. Pure A on the left, a
#   mixed zone where each cell is A or C at random, targets B on the right.
#   C transfers its profile - which carries a gene that A owns - in
#   proportion to the number of C cells among the 15 nearest; A transfers a
#   small amount of its own profile in proportion to A exposure. A exposure
#   co-varies with C's transfer across columns, so the shared gene's content
#   rises with A exposure beyond what A delivered. The question: how much
#   contamination each pair is credited with against its planted truth.
#
# Outputs: results/dose_truth_A.csv, results/dose_truth_B.csv,
# figures/dose_truth.png.
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages({library(cellAdmixCore); library(Matrix); library(ggplot2)})
gm_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/generative_model"
scratch <- Sys.getenv("GM_SCRATCH", tempdir())
ladder <- c(30L, 60L, 120L, 240L)
arms <- list(
  marginal = list(dose_mode = "marginal"),
  joint = list(dose_mode = "joint"),
  marginal_ladder = list(dose_mode = "marginal", dose_ladder = ladder),
  joint_ladder = list(dose_mode = "joint", dose_ladder = ladder))

knn_count <- function(x, y, is_type, k = 15L) {
  n <- length(x)
  out <- integer(n)
  for (i in seq_len(n)) {
    d <- (x - x[i])^2 + (y - y[i])^2
    d[i] <- Inf
    nn <- order(d)[seq_len(k)]
    out[i] <- sum(is_type[nn])
  }
  out
}

build_dataset <- function(tx, cells, tag) {
  csv_path <- file.path(scratch, paste0("dose_truth_", tag, ".csv"))
  utils::write.csv(tx, csv_path, row.names = FALSE)
  cellAdmix(csv_path, output_dir = file.path(scratch, paste0("dose_truth_out_", tag)),
    schema = celladmix_schema(x = "x", y = "y", z = "z", gene = "gene",
      cell = "cell_id", cell_type = "cell_type"),
    annotation = stats::setNames(cells$cell_type, cells$cell_id),
    annotation_name = "manual")
}

fit_arms <- function(audit, extra = list()) {
  lapply(names(arms), function(a) {
    opts <- c(list(init = "clusters", num_threads = 4L), arms[[a]], extra)
    t0 <- Sys.time()
    model <- do.call(audit$fit_generative, opts)
    list(arm = a, model = model, after = model$correct()$counts(),
      seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  })
}

# ---- Scenario A --------------------------------------------------------------
set.seed(21)
n_x <- 60; n_y <- 40; spacing <- 30
cells <- expand.grid(ix = seq_len(n_x), iy = seq_len(n_y))
cells$cell_id <- sprintf("c%04d", seq_len(nrow(cells)))
cells$cell_type <- ifelse(cells$ix <= 10, "A", "B")
cells$x <- cells$ix * spacing; cells$y <- cells$iy * spacing
genes_a <- paste0("a", 1:4); genes_b <- paste0("b", 1:4)
strict_a <- c("a1", "a2"); shared_a <- c("a3", "a4")   # B expresses a3, a4 natively
lam_ambient <- 0.15; hidden_scale <- 6; hidden_range <- 6
cells$contact <- 0L; cells$hidden <- 0L; cells$ambient <- 0L
cells$cont_strict <- 0L; cells$cont_shared <- 0L
tx <- vector("list", nrow(cells))
for (i in seq_len(nrow(cells))) {
  ci <- cells[i, ]
  g <- if (ci$cell_type == "A") sample(genes_a, 60L, replace = TRUE) else
    c(sample(genes_b, 56L, replace = TRUE), rep(shared_a, each = 2L))
  g <- c(g, rep("h1", 12L))
  if (ci$cell_type == "B") {
    dist_cols <- ci$ix - 10
    n_contact <- if (dist_cols <= 1) 14L else 0L
    n_hidden <- rpois(1, hidden_scale * exp(-(dist_cols - 1) / hidden_range) *
      (dist_cols > 1))
    n_amb <- rpois(1, lam_ambient * length(genes_a))
    cells$contact[i] <- n_contact; cells$hidden[i] <- n_hidden
    cells$ambient[i] <- n_amb
    add <- sample(genes_a, n_contact + n_hidden + n_amb, replace = TRUE)
    cells$cont_strict[i] <- sum(add %in% strict_a)
    cells$cont_shared[i] <- sum(add %in% shared_a)
    g <- c(g, add)
  }
  tx[[i]] <- data.frame(gene = g,
    x = ci$x + stats::runif(length(g), -10, 10),
    y = ci$y + stats::runif(length(g), -10, 10),
    z = 0, cell_id = ci$cell_id, cell_type = ci$cell_type,
    stringsAsFactors = FALSE)
}
tx <- do.call(rbind, tx)
ds <- build_dataset(tx, cells, "A")
audit <- ds$audit_admixture(min_target_cells = 100L, min_reference_cells = 50L,
  min_excess = 50)
before <- ds$counts()
b_cells <- cells$cell_id[cells$cell_type == "B"]
cells$exposure <- knn_count(cells$x, cells$y, cells$cell_type == "A")
res_a <- fit_arms(audit)
ctrl <- do.call(audit$fit_generative, c(list(init = "clusters", num_threads = 4L),
  arms$joint_ladder, list(use_ambient = FALSE, permute_covariates_seed = 11L)))
res_a[[length(res_a) + 1]] <- list(arm = "joint_ladder_control", model = ctrl,
  after = ctrl$correct()$counts(), seconds = NA_real_)

rows_a <- list()
for (r in res_a) {
  bc <- cells[match(b_cells, cells$cell_id), ]
  bc$removed_strict <- as.numeric(Matrix::colSums(before[strict_a, b_cells]) -
    Matrix::colSums(r$after[strict_a, b_cells]))
  bc$removed_shared <- as.numeric(Matrix::colSums(before[shared_a, b_cells]) -
    Matrix::colSums(r$after[shared_a, b_cells]))
  own_rem <- sum(before[genes_b, b_cells]) - sum(r$after[genes_b, b_cells])
  bc$dist <- bc$ix - 10
  bc$native_shared <- 4L
  by_col <- aggregate(cbind(removed_strict, removed_shared, cont_strict, cont_shared,
    native_shared, contact, hidden, ambient) ~ dist, bc, sum)
  pr <- r$model$pairs
  rows_a[[r$arm]] <- data.frame(arm = r$arm, dist = by_col$dist,
    removed_strict = by_col$removed_strict, removed_shared = by_col$removed_shared,
    cont_strict = by_col$cont_strict, cont_shared = by_col$cont_shared,
    native_shared = by_col$native_shared,
    contact = by_col$contact, hidden = by_col$hidden, ambient = by_col$ambient,
    exposed_cells = tapply(bc$exposure > 0, bc$dist, sum)[as.character(by_col$dist)],
    own_removed = own_rem,
    prior_unexposed = pr$prior_unexposed_molecules[1],
    seconds = r$seconds)
}
res_tab_a <- do.call(rbind, rows_a)
write.csv(res_tab_a, file.path(gm_dir, "results", "dose_truth_A.csv"), row.names = FALSE)
sum_a <- do.call(rbind, lapply(split(res_tab_a, res_tab_a$arm), function(d) {
  near <- d$dist <= 2; far <- d$dist > 2 & d$dist <= 9; deep <- d$dist > 9
  data.frame(arm = d$arm[1],
    strict_removed_exposed = sum(d$removed_strict[near]) / sum(d$cont_strict[near]),
    strict_removed_hidden_zone = sum(d$removed_strict[far]) / sum(d$cont_strict[far]),
    shared_removed_exposed = sum(d$removed_shared[near]) / sum(d$cont_shared[near]),
    shared_removed_hidden_zone = sum(d$removed_shared[far]) / sum(d$cont_shared[far]),
    shared_native_removed_deep = sum(d$removed_shared[deep]) / sum(d$native_shared[deep]),
    hidden_planted = sum(d$hidden[far]),
    own_removed = d$own_removed[1], prior_unexposed = d$prior_unexposed[1],
    seconds = d$seconds[1])
}))
print(sum_a, row.names = FALSE)
write.csv(sum_a, file.path(gm_dir, "results", "dose_truth_A_summary.csv"), row.names = FALSE)

# ---- Scenario B --------------------------------------------------------------
# Two variants of what C shares with A: a single gene (a4, which the audit's
# screen flags as disproportionate to A's profile and drops from A's guide
# panel), or a proportional slice of A's whole profile (a1..a4 in equal
# measure - nothing for the screen to flag, so the marginal fit credits A
# with content C delivered).
scenario_b <- function(variant, seed) {
  set.seed(seed)
  cells2 <- expand.grid(ix = seq_len(n_x), iy = seq_len(n_y))
  cells2$cell_id <- sprintf("d%04d", seq_len(nrow(cells2)))
  cells2$cell_type <- ifelse(cells2$ix <= 8, "A",
    ifelse(cells2$ix <= 12, ifelse(stats::runif(nrow(cells2)) < 0.5, "A", "C"), "B"))
  cells2$x <- cells2$ix * spacing; cells2$y <- cells2$iy * spacing
  genes_c <- paste0("c", 1:3)
  cells2$c_exposure <- knn_count(cells2$x, cells2$y, cells2$cell_type == "C")
  cells2$a_exposure <- knn_count(cells2$x, cells2$y, cells2$cell_type == "A")
  cells2$planted_c <- 0L; cells2$planted_a <- 0L; cells2$planted_a_from_c <- 0L
  c_own <- if (variant == "gene") c(rep(genes_c, each = 16L), rep("a4", 12L)) else
    c(rep(genes_c, each = 12L), rep(genes_a, each = 6L))
  c_genes <- unique(c_own)
  c_prob <- as.numeric(table(factor(c_own, levels = c_genes)))
  tx2 <- vector("list", nrow(cells2))
  for (i in seq_len(nrow(cells2))) {
    ci <- cells2[i, ]
    g <- switch(ci$cell_type,
      A = sample(genes_a, 60L, replace = TRUE),
      B = sample(genes_b, 60L, replace = TRUE),
      C = sample(c_genes, 60L, replace = TRUE, prob = c_prob))
    g <- c(g, rep("h1", 12L))
    if (ci$cell_type == "B") {
      # transfer from C samples C's profile, scaled by the number of C cells
      # among the 15 nearest; a smaller genuine transfer from A samples A's
      # profile, scaled by A exposure
      n_c <- rpois(1, 3.0 * ci$c_exposure)
      n_a <- rpois(1, 0.6 * ci$a_exposure)
      add_c <- sample(c_genes, n_c, replace = TRUE, prob = c_prob)
      add_a <- sample(genes_a, n_a, replace = TRUE)
      cells2$planted_c[i] <- n_c; cells2$planted_a[i] <- n_a
      cells2$planted_a_from_c[i] <- sum(add_c %in% genes_a)
      g <- c(g, add_c, add_a)
    }
    tx2[[i]] <- data.frame(gene = g,
      x = ci$x + stats::runif(length(g), -10, 10),
      y = ci$y + stats::runif(length(g), -10, 10),
      z = 0, cell_id = ci$cell_id, cell_type = ci$cell_type,
      stringsAsFactors = FALSE)
  }
  tx2 <- do.call(rbind, tx2)
  ds2 <- build_dataset(tx2, cells2, paste0("B_", variant))
  audit2 <- ds2$audit_admixture(min_target_cells = 100L, min_reference_cells = 50L,
    min_excess = 50)
  cat("variant:", variant, "\n")
  print(audit2$pairs(detected_only = TRUE)[, c("source", "target", "admixed_molecules", "rate")],
    row.names = FALSE)
  cat("A -> B guide panel:", paste(audit2$markers("A", "B")$pool, collapse = " "), "\n")
  before2 <- ds2$counts()
  b2 <- cells2$cell_id[cells2$cell_type == "B"]
  res_b <- fit_arms(audit2)
  rows_b <- list()
  for (r in res_b) {
    pr <- r$model$pairs
    a_row <- pr[pr$source == "A" & pr$target == "B", ]
    c_row <- pr[pr$source == "C" & pr$target == "B", ]
    rows_b[[r$arm]] <- data.frame(variant = variant, arm = r$arm,
      A_prior = if (nrow(a_row)) a_row$estimated_dose_molecules else 0,
      A_posterior = if (nrow(a_row)) a_row$removed_molecules_expected else 0,
      A_cross_attributed = if (nrow(a_row)) a_row$cross_attributed_molecules else 0,
      A_planted = sum(cells2$planted_a),
      C_prior = if (nrow(c_row)) c_row$estimated_dose_molecules else 0,
      C_posterior = if (nrow(c_row)) c_row$removed_molecules_expected else 0,
      C_planted = sum(cells2$planted_c),
      a_genes_from_c = sum(cells2$planted_a_from_c),
      removed_c_genes = sum(before2[genes_c, b2]) - sum(r$after[genes_c, b2]),
      removed_a_genes = sum(before2[genes_a, b2]) - sum(r$after[genes_a, b2]),
      own_removed = sum(before2[genes_b, b2]) - sum(r$after[genes_b, b2]),
      seconds = r$seconds)
  }
  do.call(rbind, rows_b)
}
res_tab_b <- rbind(scenario_b("gene", 22), scenario_b("profile", 23))
print(res_tab_b, row.names = FALSE)
write.csv(res_tab_b, file.path(gm_dir, "results", "dose_truth_B.csv"), row.names = FALSE)

# ---- Figure ------------------------------------------------------------------
pa <- res_tab_a[res_tab_a$arm != "joint_ladder_control" & res_tab_a$dist <= 14, ]
pa <- rbind(
  data.frame(arm = pa$arm, dist = pa$dist, genes = "strict genes (a1, a2)",
    fraction = pa$removed_strict / pmax(pa$cont_strict, 1)),
  data.frame(arm = pa$arm, dist = pa$dist, genes = "shared genes (a3, a4)",
    fraction = pa$removed_shared / pmax(pa$cont_shared, 1)))
pa$arm <- factor(pa$arm, levels = names(arms))
g1 <- ggplot(pa, aes(dist, fraction, color = arm)) +
  geom_vline(xintercept = 2.5, linetype = "22", color = "grey60") +
  geom_hline(yintercept = 1, linetype = "13", color = "grey60") +
  geom_line() + geom_point(size = 1.2) +
  facet_wrap(~genes, ncol = 1) +
  scale_color_manual(values = c(marginal = "#7f8c8d", joint = "#2980b9",
    marginal_ladder = "#e67e22", joint_ladder = "#1e8449"), name = NULL) +
  labs(x = "target column beyond the source boundary",
    y = "removed / planted contamination",
    title = "(a) hidden contamination beyond the exposure window") +
  theme_classic(base_size = 10) + theme(legend.position = "bottom")
pb <- res_tab_b
pb$arm <- factor(pb$arm, levels = names(arms))
pb$variant <- factor(ifelse(pb$variant == "gene", "C shares one A gene",
  "C shares a slice of A's profile"), levels = c("C shares one A gene",
  "C shares a slice of A's profile"))
pb_long <- rbind(data.frame(arm = pb$arm, variant = pb$variant, pair = "A -> B",
  ratio = pb$A_posterior / pb$A_planted),
  data.frame(arm = pb$arm, variant = pb$variant, pair = "C -> B",
    ratio = pb$C_posterior / pb$C_planted))
g2 <- ggplot(pb_long, aes(arm, ratio, fill = pair)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = 1, linetype = "13", color = "grey40") +
  facet_wrap(~variant, ncol = 1) +
  scale_fill_manual(values = c("#c0392b", "#2980b9"), name = NULL) +
  labs(x = NULL, y = "attributed contamination / planted transfer",
    title = "(b) a co-located source sharing genes with A") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(gm_dir, "figures", "dose_truth.png"),
  cowplot::plot_grid(g1, g2, ncol = 2, rel_widths = c(1.1, 1)),
  width = 10.5, height = 6.2, dpi = 150)
cat("DOSE TRUTH DONE\n")
