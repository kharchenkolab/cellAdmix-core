# Reusable analysis and plotting helpers for the generative admixture model.
#
# These functions take the package's public objects (the dataset, the audit,
# a fitted generative model, a correction, count matrices) and return
# tables or ggplot objects; they are the figure designs of docs/generative.md
# factored out of the example notebooks, kept outside the package until the
# model's place in it settles. Source this file from a notebook:
#   source("../../analysis/generative_model/tools/generative_plots.R")
#
# Conventions: `before`/`after` are gene x cell sparse matrices with cell
# ids as column names (ds$counts(), correction$counts()); `annotation` is a
# named character vector of cell types by cell id; `cells_xy` is ds$cells()
# (cell_id, x, y); `exposure` is a named integer vector of source-type
# neighbor counts per cell id (generative_exposure()).

# Short display name of a cell type: the part before a slash, lowercased,
# unless that part is too short to stand alone (e.g. "T / NK").
generative_short_label <- function(x, lower = TRUE) {
  s <- trimws(sub("/.*", "", x))
  keep_case <- nchar(sub(" .*", "", s)) <= 2
  s <- ifelse(lower & !keep_case, tolower(s), s)
  ifelse(nchar(s) < 3, as.character(x), s)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(Matrix)
})

generative_palette <- c("own expression" = "#2980b9", contamination = "#c0392b",
  ambient = "#8a8a8a", induced = "#1e8449")

# Number of source-type cells among each cell's neighbor_k nearest cells.
generative_exposure <- function(cells_xy, annotation, source, neighbor_k = 15L) {
  counts <- cellAdmixCore:::.celladmix_source_exposure_counts(
    cells_xy, annotation, as.integer(neighbor_k))$counts
  stats::setNames(counts[, source], as.character(cells_xy$cell_id))
}

# Which cell type expresses each gene most (counts per million over the
# type's cells together); the audit's gene-ownership rule.
generative_gene_owner <- function(before, annotation) {
  profiles <- cellAdmixCore:::.celladmix_audit_profiles(before,
    annotation[colnames(before)])
  list(profiles = profiles,
    owner = stats::setNames(colnames(profiles)[max.col(profiles, ties.method = "first")],
      rownames(profiles)))
}

# Expression profile of the source cells bordering the target type.
generative_interface_profile <- function(before, cells_xy, annotation, source, target) {
  S_cells <- names(annotation)[!is.na(annotation) & annotation == source]
  dist_T <- cellAdmixCore:::.celladmix_source_nearest_distance(
    cells_xy, annotation)[, target]
  cellAdmixCore:::.celladmix_audit_interface_profile(before, S_cells, dist_T)
}

# Exposure bins 0, 1, 2, 3+ as a factor.
generative_bins <- function(exposure) {
  factor(pmin(exposure, 3), levels = 0:3, labels = c("0", "1", "2", "3+"))
}

# Per-bin rate of one gene (per 1,000 molecules) with exact Poisson intervals.
generative_bin_rates <- function(m, gene, cells, bins, totals) {
  mm <- tapply(as.numeric(m[gene, cells]), bins, sum, default = 0)
  MM <- tapply(totals, bins, sum, default = 0)
  data.frame(bin = levels(bins), rate = as.numeric(mm / pmax(MM, 1)) * 1e3,
    lo = stats::qgamma(0.025, pmax(as.numeric(mm), 1e-9)) / pmax(as.numeric(MM), 1) * 1e3,
    hi = stats::qgamma(0.975, as.numeric(mm) + 1) / pmax(as.numeric(MM), 1) * 1e3,
    stringsAsFactors = FALSE)
}

# ---- admixture pattern ---------------------------------------------------------

# Share of each type's molecules removed by a correction.
generative_plot_burden <- function(before, after, annotation) {
  types <- annotation[colnames(before)]
  removed <- tapply(Matrix::colSums(before) - Matrix::colSums(after), types, sum)
  total <- tapply(Matrix::colSums(before), types, sum)
  burden <- data.frame(type = names(removed), share = as.numeric(removed / total))
  ggplot(burden, aes(reorder(type, share), share)) +
    geom_col(fill = generative_palette[["contamination"]], alpha = 0.85) +
    scale_y_continuous(labels = scales::percent) +
    coord_flip() +
    labs(x = NULL, y = "share of the type's molecules removed as admixture") +
    theme_classic(base_size = 11)
}

# Mean composition of the target type's content by exposure bin (own
# expression, contamination, ambient, induced): docs/generative.md Figure 1b.
generative_plot_composition <- function(model, target, exposure, x_label) {
  comp_all <- model$composition(target = target)
  cont <- tapply(comp_all$contamination, comp_all$cell_id, sum)
  indv <- tapply(comp_all$induced, comp_all$cell_id, sum)
  amb <- tapply(comp_all$ambient, comp_all$cell_id, max)
  per_cell <- data.frame(cell_id = names(cont), contamination = as.numeric(cont),
    induced = as.numeric(indv), ambient = as.numeric(amb), stringsAsFactors = FALSE)
  per_cell$own <- pmax(1 - per_cell$contamination - per_cell$induced - per_cell$ambient, 0)
  per_cell$bin <- generative_bins(exposure[per_cell$cell_id])
  agg <- aggregate(cbind(own, contamination, ambient, induced) ~ bin, per_cell, mean)
  shares <- data.frame(bin = rep(agg$bin, 4),
    part = rep(c("own expression", "contamination", "ambient", "induced"), each = nrow(agg)),
    share = c(agg$own, agg$contamination, agg$ambient, agg$induced))
  shares$part <- factor(shares$part, levels = rev(names(generative_palette)))
  ggplot(shares, aes(bin, share, fill = part)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = generative_palette, name = NULL,
      breaks = names(generative_palette)) +
    labs(x = x_label, y = sprintf("mean share of %s content", generative_short_label(target))) +
    theme_classic(base_size = 11)
}

# Fitted contamination fraction of one pair drawn over the tissue.
generative_plot_composition_map <- function(model, source, target, cells_xy,
                                            cap = 0.6, title = NULL, legend = NULL) {
  comp <- model$composition(source, target)
  xy <- cells_xy[match(comp$cell_id, as.character(cells_xy$cell_id)), ]
  ggplot(cbind(comp, x = xy$x, y = xy$y), aes(x, y, color = pmin(contamination, cap))) +
    geom_point(size = 0.25) +
    scale_color_viridis_c(option = "inferno",
      name = legend %||% sprintf("%s share\nof %s\ncontent", generative_short_label(source),
        generative_short_label(target))) +
    coord_equal() +
    labs(x = NULL, y = NULL, title = title %||% sprintf("%s contamination of %s cells",
      generative_short_label(source, lower = FALSE), generative_short_label(target))) +
    theme_void(base_size = 11)
}

# Three target cells spanning the contamination range of one pair: one
# unexposed with the least contamination, one near a given fraction, one
# with the most. Returns a named character vector of cell ids.
generative_pick_example_cells <- function(model, source, target, before, exposure,
                                          min_molecules = 100, moderate = 0.25) {
  comp <- model$composition(source, target)
  tot <- Matrix::colSums(before[, comp$cell_id])
  elig <- comp$cell_id[tot[comp$cell_id] >= min_molecules]
  cc <- stats::setNames(comp$contamination, comp$cell_id)[elig]
  picks <- c(`unexposed` = names(which.min(cc + (exposure[elig] > 0))),
    `interface, moderate` = names(which.min(abs(cc - moderate))),
    `interface, heavy` = names(which.max(cc)))
  stopifnot(!anyDuplicated(picks))
  picks
}

# Stacked composition of individual target cells (all sources, ambient, own).
generative_plot_cell_composition <- function(model, target, picks) {
  comp_all <- model$composition(target = target)
  anat <- comp_all[comp_all$cell_id %in% picks, ]
  anat$which <- names(picks)[match(anat$cell_id, picks)]
  amb <- anat[!duplicated(anat$cell_id), c("which", "ambient")]
  stack <- rbind(data.frame(which = anat$which, part = anat$source, share = anat$contamination),
    data.frame(which = amb$which, part = "ambient", share = amb$ambient))
  own <- 1 - tapply(stack$share, stack$which, sum)
  stack <- rbind(stack, data.frame(which = names(own), part = "own expression",
    share = as.numeric(own)))
  sources <- rev(sort(unique(anat$source)))
  stack$part <- factor(stack$part, levels = c("own expression", "ambient", sources))
  stack$which <- factor(stack$which, levels = names(picks))
  ggplot(stack, aes(which, share, fill = part)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c("own expression" = "#a8c6df", ambient = "#8a8a8a",
      stats::setNames(grDevices::hcl.colors(length(sources), "Reds 3"), sources))) +
    labs(x = NULL, y = "share of the cell's molecules", fill = NULL) +
    theme_classic(base_size = 11)
}

# A target cell in the package's standard molecule-level rendering, with
# molecules colored by gene origin: the target's own marker genes, the
# source's marker panel, and the rest.
generative_plot_gene_origin_cell <- function(nmf_fit, cell_id, target, own_genes,
                                             source_markers, classes, title) {
  prep <- celladmix_prepare_cell_example(nmf_fit,
    data.frame(target_cell = cell_id, target_cell_type = target,
      top_admix_factor = NA, stringsAsFactors = FALSE),
    markers = source_markers)
  g <- prep$molecules$gene
  prep$molecules$role <- factor(ifelse(g %in% source_markers, classes[[2]],
    ifelse(g %in% own_genes, classes[[1]], classes[[3]])), levels = classes)
  prep$role_levels <- classes
  celladmix_plot_cell_example(prep, title = title)
}

# Fitted contamination fraction against the exposure-derived prior.
generative_plot_dose_vs_fitted <- function(model, source, target, exposure, legend) {
  comp <- model$composition(source, target)
  ggplot(comp, aes(dose, contamination, color = factor(pmin(exposure[cell_id], 3)))) +
    geom_point(size = 0.4, alpha = 0.5) +
    geom_abline(color = "grey40", linewidth = 0.4) +
    scale_color_viridis_d(name = legend, end = 0.9) +
    labs(x = "dose: prior contamination fraction from exposure",
      y = "fitted contamination fraction") +
    theme_classic(base_size = 11)
}

# ---- induced expression ---------------------------------------------------------

# Neighborhood-associated expression changes on the corrected counts: every
# gene of the target type fitted with a Poisson regression on all sources'
# exposures jointly (quasi errors), with the two guards of
# docs/generative.md Figure 4 - the fraction of the gene's content removed
# as transfer, and its share of each source's interface-local transfer
# profile. Returns one row per gene and source.
generative_induced_changes <- function(before, after, model, cells_xy, annotation,
                                       target, neighbor_k = 15L, exposure_cap = 4,
                                       min_counts = 100) {
  pairs <- model$pairs[model$pairs$target == target, ]
  cells <- model$composition(pairs$source[[1]], target)$cell_id
  E <- sapply(pairs$source, function(src)
    pmin(generative_exposure(cells_xy, annotation, src, neighbor_k)[cells], exposure_cap))
  removed_frac <- 1 - Matrix::rowSums(after[, cells]) / pmax(Matrix::rowSums(before[, cells]), 1)
  K <- after[, cells]
  t_kept <- Matrix::colSums(K)
  ok <- t_kept > 10
  K <- K[, ok]; E <- E[ok, , drop = FALSE]; t_kept <- t_kept[ok]
  X <- cbind(1, E)
  fit_gene <- function(y) {
    beta <- c(log(max(sum(y), 0.5) / sum(t_kept)), numeric(ncol(E)))
    for (it in 1:25) {
      mu <- t_kept * exp(pmin(pmax(X %*% beta, -30), 5))
      H <- crossprod(X, X * as.numeric(mu))
      step <- tryCatch(solve(H + diag(1e-8, ncol(X)), crossprod(X, y - mu)),
        error = function(e) NULL)
      if (is.null(step)) return(NULL)
      beta <- beta + step
      if (max(abs(step)) < 1e-6) break
    }
    mu <- t_kept * exp(pmin(pmax(X %*% beta, -30), 5))
    disp <- max(sum((y - mu)^2 / pmax(mu, 1e-9)) / max(length(y) - ncol(X), 1), 1)
    se <- sqrt(pmax(diag(solve(crossprod(X, X * as.numeric(mu)) + diag(1e-8, ncol(X)))),
      1e-300) * disp)
    cbind(beta = beta[-1], z = beta[-1] / se[-1])
  }
  owner <- generative_gene_owner(before, annotation)$owner
  psi <- lapply(pairs$source, function(src) {
    p <- generative_interface_profile(before, cells_xy, annotation, src, target)
    p[names(owner)[owner == target]] <- 0
    p
  })
  rows <- lapply(rownames(K), function(g) {
    y <- as.numeric(K[g, ])
    if (sum(y) < min_counts) return(NULL)
    f <- fit_gene(y)
    if (is.null(f)) return(NULL)
    data.frame(gene = g, source = pairs$source,
      log_fold_per_neighbor = f[, "beta"], z = f[, "z"],
      removed_frac = removed_frac[[g]],
      psi_share_e4 = vapply(psi, function(p) p[[g]] * 1e4, 0),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$clean <- out$removed_frac < 0.10 & out$psi_share_e4 < 20
  rownames(out) <- NULL
  out
}

# Volcano of the changes for one source: docs/generative.md Figure 4.
generative_plot_volcano <- function(changes, source, z_cap = 60, label_n = 7,
                                    source_label = NULL) {
  d <- changes[changes$source == source & is.finite(changes$z), ]
  src <- source_label %||% generative_short_label(source)
  d$lfc <- d$log_fold_per_neighbor / log(2)
  d$zc <- pmin(abs(d$z), z_cap)
  d$tier <- ifelse(!d$clean, "transfer-entangled",
    ifelse(d$z > 0, paste("increased near", src), paste("decreased near", src)))
  lab <- d[d$clean, ]
  lab <- lab[order(-lab$zc), ][seq_len(min(label_n, sum(d$clean))), ]
  ggplot(d, aes(lfc, zc)) +
    geom_point(data = d[!d$clean, ], shape = 1, size = 1.3, color = "#b6b6b6") +
    geom_point(data = d[d$clean, ], aes(color = tier), size = 1.5) +
    geom_hline(yintercept = 6, color = "grey50", linetype = "13") +
    ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3) +
    scale_color_manual(values = stats::setNames(c("#1e8449", "#7d3c98"),
      c(paste("increased near", src), paste("decreased near", src))), name = NULL) +
    labs(x = sprintf("fold change per %s neighbor (log2)", src),
      y = sprintf("significance |z| (capped at %d)", z_cap)) +
    theme_classic(base_size = 11)
}

# The proportionality screen for one pair: measured exposure-linked excess
# of each source-owned gene against the excess expected from transfer, with
# the significance bound and the flagged genes drawn at both expectations
# (open ring: average source profile; filled: bordering cells' profile):
# docs/generative.md Figure 2a.
generative_plot_transfer_screen <- function(before, model, cells_xy, annotation,
                                            source, target, exposure, z = 8,
                                            profile_cv = 0.15, label_excess = 5e4) {
  own <- generative_gene_owner(before, annotation)
  cells <- model$composition(source, target)$cell_id
  exp_cells <- cells[exposure[cells] > 0]
  un_cells <- cells[exposure[cells] == 0]
  t_exp <- sum(Matrix::colSums(before[, exp_cells]))
  t_un <- sum(Matrix::colSums(before[, un_cells]))
  gset <- names(own$owner)[own$owner == source]
  m_exp <- Matrix::rowSums(before[gset, exp_cells])
  m_un <- Matrix::rowSums(before[gset, un_cells])
  excess <- m_exp - m_un / t_un * t_exp
  v <- m_exp + (t_exp / t_un)^2 * m_un + 1
  expected_of <- function(prof) {
    p <- prof
    p[names(p) %in% names(own$owner)[own$owner == target]] <- 0
    p <- p[gset]
    slope <- max(sum(excess * p / v) / sum(p^2 / v), 0)
    slope * p
  }
  exp_int <- expected_of(generative_interface_profile(before, cells_xy, annotation,
    source, target))
  exp_glob <- expected_of(own$profiles[, source] / sum(own$profiles[, source]))
  scr <- data.frame(gene = gset, excess = excess, expected = exp_int,
    expected_global = exp_glob, stringsAsFactors = FALSE)
  induced <- model$induced
  scr$flagged <- scr$gene %in% induced$gene[induced$source == source & induced$target == target]
  scr <- scr[scr$excess > 10 & scr$expected > 1, ]
  vfit <- stats::lm(ly ~ lx, data.frame(lx = log(scr$expected), ly = log(v[scr$gene])))
  xg <- exp(seq(log(min(scr$expected)), log(max(scr$expected) * 1.5), length.out = 200))
  bound <- data.frame(expected = xg,
    excess = xg + z * sqrt(exp(predict(vfit, data.frame(lx = log(xg)))) + (profile_cv * xg)^2))
  flag <- scr[scr$flagged, ]
  src <- generative_short_label(source); tgt <- generative_short_label(target)
  ggplot(scr, aes(expected, excess)) +
    geom_abline(color = "grey40", linewidth = 0.4) +
    geom_line(data = bound, linetype = "22", color = "grey50", linewidth = 0.4) +
    geom_point(data = scr[!scr$flagged, ], size = 1.3, color = "#9bb5c9") +
    geom_segment(data = flag, aes(x = expected_global, xend = expected, y = excess,
      yend = excess), color = "#1e8449", linewidth = 0.5) +
    geom_point(data = flag, aes(x = expected_global), size = 3.6, shape = 1, stroke = 1.1,
      color = "#c0392b") +
    geom_point(data = flag, size = 2, color = "#1e8449") +
    ggrepel::geom_text_repel(data = scr[scr$flagged | scr$excess > label_excess, ],
      aes(label = gene), size = 3) +
    scale_x_log10() + scale_y_log10() +
    labs(x = sprintf("excess expected in %s cells from %s transfer (molecules)", tgt, src),
      y = sprintf("measured excess in %s cells (molecules)", tgt)) +
    theme_classic(base_size = 11)
}

# Per-bin rates of chosen genes in the target cells, observed and retained,
# with exact intervals, the removed band, and an optional third series (for
# example the correction that also removes the induced share):
# docs/generative.md Figure 1c/d.
generative_plot_gene_contrast <- function(before, after, cells, exposure, genes,
                                          x_label, y_label, after_alt = NULL,
                                          alt_label = "induced share also removed") {
  tot <- Matrix::colSums(before[, cells])
  bins <- generative_bins(exposure[cells])
  series <- list(list(m = before, lab = "observed"), list(m = after, lab = "retained"))
  if (!is.null(after_alt)) series[[3]] <- list(m = after_alt, lab = alt_label)
  df <- do.call(rbind, lapply(genes, function(g) do.call(rbind, lapply(series, function(s)
    cbind(gene = g, counts = s$lab, generative_bin_rates(s$m, g, cells, bins, tot))))))
  levels_c <- vapply(series, `[[`, "", "lab")
  df$counts <- factor(df$counts, levels = levels_c)
  df$gene <- factor(df$gene, levels = genes)
  band <- merge(df[df$counts == "retained", c("gene", "bin", "rate")],
    df[df$counts == "observed", c("gene", "bin", "rate")],
    by = c("gene", "bin"), suffixes = c("_ret", "_obs"))
  base0 <- df[df$counts == "observed" & df$bin == "0", c("gene", "rate")]
  cols <- stats::setNames(c("#c0392b", "#2980b9", "#7d3c98")[seq_along(levels_c)], levels_c)
  shapes <- stats::setNames(c(15, 16, 16)[seq_along(levels_c)], levels_c)
  lines <- stats::setNames(c("solid", "solid", "22")[seq_along(levels_c)], levels_c)
  ggplot(df, aes(bin, rate, group = counts)) +
    geom_ribbon(data = band, aes(bin, ymin = rate_ret, ymax = rate_obs, group = gene),
      fill = "#c0392b", alpha = 0.12, inherit.aes = FALSE) +
    geom_hline(data = base0, aes(yintercept = rate), color = "grey50", linetype = "13",
      linewidth = 0.5) +
    geom_line(aes(color = counts, linetype = counts)) +
    geom_errorbar(aes(ymin = lo, ymax = hi, color = counts), width = 0.1, linewidth = 0.4,
      show.legend = FALSE) +
    geom_point(aes(color = counts, shape = counts), size = 1.9) +
    facet_wrap(~gene, scales = "free_y") +
    scale_color_manual(values = cols, name = NULL) +
    scale_shape_manual(values = shapes, name = NULL) +
    scale_linetype_manual(values = lines, name = NULL) +
    labs(x = x_label, y = y_label) +
    theme_classic(base_size = 11)
}

# Exposure-linked excess of a gene in the target cells of a pair.
generative_gene_excess <- function(m, before, gene, source, target, cells_xy, annotation) {
  t_cells <- names(annotation)[!is.na(annotation) & annotation == target]
  t_cells <- intersect(t_cells, colnames(before))
  e <- generative_exposure(cells_xy, annotation, source)[t_cells]
  tt <- Matrix::colSums(before[, t_cells])
  sum(m[gene, t_cells[e > 0]]) -
    sum(m[gene, t_cells[e == 0]]) / max(sum(tt[e == 0]), 1) * sum(tt[e > 0])
}

# Observed and retained exposure-linked excess of the largest retained genes.
generative_plot_induced_retention <- function(before, after, model, cells_xy, annotation, n = 8) {
  top <- head(model$induced[order(-model$induced$excess), ], n)
  ret <- do.call(rbind, lapply(seq_len(nrow(top)), function(i) {
    r <- top[i, ]
    data.frame(gene = r$gene, counts = c("observed excess", "excess after correction"),
      excess = c(generative_gene_excess(before, before, r$gene, r$source, r$target, cells_xy, annotation),
        generative_gene_excess(after, before, r$gene, r$source, r$target, cells_xy, annotation)))
  }))
  ret$counts <- factor(ret$counts, levels = c("observed excess", "excess after correction"))
  ggplot(ret, aes(reorder(gene, -excess), excess / 1e3, fill = counts)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("observed excess" = "#c0392b",
      "excess after correction" = "#1e8449")) +
    labs(x = NULL, y = "exposure-linked excess (thousand molecules)", fill = NULL) +
    theme_classic(base_size = 11)
}

# The largest retained gene of several pairs, each in its target cells
# against its source's exposure, observed and corrected with intervals.
generative_plot_induced_pairs <- function(before, after, model, cells_xy, annotation,
                                          n = 4, exclude = character()) {
  induced <- model$induced
  top <- do.call(rbind, lapply(split(induced, paste(induced$source, induced$target)),
    function(d) d[which.max(d$excess), ]))
  top <- top[order(-top$excess), ]
  top <- top[!top$gene %in% exclude, ][seq_len(min(n, nrow(top))), ]
  df <- do.call(rbind, lapply(seq_len(nrow(top)), function(i) {
    r <- top[i, ]
    tc <- model$composition(r$source, r$target)$cell_id
    bins <- generative_bins(generative_exposure(cells_xy, annotation, r$source)[tc])
    tot <- Matrix::colSums(before[, tc])
    lab <- sprintf("%s\n(%s → %s)", r$gene, generative_short_label(r$source, lower = FALSE), generative_short_label(r$target, lower = FALSE))
    rbind(cbind(panel = lab, counts = "observed", generative_bin_rates(before, r$gene, tc, bins, tot)),
      cbind(panel = lab, counts = "corrected", generative_bin_rates(after, r$gene, tc, bins, tot)))
  }))
  ggplot(df, aes(bin, rate, color = counts, group = counts)) +
    geom_line() + geom_point(size = 1.4) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.4) +
    facet_wrap(~panel, scales = "free_y", nrow = 2) +
    scale_color_manual(values = c(observed = "#c0392b", corrected = "#2980b9")) +
    labs(x = "source-type cells among 15 nearest neighbors",
      y = "rate in target cells (per 1,000 molecules)", color = NULL) +
    theme_classic(base_size = 10)
}

# Per-cell induced activity of a pair by exposure bin.
generative_plot_induced_activity <- function(model, source, target, exposure, x_label) {
  comp <- model$composition(source, target)
  ggplot(comp, aes(factor(pmin(exposure[cell_id], 3)), induced_activity)) +
    geom_boxplot(outlier.size = 0.3, fill = "#a9dfbf") +
    labs(x = x_label,
      y = "fitted induced activity of the cell\n(multiplier of the pair's induced term)") +
    theme_classic(base_size = 11)
}

# One gene's corrected rate drawn over the cells of a pair's target type.
generative_plot_gene_map <- function(after, model, source, target, cells_xy, gene,
                                     cap = 250, title = NULL) {
  comp <- model$composition(source, target)
  cells <- comp$cell_id
  rate <- after[gene, cells] / pmax(Matrix::colSums(after[, cells]), 1) * 1e3
  xy <- cells_xy[match(cells, as.character(cells_xy$cell_id)), ]
  ggplot(data.frame(x = xy$x, y = xy$y, rate = pmin(rate, cap)), aes(x, y, color = rate)) +
    geom_point(size = 0.25) +
    scale_color_viridis_c(name = sprintf("corrected %s\n(per 1,000\nmolecules)", gene)) +
    coord_equal() +
    labs(x = NULL, y = NULL, title = title %||%
      sprintf("The retained %s program in %s cells", gene, generative_short_label(target))) +
    theme_void(base_size = 11)
}

# Side-by-side cell-state UMAPs (ds$cell_state_umap(), correction$cell_state_umap()).
generative_plot_state_umaps <- function(original_state, corrected_state, annotation,
                                        titles = c("Original counts", "Generative correction")) {
  type_levels <- sort(unique(annotation))
  palette <- stats::setNames(grDevices::hcl.colors(length(type_levels), "Dark 3"), type_levels)
  panel <- function(df, title) {
    ggplot(df, aes(umap_1, umap_2, color = factor(cell_type, levels = type_levels))) +
      geom_point(size = 0.4, alpha = 0.7) +
      scale_color_manual(values = palette, drop = FALSE, name = NULL) +
      labs(title = title, x = "UMAP 1", y = "UMAP 2") +
      theme_classic(base_size = 10) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  }
  p1 <- panel(original_state, titles[[1]])
  p2 <- panel(corrected_state, titles[[2]])
  cowplot::plot_grid(
    cowplot::plot_grid(p1 + theme(legend.position = "none"),
      p2 + theme(legend.position = "none"), ncol = 2, align = "hv"),
    cowplot::get_legend(p1), ncol = 2, rel_widths = c(1, 0.22))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
