# Composite corrector: exposure-regression removal with induced-gene
# retention.
#
# The exposure-regression corrector v3 (08_phase_a.R) removes the measured
# exposure-linked excess but, by construction, removes induced expression
# along with transferred material. The generative-model evaluation
# (analysis/generative_model/FINDINGS.md) showed that the induced share is
# identifiable exactly where a gene's excess is disproportionate to the
# interface-local source profile, and that retaining it costs little
# removal power elsewhere. This script combines the two mechanisms:
#
#  1. the interface-local overdispersed screen (now the package's audit
#     screen) flags (pair, gene) combinations with disproportionate excess
#     and yields each one's induced share f = (excess - expected) / excess;
#  2. the regression corrector's dose is measured on unflagged guide genes
#     only, and its delivery caps each flagged gene's removable content at
#     the proportional (transferred) share 1 - f, waterfilling the rest to
#     the other source-owned genes; the ambient tier likewise retains the
#     f share of flagged strict genes.
#
# Arms: validation (A-guided, ambient off), production (pool-guided,
# ambient on), shuffled-exposure control. Scored with the gene-split
# harness plus an admixture-only yardstick (B halves excluding flagged
# genes), a per-gene retention ledger, and downstream kNN purity.
.libPaths(c("/tmp/celladmix_r_lib", .libPaths()))
suppressMessages(library(cellAdmixCore))
ag_dir <- "/home/pkharchenko/cellAdmix/cellAdmix-core/analysis/audit_guided"
source(file.path(ag_dir, "01_split_eval.R"))

dataset <- Sys.getenv("COMPOSITE_DATASET", "pancreas")
root <- "/home/pkharchenko/cellAdmix/cellAdmix-core"
if (dataset == "pancreas") {
  setwd(file.path(root, "examples/xenium_pancreas_membrane_377_full"))
  annotation <- read.csv("annotations/annotation.csv.gz", stringsAsFactors = FALSE)
  cell_annotation <- setNames(annotation$merged_annotation, annotation$cell_id)
  ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)
  fit <- ds$fit(nmf_variant = "invsqrt_kl", verbose = FALSE)
} else if (dataset == "nsclc") {
  setwd(file.path(root, "examples/cosmx_nsclc_giotto"))
  meta <- read.csv("prepared/cell_metadata_all.csv.gz", stringsAsFactors = FALSE)
  cell_annotation <- setNames(meta$cell_type_coarse, meta$cell)
  ds <- cellAdmix("prepared/molecules_all.csv.gz", output_dir = "out",
    annotation = cell_annotation)
  fit <- ds$fit(verbose = FALSE)
} else stop("unknown dataset")
res_name <- function(x) file.path(ag_dir, "results",
  paste0(if (dataset == "pancreas") "" else paste0(dataset, "_"), x))

defs <- build_split_defs(fit, cell_annotation)
counts0 <- defs$counts_before
profiles <- defs$profiles
top_type <- colnames(profiles)[max.col(profiles, ties.method = "first")]
totals0 <- defs$totals_before

# The screen, applied to every source-owned gene of each pair (not just
# the audit's marker panel): the corrector can remove any source-owned
# gene, so retention must cover them all.
cells_df <- fit$cell_factors()
near_dist <- cellAdmixCore:::.celladmix_source_nearest_distance(cells_df,
  cell_annotation)
flags <- list()
for (nm in names(defs$pairs)) {
  d <- defs$pairs[[nm]]
  gset <- rownames(profiles)[top_type == d$S]
  S_cells <- names(defs$cell_types)[!is.na(defs$cell_types) &
    defs$cell_types == d$S]
  psi_int <- cellAdmixCore:::.celladmix_audit_interface_profile(counts0,
    S_cells, near_dist[, d$T_type])
  sc <- cellAdmixCore:::.celladmix_audit_screen_panel(counts0, gset,
    d$T_cells, d$e[d$T_cells], totals0, psi_int, n_pool = length(gset))
  st <- sc$induced_stats
  if (nrow(st)) {
    st$induced_share <- pmax(st$excess - st$expected, 0) / st$excess
    flags[[nm]] <- st
  }
}
message(sprintf("flagged %d pair-gene combinations across %d pairs",
  sum(vapply(flags, nrow, 0L)), length(flags)))

pava <- function(y, w) {
  n <- length(y)
  val <- y; wt <- w; idx <- as.list(seq_len(n))
  i <- 1
  while (i < length(val)) {
    if (val[i] > val[i + 1] + 1e-15) {
      merged_w <- wt[i] + wt[i + 1]
      merged_v <- (val[i] * wt[i] + val[i + 1] * wt[i + 1]) / merged_w
      val[i] <- merged_v; wt[i] <- merged_w; idx[[i]] <- c(idx[[i]], idx[[i + 1]])
      val <- val[-(i + 1)]; wt <- wt[-(i + 1)]; idx[[i + 1]] <- NULL
      i <- max(1, i - 1)
    } else i <- i + 1
  }
  out <- numeric(n)
  for (j in seq_along(val)) out[idx[[j]]] <- val[j]
  out
}

stratum_budgets <- function(cur, d, guide, e) {
  strata <- sort(unique(e))
  mk <- ag_mcount(cur, guide, d$T_cells)
  tt <- totals0[d$T_cells]
  m_s <- tapply(mk, e, sum)[as.character(strata)]
  t_s <- tapply(tt, e, sum)[as.character(strata)]
  rate <- pava(as.numeric(m_s) / pmax(as.numeric(t_s), 1), as.numeric(t_s))
  r0 <- rate[strata == 0]
  if (!length(r0)) return(NULL)
  excess <- pmax(rate - r0, 0)
  data.frame(stratum = strata, excess_rate = excess, totals = as.numeric(t_s))
}

composite_correct <- function(guide_half = c("A", "pool"), ambient = FALSE,
                              retention = TRUE, shuffle = FALSE,
                              passes = 2L, seed = 1L) {
  guide_half <- match.arg(guide_half)
  set.seed(seed)
  cur <- counts0
  removed_ambient <- 0
  flagged_of <- function(nm) if (retention && !is.null(flags[[nm]]))
    flags[[nm]] else data.frame(gene = character(0), induced_share = numeric(0))
  if (ambient) {
    for (nm in names(defs$pairs)) {
      d <- defs$pairs[[nm]]
      fl <- flagged_of(nm)
      gs <- intersect(d$strict, rownames(cur))
      if (!length(gs)) next
      # Flagged strict genes keep their induced share; the rest of the
      # strict content is contamination at any exposure and is removed.
      keep <- setNames(rep(0, length(gs)), gs)
      keep[intersect(gs, fl$gene)] <-
        fl$induced_share[match(intersect(gs, fl$gene), fl$gene)]
      removed_ambient <- removed_ambient +
        sum(cur[gs, d$T_cells] * (1 - keep[gs]))
      cur[gs, d$T_cells] <- cur[gs, d$T_cells] * keep[gs]
    }
  }
  for (pass in seq_len(passes)) {
    total_budget <- 0; total_delivered <- 0
    for (nm in names(defs$pairs)) {
      d <- defs$pairs[[nm]]
      fl <- flagged_of(nm)
      e <- d$e[d$T_cells]
      if (shuffle) e <- setNames(sample(as.numeric(e)), names(e))
      guide <- if (guide_half == "A") d$A else d$pool
      if (ambient) guide <- setdiff(guide, d$strict)
      # The dose is measured on unflagged genes only: an induced gradient
      # is not transfer and must not inflate the removal budget.
      guide <- setdiff(guide, fl$gene)
      if (length(guide) < 3) next
      sb <- stratum_budgets(cur, d, guide, e)
      if (is.null(sb) || sum(sb$excess_rate) <= 0) next
      gset <- rownames(profiles)[top_type == d$S]
      if (ambient) gset <- setdiff(gset, d$strict)
      gset <- intersect(gset, rownames(cur))
      if (!length(gset)) next
      share_guide <- sum(profiles[intersect(guide, rownames(profiles)), d$S]) /
        max(sum(profiles[, d$S]), 1)
      share_gset <- sum(profiles[gset, d$S]) / max(sum(profiles[, d$S]), 1)
      baseline_T0 <- bench_pseudobulk(counts0, d$T_cells[d$e[d$T_cells] == 0])
      psi <- profiles[gset, d$S]
      w <- psi * psi / (psi + pmax(baseline_T0[gset], 0))
      w <- w / max(sum(w), 1e-12)
      # Cap each flagged gene's removable content at its proportional
      # (transferred) share; the induced share stays with the cell.
      cap <- setNames(rep(1, length(gset)), gset)
      shared <- intersect(gset, fl$gene)
      cap[shared] <- 1 - fl$induced_share[match(shared, fl$gene)]
      pair_mat <- as.matrix(cur[gset, d$T_cells, drop = FALSE])
      for (si in which(sb$excess_rate > 0)) {
        s <- sb$stratum[si]
        in_s <- e[d$T_cells] == s
        budget <- sb$excess_rate[si] * sb$totals[si] * share_gset / max(share_guide, 1e-6)
        total_budget <- total_budget + budget
        sub <- pair_mat[, in_s, drop = FALSE]
        avail_true <- rowSums(sub)
        avail_g <- avail_true * cap[gset]
        rem_g <- numeric(length(gset)); resid <- budget
        for (round in 1:4) {
          if (resid <= 0.01 * budget) break
          hw <- w * (avail_g - rem_g > 0.5)
          if (sum(hw) <= 0) break
          step <- pmin(avail_g - rem_g, hw / sum(hw) * resid)
          rem_g <- rem_g + step
          resid <- resid - sum(step)
        }
        total_delivered <- total_delivered + sum(rem_g)
        frac <- ifelse(avail_true > 1e-12,
          pmin(1, rem_g / pmax(avail_true, 1e-12)), 0)
        pair_mat[, in_s] <- sub * (1 - frac)
      }
      cur[gset, d$T_cells] <- as(pair_mat, "CsparseMatrix")
    }
    message(sprintf("  pass %d: budget %.0f, delivered %.0f (%.0f%%)",
      pass, total_budget, total_delivered,
      100 * total_delivered / max(total_budget, 1)))
    if (total_delivered < 0.05 * sum(totals0) / 1e3) break
  }
  if (ambient) message(sprintf("  ambient tier removed: %.0f molecules",
    removed_ambient))
  cur
}

# Admixture-only yardstick: power on B halves excluding flagged genes.
power_B_admixture_only <- function(counts_after) {
  rows <- lapply(names(defs$pairs), function(nm) {
    d <- defs$pairs[[nm]]
    fl <- if (!is.null(flags[[nm]])) flags[[nm]]$gene else character(0)
    Bx <- setdiff(d$B, fl)
    if (length(Bx) < 2) return(NULL)
    rates_of <- function(m) bench_bin_rates(ag_mcount(m, Bx, d$T_cells),
      totals0[d$T_cells], d$bins)
    rb <- rates_of(counts0)
    ex <- bench_detect(rb)$excess_molecules
    data.frame(pair = nm, excess_Bx = ex,
      power_Bx = bench_power(rb, rates_of(counts_after)))
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

# Retention ledger: exposure-linked excess of each flagged gene before and
# after correction.
retention_ledger <- function(counts_after) {
  rows <- list()
  for (nm in names(flags)) {
    d <- defs$pairs[[nm]]
    exp_cells <- d$T_cells[d$e[d$T_cells] > 0]
    un_cells <- d$T_cells[d$e[d$T_cells] == 0]
    t_exp <- sum(totals0[exp_cells]); t_un <- sum(totals0[un_cells])
    for (i in seq_len(nrow(flags[[nm]]))) {
      g <- flags[[nm]]$gene[[i]]
      if (!g %in% rownames(counts0)) next
      ex_of <- function(m) sum(m[g, exp_cells]) -
        sum(m[g, un_cells]) / max(t_un, 1) * t_exp
      eb <- ex_of(counts0); ea <- ex_of(counts_after)
      rows[[length(rows) + 1]] <- data.frame(pair = nm, gene = g,
        excess_before = eb, excess_after = ea,
        retained = ifelse(eb > 0, ea / eb, NA),
        induced_share = flags[[nm]]$induced_share[[i]])
    }
  }
  do.call(rbind, rows)
}

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

message("== validation arm (A-guided, no ambient, retention on) ==")
cv <- composite_correct("A", ambient = FALSE)
rv <- eval_split(defs, cv, "composite_validation")
print(rv$summary, row.names = FALSE)
write.csv(rv$pairs, res_name("composite_validation_pairs.csv"),
  row.names = FALSE)

message("\n== production arm (pool-guided, ambient on, retention on) ==")
cp <- composite_correct("pool", ambient = TRUE)
rp <- eval_split(defs, cp, "composite_production")
print(rp$summary, row.names = FALSE)
write.csv(rp$pairs, res_name("composite_production_pairs.csv"),
  row.names = FALSE)
write.csv(rp$types, res_name("composite_production_types.csv"),
  row.names = FALSE)

adm <- power_B_admixture_only(cp)
wadm <- sum(adm$power_Bx * adm$excess_Bx, na.rm = TRUE) /
  sum(adm$excess_Bx[!is.na(adm$power_Bx)])
message(sprintf(
  "admixture-only yardstick (B minus flagged): weighted %.3f, min %.3f",
  wadm, min(adm$power_Bx, na.rm = TRUE)))
low <- adm[!is.na(adm$power_Bx) & adm$power_Bx < 0.8, ]
if (nrow(low)) print(low, row.names = FALSE) else
  message("  every pair at or above 0.8")
write.csv(adm, res_name("composite_admixture_only.csv"),
  row.names = FALSE)

led <- retention_ledger(cp)
wret <- sum(led$excess_after) / max(sum(led$excess_before), 1)
message(sprintf(
  "retention: flagged excess kept %.1f%% in aggregate (%d pair-gene rows)",
  100 * wret, nrow(led)))
write.csv(led, res_name("composite_retention.csv"),
  row.names = FALSE)

message("\n== retention off (ablation): production machinery, no filter ==")
cn <- composite_correct("pool", ambient = TRUE, retention = FALSE)
rn <- eval_split(defs, cn, "composite_no_retention")
print(rn$summary, row.names = FALSE)
led_n <- retention_ledger(cn)
message(sprintf("  flagged excess kept without the filter: %.1f%%",
  100 * sum(led_n$excess_after) / max(sum(led_n$excess_before), 1)))

message("\n== shuffled-exposure control (retention on, ambient off) ==")
cs <- composite_correct("pool", ambient = FALSE, shuffle = TRUE, passes = 1L)
rs <- eval_split(defs, cs, "composite_shuffled")
print(rs$summary, row.names = FALSE)

message("\n== downstream kNN purity ==")
p0 <- knn_purity(counts0)
pp <- knn_purity(cp)
message(sprintf("original %.3f | composite %.3f", p0, pp))
message("COMPOSITE DONE")
