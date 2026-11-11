# Factor-space comparison on a cohort of samples: the pooled-clustering
# consensus of the package against the earlier star construction (every
# sample's primary factors matched one to one to one reference sample's and
# averaged, unmatched factors appended). Reports, per construction, the
# number of shared factors, duplicate factors (pairs of shared factors whose
# gene-ownership correlation reaches the matching threshold), sample
# coverage, and every sample's loss gap against its restart spread. Uses the
# cached per-sample fits; no scoring.
#
# Usage: Rscript cohort_factor_space.R <out_dir> <threads> <nmf_variant>
#          <sample>=<bundle>=<annotation>=<separate_dir> ...
.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/home/pkharchenko/cellAdmix/wt-dev_multisample/.Rlib_dev"), .libPaths()))
suppressMessages(library(cellAdmixCore))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[[1]]; threads <- as.integer(args[[2]]); variant <- args[[3]]; rest <- args[-(1:3)]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
read_annotation <- function(path) {
  ann <- read.csv(path, stringsAsFactors = FALSE)
  ann <- ann[!is.na(ann$cell_type) & ann$cell_type != "Unlabeled", ]
  stats::setNames(as.character(ann$cell_type), as.character(ann$cell_id))
}
specs <- lapply(strsplit(rest, "=", fixed = TRUE), function(p) list(sample = p[[1]], bundle = p[[2]], annotation = p[[3]], dir = p[[4]]))
datasets <- lapply(specs, function(spec) cellAdmix(spec$bundle, output_dir = spec$dir,
  annotation = read_annotation(spec$annotation), num_threads = threads))
names(datasets) <- vapply(specs, `[[`, character(1), "sample")
collection <- cellAdmixCollection(datasets)
collection$fit(nmf_variant = variant, nmf_n_runs = 10L)
fits <- collection$fits
samples <- names(fits)
molecules <- vapply(fits, function(fit) as.numeric(fit$run$n_transcripts), numeric(1))
weights <- molecules / sum(molecules)
own <- lapply(fits, cellAdmixCore:::.celladmix_fit_profiles_raw)
genes <- cellAdmixCore:::.celladmix_common_genes(own)
own <- lapply(own, function(h) h[, genes, drop = FALSE])
pools <- lapply(samples, function(s) lapply(cellAdmixCore:::.celladmix_fit_member_profiles(fits[[s]]), function(h) h[, genes, drop = FALSE]))
names(pools) <- samples
selected <- vapply(fits, function(fit) as.integer(fit$run$nmf_diagnostics$selected_run), integer(1))
own_loss <- vapply(samples, function(s) cellAdmixCore:::.celladmix_training_loss(fits[[s]]$run$path, own[[s]], threads), numeric(1))
spread <- vapply(samples, function(s) { o <- fits[[s]]$run$nmf_candidate_final_objectives; (max(o) - min(o)) / min(o) }, numeric(1))
normalize <- cellAdmixCore:::.celladmix_normalize_profiles

# Star construction: one-to-one matching of every other sample's primary
# factors to the reference sample's (the sample with most molecules), matched
# factors averaged, unmatched ones appended.
star <- function(match_min = 0.3) {
  reference <- samples[[which.max(molecules)]]
  ref_norm <- normalize(own[[reference]]); k <- nrow(ref_norm)
  sums <- ref_norm * weights[[reference]]; wt <- rep(weights[[reference]], k)
  scale <- rowSums(own[[reference]]) * weights[[reference]]
  appended <- list(); covered <- rep(1L, k)
  for (s in setdiff(samples, reference)) {
    h <- normalize(own[[s]]); m <- cellAdmixCore:::.celladmix_match_factors(h, ref_norm)
    for (i in seq_len(nrow(h))) {
      j <- m$match[[i]]; r <- m$correlation[[i]]
      if (!is.na(j) && is.finite(r) && r >= match_min) {
        sums[j, ] <- sums[j, ] + weights[[s]] * h[i, ]; wt[j] <- wt[j] + weights[[s]]
        scale[j] <- scale[j] + weights[[s]] * sum(own[[s]][i, ]); covered[j] <- covered[j] + 1L
      } else appended[[length(appended) + 1L]] <- own[[s]][i, ]
    }
  }
  profiles <- (sums / wt) * (scale / wt)
  for (p in appended) { profiles <- rbind(profiles, p); covered <- c(covered, 1L) }
  dimnames(profiles) <- list(paste0("F", seq_len(nrow(profiles))), genes)
  list(profiles = profiles, coverage = covered / length(samples))
}

describe <- function(label, profiles, coverage) {
  corr <- cellAdmixCore:::.celladmix_profile_correlation(profiles, profiles)
  diag(corr) <- NA
  duplicates <- sum(corr >= 0.3, na.rm = TRUE) / 2
  gaps <- vapply(samples, function(s) cellAdmixCore:::.celladmix_training_loss(fits[[s]]$run$path, profiles, threads) / own_loss[[s]] - 1, numeric(1))
  cat(sprintf("%-18s factors %3d  duplicate pairs (ownership r >= 0.3) %3d  max pair r %.2f  coverage: full %d, single-sample %d  loss gap above spread: %d of %d samples (max gap %+.3f)\n",
    label, nrow(profiles), duplicates, max(corr, na.rm = TRUE), sum(coverage >= 1), sum(coverage <= 1 / length(samples) + 1e-9),
    sum(gaps > spread), length(samples), max(gaps)))
  data.frame(construction = label, factors = nrow(profiles), duplicate_pairs = duplicates, max_pair_correlation = max(corr, na.rm = TRUE),
    full_coverage = sum(coverage >= 1), single_sample = sum(coverage <= 1 / length(samples) + 1e-9),
    samples_flagged = sum(gaps > spread), max_gap = max(gaps), stringsAsFactors = FALSE)
}
cat(sprintf("%d samples; restart spreads %s\n", length(samples), paste(sprintf("%.3f", spread), collapse = " ")))
rows <- list()
st <- star(); rows[[1]] <- describe("star (pairwise)", st$profiles, st$coverage)
cons <- cellAdmixCore:::.celladmix_consensus_profiles(pools, weights, selected, match_min = 0.3, min_restarts = 2L)
rows[[2]] <- describe("pooled clustering", cons$profiles, cons$factors$recurrence)
write.csv(do.call(rbind, rows), file.path(out_dir, "factor_space.csv"), row.names = FALSE)
write.csv(cons$factors, file.path(out_dir, "pooled_factors.csv"), row.names = FALSE)
gaps <- data.frame(sample = samples, restart_spread = spread,
  star_gap = vapply(samples, function(s) cellAdmixCore:::.celladmix_training_loss(fits[[s]]$run$path, st$profiles, threads) / own_loss[[s]] - 1, numeric(1)),
  pooled_gap = vapply(samples, function(s) cellAdmixCore:::.celladmix_training_loss(fits[[s]]$run$path, cons$profiles, threads) / own_loss[[s]] - 1, numeric(1)))
write.csv(gaps, file.path(out_dir, "loss_gaps.csv"), row.names = FALSE)
print(gaps)
