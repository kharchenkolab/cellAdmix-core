# How the shared profile of a cluster should be formed. For each cut of the
# pooled clustering, builds the shared profiles four ways (mean of unit-sum
# components, per-gene median, the medoid component, agreement-weighted mean)
# and reports each sample's training loss under them relative to its own fit,
# next to its restart spread. Uses cached per-sample fits only.
#
# Usage: Rscript consensus_variants.R <nmf_variant> <threads>
#          <sample>=<bundle>=<annotation>=<separate_dir> ...
.libPaths(c(Sys.getenv("CELLADMIX_R_LIB", "/home/pkharchenko/cellAdmix/wt-dev_multisample/.Rlib_dev"), .libPaths()))
suppressMessages(library(cellAdmixCore))
args <- commandArgs(trailingOnly = TRUE)
variant <- args[[1]]
threads <- as.integer(args[[2]])
rest <- args[-(1:2)]
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
own_profiles <- lapply(fits, cellAdmixCore:::.celladmix_fit_profiles_raw)
genes <- cellAdmixCore:::.celladmix_common_genes(own_profiles)
pools <- lapply(samples, function(s) lapply(cellAdmixCore:::.celladmix_fit_member_profiles(fits[[s]]), function(h) h[, genes, drop = FALSE]))
names(pools) <- samples
selected <- vapply(fits, function(fit) as.integer(fit$run$nmf_diagnostics$selected_run), integer(1))
own_loss <- vapply(samples, function(s) cellAdmixCore:::.celladmix_training_loss(fits[[s]]$run$path, own_profiles[[s]][, genes, drop = FALSE], threads), numeric(1))
spread <- vapply(samples, function(s) { o <- fits[[s]]$run$nmf_candidate_final_objectives; (max(o) - min(o)) / min(o) }, numeric(1))
cat(sprintf("%s: restart spread %s\n", variant, paste(sprintf("%s %.3f", samples, spread), collapse = ", ")))
stacked <- do.call(rbind, unlist(pools, recursive = FALSE))
normalized <- cellAdmixCore:::.celladmix_normalize_profiles(stacked)
row_sums <- rowSums(stacked)
gaps <- function(profiles) vapply(samples, function(s) cellAdmixCore:::.celladmix_training_loss(fits[[s]]$run$path, profiles, threads) / own_loss[[s]] - 1, numeric(1))
rows <- list()
for (cut in c(0.3, 0.5, 0.7)) {
  cons <- cellAdmixCore:::.celladmix_consensus_profiles(pools, weights, selected, match_min = cut, min_restarts = 2L)
  comp <- cons$components
  k <- nrow(cons$profiles)
  build <- function(how) {
    out <- cons$profiles
    for (i in seq_len(k)) {
      idx <- which(comp$consensus_factor == i)
      scale <- mean(row_sums[idx])
      out[i, ] <- switch(how,
        mean = cons$profiles[i, ],
        median = apply(normalized[idx, , drop = FALSE], 2, stats::median) * scale,
        medoid = stacked[idx[which.max(replace(comp$agreement[idx], is.na(comp$agreement[idx]), 1))], ],
        weighted = {
          w <- pmax(replace(comp$agreement[idx], is.na(comp$agreement[idx]), 1) - cut, 0) + 1e-6
          colSums(normalized[idx, , drop = FALSE] * w) / sum(w) * scale
        })
    }
    out
  }
  for (how in c("mean", "median", "medoid", "weighted")) {
    g <- gaps(build(how))
    rows[[length(rows) + 1L]] <- data.frame(cut = cut, construction = how, factors = k,
      unreproduced = nrow(cons$unreproduced), t(round(g, 4)), check.names = FALSE)
    cat(sprintf("cut %.1f %-9s factors %2d unreproduced %3d  gap %s\n", cut, how, k, nrow(cons$unreproduced),
      paste(sprintf("%s %+.4f", samples, g), collapse = "  ")))
  }
}
# Reference: each sample's loss under the other sample's primary profiles.
for (s in samples) for (o in setdiff(samples, s)) {
  g <- cellAdmixCore:::.celladmix_training_loss(fits[[s]]$run$path, own_profiles[[o]][, genes, drop = FALSE], threads) / own_loss[[s]] - 1
  cat(sprintf("%s under %s's own profiles: gap %+.4f\n", s, o, g))
}
