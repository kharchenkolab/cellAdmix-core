#' Construct a Multi-Sample cellAdmix Collection
#'
#' Groups per-sample datasets that share one cell-type annotation scheme so
#' that factorization, scoring, and correction treat every sample by the same
#' procedure: one shared factor space (a consensus of the per-sample
#' factorizations), one set of removal rules derived from the evidence pooled
#' across samples, and sample-specific molecule labels and removal amounts.
#'
#' @param datasets Named list of [CellAdmixDataset] objects, one per sample,
#'   each with an active annotation. The names are the sample names.
#' @param sample_allocation How quantities with a per-dataset budget (profile
#'   averaging weights, scoring candidate caps) are split across samples:
#'   `"proportional"` to each sample's molecule count, or `"equal"`.
#'
#' @return A `CellAdmixCollection` object.
#' @export
cellAdmixCollection <- function(datasets, sample_allocation = c("proportional", "equal")) {
  CellAdmixCollection$new(datasets, sample_allocation = match.arg(sample_allocation))
}

#' Multi-Sample Collection of cellAdmix Datasets
#'
#' R6 object holding the per-sample datasets of a multi-sample analysis and
#' the fits made on them. Construct with [cellAdmixCollection()]. The
#' workflow is `fit()` (per-sample factorizations, cached like any fit),
#' `fit_consensus()` (match factors across samples, average them into
#' shared profiles, and label every sample against those profiles), then
#' scoring, rules, and correction on the returned [CellAdmixConsensusFit].
#' @export
CellAdmixCollection <- R6::R6Class(
  "CellAdmixCollection",
  public = list(
    datasets = NULL,
    sample_allocation = NULL,
    fits = NULL,

    initialize = function(datasets, sample_allocation = "proportional") {
      if (!is.list(datasets) || !length(datasets)) {
        stop("datasets must be a non-empty named list of CellAdmixDataset objects")
      }
      if (is.null(names(datasets)) || any(!nzchar(names(datasets))) ||
          anyDuplicated(names(datasets))) {
        stop("datasets must carry unique, non-empty sample names")
      }
      for (name in names(datasets)) {
        if (!inherits(datasets[[name]], "CellAdmixDataset")) {
          stop("datasets[['", name, "']] is not a CellAdmixDataset")
        }
        if (is.null(datasets[[name]]$active_annotation)) {
          stop("dataset '", name, "' has no active annotation")
        }
      }
      self$datasets <- datasets
      self$sample_allocation <- match.arg(sample_allocation, c("proportional", "equal"))
      self$fits <- list()
    },

    print = function(...) {
      cat("cellAdmix collection\n")
      cat("  samples:", paste(names(self$datasets), collapse = ", "), "\n")
      cat("  allocation:", self$sample_allocation, "\n")
      cat("  per-sample fits:", if (length(self$fits)) paste(names(self$fits), collapse = ", ") else "none", "\n")
      invisible(self)
    },

    samples = function() {
      out <- do.call(rbind, lapply(names(self$datasets), function(name) {
        ds <- self$datasets[[name]]
        ann <- ds$annotation(as_vector = TRUE)
        data.frame(sample = name, annotated_cells = length(ann),
          cell_types = length(unique(ann)), output_dir = ds$output_dir,
          stringsAsFactors = FALSE)
      }))
      rownames(out) <- NULL
      out
    },

    fit = function(...) {
      # One ordinary fit per sample; cached runs are reused like any fit.
      for (name in names(self$datasets)) {
        self$fits[[name]] <- self$datasets[[name]]$fit(...)
      }
      invisible(self)
    },

    fit_consensus = function(..., reference = NULL, match_min = 0.3,
                             unmatched = c("append", "drop"), name = "consensus",
                             num_threads = NULL) {
      unmatched <- match.arg(unmatched)
      state <- private$prepare_profiles(...)
      reference <- reference %||% state$samples[[which.max(state$molecules)]]
      if (!(reference %in% state$samples)) {
        stop("reference must name one of the samples")
      }
      consensus <- .celladmix_consensus_profiles(state$profiles, state$weights, reference,
        match_min = match_min, unmatched = unmatched)
      # Each factor's restart stability within its own sample, next to its
      # cross-sample match: an unmatched factor that is also unstable within
      # its sample is restart noise rather than sample-specific structure.
      consensus$matching$own_stability <- .celladmix_own_factor_stability(
        state$fits, consensus$matching$sample, consensus$matching$sample_factor)
      consensus$sample_specific$own_stability <- .celladmix_own_factor_stability(
        state$fits, consensus$sample_specific$sample, consensus$sample_specific$sample_factor)
      # Each sample's restart profiles, aligned to the consensus factor order,
      # become the member pool of its fixed-profile run.
      members <- lapply(state$samples, function(sample) {
        pool <- .celladmix_fit_member_profiles(state$fits[[sample]])
        pool <- lapply(pool, function(h) h[, state$genes, drop = FALSE])
        .celladmix_align_members(pool, consensus$profiles, match_min = match_min)
      })
      names(members) <- state$samples
      run_id <- paste0(.celladmix_clean_name(name, "run"), "_",
        .celladmix_matrix_hash(consensus$profiles))
      CellAdmixConsensusFit$new(
        collection = self,
        fits = private$project(consensus$profiles, members, run_id, num_threads),
        own_fits = state$fits,
        profiles = consensus$profiles,
        matching = consensus$matching,
        sample_specific = consensus$sample_specific,
        loss = private$loss_table(state, consensus$profiles),
        params = list(reference = reference, match_min = match_min,
          unmatched = unmatched, weights = state$weights, name = name, run_id = run_id,
          source = "consensus")
      )
    },

    fit_profiles = function(profiles, members = NULL, name = "profiles",
                            match_min = 0.3, num_threads = NULL, ...) {
      # Project every sample onto given profiles (factor x gene matrix with
      # gene column names): profiles of a joint factorization, or one sample's
      # own, so that alternatives to the consensus run through the same
      # scoring, rules, and diagnostics. `members`, when given, is a list of
      # matrices shaped like `profiles` that becomes every sample's member
      # pool after alignment to the profiles' factor order.
      state <- private$prepare_profiles(...)
      profiles <- as.matrix(profiles)
      if (is.null(colnames(profiles))) {
        stop("profiles must carry gene names as column names")
      }
      genes <- .celladmix_common_genes(c(list(profiles), unname(state$profiles)))
      profiles <- profiles[, genes, drop = FALSE]
      dimnames(profiles) <- list(paste0("F", seq_len(nrow(profiles))), genes)
      profiles_norm <- .celladmix_normalize_profiles(profiles)
      matching <- do.call(rbind, lapply(state$samples, function(sample) {
        own <- .celladmix_normalize_profiles(state$profiles[[sample]])
        match <- .celladmix_match_factors(own, profiles_norm)
        matched <- !is.na(match$match) & is.finite(match$correlation) &
          match$correlation >= match_min
        data.frame(sample = sample, sample_factor = seq_len(nrow(own)),
          consensus_factor = ifelse(matched, match$match, NA_integer_),
          correlation = match$correlation, matched = matched, stringsAsFactors = FALSE)
      }))
      rownames(matching) <- NULL
      matching$own_stability <- .celladmix_own_factor_stability(
        state$fits, matching$sample, matching$sample_factor)
      aligned <- if (length(members)) {
        .celladmix_align_members(lapply(members, function(h) {
          h <- as.matrix(h)
          if (is.null(colnames(h))) stop("members must carry gene names as column names")
          h
        }), profiles, match_min = match_min)
      } else {
        list()
      }
      member_pools <- stats::setNames(rep(list(aligned), length(state$samples)), state$samples)
      run_id <- paste0(.celladmix_clean_name(name, "run"), "_", .celladmix_matrix_hash(profiles))
      CellAdmixConsensusFit$new(
        collection = self,
        fits = private$project(profiles, member_pools, run_id, num_threads),
        own_fits = state$fits,
        profiles = profiles,
        matching = matching,
        sample_specific = data.frame(consensus_factor = integer(), sample = character(),
          sample_factor = integer(), stringsAsFactors = FALSE),
        loss = private$loss_table(state, profiles),
        params = list(reference = NA_character_, match_min = match_min, name = name,
          run_id = run_id, source = "profiles")
      )
    }
  ),
  private = list(
    # Per-sample fits with their raw-space profiles over the common gene panel
    # and the averaging weights.
    prepare_profiles = function(...) {
      if (!length(self$fits)) {
        self$fit(...)
      }
      fits <- self$fits
      samples <- names(fits)
      if (length(samples) < 2L) {
        stop("a collection fit needs at least two samples")
      }
      molecules <- vapply(fits, function(fit) as.numeric(fit$run$n_transcripts), numeric(1))
      weights <- if (identical(self$sample_allocation, "equal")) {
        stats::setNames(rep(1, length(samples)), samples)
      } else {
        molecules / sum(molecules)
      }
      profiles <- lapply(fits, .celladmix_fit_profiles_raw)
      genes <- .celladmix_common_genes(profiles)
      profiles <- lapply(profiles, function(h) h[, genes, drop = FALSE])
      list(fits = fits, samples = samples, molecules = molecules, weights = weights,
        profiles = profiles, genes = genes)
    },

    # How well shared profiles represent each sample's own training rows.
    loss_table = function(state, shared) {
      out <- do.call(rbind, lapply(state$samples, function(sample) {
        run <- state$fits[[sample]]$run
        own <- .celladmix_training_loss(run$path, state$profiles[[sample]])
        under_shared <- .celladmix_training_loss(run$path, shared)
        objectives <- run$nmf_candidate_final_objectives
        spread <- if (length(objectives) > 1L && min(objectives) > 0) {
          (max(objectives) - min(objectives)) / min(objectives)
        } else {
          NA_real_
        }
        data.frame(sample = sample, own_loss = own, consensus_loss = under_shared,
          loss_gap = under_shared / own - 1, restart_spread = spread,
          stringsAsFactors = FALSE)
      }))
      rownames(out) <- NULL
      out
    },

    # Fixed-profile projection of every sample under one run id.
    project = function(shared, member_pools, run_id, num_threads) {
      projected <- lapply(names(self$fits), function(sample) {
        own <- self$fits[[sample]]
        ds <- self$datasets[[sample]]
        variant <- own$run$pipeline_options$nmf_variant %||% "ls_nmf"
        pool <- member_pools[[sample]]
        ds$fit(
          nmf_fixed_h = shared,
          nmf_fixed_h_members = if (length(pool)) pool else NULL,
          run_id = run_id,
          nmf_variant = variant,
          nmf_n_runs = 1L,
          num_threads = num_threads %||% ds$num_threads
        )
      })
      names(projected) <- names(self$fits)
      projected
    }
  )
)

#' Consensus Factorization Across Samples
#'
#' Result of [CellAdmixCollection]`$fit_consensus()`: the shared gene
#' profiles, the factor matching that produced them, one fixed-profile fit
#' per sample (whose molecule labels and ensemble members share the consensus
#' factor index), and diagnostics of how well the consensus represents each
#' sample. Scoring methods return a [CellAdmixCollectionScore] whose rules
#' pool the evidence of all samples.
#' @export
CellAdmixConsensusFit <- R6::R6Class(
  "CellAdmixConsensusFit",
  public = list(
    collection = NULL,
    fits = NULL,
    own_fits = NULL,
    profiles = NULL,
    matching = NULL,
    sample_specific = NULL,
    loss = NULL,
    params = NULL,

    initialize = function(collection, fits, own_fits, profiles, matching,
                          sample_specific, loss, params = list()) {
      self$collection <- collection
      self$fits <- fits
      self$own_fits <- own_fits
      self$profiles <- profiles
      self$matching <- matching
      self$sample_specific <- sample_specific
      self$loss <- loss
      self$params <- params
    },

    print = function(...) {
      cat("cellAdmix consensus fit\n")
      cat("  samples:", paste(names(self$fits), collapse = ", "), "\n")
      cat("  factors:", nrow(self$profiles),
        if (nrow(self$sample_specific)) paste0("(", nrow(self$sample_specific), " sample-specific)") else "", "\n")
      cat("  reference sample:", self$params$reference, "\n")
      matched <- self$matching[self$matching$sample != self$params$reference, , drop = FALSE]
      if (nrow(matched)) {
        cat("  matched ownership correlation (non-reference factors): median",
          round(stats::median(matched$correlation, na.rm = TRUE), 3), "\n")
      }
      invisible(self)
    },

    samples = function() names(self$fits),
    loadings = function() self$profiles,

    stability = function() {
      # Within-sample restart stability of each own fit next to the agreement
      # of that sample's aligned restarts with the consensus profiles.
      out <- do.call(rbind, lapply(names(self$fits), function(sample) {
        own <- self$own_fits[[sample]]$run
        shared <- self$fits[[sample]]$run
        data.frame(
          sample = sample,
          within_sample_stability = own$nmf_candidate_best_match_correlation_mean %||% NA_real_,
          consensus_agreement = shared$nmf_candidate_best_match_correlation_mean %||% NA_real_,
          members = shared$fixed_profile_members %||% 0L,
          stringsAsFactors = FALSE)
      }))
      rownames(out) <- NULL
      out
    },

    sample_usage = function() {
      # Share of each cell type's molecules labeled with each factor, per
      # sample, from the fixed-profile fits' cell factor fractions.
      out <- do.call(rbind, lapply(names(self$fits), function(sample) {
        fit <- self$fits[[sample]]
        cells <- fit$cell_factors()
        annotation <- fit$dataset$annotation(fit$annotation_name, as_vector = TRUE)
        cells$cell_type <- as.character(annotation[as.character(cells$cell_id)])
        cells <- cells[!is.na(cells$cell_type), , drop = FALSE]
        fraction_cols <- grep("^factor_[0-9]+_fraction$", names(cells), value = TRUE)
        fraction_cols <- fraction_cols[order(as.integer(sub("^factor_([0-9]+)_fraction$", "\\1", fraction_cols)))]
        molecules <- as.matrix(cells[, fraction_cols, drop = FALSE]) * cells$transcript_count
        by_type <- rowsum(molecules, cells$cell_type)
        share <- by_type / pmax(rowSums(by_type), 1)
        data.frame(
          sample = sample,
          cell_type = rep(rownames(by_type), ncol(by_type)),
          factor = rep(seq_len(ncol(by_type)), each = nrow(by_type)),
          molecules = as.vector(by_type),
          share = as.vector(share),
          stringsAsFactors = FALSE)
      }))
      rownames(out) <- NULL
      out
    },

    usage_divergence = function() {
      # Per cell type: half the L1 distance between a sample's factor usage
      # and the pooled usage, maximized over samples. Zero when every sample
      # labels the type's molecules alike.
      usage <- self$sample_usage()
      pooled <- stats::aggregate(molecules ~ cell_type + factor, usage, sum)
      totals <- stats::aggregate(molecules ~ cell_type, pooled, sum)
      pooled$share <- pooled$molecules / totals$molecules[match(pooled$cell_type, totals$cell_type)]
      key <- function(df) paste(df$cell_type, df$factor, sep = "\r")
      usage$pooled_share <- pooled$share[match(key(usage), key(pooled))]
      usage$deviation <- abs(usage$share - usage$pooled_share)
      per_sample <- stats::aggregate(deviation ~ cell_type + sample, usage, function(x) sum(x) / 2)
      out <- stats::aggregate(deviation ~ cell_type, per_sample, max)
      names(out)[names(out) == "deviation"] <- "divergence"
      out$molecules <- totals$molecules[match(out$cell_type, totals$cell_type)]
      out[order(-out$divergence), , drop = FALSE]
    },

    cells = function() {
      out <- do.call(rbind, lapply(names(self$fits), function(sample) {
        cells <- self$fits[[sample]]$cell_factors()
        cells$sample <- sample
        cells$orig_cell_id <- as.character(cells$cell_id)
        cells$cell_id <- paste0(sample, "_", cells$orig_cell_id)
        cells
      }))
      rownames(out) <- NULL
      out
    },

    counts = function(sample = NULL, ...) {
      if (!is.null(sample)) {
        return(self$fits[[sample]]$counts(...))
      }
      .celladmix_bind_sample_counts(lapply(self$fits, function(fit) fit$counts(...)))
    },

    audit_admixture = function(...) {
      audits <- lapply(self$fits, function(fit) fit$audit_admixture(...))
      CellAdmixCollectionAudit$new(self, audits)
    },

    score = function(method = c("bridge", "membrane"), name = NULL, ...) {
      method <- match.arg(method)
      private$score_samples(method, name = name %||% method, ...)
    },
    score_bridge = function(name = "bridge", ...) private$score_samples("bridge", name = name, ...),
    score_membrane = function(name = "membrane", ...) private$score_samples("membrane", name = name, ...),

    correct = function(score, rules = NULL, ...) {
      if (!inherits(score, "CellAdmixCollectionScore")) {
        stop("score must be a CellAdmixCollectionScore")
      }
      score$correct(rules = rules, ...)
    }
  ),
  private = list(
    score_samples = function(method, name, ...) {
      dots <- list(...)
      shares <- private$allocation_shares()
      defaults <- if (identical(method, "bridge")) {
        c(candidate_pairs_per_type_pair = 1000L, max_cells_per_type_pair = 400L)
      } else {
        c(candidate_pairs_per_type_pair = 400L, max_cells_per_type_pair = 400L)
      }
      scores <- lapply(names(self$fits), function(sample) {
        args <- dots
        # The candidate budget of k single-sample analyses is split across the
        # samples by the collection's allocation rule.
        for (cap in names(defaults)) {
          base <- dots[[cap]] %||% defaults[[cap]]
          args[[cap]] <- max(1L, as.integer(round(base * shares[[sample]])))
        }
        fn <- if (identical(method, "membrane")) self$fits[[sample]]$score_membrane else self$fits[[sample]]$score_bridge
        do.call(fn, c(list(name = name), args))
      })
      names(scores) <- names(self$fits)
      CellAdmixCollectionScore$new(self, method, name, scores, params = dots)
    },

    allocation_shares = function() {
      samples <- names(self$fits)
      if (identical(self$collection$sample_allocation, "equal")) {
        return(stats::setNames(rep(1, length(samples)), samples))
      }
      molecules <- vapply(self$fits, function(fit) as.numeric(fit$run$n_transcripts), numeric(1))
      stats::setNames(length(samples) * molecules / sum(molecules), samples)
    }
  )
)

#' Pooled Admixture Score Across Samples
#'
#' Per-sample scores of a consensus fit together with the pooled summary:
#' per target/source/factor, pair counts are summed, mean scores are
#' pair-weighted, and the one-sided p-values are combined by Stouffer's
#' method with weights proportional to the square root of each sample's pair
#' count. Rules derived from the pooled summary carry the fraction of samples
#' whose own evidence supports them (`sample_support`), as a report, never a
#' gate.
#' @export
CellAdmixCollectionScore <- R6::R6Class(
  "CellAdmixCollectionScore",
  public = list(
    consensus = NULL,
    method = NULL,
    name = NULL,
    sample_scores = NULL,
    params = NULL,

    initialize = function(consensus, method, name, sample_scores, params = list()) {
      self$consensus <- consensus
      self$method <- method
      self$name <- .celladmix_clean_name(name, "score")
      self$sample_scores <- sample_scores
      self$params <- params
    },

    print = function(...) {
      cat("cellAdmix collection score\n")
      cat("  method:", self$method, "\n")
      cat("  samples:", paste(names(self$sample_scores), collapse = ", "), "\n")
      cat("  pooled summary rows:", nrow(self$summary()), "\n")
      invisible(self)
    },

    samples = function() names(self$sample_scores),

    sample_summaries = function() {
      .celladmix_stack_samples(
        lapply(self$sample_scores, function(score) score$summary()),
        empty_columns = c("target_cell_type", "source_cell_type", "factor", "n_pairs",
          "mean_score", "p_value", "neg_log10_p"))
    },

    pairs = function() {
      .celladmix_stack_samples(lapply(self$sample_scores, function(score) score$pairs()))
    },

    summary = function() {
      if (is.null(private$.pooled)) {
        private$.pooled <- .celladmix_pool_score_summaries(self$sample_summaries())
      }
      private$.pooled
    },

    annotation = function(p_thresh = 0.1, adjust_p = FALSE) {
      celladmix_bridge_annotation(self$summary(), p_thresh = p_thresh, adjust_p = adjust_p)
    },

    rules = function(p_thresh = 0.1, adjust_p = FALSE, targets = NULL,
                     native_check = TRUE, native_median_thresh = 0.1,
                     native_expr_thresh = 0.05, native_outlier_min_frac = 0.1,
                     neighbor_k = 15L) {
      rules <- celladmix_bridge_rules(self$summary(), p_thresh = p_thresh,
        adjust_p = adjust_p, target_cell_types = targets)
      rules <- private$add_sample_support(rules, p_thresh = p_thresh)
      if (isTRUE(native_check)) {
        rules <- private$pooled_native_check(rules,
          median_thresh = native_median_thresh, expr_thresh = native_expr_thresh,
          outlier_min_frac = native_outlier_min_frac, neighbor_k = neighbor_k)
      }
      rules
    },

    correct = function(rules = NULL, name = NULL, p_thresh = 0.1, adjust_p = FALSE,
                       targets = NULL, ensemble = NULL, vote = 0.3, ...) {
      if (is.null(rules)) {
        rules <- self$rules(p_thresh = p_thresh, adjust_p = adjust_p, targets = targets)
      }
      kept <- rules
      if ("keep" %in% names(kept)) {
        kept <- kept[is.na(kept$keep) | kept$keep, , drop = FALSE]
      }
      name <- .celladmix_clean_name(name %||% paste0(self$name, "_clean"), "correction")
      corrections <- lapply(names(self$sample_scores), function(sample) {
        self$consensus$fits[[sample]]$correct(
          self$sample_scores[[sample]], rules = kept, name = name,
          ensemble = ensemble, vote = vote, ...)
      })
      names(corrections) <- names(self$sample_scores)
      CellAdmixCollectionCorrection$new(self, name, corrections, rules = rules)
    }
  ),
  private = list(
    .pooled = NULL,

    add_sample_support = function(rules, p_thresh) {
      # Fraction of samples whose own summary reaches p_thresh for the rule's
      # factor and target (any source), plus their names. Report only.
      summaries <- self$sample_summaries()
      n_samples <- length(self$sample_scores)
      support <- vapply(seq_len(nrow(rules)), function(i) {
        rows <- summaries[summaries$factor == rules$factor[[i]] &
          summaries$target_cell_type == rules$target_cell_type[[i]], , drop = FALSE]
        if (!nrow(rows)) return(0)
        hits <- stats::aggregate(p_value ~ sample, rows, min)
        sum(hits$p_value < p_thresh, na.rm = TRUE) / n_samples
      }, numeric(1))
      rules$sample_support <- support
      rules$samples_significant <- vapply(seq_len(nrow(rules)), function(i) {
        rows <- summaries[summaries$factor == rules$factor[[i]] &
          summaries$target_cell_type == rules$target_cell_type[[i]], , drop = FALSE]
        if (!nrow(rows)) return("")
        hits <- stats::aggregate(p_value ~ sample, rows, min)
        paste(sort(hits$sample[!is.na(hits$p_value) & hits$p_value < p_thresh]), collapse = ",")
      }, character(1))
      rules
    },

    pooled_native_check = function(rules, median_thresh, expr_thresh,
                                   outlier_min_frac, neighbor_k) {
      # Source exposure is computed within each sample (its own cell
      # neighborhoods) and the per-cell tables are stacked, so the check's
      # statistics pool the samples without neighborhoods crossing them.
      fits <- self$consensus$fits
      annotations <- lapply(fits, function(fit) fit$dataset$annotation(fit$annotation_name, as_vector = TRUE))
      types <- sort(unique(unlist(lapply(annotations, function(ann) as.character(ann[!is.na(ann)])), use.names = FALSE)))
      pieces <- lapply(names(fits), function(sample) {
        cells <- fits[[sample]]$cell_factors()
        annotation <- annotations[[sample]]
        exposure <- .celladmix_source_exposure_counts(cells, annotation, neighbor_k, types = types)
        cells$cell_id <- paste0(sample, "_", as.character(cells$cell_id))
        names(annotation) <- paste0(sample, "_", names(annotation))
        list(cells = cells, counts = exposure$counts, annotation = annotation)
      })
      cells <- do.call(rbind, lapply(pieces, `[[`, "cells"))
      counts <- do.call(rbind, lapply(pieces, `[[`, "counts"))
      annotation <- unlist(unname(lapply(pieces, `[[`, "annotation")))
      .celladmix_apply_native_check(rules, fit = NULL,
        median_thresh = median_thresh, expr_thresh = expr_thresh,
        outlier_min_frac = outlier_min_frac, neighbor_k = neighbor_k,
        cell_factors = cells, annotation = annotation,
        exposure = list(counts = counts, types = types))
    }
  )
)

#' Per-Sample Corrections Under Shared Rules
#'
#' The corrections of every sample of a consensus fit under one pooled rule
#' set, with per-sample and pooled removal summaries and combined or
#' per-sample corrected counts.
#' @export
CellAdmixCollectionCorrection <- R6::R6Class(
  "CellAdmixCollectionCorrection",
  public = list(
    score = NULL,
    name = NULL,
    corrections = NULL,
    rules = NULL,

    initialize = function(score, name, corrections, rules = NULL) {
      self$score <- score
      self$name <- name
      self$corrections <- corrections
      self$rules <- rules
    },

    print = function(...) {
      cat("cellAdmix collection correction\n")
      cat("  samples:", paste(names(self$corrections), collapse = ", "), "\n")
      removed <- vapply(self$corrections, function(x) as.numeric(x$run$n_removed %||% NA_real_), numeric(1))
      cat("  molecules removed:", paste(format(removed, big.mark = ","), collapse = ", "), "\n")
      invisible(self)
    },

    samples = function() names(self$corrections),

    summary = function(by_sample = TRUE) {
      per_sample <- .celladmix_stack_samples(
        lapply(self$corrections, function(correction) correction$summary()))
      if (isTRUE(by_sample)) {
        return(per_sample)
      }
      pooled <- stats::aggregate(
        cbind(n_cells, n_modified_cells, molecules_before, molecules_after, molecules_removed) ~ cell_type,
        per_sample, sum)
      pooled$fraction_cells_modified <- pooled$n_modified_cells / pmax(pooled$n_cells, 1)
      pooled$fraction_molecules_removed <- pooled$molecules_removed / pmax(pooled$molecules_before, 1)
      pooled
    },

    ensemble = function() lapply(self$corrections, function(x) x$ensemble()),

    cells = function() {
      out <- do.call(rbind, lapply(names(self$corrections), function(sample) {
        cells <- self$corrections[[sample]]$cells()
        cells$sample <- sample
        cells$orig_cell_id <- as.character(cells$cell_id)
        cells$cell_id <- paste0(sample, "_", cells$orig_cell_id)
        cells
      }))
      rownames(out) <- NULL
      out
    },

    counts = function(sample = NULL, ...) {
      if (!is.null(sample)) {
        return(self$corrections[[sample]]$counts(...))
      }
      .celladmix_bind_sample_counts(lapply(self$corrections, function(x) x$counts(...)))
    }
  )
)

#' Per-Sample Admixture Audits of a Consensus Fit
#'
#' One [CellAdmixAudit] per sample, with pooled per-pair estimates and
#' per-sample cleanup verification.
#' @export
CellAdmixCollectionAudit <- R6::R6Class(
  "CellAdmixCollectionAudit",
  public = list(
    consensus = NULL,
    audits = NULL,

    initialize = function(consensus, audits) {
      self$consensus <- consensus
      self$audits <- audits
    },

    samples = function() names(self$audits),

    pairs = function(by_sample = FALSE, detected_only = FALSE) {
      per_sample <- .celladmix_stack_samples(
        lapply(self$audits, function(audit) audit$pairs(detected_only = detected_only)),
        empty_columns = c("source", "target", "rate", "admixed_molecules", "detected"))
      if (isTRUE(by_sample) || !nrow(per_sample)) {
        if (isTRUE(by_sample)) return(per_sample)
        return(data.frame(source = character(), target = character(),
          admixed_molecules = numeric(), rate = numeric(), detected = logical(),
          n_samples = integer(), n_samples_detected = integer(), stringsAsFactors = FALSE))
      }
      # Pooled estimate per pair: admixed molecules summed, the rate taken
      # over the samples' summed target molecule counts, detection when any
      # sample detects the pair.
      per_sample$target_molecules <- ifelse(per_sample$rate > 0,
        per_sample$admixed_molecules / per_sample$rate, NA_real_)
      key <- paste(per_sample$source, per_sample$target, sep = "\r")
      pooled <- do.call(rbind, lapply(split(per_sample, key), function(rows) {
        target_molecules <- sum(rows$target_molecules, na.rm = TRUE)
        data.frame(
          source = rows$source[[1]], target = rows$target[[1]],
          admixed_molecules = sum(rows$admixed_molecules),
          rate = if (target_molecules > 0) sum(rows$admixed_molecules) / target_molecules else 0,
          detected = any(rows$detected),
          n_samples = nrow(rows),
          n_samples_detected = sum(rows$detected),
          stringsAsFactors = FALSE)
      }))
      rownames(pooled) <- NULL
      if (isTRUE(detected_only)) pooled <- pooled[pooled$detected, , drop = FALSE]
      pooled[order(-pooled$admixed_molecules), , drop = FALSE]
    },

    evaluate = function(correction, ...) {
      corrections <- if (inherits(correction, "CellAdmixCollectionCorrection")) {
        correction$corrections
      } else {
        correction
      }
      reports <- lapply(names(self$audits), function(sample) {
        self$audits[[sample]]$evaluate(corrections[[sample]], ...)
      })
      names(reports) <- names(self$audits)
      CellAdmixCollectionCleanupReport$new(reports)
    }
  )
)

#' Per-Sample Cleanup Reports
#'
#' The [CellAdmixCleanupReport] of every sample, with stacked summaries.
#' @export
CellAdmixCollectionCleanupReport <- R6::R6Class(
  "CellAdmixCollectionCleanupReport",
  public = list(
    reports = NULL,
    initialize = function(reports) self$reports <- reports,
    samples = function() names(self$reports),
    summary = function() {
      .celladmix_stack_samples(lapply(self$reports, function(report) report$summary()))
    },
    pairs = function() {
      .celladmix_stack_samples(lapply(self$reports, function(report) report$pairs()))
    }
  )
)

# ---- consensus helpers -------------------------------------------------------

# Stack per-sample data frames with a leading `sample` column; empty or
# missing tables are skipped, and the result is an empty data frame with the
# requested columns when nothing remains.
.celladmix_stack_samples <- function(tables, empty_columns = character()) {
  pieces <- lapply(names(tables), function(sample) {
    df <- tables[[sample]]
    if (is.null(df)) {
      return(NULL)
    }
    if (!is.data.frame(df)) {
      # One-row summaries arrive as named lists.
      df <- as.data.frame(df, stringsAsFactors = FALSE)
    }
    if (!nrow(df)) {
      return(NULL)
    }
    df$sample <- rep(sample, nrow(df))
    df[, c("sample", setdiff(names(df), "sample")), drop = FALSE]
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (!length(pieces)) {
    out <- data.frame(sample = character(), stringsAsFactors = FALSE)
    for (column in empty_columns) out[[column]] <- logical(0)
    return(out)
  }
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}


# A fit's gene profiles (factor x gene, gene column names) in raw gene space.
# Weighted least-squares fits store raw-space profiles. The inverse-square-root
# KL variant fits in a feature space where every gene is scaled by a weight
# recorded in the run; dividing by that weight returns the profiles to raw
# space, where profiles of different samples can be compared and averaged.
.celladmix_fit_profiles_raw <- function(fit, h = NULL) {
  run <- fit$run
  h <- as.matrix(h %||% run$h)
  if (is.null(colnames(h))) {
    colnames(h) <- run$genes
  }
  variant <- run$pipeline_options$nmf_variant %||% "ls_nmf"
  weights <- run$nmf_gene_weights
  if (identical(variant, "invsqrt_kl") && length(weights) == ncol(h)) {
    scale <- ifelse(weights > 0, 1 / weights, 0)
    h <- sweep(h, 2, scale, "*")
  }
  rownames(h) <- paste0("F", seq_len(nrow(h)))
  h
}

# Restart stability of given factors within their own sample's fit (the
# fit's matched ownership stability per factor), NA when not recorded.
.celladmix_own_factor_stability <- function(fits, samples, factors) {
  vapply(seq_along(samples), function(i) {
    stability <- fits[[samples[[i]]]]$run$nmf_factor_stability
    k <- as.integer(factors[[i]])
    if (is.null(stability) || is.na(k) || k < 1L || k > length(stability)) NA_real_ else as.numeric(stability[[k]])
  }, numeric(1))
}

# Every persisted restart profile matrix of a fit, in raw gene space.
.celladmix_fit_member_profiles <- function(fit) {
  pool <- .celladmix_ensemble_h(fit$run$path)
  lapply(pool, function(h) .celladmix_fit_profiles_raw(fit, h))
}

# Gene panel shared by all profile matrices; the panels must be identical.
.celladmix_common_genes <- function(profiles) {
  genes <- colnames(profiles[[1]])
  for (h in profiles[-1]) {
    if (!setequal(colnames(h), genes)) {
      missing <- union(setdiff(genes, colnames(h)), setdiff(colnames(h), genes))
      stop("samples do not share one gene panel (", length(missing),
        " genes differ, e.g. ", paste(utils::head(missing, 3), collapse = ", "), ")")
    }
  }
  genes
}

# Rows rescaled to sum to one.
.celladmix_normalize_profiles <- function(h) {
  sums <- rowSums(h)
  sums[sums <= 0] <- 1
  h / sums
}

# Shared profiles from per-sample profile matrices (raw gene space, common
# gene columns). The other samples' factors are matched one-to-one to the
# reference sample's factors by gene-ownership correlation; matches at or
# above `match_min` are averaged (each profile normalized to unit sum, the
# average rescaled to the weighted mean of the original row sums so that
# labeling sees the usual profile scale); unmatched factors are appended as
# sample-specific factors or dropped.
.celladmix_consensus_profiles <- function(profiles, weights, reference, match_min, unmatched) {
  genes <- colnames(profiles[[reference]])
  reference_raw <- profiles[[reference]]
  reference_norm <- .celladmix_normalize_profiles(reference_raw)
  n_factors <- nrow(reference_norm)
  w_ref <- weights[[reference]]
  sums <- reference_norm * w_ref
  weight_total <- rep(w_ref, n_factors)
  scale_total <- rowSums(reference_raw) * w_ref
  matching <- data.frame(sample = reference, sample_factor = seq_len(n_factors),
    consensus_factor = seq_len(n_factors), correlation = 1, matched = TRUE,
    stringsAsFactors = FALSE)
  appended <- list()
  for (sample in setdiff(names(profiles), reference)) {
    h_raw <- profiles[[sample]]
    h <- .celladmix_normalize_profiles(h_raw)
    match <- .celladmix_match_factors(h, reference_norm)
    for (i in seq_len(nrow(h))) {
      j <- match$match[[i]]
      r <- match$correlation[[i]]
      ok <- !is.na(j) && is.finite(r) && r >= match_min
      if (ok) {
        sums[j, ] <- sums[j, ] + weights[[sample]] * h[i, ]
        weight_total[j] <- weight_total[j] + weights[[sample]]
        scale_total[j] <- scale_total[j] + weights[[sample]] * sum(h_raw[i, ])
      } else if (identical(unmatched, "append")) {
        appended[[length(appended) + 1L]] <- list(sample = sample, factor = i,
          profile = h_raw[i, ], correlation = r)
      }
      matching <- rbind(matching, data.frame(sample = sample, sample_factor = i,
        consensus_factor = if (ok) j else NA_integer_, correlation = r, matched = ok,
        stringsAsFactors = FALSE))
    }
  }
  consensus <- (sums / weight_total) * (scale_total / weight_total)
  sample_specific <- data.frame(consensus_factor = integer(), sample = character(),
    sample_factor = integer(), stringsAsFactors = FALSE)
  for (item in appended) {
    consensus <- rbind(consensus, item$profile)
    k <- nrow(consensus)
    matching$consensus_factor[matching$sample == item$sample &
      matching$sample_factor == item$factor] <- k
    sample_specific <- rbind(sample_specific, data.frame(consensus_factor = k,
      sample = item$sample, sample_factor = item$factor, stringsAsFactors = FALSE))
  }
  dimnames(consensus) <- list(paste0("F", seq_len(nrow(consensus))), genes)
  rownames(matching) <- NULL
  list(profiles = consensus, matching = matching, sample_specific = sample_specific)
}

# Restart profiles of one sample reordered to the consensus factor index.
# Consensus factors without a counterpart in a restart take the consensus
# profile, so every member is a complete profile set in the shared index.
.celladmix_align_members <- function(pool, consensus, match_min) {
  consensus_norm <- .celladmix_normalize_profiles(consensus)
  lapply(pool, function(h) {
    h <- h[, colnames(consensus), drop = FALSE]
    match <- .celladmix_match_factors(.celladmix_normalize_profiles(h), consensus_norm)
    aligned <- consensus
    for (i in seq_len(nrow(h))) {
      j <- match$match[[i]]
      if (!is.na(j) && is.finite(match$correlation[[i]]) && match$correlation[[i]] >= match_min) {
        aligned[j, ] <- h[i, ]
      }
    }
    aligned
  })
}

# Short content hash of a numeric matrix, for run ids tied to the profiles.
.celladmix_matrix_hash <- function(x) {
  bytes <- as.numeric(serialize(signif(as.numeric(x), 8), NULL, ascii = FALSE))
  sprintf("%08x", as.integer(sum(bytes * ((seq_along(bytes) %% 251) + 1)) %% 2147483647))
}

# Pool per-sample score summaries into one summary per target/source/factor:
# pair counts summed, mean scores pair-weighted, one-sided p-values combined by
# Stouffer's method with weights proportional to sqrt(pairs).
.celladmix_pool_score_summaries <- function(summaries) {
  if (is.null(summaries) || !nrow(summaries)) {
    return(data.frame(target_cell_type = character(), source_cell_type = character(),
      factor = integer(), n_pairs = integer(), mean_score = numeric(),
      p_value = numeric(), neg_log10_p = numeric(), n_samples = integer(),
      stringsAsFactors = FALSE))
  }
  key <- paste(summaries$target_cell_type, summaries$source_cell_type, summaries$factor, sep = "\r")
  has_null <- "mean_null_score" %in% names(summaries)
  pooled <- do.call(rbind, lapply(split(summaries, key), function(rows) {
    n <- as.numeric(rows$n_pairs)
    n[!is.finite(n) | n <= 0] <- 1
    p <- pmin(pmax(as.numeric(rows$p_value), 1e-300), 1 - 1e-16)
    z <- stats::qnorm(1 - p)
    z_pooled <- sum(sqrt(n) * z) / sqrt(sum(n))
    p_pooled <- stats::pnorm(-z_pooled)
    out <- data.frame(
      target_cell_type = rows$target_cell_type[[1]],
      source_cell_type = rows$source_cell_type[[1]],
      factor = as.integer(rows$factor[[1]]),
      n_pairs = sum(as.numeric(rows$n_pairs)),
      mean_score = sum(n * rows$mean_score) / sum(n),
      stringsAsFactors = FALSE)
    if (has_null) {
      out$mean_null_score <- sum(n * rows$mean_null_score) / sum(n)
    }
    out$p_value <- p_pooled
    out$neg_log10_p <- -log10(max(p_pooled, 1e-300))
    out$n_samples <- nrow(rows)
    out
  }))
  rownames(pooled) <- NULL
  pooled[order(pooled$target_cell_type, pooled$source_cell_type, pooled$factor), , drop = FALSE]
}

# Column-bind per-sample sparse count matrices (genes x cells) with sample-
# prefixed cell ids, aligning gene rows by name.
.celladmix_bind_sample_counts <- function(counts_by_sample) {
  genes <- rownames(counts_by_sample[[1]])
  mats <- lapply(names(counts_by_sample), function(sample) {
    m <- counts_by_sample[[sample]]
    if (!setequal(rownames(m), genes)) {
      stop("samples do not share one gene panel")
    }
    m <- m[genes, , drop = FALSE]
    colnames(m) <- paste0(sample, "_", colnames(m))
    m
  })
  do.call(cbind, mats)
}
