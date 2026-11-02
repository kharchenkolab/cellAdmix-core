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
#' `fit_consensus()` (cluster the pooled restart factors of all samples into
#' shared profiles and project every sample onto them), then
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

    fit_consensus = function(..., match_min = 0.3, min_restarts = 2L, min_recurrence = 0.5,
                             name = "consensus", num_threads = NULL) {
      state <- private$prepare_profiles(...)
      min_restarts <- max(1L, as.integer(min_restarts))
      min_recurrence <- min(1, max(0, as.numeric(min_recurrence)))
      # Every restart of every sample enters the pool; each sample's primary
      # restart is marked so the matching table can describe the sample's own
      # factors.
      pools <- lapply(state$samples, function(sample) {
        pool <- .celladmix_fit_member_profiles(state$fits[[sample]])
        lapply(pool, function(h) h[, state$genes, drop = FALSE])
      })
      names(pools) <- state$samples
      selected <- vapply(state$fits, function(fit) {
        as.integer(fit$run$nmf_diagnostics$selected_run)
      }, integer(1))
      consensus <- .celladmix_consensus_profiles(pools, state$weights, selected,
        match_min = match_min, min_restarts = min_restarts, min_recurrence = min_recurrence)
      if (!nrow(consensus$profiles)) {
        stop("no factor recurs across restarts at match_min = ", match_min,
          "; lower match_min, min_restarts, or min_recurrence")
      }
      # The matching table describes each sample's primary factors: the
      # consensus factor whose cluster they joined, their ownership correlation
      # with that consensus profile, their within-cluster agreement, and their
      # restart stability within the sample.
      matching <- do.call(rbind, lapply(state$samples, function(sample) {
        comp <- consensus$components[consensus$components$sample == sample &
          consensus$components$selected, , drop = FALSE]
        comp <- comp[order(comp$factor), , drop = FALSE]
        corr <- .celladmix_profile_correlation(pools[[sample]][[selected[[sample]]]], consensus$profiles)
        correlation <- vapply(seq_len(nrow(comp)), function(i) {
          k <- comp$consensus_factor[[i]]
          if (is.na(k)) NA_real_ else corr[comp$factor[[i]], k]
        }, numeric(1))
        data.frame(sample = sample, sample_factor = comp$factor,
          consensus_factor = comp$consensus_factor, correlation = correlation,
          matched = !is.na(comp$consensus_factor), agreement = comp$agreement,
          stringsAsFactors = FALSE)
      }))
      rownames(matching) <- NULL
      matching$own_stability <- .celladmix_own_factor_stability(
        state$fits, matching$sample, matching$sample_factor)
      # One member pool shared by every sample: the restart profiles of all
      # samples aligned to the consensus factor order and interleaved, so an
      # ensemble of any size draws from each sample alike and every member
      # votes with the same decisions in every sample. The consensus itself
      # leads the pool as member 0.
      aligned <- lapply(pools, function(pool) {
        .celladmix_align_members(pool, consensus$profiles, match_min = match_min)
      })
      shared <- .celladmix_interleave(aligned)
      members <- stats::setNames(rep(list(shared), length(state$samples)), state$samples)
      run_id <- paste0(.celladmix_clean_name(name, "run"), "_",
        .celladmix_profiles_hash(consensus$profiles, shared))
      CellAdmixConsensusFit$new(
        collection = self,
        fits = private$project(consensus$profiles, members, run_id, num_threads),
        own_fits = state$fits,
        profiles = consensus$profiles,
        matching = matching,
        factors = consensus$factors,
        components = consensus$components,
        unreproduced = consensus$unreproduced,
        loss = private$loss_table(state, consensus$profiles),
        params = list(match_min = match_min, min_restarts = min_restarts,
          min_recurrence = min_recurrence,
          weights = state$weights, restarts_total = consensus$restarts_total,
          name = name, run_id = run_id, source = "consensus")
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
      run_id <- paste0(.celladmix_clean_name(name, "run"), "_", .celladmix_profiles_hash(profiles, aligned))
      # With given profiles the per-factor table counts the samples whose
      # primary factors matched each profile.
      matched <- matching[matching$matched, , drop = FALSE]
      samples_per_factor <- tabulate(matched$consensus_factor, nbins = nrow(profiles))
      factors <- data.frame(consensus_factor = seq_len(nrow(profiles)), cluster = NA_integer_,
        samples = samples_per_factor, coverage = samples_per_factor / length(state$samples),
        restarts = NA_integer_, components = NA_integer_, max_recurrence = NA_real_,
        agreement = NA_real_, stringsAsFactors = FALSE)
      components <- data.frame(sample = matching$sample, restart = NA_integer_,
        factor = matching$sample_factor, group = NA_integer_, selected = TRUE,
        cluster = NA_integer_, agreement = NA_real_, consensus_factor = matching$consensus_factor,
        stringsAsFactors = FALSE)
      CellAdmixConsensusFit$new(
        collection = self,
        fits = private$project(profiles, member_pools, run_id, num_threads),
        own_fits = state$fits,
        profiles = profiles,
        matching = matching,
        factors = factors,
        components = components,
        unreproduced = components[0, , drop = FALSE],
        loss = private$loss_table(state, profiles),
        params = list(match_min = match_min, name = name, run_id = run_id, source = "profiles")
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
#' profiles, the pooled clustering that produced them, one fixed-profile fit
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
    factors = NULL,
    components = NULL,
    unreproduced = NULL,
    loss = NULL,
    params = NULL,

    initialize = function(collection, fits, own_fits, profiles, matching, factors,
                          components, unreproduced, loss, params = list()) {
      self$collection <- collection
      self$fits <- fits
      self$own_fits <- own_fits
      self$profiles <- profiles
      self$matching <- matching
      self$factors <- factors
      self$components <- components
      self$unreproduced <- unreproduced
      self$loss <- loss
      self$params <- params
    },

    print = function(...) {
      cat("cellAdmix consensus fit\n")
      cat("  samples:", paste(names(self$fits), collapse = ", "), "\n")
      single <- sum(self$factors$samples == 1L, na.rm = TRUE)
      cat("  factors:", nrow(self$profiles),
        if (single) paste0("(", single, " from one sample only)") else "", "\n")
      if (identical(self$params$source, "consensus")) {
        cat("  sample coverage per factor: median",
          round(stats::median(self$factors$coverage), 2),
          "; unreproduced components:", nrow(self$unreproduced), "\n")
      }
      if (nrow(self$matching)) {
        cat("  primary factors in the shared space:",
          sum(self$matching$matched), "of", nrow(self$matching),
          "; ownership correlation median",
          round(stats::median(self$matching$correlation, na.rm = TRUE), 3), "\n")
      }
      invisible(self)
    },

    samples = function() names(self$fits),
    loadings = function() self$profiles,

    recurrence = function() {
      # Per consensus factor and sample: the number of that sample's restarts
      # with a component in the factor's cluster, and that number as a
      # fraction of the sample's restarts. Zero rows mark samples that do not
      # carry the factor at all.
      comp <- self$components[!is.na(self$components$consensus_factor), , drop = FALSE]
      totals <- self$params$restarts_total
      grid <- expand.grid(consensus_factor = seq_len(nrow(self$profiles)),
        sample = names(self$fits), stringsAsFactors = FALSE)
      grid$restarts <- vapply(seq_len(nrow(grid)), function(i) {
        rows <- comp$consensus_factor == grid$consensus_factor[[i]] & comp$sample == grid$sample[[i]]
        length(unique(comp$restart[rows]))
      }, integer(1))
      grid$restarts_total <- if (is.null(totals)) NA_integer_ else as.integer(totals[grid$sample])
      grid$recurrence <- grid$restarts / grid$restarts_total
      grid[order(grid$consensus_factor, grid$sample), , drop = FALSE]
    },

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

    usage_divergence = function(by_sample = FALSE) {
      # Per cell type: half the L1 distance between a sample's factor usage
      # and the pooled usage, maximized over samples (or, with `by_sample`,
      # one row per cell type and sample). Zero when every sample labels the
      # type's molecules alike.
      usage <- self$sample_usage()
      pooled <- stats::aggregate(molecules ~ cell_type + factor, usage, sum)
      totals <- stats::aggregate(molecules ~ cell_type, pooled, sum)
      pooled$share <- pooled$molecules / totals$molecules[match(pooled$cell_type, totals$cell_type)]
      key <- function(df) paste(df$cell_type, df$factor, sep = "\r")
      usage$pooled_share <- pooled$share[match(key(usage), key(pooled))]
      usage$deviation <- abs(usage$share - usage$pooled_share)
      per_sample <- stats::aggregate(deviation ~ cell_type + sample, usage, function(x) sum(x) / 2)
      names(per_sample)[names(per_sample) == "deviation"] <- "divergence"
      per_sample$molecules <- totals$molecules[match(per_sample$cell_type, totals$cell_type)]
      if (isTRUE(by_sample)) {
        out <- per_sample[order(per_sample$sample, -per_sample$divergence), , drop = FALSE]
        rownames(out) <- NULL
        return(out)
      }
      out <- stats::aggregate(divergence ~ cell_type, per_sample, max)
      out$molecules <- totals$molecules[match(out$cell_type, totals$cell_type)]
      out[order(-out$divergence), , drop = FALSE]
    },

    diagnostics = function() {
      # One row per sample: how well the shared factor space represents it
      # (the training loss under the shared profiles relative to its own fit,
      # with the sample's restart spread as the yardstick), how many of its
      # primary factors entered the shared space, its restart stability next
      # to its agreement with the consensus, and the largest factor-usage
      # divergence of any of its cell types. Flags are reports:
      # `representation_flag` when the loss gap exceeds the restart spread,
      # `shared_factors_flag` when fewer than half of the sample's primary
      # factors are in the shared space.
      first <- function(x) if (length(x)) x[[1]] else NA
      loss <- self$loss
      stability <- self$stability()
      deviation <- self$usage_divergence(by_sample = TRUE)
      out <- do.call(rbind, lapply(names(self$fits), function(sample) {
        m <- self$matching[self$matching$sample == sample, , drop = FALSE]
        l <- loss[loss$sample == sample, , drop = FALSE]
        st <- stability[stability$sample == sample, , drop = FALSE]
        d <- deviation[deviation$sample == sample, , drop = FALSE]
        worst <- if (nrow(d)) d[which.max(d$divergence), , drop = FALSE] else NULL
        gap <- as.numeric(first(l$loss_gap))
        spread <- as.numeric(first(l$restart_spread))
        shared <- sum(m$matched)
        data.frame(
          sample = sample,
          loss_gap = gap,
          restart_spread = spread,
          representation_flag = isTRUE(gap > spread),
          primary_factors = nrow(m),
          primary_factors_shared = shared,
          shared_fraction = if (nrow(m)) shared / nrow(m) else NA_real_,
          shared_factors_flag = nrow(m) > 0 && shared / nrow(m) < 0.5,
          within_sample_stability = as.numeric(first(st$within_sample_stability)),
          consensus_agreement = as.numeric(first(st$consensus_agreement)),
          max_usage_divergence = if (is.null(worst)) NA_real_ else worst$divergence,
          divergent_cell_type = if (is.null(worst)) NA_character_ else worst$cell_type,
          stringsAsFactors = FALSE)
      }))
      rownames(out) <- NULL
      out
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
#'
#' `correct()` applies the rules with a molecule-vote ensemble shared across
#' samples: every ensemble member is scored in every sample, its evidence is
#' pooled like the primary evidence, and one member rule set votes in every
#' sample. The member scorings run concurrently, one (sample, member) job per
#' worker process, with as many workers as fit the thread budget given each
#' dataset's `num_threads`; `num_threads` sets the budget (default: the number
#' of available cores).
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
                       targets = NULL, ensemble = NULL, vote = 0.3, num_threads = NULL, ...) {
      if (is.null(rules)) {
        rules <- self$rules(p_thresh = p_thresh, adjust_p = adjust_p, targets = targets)
      }
      kept <- rules
      if ("keep" %in% names(kept)) {
        kept <- kept[is.na(kept$keep) | kept$keep, , drop = FALSE]
      }
      if (!nrow(kept)) {
        stop("no rules to apply")
      }
      name <- .celladmix_clean_name(name %||% paste0(self$name, "_clean"), "correction")
      fits <- self$consensus$fits
      samples <- names(fits)

      # The molecule-vote ensemble is shared across samples: the member pool
      # is the same in every sample's run, each member is scored in every
      # sample, its evidence pooled like the primary evidence, and one member
      # rule set votes in every sample.
      ensemble <- if (is.null(ensemble)) 10L else max(1L, as.integer(ensemble))
      member_rules <- NULL
      members <- integer(0)
      if (ensemble > 1L) {
        available <- min(vapply(samples, function(sample) {
          .celladmix_ensemble_prepare(fits[[sample]]$run$path, fits[[sample]]$dataset$num_threads)
        }, numeric(1)))
        if (available < 2L) {
          message("Runs carry no shared ensemble member pool; applying single-fit corrections.")
        } else {
          selected <- as.integer(fits[[1]]$run$nmf_diagnostics$selected_run) - 1L
          members <- unique(c(selected, setdiff(seq_len(available) - 1L, selected)))
          members <- members[seq_len(min(ensemble, length(members)))]
          member_rules <- private$pooled_member_rules(members, selected, kept,
            p_thresh = p_thresh, adjust_p = adjust_p, targets = targets,
            num_threads = num_threads)
        }
      }
      min_votes <- max(1L, as.integer(ceiling(vote * length(members))))
      corrections <- lapply(samples, function(sample) {
        fit <- fits[[sample]]
        score_sample <- self$sample_scores[[sample]]
        out_dir <- file.path(fit$run$paths$corrected_dir, name)
        corrected <- if (is.null(member_rules)) {
          celladmix_correct(fit$run, rules = kept, out_dir = out_dir,
            annotation = score_sample$annotation_vector())
        } else {
          celladmix_correct(fit$run, rules = member_rules$rules, out_dir = out_dir,
            annotation = score_sample$annotation_vector(),
            rule_member = member_rules$member, min_votes = min_votes)
        }
        correction <- CellAdmixCorrection$new(name, fit, score_sample, corrected,
          params = list(p_thresh = p_thresh, adjust_p = adjust_p, targets = targets,
            ensemble = max(1L, length(members)),
            vote = if (length(members) > 1L) vote else NA_real_,
            pooled_members = !is.null(member_rules)),
          rules = rules)
        fit$corrections_registry[[name]] <- correction
        correction$save_metadata()
        correction
      })
      names(corrections) <- samples
      if (!is.null(member_rules)) {
        rules$support <- .celladmix_rule_support(rules, member_rules$pair_sets)
        removed <- vapply(corrections, function(x) as.numeric(x$run$n_removed %||% NA_real_), numeric(1))
        message(sprintf(
          "Pooled ensemble correction over %d members shared by %d samples: molecules removed by >= %d members (vote >= %.2f) are dropped (%s).",
          length(members), length(samples), min_votes, vote,
          paste(sprintf("%s: %s", samples, format(removed, big.mark = ",")), collapse = ", ")))
      }
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

    # Score the given ensemble members in every sample. The (sample, member)
    # jobs are independent and run in forked worker processes, as many as the
    # thread budget allows given each dataset's own thread count; each worker
    # writes its summary to the run's cache, which the parent then reads.
    score_members = function(members, num_threads = NULL) {
      samples <- names(self$sample_scores)
      jobs <- expand.grid(sample = samples, member = members, stringsAsFactors = FALSE)
      if (!nrow(jobs)) {
        return(invisible(NULL))
      }
      per_job <- max(vapply(self$consensus$fits, function(fit) {
        as.integer(fit$dataset$num_threads %||% 1L)
      }, integer(1)))
      workers <- .celladmix_parallel_workers(num_threads, per_job, nrow(jobs))
      run_job <- function(i) {
        private$member_summary(jobs$sample[[i]], jobs$member[[i]])
        NULL
      }
      if (workers <= 1L) {
        for (i in seq_len(nrow(jobs))) run_job(i)
        return(invisible(NULL))
      }
      results <- parallel::mclapply(seq_len(nrow(jobs)), run_job,
        mc.cores = workers, mc.preschedule = FALSE)
      failed <- vapply(results, inherits, logical(1), "try-error")
      if (any(failed)) {
        first <- results[[which(failed)[[1]]]]
        stop("scoring of ensemble members failed: ",
          conditionMessage(attr(first, "condition")), call. = FALSE)
      }
      invisible(NULL)
    },

    # One ensemble member's score summary in one sample, with the same method
    # and parameters as the primary score; persisted in the sample's run so
    # re-corrections reuse it.
    member_summary = function(sample, member) {
      fit <- self$consensus$fits[[sample]]
      score_sample <- self$sample_scores[[sample]]
      key <- paste(self$method, member, score_sample$annotation_hash,
        paste(deparse(score_sample$params), collapse = ""), sep = "|")
      path <- file.path(fit$run$paths$scores_dir,
        sprintf("collection_member_summary_%s_m%d.rds", self$method, member))
      if (file.exists(path)) {
        stored <- tryCatch(readRDS(path), error = function(e) NULL)
        if (!is.null(stored) && identical(stored$key, key)) {
          return(stored$summary)
        }
      }
      summary <- .celladmix_score_member(fit, score_sample, member)$summary
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(list(key = key, summary = summary), path)
      summary
    },

    # Per-member rule sets for the shared molecule-vote ensemble: each
    # member's summaries from every sample are pooled like the primary
    # summaries, thresholded, vetted by the pooled native-factor check on
    # the member's own cell factor fractions, and restricted to the primary
    # rules' source -> target pairs. Returns the stacked (factor, target)
    # rows with their member indices and each member's kept pair set.
    pooled_member_rules = function(members, selected, primary_rules, p_thresh,
                                   adjust_p, targets, num_threads = NULL) {
      private$score_members(setdiff(members, selected), num_threads)
      pair_key <- function(df) paste(df$source_cell_type, df$target_cell_type, sep = "\r")
      restrict <- unique(pair_key(primary_rules))
      frames <- list()
      frame_members <- integer(0)
      pair_sets <- list()
      for (member in members) {
        rules_member <- if (member == selected) {
          primary_rules
        } else {
          summaries <- .celladmix_stack_samples(
            lapply(stats::setNames(nm = names(self$sample_scores)), function(sample) {
              private$member_summary(sample, member)
            }))
          pooled <- .celladmix_pool_score_summaries(summaries)
          rules_member <- celladmix_bridge_rules(pooled, p_thresh = p_thresh,
            adjust_p = adjust_p, target_cell_types = targets)
          rules_member <- private$pooled_native_check(rules_member,
            cell_factors_fn = function(fit) .celladmix_member_cell_factors(fit, member))
          rules_member[is.na(rules_member$keep) | rules_member$keep, , drop = FALSE]
        }
        if (nrow(rules_member)) {
          rules_member <- rules_member[pair_key(rules_member) %in% restrict, , drop = FALSE]
        }
        pair_sets[[as.character(member)]] <- unique(pair_key(rules_member))
        if (nrow(rules_member)) {
          frames[[length(frames) + 1L]] <- data.frame(
            factor = as.integer(rules_member$factor),
            target_cell_type = as.character(rules_member$target_cell_type),
            stringsAsFactors = FALSE)
          frame_members <- c(frame_members, rep(as.integer(member), nrow(rules_member)))
        }
      }
      list(
        rules = if (length(frames)) do.call(rbind, frames) else
          data.frame(factor = integer(0), target_cell_type = character(0), stringsAsFactors = FALSE),
        member = frame_members,
        pair_sets = pair_sets)
    },

    pooled_native_check = function(rules, median_thresh = 0.1, expr_thresh = 0.05,
                                   outlier_min_frac = 0.1, neighbor_k = 15L,
                                   cell_factors_fn = function(fit) fit$cell_factors()) {
      # Source exposure is computed within each sample (its own cell
      # neighborhoods) and the per-cell tables are stacked, so the check's
      # statistics pool the samples without neighborhoods crossing them.
      # `cell_factors_fn` supplies the per-sample cell factor fractions: the
      # run's own, or one ensemble member's.
      fits <- self$consensus$fits
      annotations <- lapply(fits, function(fit) fit$dataset$annotation(fit$annotation_name, as_vector = TRUE))
      types <- sort(unique(unlist(lapply(annotations, function(ann) as.character(ann[!is.na(ann)])), use.names = FALSE)))
      pieces <- lapply(names(fits), function(sample) {
        cells <- cell_factors_fn(fits[[sample]])
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
#' set, with per-sample and pooled removal summaries, combined or per-sample
#' corrected counts, and per-sample diagnostics (`diagnostics()`): the rules
#' applied to a sample without evidence from that sample, the share of its
#' removed molecules they account for, and, given a cleanup report, the
#' audit's leakage removed and own-marker false removal.
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

    diagnostics = function(report = NULL) {
      # One row per sample: the number of applied rules that lack local
      # evidence (the sample is absent from the rule's `samples_significant`),
      # the fraction of the sample's removed molecules attributable to those
      # rules (a removed molecule is attributed to the rule of its cell's type
      # and its primary factor label), and, when a cleanup report is given,
      # the audit's leakage removed, median pair sensitivity, and own-marker
      # false removal, with `false_removal_flag` when a sample's false removal
      # exceeds twice the median across samples. Flags are reports.
      rules <- self$rules
      if (!is.null(rules) && "keep" %in% names(rules)) {
        rules <- rules[is.na(rules$keep) | rules$keep, , drop = FALSE]
      }
      fits <- self$score$consensus$fits
      out <- do.call(rbind, lapply(names(self$corrections), function(sample) {
        removed <- .celladmix_removed_by_type_factor(fits[[sample]], self$corrections[[sample]])
        total <- sum(removed$removed)
        if (is.null(rules) || !nrow(rules)) {
          return(data.frame(sample = sample, rules_applied = 0L, rules_without_local_evidence = 0L,
            molecules_removed = total, removed_without_local_evidence = NA_real_,
            stringsAsFactors = FALSE))
        }
        significant <- as.character(rules$samples_significant %||% rep("", nrow(rules)))
        significant[is.na(significant)] <- ""
        has_local <- vapply(strsplit(significant, ",", fixed = TRUE), function(x) sample %in% x, logical(1))
        without <- rules[!has_local, , drop = FALSE]
        attributed <- removed$removed[match(
          paste(without$target_cell_type, as.integer(without$factor), sep = "\r"),
          paste(removed$cell_type, removed$factor, sep = "\r"))]
        data.frame(sample = sample, rules_applied = nrow(rules),
          rules_without_local_evidence = nrow(without), molecules_removed = total,
          removed_without_local_evidence = sum(attributed, na.rm = TRUE) / max(total, 1),
          stringsAsFactors = FALSE)
      }))
      if (!is.null(report)) {
        s <- report$summary()
        cols <- intersect(c("leakage_removed_overall", "median_pair_sensitivity",
          "own_marker_false_removal"), names(s))
        for (col in cols) out[[col]] <- s[[col]][match(out$sample, s$sample)]
        if ("own_marker_false_removal" %in% cols) {
          reference <- stats::median(out$own_marker_false_removal, na.rm = TRUE)
          out$false_removal_flag <- nrow(out) > 1L & !is.na(out$own_marker_false_removal) &
            out$own_marker_false_removal > 2 * reference
        }
      }
      rownames(out) <- NULL
      out
    },

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

# Molecules removed from a sample per annotated cell type and primary factor
# label: each cell's per-factor molecule counts before the correction (factor
# fractions times molecule count in the fit) minus those after it.
.celladmix_removed_by_type_factor <- function(fit, correction) {
  before <- fit$cell_factors()
  after <- correction$cells()
  annotation <- fit$dataset$annotation(fit$annotation_name, as_vector = TRUE)
  fraction_cols <- grep("^factor_[0-9]+_fraction$", names(before), value = TRUE)
  fraction_cols <- fraction_cols[order(as.integer(sub("^factor_([0-9]+)_fraction$", "\\1", fraction_cols)))]
  mass_before <- as.matrix(before[, fraction_cols, drop = FALSE]) * before$transcript_count
  idx <- match(as.character(before$cell_id), as.character(after$cell_id))
  mass_after <- as.matrix(after[idx, fraction_cols, drop = FALSE]) * after$transcript_count[idx]
  mass_after[is.na(mass_after)] <- 0
  removed <- pmax(mass_before - mass_after, 0)
  type <- as.character(annotation[as.character(before$cell_id)])
  keep <- !is.na(type)
  if (!any(keep)) {
    return(data.frame(cell_type = character(), factor = integer(), removed = numeric(),
      stringsAsFactors = FALSE))
  }
  by_type <- rowsum(removed[keep, , drop = FALSE], type[keep])
  data.frame(
    cell_type = rep(rownames(by_type), ncol(by_type)),
    factor = rep(seq_len(ncol(by_type)), each = nrow(by_type)),
    removed = as.vector(by_type),
    stringsAsFactors = FALSE)
}

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

# Weighted median: the midpoint of the lower weighted median (the smallest
# value whose cumulative weight reaches half the total) and the upper one
# (the largest value whose cumulative weight from above reaches half), which
# reduces to the ordinary median for equal weights. The Python binding uses
# the same definition.
.celladmix_weighted_median <- function(x, w) {
  o <- order(x)
  x <- x[o]
  w <- w[o]
  half <- sum(w) / 2
  lower <- x[[which(cumsum(w) >= half)[[1]]]]
  upper <- x[[rev(which(rev(cumsum(rev(w))) >= half))[[1]]]]
  (lower + upper) / 2
}

# Rows rescaled to sum to one.
.celladmix_normalize_profiles <- function(h) {
  sums <- rowSums(h)
  sums[sums <= 0] <- 1
  h / sums
}

# Shared profiles by pooled clustering of every restart of every sample.
# `pools` is a named list (sample) of lists (restart) of raw-space profile
# matrices over common gene columns; `selected` gives each sample's primary
# restart (1-based). The pooled components are clustered by average linkage
# on gene-ownership correlation, cut at `match_min`; a cluster whose
# components come from at least `min_restarts` distinct restarts, of any
# samples, and from at least a fraction `min_recurrence` of the restarts of
# some sample, becomes a consensus factor, and the remaining components are
# reported as unreproduced. The consensus profile of a cluster is the
# weighted per-gene median of its unit-sum component profiles (every
# component weighted by its sample's weight, so a program found by more of a
# sample's restarts carries more of that sample's evidence), renormalized to
# unit sum and rescaled to the weighted mean of the original row sums so
# that labeling sees the usual profile scale; the median keeps the profile
# of a loosely agreeing cluster from blurring. Factors are ordered by the number
# of samples covered, then restarts, then profile mass.
.celladmix_consensus_profiles <- function(pools, weights, selected, match_min, min_restarts,
                                          min_recurrence = 0) {
  genes <- colnames(pools[[1]][[1]])
  blocks <- list()
  meta <- list()
  group <- 0L
  for (sample in names(pools)) {
    for (r in seq_along(pools[[sample]])) {
      h <- pools[[sample]][[r]]
      group <- group + 1L
      blocks[[group]] <- h
      meta[[group]] <- data.frame(sample = sample, restart = r, factor = seq_len(nrow(h)),
        group = group, selected = identical(as.integer(r), as.integer(selected[[sample]])),
        stringsAsFactors = FALSE)
    }
  }
  stacked <- do.call(rbind, blocks)
  components <- do.call(rbind, meta)
  rownames(components) <- NULL
  clustering <- .celladmix_cluster_factors(stacked, components$group, match_min)
  components$cluster <- as.integer(clustering$cluster)
  components$agreement <- as.numeric(clustering$within_correlation)
  row_sums <- rowSums(stacked)
  restarts_total <- vapply(pools, length, integer(1))
  summary <- do.call(rbind, lapply(sort(unique(components$cluster)), function(k) {
    idx <- which(components$cluster == k)
    per_sample <- tapply(components$group[idx], components$sample[idx], function(g) length(unique(g)))
    data.frame(cluster = k,
      samples = length(unique(components$sample[idx])),
      restarts = length(unique(components$group[idx])),
      components = length(idx),
      max_recurrence = max(as.numeric(per_sample) / as.numeric(restarts_total[names(per_sample)])),
      agreement = if (length(idx) > 1L) mean(components$agreement[idx], na.rm = TRUE) else NA_real_,
      mass = mean(row_sums[idx]),
      stringsAsFactors = FALSE)
  }))
  summary <- summary[order(-summary$samples, -summary$restarts, -summary$mass), , drop = FALSE]
  reproduced <- summary[summary$restarts >= min_restarts & summary$max_recurrence >= min_recurrence, , drop = FALSE]
  consensus <- matrix(0, nrow(reproduced), length(genes),
    dimnames = list(if (nrow(reproduced)) paste0("F", seq_len(nrow(reproduced))) else NULL, genes))
  normalized <- .celladmix_normalize_profiles(stacked)
  for (i in seq_len(nrow(reproduced))) {
    idx <- which(components$cluster == reproduced$cluster[[i]])
    w <- as.numeric(weights[components$sample[idx]])
    profile <- apply(normalized[idx, , drop = FALSE], 2, .celladmix_weighted_median, w = w)
    profile <- profile / max(sum(profile), .Machine$double.eps)
    consensus[i, ] <- profile * (sum(w * row_sums[idx]) / sum(w))
  }
  lookup <- stats::setNames(seq_len(nrow(reproduced)), reproduced$cluster)
  components$consensus_factor <- unname(lookup[as.character(components$cluster)])
  factors <- data.frame(consensus_factor = seq_len(nrow(reproduced)), cluster = reproduced$cluster,
    samples = reproduced$samples, coverage = reproduced$samples / length(pools),
    restarts = reproduced$restarts, components = reproduced$components,
    max_recurrence = reproduced$max_recurrence, agreement = reproduced$agreement,
    stringsAsFactors = FALSE)
  rownames(factors) <- NULL
  unreproduced <- components[is.na(components$consensus_factor), , drop = FALSE]
  rownames(unreproduced) <- NULL
  list(profiles = consensus, components = components, factors = factors,
    unreproduced = unreproduced, restarts_total = restarts_total)
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

# Number of concurrent worker processes for independent jobs that each use
# `per_job` threads, within a total thread budget (default: all cores).
.celladmix_parallel_workers <- function(num_threads, per_job, n_jobs) {
  if (identical(.Platform$OS.type, "windows")) {
    return(1L)
  }
  budget <- if (is.null(num_threads)) parallel::detectCores() else as.integer(num_threads)
  if (is.na(budget) || budget < 1L) {
    budget <- 1L
  }
  max(1L, min(as.integer(n_jobs), budget %/% max(1L, as.integer(per_job))))
}

# Round-robin interleaving of several lists (a list per sample) into one.
.celladmix_interleave <- function(lists) {
  lists <- lists[vapply(lists, length, integer(1)) > 0]
  if (!length(lists)) return(list())
  out <- list()
  for (i in seq_len(max(vapply(lists, length, integer(1))))) {
    for (items in lists) {
      if (i <= length(items)) out[[length(out) + 1L]] <- items[[i]]
    }
  }
  out
}

# Content hash of shared profiles together with their member pool, so a
# cached fixed-profile run is reused only for the same profiles and members.
.celladmix_profiles_hash <- function(profiles, members = list()) {
  stacked <- if (length(members)) rbind(profiles, do.call(rbind, members)) else profiles
  .celladmix_matrix_hash(stacked)
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
