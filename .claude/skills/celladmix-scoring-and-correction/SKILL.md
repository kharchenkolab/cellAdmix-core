---
name: celladmix-scoring-and-correction
description: How cellAdmix scores factors for admixture evidence (bridge, membrane, coherence, factor sources), derives and thresholds removal rules, applies the native-factor check, and corrects molecules by the molecule-vote ensemble; how to read the resulting tables. Use when tuning, interpreting, or modifying scoring or correction.
---

Use this skill when choosing a scoring method, interpreting a score summary or a rule table, deciding why a rule was vetoed, tuning the ensemble vote, or changing the scoring and correction code. The authoritative definitions with formulas are in `docs/scoring_methods.md`; the code is `src/bridge.cpp`, `src/membrane.cpp`, `src/coherence.cpp`, `src/pipeline.cpp` (`apply_removal_rule`), `src/run_store.cpp` (ensemble labels and voting), `r/R/api.R` (`celladmix_score_*`, `celladmix_*_rules`), `r/R/score.R` (`.celladmix_apply_native_check`, ensemble member rules), `python/src/celladmix/score.py`, `python/src/celladmix/_score_utils.py`.

## What scoring asks

The NMF factors are molecular programs; some are native cell states, others are molecules that entered a target cell from a neighboring source cell. Scoring works on the fitted molecule labels and asks, for each ordered triple (target cell type T, source cell type S, factor f): does factor f in T cells show the spatial pattern of admixture from S? Every method writes a summary table with one row per triple holding `mean_score`, a one-sided `p_value`, and `neg_log10_p`, plus a per-pair (or per-cell) table behind it (`score$pairs()` / `score.pairs`). Two separate questions are answered: the **source identity** of a factor (a factor annotation) and the **target evidence** for cleanup (a target-specific call).

## Methods

**Bridge** (`score_bridge`; no image needed). Candidate ordered cell pairs come from centroid neighbors (`candidate_mode = "cell_center"`, the binding default; `"molecule_global"` uses molecule-level nearest neighbors). For each pair, each target molecule's `crossing_k` (20) nearest molecules within the two-cell set are queried; the score is `|p_bf - p_af| * q_abf`, the factor-fraction difference times the fraction of factor molecules whose neighborhood crosses into the source cell. Groups need `min_pairs` (5) pairs with at least `min_factor_molecules` (5) factor molecules each and are capped at `max_cells_per_type_pair` (400). The p-value is a paired one-sided Wilcoxon test against a matched null that reassigns factors on the same pairs (`compute_null`, `null_iterations = 3`). Extra columns: `total_crossing_count`, `crossing_fraction`, `mean_null_score`. Conservative: needs a visible, segmented source cell and real molecular contact.

**Membrane** (`score_membrane`; needs a membrane-like stain image). For each target molecule of factor f, the maximum normalized stain intensity along the segment to the target centroid, weighted by a source-facing directional weight, is compared in log ratio with a control molecule of the same cell at a similar radius (`control_distance_fraction` 0.025 of the cell radius). The pair score is the mean per-molecule log ratio; the summary adds `q75_score` and `fraction_positive_pairs`, and the p-value is a one-sided Wilcoxon signed-rank test against zero. Candidates are centroid neighbors (`cell_candidate_k` 50) optionally filtered by a halo around estimated radii (`cell_candidate_halo`, -1 disables), capped per type pair (`candidate_pairs_per_type_pair` 400). Image: `image_path`, `pixel_size`, `x_offset`, `y_offset` (discovered from the Xenium manifest when omitted), intensities rescaled by low/high quantiles, `line_samples` 16.

**Coherence** (`score_coherence`; R only; no source cell needed). Combines a source/factor enrichment prior over annotated types (`source_log_enrichment`, `source_probability`) with a within-cell molecular graph (`k_neighbors` 10) measuring whether source-like factor molecules form coherent patches: mean score, patch score, `active_fraction`, `largest_patch_fraction`, each against matched null subsets (`null_iterations` 20). `stain_mode = "membrane_barrier"` downweights graph edges crossing strong stain (`membrane_alpha` 3). Rules use `score_mode = "patch"` by default (`"mean"`, `"combined"`) and gates `min_delta_score`, `min_delta_patch_score`, `min_active_fraction`, `min_largest_patch_fraction`, `min_cells`, `min_molecules`.

**Factor source score** (`score_factor_sources`). Not a cleanup score: cosine similarity between a factor's gene loadings and each type's one-vs-rest marker weight vector (`min_marker_logfc`, `top_markers_per_type`, `specificity_power`, `factor_power`); reports `score`, `rank`, `best_score`, `second_score`, `margin`, `is_best` per factor x type. Use it as a source prior and to name factors.

Choosing: on stained Xenium data, membrane scoring with `nmf_variant = "invsqrt_kl"` cleans best; without a usable stain, bridge scoring with `ls_nmf` (the default fit) does (`docs/audit.md`, Table 1). `--score auto` in the batch scripts applies this choice.

## From scores to rules

`score$rules(p_thresh = 0.1, adjust_p = FALSE, targets = NULL, ...)` (R) and `score.rules(p_thresh=0.1, adjust_p=False, targets=None, ...)` (Python):

1. Source calls: per factor, the best-supported source type (bridge/membrane: from the summary collapse; coherence: highest `source_probability`).
2. Removal calls: target types whose row `(T, S*, f)` passes `p_value < p_thresh` (Benjamini-Hochberg adjusted when `adjust_p`), excluding the source type itself and any method gates. `targets` restricts to named target types.
3. One rule per removal call: `factor`, `source_cell_type`, `target_cell_type`, `p_value`, `neg_log10_p`, `rule_id` (`<factor>_<target>`).
4. The native-factor check annotates each rule (below).

## Native-factor check

A false-positive control for rules whose factor is genuinely expressed by the target type (for example a malignant-cell factor in epithelial cells). Target cells of type T are split by source exposure, counted among their `neighbor_k` (15) nearest cells: **source-distant** (no S cells) versus **source-exposed**. With `x_af` the cell-level fraction of factor f, a rule is flagged `keep = FALSE` with the reason in `native_check` when any of:

- `native_median`: the median of `x_af` over distant cells exceeds `native_median_thresh` (0.1); the factor persists without source contact.
- `native_outlier`: the fraction of distant T cells with `x_af` > `native_expr_thresh` (0.05) exceeds Q3 + 1.5 IQR of that rate over all non-source types and exceeds `native_outlier_min_frac` (0.1); the target expresses the factor distinctly in isolation.
- `no_gradient`: mean `x_af` in exposed cells minus distant cells is below 0; source proximity does not explain the signal.
- Unverifiable cases, excluded to be safe: `source_not_in_annotation`, `target_not_in_annotation`, `missing_factor_fractions` (R only), `no_distant_cells` or `no_exposed_cells` (fewer than two cells in a group), `no_expression` (the factor is absent from both groups). Passing rules read `pass`.

Evidence columns: `native_distant_n`, `native_exposed_n`, `native_distant_median`, `native_exposure_gradient`, `native_distant_expr_frac`. Flagged rules stay in the table; `correct()` drops them; heatmaps mark the call `(*)`; example selection skips them. `native_check = FALSE` disables the check. The same cell-neighborhood construction (`cell_neighbor_type_counts` in `src/clustering.cpp`) serves the audit's exposure bins.

## Molecule-vote ensemble correction

`score$correct(rules = NULL, ensemble = NULL, vote = 0.3, name = NULL, ...)`:

1. Members are the fit's NMF restarts (the `ensemble_h.parquet` pool), the selected restart first, up to `ensemble` (default 10). Each member's molecule labels are computed once by projection and within-cell smoothing and cached as `ensemble_labels_m<k>.parquet` (`ensure_ensemble_labels` in the core).
2. Each member is scored with the same method and parameters (`ensemble_member = k`), its rules derived at the same `p_thresh`, and vetted by its own native-factor check on its own cell factor fractions. When `rules` is passed explicitly, every member is restricted to those source/target pairs.
3. A molecule is removed when at least `min_votes = ceil(vote * n_members)` members remove it (default 3 of 10). `ensemble = 1`, a single-restart fit, or a run without a stored pool gives the single-fit correction with a message.
4. The returned rules carry `support`, the fraction of members that kept each source/target pair; `correction$ensemble()` reports `members`, `min_votes`, `n_removed`, and the vote histogram.

Single-fit corrections are seed-dependent (two seeds share about half their removed molecules); the vote removes that dependence and its threshold is a sensitivity/specificity dial: lower `vote` removes more, higher `vote` protects more (`docs/audit.md`, Figure 3). Majority voting discards genuine cleanup because the factor structure that cleans a given pair arises in a minority of restarts.

## Reading the outputs

- `score$summary()`: one row per (target, source, factor) with the statistics above; filter `p_value < 0.1` and inspect `mean_score` and the pair count to judge support. `score$plot_heatmap()` shows the factor-by-target matrix with `*` for removal calls and `(*)` for vetoed ones; `score$plot_pairs()` shows the top-scoring cell pairs; `score$plot_examples()` overlays molecules on stains and boundaries for selected target cells.
- Rule table: a `keep = TRUE` row is applied; read `native_check` for the veto reason and `support` for ensemble agreement.
- `correction$summary()`: per cell type, molecules before and after, fraction removed, fraction of cells modified; a removal above 60% of a type's molecules triggers a warning and usually means the factorization has no native factor for that type (raise the rank or check the annotation level).
- Verify with the audit, which measures cleanup independently of the scoring (`celladmix-audit-and-stability`).
