---
name: celladmix-audit-and-stability
description: The admixture audit of cellAdmix (what it estimates per cell-type pair, how it verifies a correction, how to read its summary) and the NMF restart-stability diagnostic (definition, fields, use for factor trust and rank). Use when measuring admixture, judging a cleanup, or deciding whether factors and rank are trustworthy.
---

Use this skill when you need to quantify how much admixture a dataset carries, verify whether a correction removed it without destroying native expression, interpret `audit$pairs()` or `report$summary()`, or read the per-factor stability values of a fit. The authoritative pages are `docs/audit.md` and `docs/nmf_stability.md`; the code is `r/R/audit.R`, `python/src/celladmix/audit.py`, and `src/nmf_stability.cpp` with `include/celladmix/nmf_stability.hpp`.

## The admixture audit

### Principle

No ground truth marks which molecules leaked between cells, so the audit uses the spatial structure of admixture itself. A target cell can only be contaminated by a cell type it physically borders, so the content of source-specific genes in target cells must rise with the number of source-type neighbors, while target cells with no such neighbors are an internal negative control. The audit needs only counts, cell positions, and the annotation; it never uses the factorization, so it can judge any correction, including the generative model's.

### Measurement per ordered pair (source S -> target T)

- **Exposure.** Each T cell's exposure is the number of S cells among its `neighbor_k` (15) nearest cells; T cells are binned at 0, 1, 2, 3+ neighbors.
- **Marker panel.** Up to `n_pool` (20) genes whose top expresser is S, ranked by expression in S relative to unexposed T cells. Candidates whose exposure-linked excess is disproportionate to the source expression profile (measured from S cells bordering T, with an allowance for counting noise and profile uncertainty) are likely induced by proximity rather than transferred; they are excluded, replaced by the next candidates, and listed in `markers(source, target)$induced` with `induced_stats` (`excess`, `expected`, `fold`, `z`). The **strict tier** holds genes essentially absent from reference T cells (baseline under 5% of the source level).
- **Dose-response.** Per bin, the marker rate is the marker molecules over all molecules of the bin's cells. The excess of exposed bins over the reference, times the bin totals, is the demonstrated leakage L on the panel.
- **Reference level.** Zero observed neighbors does not mean contamination-free (cells just above or below the section plane). The reference is therefore taken from T cells with zero S cells among a progressively larger neighborhood (30, 60, 120, 240 nearest cells; a monotone fit across the ladder read at its deepest well-populated rung). `reference_kind` records the neighborhood used, `reference_inflation` the ratio of the unexposed rate to that reference, `reference_trend` the decline over the ladder's final step (values well above 0 mean the estimate is still conservative).
- **Extrapolation.** The panel carries a share `coverage` of the source transcriptome; dividing by it gives the pair's admixed molecules and `rate` (fraction of all T molecules that came from S).
- **Detection.** A pair is detected when the exposed counts exceed the reference expectation by a one-sided Poisson test at Benjamini-Hochberg `q_thresh` (0.01) with at least `min_excess` (200) leaked panel molecules, given `min_target_cells` (200) and `min_reference_cells` (100).

### Running it

```r
audit <- fit$audit_admixture()     # or ds$audit_admixture() with no fit
audit$pairs(detected_only = TRUE)  # source, target, rate, admixed_molecules, excess, excess_strict,
                                   # coverage, q_value, detected, reference_kind, reference_inflation,
                                   # reference_trend, n_exposed, n_reference, n_markers, n_strict, n_induced
audit$markers("Exocrine epithelial", "Endothelial")
audit$plot_map()                                     # per-pair rate (% of target molecules), value = "molecules" for counts
audit$plot_exposure()                                # pooled excess-exposure profile with intervals
audit$plot_exposure("Exocrine epithelial", "Endothelial", correction = correction)
audit$plot_reference("Exocrine epithelial", "Endothelial")   # the reference ladder
audit$plot_remaining(list(membrane = correction))    # estimated admixture left per correction
```

Python mirrors every call (`fit.audit_admixture()`, `audit.pairs()`, `audit.plot_map()`, ...).

### Evaluating a correction

`report <- audit$evaluate(correction, warn_uncovered = TRUE)` recomputes the dose-response on the corrected counts for every detected pair, keeping bin totals at their pre-correction values and re-estimating the reference. `report$pairs()` holds `source`, `target`, `admixed_molecules`, `excess`, `excess_strict`, `sensitivity` (fraction of the exposure-linked excess removed, all panel genes) and `sensitivity_strict` (strict tier). `report$false_removal()` holds, per cell type, `own_marker_molecules` and `false_removal`, the fraction of the type's own-marker molecules in its own cells that the correction removed (an estimate of the false-positive rate on near-surely-genuine molecules). `report$summary()` returns `detected_pairs`, `estimated_admixed_molecules`, `leakage_removed_overall` (molecule-weighted sensitivity), `median_pair_sensitivity`, `own_marker_false_removal` (pooled), and in R `worst_false_removal` (the type with the largest rate). `report$plot_cleanup()` draws per-pair sensitivity bars with admixed-molecule counts overlaid.

Warnings are raised from the corrected counts themselves: a detected pair whose molecules were measurably not removed (sensitivity below 0.2 with excess at least five times `min_excess`), and severe over-removal when a type loses more than 25% of its own-marker molecules (gate `options(celladmix.own_marker_warn_min = 1000)`; Python `OWN_MARKER_WARN_MIN`).

Reading the two headline numbers: sensitivity samples leakage carried by source-specific genes, the easiest to attribute, so it is optimistic; own-marker false removal samples the molecules easiest to keep, so it is optimistic the other way. Shared genes between source and target lie in the unmeasured middle. A correction that deletes every molecule of a marker gene still reaches sensitivity 1; the damage shows up in false removal, not in sensitivity.

### What the audit established

Single-fit corrections are irreproducible across seeds (removed-molecule sets overlap by about half), the molecule-vote ensemble at 3 of 10 votes is uniformly in the upper range of single fits at 0.03-3.2% own-marker false removal, and the factorization variant should follow the scoring method (`invsqrt_kl` with membrane, `ls_nmf` with bridge). The audit also fits the generative model (`audit$fit_generative()`, `celladmix-generative`). The evaluation harness is `analysis/cleanup_benchmark/`.

## NMF restart stability

### Definition

With `nmf_n_runs > 1`, each factor of the selected restart gets a stability value. Factors are described by **gene-ownership profiles**: the loading matrix H renormalized so each gene's loadings sum to one across factors, which cancels gene abundance and makes values comparable across variants and ranks. Factors of the selected run are paired one-to-one with those of each other restart by a maximum-similarity (Hungarian) assignment on Pearson correlations of ownership profiles; a factor's stability is its mean matched correlation over the independent random restarts. Values near 1: re-found with the same gene ownership every time; near 0: nothing specific reproduces; the threshold for counting a factor as stable is 0.3 (`kStableFactorThreshold`), with values above about 0.6 solidly re-found. Two genuinely similar factors do not lower each other's score.

Restarts differ only in seed. With cluster initialization (`nmf_init = "auto"` with an annotation, or `"cluster"`), only the first restart is cluster-initialized; stability averages over the random restarts only, so it reads as "do independent starts re-find the selected solution?". The best final objective still selects the run. With `nmf_n_runs = 1` the value is a placeholder 1.0 and `stability_comparison_runs` is 0.

### Where it is recorded

The run manifest's `nmf_diagnostics` block: `selected_factor_stability` (per factor), `stable_factor_count`, `stability_threshold` (0.3), `stability_comparison_runs`, `candidate_best_match_correlations` (per restart against the selected run), `stability_metric` (`"ownership_matched"`; manifests lacking it load as `"best_match_legacy"`, whose values are not comparable and call for a refit), plus `final_objective`, `selected_seed`, `selected_run`, `candidate_final_objectives`.

```r
fit <- ds$fit(nmf_n_runs = 10)
fit$stability()            # per factor: stability and an importance measure (molecule_fraction by default)
fit$plot_stability()       # stability versus factor importance
```

```python
fit.plot_stability(min_stability=0.3)
fit.manifest["nmf_diagnostics"]["selected_factor_stability"]
```

The verbose fit log prints `stable_factors=<k>/<rank>` with the restart objective spread.

### How to use it

- Stability is a diagnostic only. It never gates scoring or correction; the ensemble correction votes across the restarts instead of trusting any one of them. Do not add filtering on stability.
- Factor trust: interpret factors below 0.3 with caution when naming factors or reading rules.
- Rank: the number of stable factors levels off at what the data resolve; a large gap between rank and stable count means the rank is too high.
- Low stability across the board on a large panel usually means undersized neighborhoods; keep the automatic `ncv_k` (it grows with the panel) rather than forcing a small value. Stability certifies reproducibility across restarts on one training subsample, not biological validity.
- `ls_nmf` gives far more reproducible factors than the KL variants (ownership correlation about 0.92 versus 0.48 between restarts) but seed dependence still enters through molecule labels, which is why the ensemble is the default for every variant.
