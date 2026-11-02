---
name: celladmix-multisample
description: The multi-sample collection API of cellAdmix - processing several samples that share one annotation scheme with one shared factor space, pooled rules with per-sample support, a shared ensemble correction, per-sample audits, and the diagnostics. Use when a task involves replicates or a cohort of sections.
---

Use this skill when several samples (replicate sections, a cohort) must be corrected by the same procedure so that they do not diverge for pipeline reasons, when you meet the `CellAdmixCollection` family of objects, or when combined outputs carry sample-prefixed cell ids. `docs/multisample.md` is the authoritative description; this skill stays at the level of that page and does not describe internals (`r/R/collection.R`, `python/src/celladmix/collection.py`).

## Aim

Every sample is treated by the same procedure: one shared factor space, one set of removal rules derived from the evidence of all samples, one ensemble vote threshold, while every spatial computation (neighborhoods, candidate pairs, exposure, stain lookups) stays inside its sample and the number of molecules removed from each sample follows that sample's data. Uniform rules, not uniform removal. Samples processed separately would each arrive at their own factors, rules, and seed-dependent molecule labels, and replicates would diverge after correction for reasons unrelated to biology.

The implementation is a collection layer over ordinary per-sample datasets: each sample keeps its own input store, bundle assets, and output directory; the collection coordinates the fits and pools what has to be pooled. The samples must share one gene panel and one annotation scheme.

## Constructor

```r
datasets <- list(
  rep1 = cellAdmix("data/rep1", output_dir = "out/rep1", annotation = annotation_rep1),
  rep2 = cellAdmix("data/rep2", output_dir = "out/rep2", annotation = annotation_rep2))
collection <- cellAdmixCollection(datasets, sample_allocation = "proportional")
```

```python
collection = ca.CellAdmixCollection({"rep1": ds1, "rep2": ds2}, sample_allocation="proportional")
```

Each dataset must have an active annotation. `sample_allocation` (`"proportional"` to molecule counts, the default, or `"equal"`) controls how quantities defined per dataset (profile averaging weights, scoring candidate caps) are split across samples.

## Shared factor space

- `collection$fit(...)` / `collection.fit(**kwargs)` runs the ordinary per-sample fits (cached like any fit).
- `collection$fit_consensus(..., match_min = 0.3, min_restarts = 2L, min_recurrence = 0, name = "consensus", num_threads = NULL)` / `collection.fit_consensus(*, match_min=0.3, min_restarts=2, min_recurrence=0.0, name="consensus", num_threads=None, **fit_kwargs)` builds the consensus factorization: the factors of every restart of every sample are pooled and clustered by average linkage on gene-ownership correlation (merging while the mean correlation between clusters stays at or above `match_min`, and never merging two clusters that hold factors of the same restart); a cluster whose components come from at least `min_restarts` distinct restarts (and, with `min_recurrence`, from at least that fraction of some sample's restarts) becomes a shared factor (the weighted per-gene median of the cluster's unit-sum components, rescaled to the usual profile scale), the rest are reported as unreproduced; every sample is then refitted against the fixed shared profiles with one shared ensemble member pool (the aligned restarts of all samples, interleaved). Samples and restarts enter symmetrically, so no reference sample or merge order is involved. Returns a `CellAdmixConsensusFit`.
- `collection$fit_profiles(profiles, members = NULL, name = "profiles", match_min = 0.3, num_threads = NULL, ...)` / `collection.fit_profiles(profiles, *, members=None, ...)` projects every sample onto given profiles (a factor-by-gene matrix with gene names as columns), such as one sample's own profiles or profiles fitted elsewhere, with an optional member pool, so alternative factor spaces run through the same scoring, rules, and diagnostics.

The same fixed-profile projection is available on a single dataset: `ds$fit(nmf_fixed_h = profiles, nmf_fixed_h_members = members)` / `ds.fit(nmf_fixed_h=..., nmf_fixed_h_members=...)`.

## Consensus fit object

`CellAdmixConsensusFit` fields: `fits` (one fixed-profile fit per sample), `own_fits` (the per-sample factorizations), `profiles` (the shared factor-by-gene matrix; also `loadings()`), `factors`, `components`, `unreproduced`, `matching`, `loss`, `params`. Methods: `samples()`, `recurrence()`, `diagnostics()`, `stability()`, `sample_usage()`, `usage_divergence(by_sample = FALSE)`, `cells()`, `counts(sample = NULL)`, `audit_admixture(...)`, `score(method = c("bridge","membrane"), name, ...)`, `score_bridge(name = "bridge", ...)`, `score_membrane(name = "membrane", ...)`, `correct(score, rules = NULL, ...)`. Coherence scoring is not available at the collection level.

## Pooled scoring and rules

```r
score <- consensus$score_bridge()      # per-sample scores with each sample's own assets; pooled summary
rules <- score$rules()                 # pooled rules with sample_support and samples_significant
```

Scoring runs per sample (membrane scoring discovers each sample's own stain image), with the scorer's candidate caps (`candidate_pairs_per_type_pair`, `max_cells_per_type_pair`) split across samples by `sample_allocation`. The pooled summary combines the per-sample summaries per target type, source type, and factor (pair counts summed, mean scores pair-weighted, one-sided p-values combined by Stouffer's method weighted by the square root of each sample's pair count), and rules are derived from it exactly as from a single-sample summary. `CellAdmixCollectionScore`: `sample_summaries()`, `pairs()`, `summary()`, `annotation(p_thresh, adjust_p)`, `rules(p_thresh = 0.1, adjust_p = FALSE, targets = NULL, native_check = TRUE, native_median_thresh = 0.1, native_expr_thresh = 0.05, native_outlier_min_frac = 0.1, neighbor_k = 15L)`, `correct(...)`.

Each rule carries `sample_support`, the fraction of samples whose own summary reaches the p-value threshold for that factor and target type, with the supporting samples listed in `samples_significant`. This is a report, not a gate: a rule resting on one sample's evidence is applied everywhere, and the column says so. The native-factor check runs on the pooled rules with source exposure computed within each sample and the per-cell statistics stacked, so neighborhoods never cross samples.

## Shared ensemble correction

```r
correction <- score$correct(rules = rules)   # ensemble = NULL (10), vote = 0.3, num_threads = NULL
correction$summary()                         # removal per cell type per sample (by_sample = TRUE)
correction$counts()                          # combined; columns "<sample>_<cell_id>"
correction$counts(sample = "rep1")           # one sample, original cell ids
correction$ensemble()                        # per-sample vote details
```

The molecule-vote ensemble is identical across samples: every member of the shared pool is scored in every sample with the primary method and parameters, its summaries pooled like the primary ones, its rule set derived once (vetted by the pooled native-factor check and restricted to the primary rules' source/target pairs), and that one member rule set votes in every sample. The vote threshold is the same as for a single sample; returned rules carry `support`. Member scorings run concurrently in worker processes, as many as fit the thread count given each dataset's `num_threads` (`num_threads` on `correct()` sets the total; default the available cores). `CellAdmixCollectionCorrection`: `samples()`, `diagnostics(report = NULL)`, `summary(by_sample = TRUE)`, `ensemble()`, `cells()`, `counts(sample = NULL)`.

## Audit

```r
audit <- consensus$audit_admixture()
audit$pairs(by_sample = FALSE, detected_only = FALSE)   # pooled per-pair table, or per sample
audit$evaluate(correction)$summary()                    # per-sample cleanup verification
```

`CellAdmixCollectionAudit` runs the audit per sample and pools the per-pair estimates; `evaluate()` returns a `CellAdmixCollectionCleanupReport` with `samples()`, `summary()`, `pairs()`.

## Diagnostics

Reported with every consensus fit (all are reports, none gate anything):

- `loss`: per sample, the loss of its training rows under the consensus profiles relative to its own profiles, in the run's objective, next to the spread of its restart objectives; a consensus within that spread fits the sample as well as one of its own restarts would.
- `factors`: one row per shared factor with the number and fraction of samples contributing to its cluster (`samples`, `coverage`), the restarts carrying it (`restarts`; one component per restart at most), and the mean within-cluster ownership correlation (`agreement`). `recurrence()`: per factor and sample, how many of the sample's restarts carry the program. `unreproduced`: components whose cluster did not reach `min_restarts`. `components`: every pooled component with its cluster.
- `matching`: each sample's primary factors with the shared factor their cluster became (or none), the ownership correlation with that shared profile, their within-cluster agreement, and their restart stability within the sample (`own_stability`).
- `stability()`: within-sample restart stability next to the agreement of each sample's aligned restarts with the consensus profiles.
- `sample_usage()`, `usage_divergence(by_sample = FALSE)`: the share of each cell type's molecules labeled with each factor per sample, and per type half the L1 distance between a sample's usage and the pooled usage, maximized over samples (or per type and sample). A factor used by one sample only, or a type whose usage differs between samples, points to a sample effect in the factorization.
- `diagnostics()` on the consensus fit: one row per sample with the loss gap and restart spread (`representation_flag` when the gap exceeds the spread), the number and fraction of the sample's primary factors in the shared space (`shared_factors_flag` below one half), stability against consensus agreement, and the largest usage divergence of any cell type with its name.
- `diagnostics(report = NULL)` on the correction: per sample, the applied rules that lack evidence from that sample (absent from `samples_significant`) and the fraction of its removed molecules they account for (a removed molecule is attributed to the rule of its cell's type and primary factor label); given the audit's cleanup report, leakage removed, median pair sensitivity, and own-marker false removal per sample with `false_removal_flag` above twice the median across samples. All flags are reports.

## Cell id namespacing

Combined outputs (`cells()`, `counts()` without `sample`) prefix cell ids with the sample name and an underscore, `<sample>_<cell_id>`, the convention of Seurat's `merge()`; per-sample accessors (`counts(sample = "rep1")`) return the original ids.

## Python equivalents

```python
collection = ca.CellAdmixCollection({"rep1": ds1, "rep2": ds2})
consensus = collection.fit_consensus()
score = consensus.score_bridge(); rules = score.rules()
correction = score.correct(rules=rules)
counts, genes, cells = correction.counts()          # cells "<sample>_<cell_id>"
audit = consensus.audit_admixture(); audit.evaluate(correction).summary()
```

The Python classes `CellAdmixCollection`, `CellAdmixConsensusFit`, `CellAdmixCollectionScore`, `CellAdmixCollectionCorrection`, `CellAdmixCollectionAudit`, and `CellAdmixCollectionCleanupReport` mirror the R methods with keyword arguments. Development scripts and the working log for this feature are in `analysis/multisample/` (not part of the packages).
