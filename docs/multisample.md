# Multi-sample processing

Several samples that share one cell-type annotation scheme, such as
replicate sections or the sections of a cohort, can be processed as one
analysis. The aim is that every sample is treated by the same procedure:
one shared factor space, one set of removal rules derived from the
evidence of all samples, and one ensemble vote threshold, while every
spatial computation stays inside its sample and the number of molecules
removed from each sample follows that sample's data. Uniform rules, not
uniform removal. Samples processed separately would each arrive at their
own factors, their own rule set, and their own seed-dependent molecule
labels, and replicates would diverge after correction for reasons that
have nothing to do with their biology.

The implementation is a collection layer over ordinary per-sample
datasets. Each sample keeps its own input store, bundle assets, and
output directory; the collection coordinates the fits and pools what has
to be pooled.

## Consensus factorization

The shared factor space is a consensus of the per-sample factorizations.

1. **Per-sample fits.** Every sample is fitted as usual, with its NMF
   restarts. Cached runs are reused.
2. **Matching.** Factors are described by their gene-ownership profiles
   (each gene's loading rescaled to sum to one across factors, the
   representation used by the restart-stability diagnostic). The other
   samples' factors are matched one to one to those of a reference sample
   (by default the sample with the most molecules) by the Hungarian
   assignment on ownership correlations. Pairs whose correlation falls
   below `match_min` (default 0.3) count as unmatched.
3. **Consensus profiles.** Matched factors are averaged after each profile
   is normalized to unit sum over genes, with weights proportional to the
   samples' molecule counts (or equal, by `sample_allocation`). The
   average is rescaled to the weighted mean of the original row sums, so
   molecule labeling sees profiles on the usual scale. Unmatched factors
   are appended to the shared set and flagged as sample-specific
   (`unmatched = "append"`, the default), which keeps every sample near
   its own optimum; `"drop"` keeps only the shared factors. Profiles are
   compared in raw gene space: the inverse-square-root KL variant fits in
   a feature space where every gene is scaled by a weight recorded in the
   run, and the profiles are divided by that weight before matching.
4. **Fixed-profile projection.** Every sample is refitted with the
   consensus profiles fixed: no factorization is run, and every molecule
   is labeled against the shared profiles in their given order. Each
   sample's restart profiles, aligned to the consensus factor order,
   become the ensemble member pool of its fixed-profile run (the
   consensus itself leads as member 0, so the primary labeling and the
   vote share one factor index). The run id carries a hash of the
   profiles, so a cached projection is reused only for the consensus that
   produced it.

The same projection can be applied to profiles from any other source
through `fit_profiles()`: the profiles of one sample (reference-sample
projection) or profiles fitted elsewhere, with an optional member pool.
This is how alternative factor spaces are compared on equal footing.

## Pooled scoring and rules

Scoring runs per sample with that sample's own assets, so membrane scoring
discovers each sample's stain image. The candidate budget of the scorer
(`candidate_pairs_per_type_pair`, `max_cells_per_type_pair`) is the budget
of one single-sample analysis per sample, split across samples by
`sample_allocation`. The pooled summary combines the per-sample summaries
per target type, source type, and factor: pair counts are summed, mean
scores are pair-weighted, and the one-sided p-values are combined by
Stouffer's method with weights proportional to the square root of each
sample's pair count. Removal rules are derived from the pooled summary
exactly as from a single-sample summary.

Each rule carries `sample_support`, the fraction of samples whose own
summary reaches the p-value threshold for the rule's factor and target
type, with the supporting samples listed in `samples_significant`. This
is a report, not a gate: a rule resting on one sample's evidence is
applied everywhere, and the column says so. The native-factor check runs
on the pooled rules with source exposure computed within each sample and
the per-cell statistics stacked, so neighborhoods never cross samples.

Correction applies the pooled rules to every sample through that sample's
own ensemble members; the molecule vote and its threshold are the same as
for a single sample. The audit runs per sample, with a pooled per-pair
table and per-sample cleanup verification.

## Diagnostics

A consensus fit reports how well the shared factor space represents each
sample:

- **Loss gap** (`loss`): the loss of each sample's training rows under
  the consensus profiles relative to its own profiles, in the run's own
  objective, alongside the relative spread of the sample's restart
  objectives. A consensus within that spread fits the sample as well as
  one of its own restarts would.
- **Matching** (`matching`): every per-sample factor with its matched
  consensus factor, the ownership correlation of the match, and the
  factor's restart stability within its own sample (`own_stability`);
  `sample_specific` lists the appended factors with the same stability. An
  unmatched factor that is also unstable within its own sample is restart
  noise rather than sample-specific structure.
- **Stability** (`stability()`): each sample's within-sample restart
  stability next to the agreement of its aligned restarts with the
  consensus profiles.
- **Factor usage** (`sample_usage()`, `usage_divergence()`): the share of
  each cell type's molecules labeled with each factor, per sample, and
  per type half the L1 distance between a sample's usage and the pooled
  usage, maximized over samples. A factor used by one sample only, or a
  type whose usage differs between samples, points to a sample effect in
  the factorization.

## Usage

R:

```r
datasets <- list(
  rep1 = cellAdmix("data/rep1", output_dir = "out/rep1", annotation = annotation_rep1),
  rep2 = cellAdmix("data/rep2", output_dir = "out/rep2", annotation = annotation_rep2))
collection <- cellAdmixCollection(datasets)          # sample_allocation = "proportional"
consensus <- collection$fit_consensus()              # per-sample fits, matching, projection
consensus$loss                                       # loss gap per sample
consensus$matching                                   # factor matching table
score <- consensus$score_bridge()                    # per-sample scores, pooled summary
rules <- score$rules()                               # pooled rules with sample_support
correction <- score$correct(rules = rules)           # per-sample ensemble corrections
correction$summary()                                 # removal per cell type per sample
correction$counts()                                  # combined; columns "<sample>_<cell_id>"
correction$counts(sample = "rep1")                   # one sample, original cell ids
audit <- consensus$audit_admixture()
audit$evaluate(correction)$summary()                 # per-sample cleanup verification
```

Python:

```python
collection = ca.CellAdmixCollection({"rep1": ds1, "rep2": ds2})
consensus = collection.fit_consensus()
score = consensus.score_bridge()
rules = score.rules()
correction = score.correct(rules=rules)
counts, genes, cells = correction.counts()          # cells "<sample>_<cell_id>"
audit = consensus.audit_admixture()
audit.evaluate(correction).summary()
```

Cell ids of combined outputs are prefixed with the sample name and an
underscore, the convention of Seurat's `merge()`; per-sample accessors
return the original ids. The samples must share one gene panel.

The fixed-profile fit is also available on a single dataset:
`ds$fit(nmf_fixed_h = profiles, nmf_fixed_h_members = members)` labels the
dataset against the given factor-by-gene matrix in its given order, with
an optional list of profile matrices of the same shape as the ensemble
member pool.
