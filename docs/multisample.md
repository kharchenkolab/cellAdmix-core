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

The shared factor space is a consensus of the per-sample factorizations,
built from every restart of every sample.

1. **Per-sample fits.** Every sample is fitted as usual, with its NMF
   restarts. Cached runs are reused.
2. **Pooled clustering.** The factors of every restart of every sample
   are pooled and described by their gene-ownership profiles (each gene's
   loading rescaled to sum to one across the factors of its own restart,
   the representation used by the restart-stability diagnostic). The
   pooled factors are clustered by average linkage on the correlation of
   these profiles, merging while the mean correlation between clusters
   stays at or above `match_min` (default 0.3) and never merging two
   clusters that both hold a factor of the same restart: one restart finds
   a program once, so a cluster is the one-to-one correspondence of a
   program across restarts and samples, and two distinct programs that a
   restart finds side by side stay separate however similar they are.
   Samples and restarts enter symmetrically, so the result does not depend
   on a reference sample or on a merge order, and samples fitted at
   different ranks pose no difficulty.
3. **Consensus profiles.** A cluster becomes a shared factor when its
   members come from at least `min_restarts` distinct restarts (default
   2), of any samples, and from at least a fraction `min_recurrence`
   (default 0.5) of the restarts of some sample: a program found by a
   single random start is not reproducible, and one found by a minority
   of every sample's restarts is an alternative decomposition rather than
   the solution any sample settles on. With the default, the shared
   factors are the programs that are the majority solution of at least
   one sample, which keeps the shared space close to the per-sample rank;
   `min_recurrence = 0` admits every program that recurs at all. The
   excluded components are reported rather than used. The shared profile is
   the per-gene median of the cluster's components after each is
   normalized to unit sum, every component weighted by its sample's weight
   (molecule count, or equal, by `sample_allocation`), so that a program
   found by more of a sample's restarts carries more of that sample's
   evidence; the median profile is renormalized to unit sum and rescaled
   to the weighted mean of the original row sums so that molecule labeling
   sees profiles on the usual scale. The median
   rather than the mean keeps the profile of a loosely agreeing cluster
   from blurring, which matters when a sample's restarts disagree among
   themselves. Factors are ordered by the number of samples they cover,
   then by the number of restarts. Profiles are compared in raw gene
   space: the inverse-square-root KL variant fits in a feature space where
   every gene is scaled by a weight recorded in the run, and the profiles
   are divided by that weight before pooling.
4. **Fixed-profile projection.** Every sample is refitted with the
   consensus profiles fixed: no factorization is run, and every molecule
   is labeled against the shared profiles in their given order. The
   restart profiles of all samples, aligned to the consensus factor
   order and interleaved sample by sample, form one ensemble member pool
   that every sample's run carries (the consensus itself leads as member
   0, so the primary labeling and the vote share one factor index). The
   run id carries a hash of the profiles and the pool, so a cached
   projection is reused only for the consensus that produced it.

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

Every molecule of every sample is labeled with the shared program that
fits it best, so a program absent from a sample still receives labels
there, and those labels mark the cells' own molecules rather than the
program. A rule on such a program, derived from the samples that do have
it, removes native expression in the samples that do not. The rules carry
`samples_without_program`, the samples none of whose own NMF restarts
carries the rule's program. By default every rule is applied in every
sample (`apply_rules = "everywhere"` on `correct()`), which treats
replicates identically and is right when every program is present in
every sample. With `apply_rules = "where_present"`, a rule is applied in a
sample only when at least one of that sample's own restarts carries the
program, so a sample is not cleaned of a program its own data never finds;
with true replicates nothing changes, and with regional sections or a
cohort the rules on regional programs are confined to the samples that
have them. The condition rests on the sample's factorization, not on the
correction outcome, and the skipped rules are counted per sample in the
correction diagnostics. The option trades cleanup for conservatism: a
program too weak for any restart of a sample to find may still contribute
admixture there, which is then left in place, and samples that happen to
carry different programs are cleaned differently. On eight regional tiles
cut from two breast replicates, confining the rules lowered own-marker
false removal from 4.1% to 2.7% but also lowered the leakage removed (0.69
to 0.65) and the median pair sensitivity (0.91 to 0.83), and the
concordance between the replicates fell from 0.991 to 0.938, below the
uncorrected 0.987, because corresponding regions of the two replicates
were no longer cleaned alike. Uniform application is therefore the
default; the option suits cohorts whose samples are not replicates of one
another and whose composition differences are real.

Correction keeps the molecule-vote ensemble identical across samples,
because the member-level scoring is the most seed-dependent step of the
pipeline and pooling only the primary rules would leave it per sample.
Every member of the shared pool is scored in every sample with the
primary method and parameters, its summaries are pooled exactly like the
primary ones, its rule set is derived once (vetted by the pooled
native-factor check on the member's own cell factor fractions and
restricted to the primary rules' source → target pairs), and that one
member rule set votes in every sample. The vote threshold is the same as
for a single sample, and the returned rules carry `support`, the fraction
of members keeping each pair.

The per-sample steps of the correction are independent jobs that run
concurrently in forked worker processes: labeling each sample's molecules
with every member of the pool (once per run, then reused), scoring every
member in every sample (one job per sample and member), and the vote in
each sample. The number of jobs running at a time is a thread budget
divided by the threads one job uses (the dataset's `num_threads` for
labeling and scoring, one for the vote), never more than the number of
jobs. `num_threads` on `correct()` sets the budget; by default it is the
`CELLADMIX_NUM_THREADS` environment variable, or else the number of
available cores. Each concurrent vote holds one sample's molecule table in
memory, so a smaller budget also bounds memory use. The audit runs per
sample, with a pooled per-pair table and per-sample cleanup verification.

## Diagnostics

A consensus fit reports how well the shared factor space represents each
sample:

- **Loss gap** (`loss`): the loss of each sample's training rows under
  the consensus profiles relative to its own profiles, in the run's own
  objective, alongside the relative spread of the sample's restart
  objectives. A consensus within that spread fits the sample as well as
  one of its own restarts would.
- **Factors** (`factors`): one row per shared factor with the number and
  fraction of samples contributing to its cluster (`samples`, `coverage`),
  the number of restarts carrying it (`restarts`, equal to the number of
  components since each restart contributes at most one), the largest
  fraction of any one sample's restarts carrying it (`max_recurrence`),
  and the mean within-cluster ownership correlation (`agreement`). A factor whose
  coverage is one sample's share is carried by that sample alone.
- **Recurrence** (`recurrence()`): per factor and sample, the number of
  that sample's restarts containing the program and the fraction of its
  restarts; zero marks a sample that does not carry the factor.
- **Unreproduced components** (`unreproduced`): the pooled components whose
  cluster did not reach `min_restarts`, with their sample, restart, factor
  index, and whether they belong to the sample's primary restart.
- **Matching** (`matching`): each sample's primary factors with the shared
  factor their cluster became (or none), the ownership correlation with
  that shared profile, their within-cluster agreement, and their restart
  stability within the sample (`own_stability`).
- **Stability** (`stability()`): each sample's within-sample restart
  stability next to the agreement of its aligned restarts with the
  consensus profiles.
- **Factor usage** (`sample_usage()`, `usage_divergence()`): the share of
  each cell type's molecules labeled with each factor, per sample, and
  per type half the L1 distance between a sample's usage and the pooled
  usage, maximized over samples (or per type and sample with
  `by_sample = TRUE`). A factor used by one sample only, or a type whose
  usage differs between samples, points to a sample effect in the
  factorization.

Two per-sample summaries gather these into one row per sample, with
flags that report rather than gate:

- **Consensus fit** (`diagnostics()`): the loss gap and restart spread
  with `representation_flag` when the gap exceeds the spread (the shared
  profiles reconstruct the sample's neighborhood compositions worse than
  one of its own restarts would); the number and fraction of the sample's
  primary factors that entered the shared space with `shared_factors_flag`
  below one half; within-sample stability next to consensus agreement; and
  the largest usage divergence of any cell type with that type's name. The
  representation flag describes reconstruction, not cleanup: a compact
  shared space of reproducible programs can reconstruct an unstable
  variant's neighborhoods poorly and still remove admixture well, so the
  per-sample audit of the correction remains the arbiter of cleanup.
- **Correction** (`diagnostics(report)`): the number of rules applied in
  the sample (`rules_applied`); among them, the number that lack evidence
  from the sample itself (the sample is absent from the rule's
  `samples_significant`) and the fraction of the sample's removed
  molecules attributable to those rules; and the number of rules whose
  program none of the sample's own restarts carries (zero recurrence in
  that sample), applied there under `apply_rules = "everywhere"` and
  skipped under `"where_present"`, with the fraction of the sample's
  removed molecules attributable to those. A removed molecule is
  attributed to the rule of its cell's type and primary factor label. A
  high fraction in the first pair means the sample's cleanup rests on
  other samples' evidence, which is the intended uniform treatment of
  replicates but a point to examine when the samples are not replicates; a
  high fraction in the second means the sample is being cleaned of a
  program it does not have, whose labels in that sample are nearest-factor
  artifacts, which happens when samples differ in composition (regional
  sections, cohorts) and is what `apply_rules = "where_present"` prevents.
  Given the cleanup report of the audit, the audit's leakage removed,
  median pair sensitivity, and own-marker false removal are added per
  sample, with `false_removal_flag` when a sample's false removal exceeds
  twice the median across samples.

A sample that trips the shared-factor flag has many primary factors that
no other restart reproduces, which happens when its own restarts disagree
among themselves (low within-sample stability); the shared space then
holds fewer programs than the sample's own fit, and admixture programs
among the lost ones are not cleaned. Lowering `min_recurrence` admits the
less reproducible programs and recovers sensitivity, at the price of a
larger shared space and more own-marker removal; on an unstable variant
neither threshold matches the balance of a stable one, and the better
remedy is a more stable NMF variant or more restarts for the per-sample
fits. On the breast replicates, `ls_nmf` with bridge scoring under the
default consensus cleaned better on every audit measure than
`invsqrt_kl` with membrane scoring under any threshold. Projecting every
sample onto one coherent decomposition through `fit_profiles()`, or a
joint factorization, are the remaining options.

## Usage

R:

```r
datasets <- list(
  rep1 = cellAdmix("data/rep1", output_dir = "out/rep1", annotation = annotation_rep1),
  rep2 = cellAdmix("data/rep2", output_dir = "out/rep2", annotation = annotation_rep2))
collection <- cellAdmixCollection(datasets)          # sample_allocation = "proportional"
consensus <- collection$fit_consensus()              # per-sample fits, matching, projection
consensus$loss                                       # loss gap per sample
consensus$factors                                    # coverage and recurrence per shared factor
consensus$recurrence()                               # per factor and sample
consensus$matching                                   # each sample's primary factors in the shared space
score <- consensus$score_bridge()                    # per-sample scores, pooled summary
rules <- score$rules()                               # pooled rules with sample_support
correction <- score$correct(rules = rules)           # per-sample ensemble corrections
score$correct(rules = rules, apply_rules = "where_present")  # skip rules on programs a sample lacks
correction$summary()                                 # removal per cell type per sample
correction$counts()                                  # combined; columns "<sample>_<cell_id>"
correction$counts(sample = "rep1")                   # one sample, original cell ids
audit <- consensus$audit_admixture()
report <- audit$evaluate(correction)
report$summary()                                     # per-sample cleanup verification
consensus$diagnostics()                              # per-sample representation flags
correction$diagnostics(report)                       # rules without local evidence, audit per sample
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
