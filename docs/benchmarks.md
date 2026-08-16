# Benchmarking admixture cleanup

Molecule-level admixture correction faces a basic evaluation problem: no one
can say, molecule by molecule, which transcripts truly leaked between
neighboring cells. This report describes a benchmark that lets the tissue
itself supply the answer at the population level, applies it to cellAdmix's
factorization variants and scoring methods across three datasets, quantifies
the run-to-run stochasticity of the correction pipeline, and evaluates
ensemble corrections that turn that stochasticity from a liability into a
calibration dial. The benchmark harness lives in `analysis/cleanup_benchmark/`.

## The neighbor benchmark

Leakage is a neighborhood phenomenon: a cell can only pick up stray molecules
from cells physically adjacent to it. For an ordered pair of cell types —
a *source* S and a *target* T — every T cell is characterized by its
*exposure*: the number of S cells among its 15 nearest cells. T cells with
zero exposure form a clean reference; whatever they contain, a T cell contains
on its own.

For each pair we select a panel of *source markers*: genes whose top
expresser is S, ranked by contrast against the zero-exposure T baseline
(rank-based selection, since an absolute baseline cutoff would itself be
skewed by contamination). Within the panel, the *strict tier* holds genes
essentially absent from reference T cells (baseline under 5% of the source
level); their excess in exposed T cells can only be leaked material. Before
correction, the pooled marker rate climbs steadily with exposure
(Figure 1, red); the climb, converted to molecule counts above the
zero-exposure baseline, is the pair's estimated leakage. Pairs enter the
benchmark when this excess is statistically solid (Poisson test,
Benjamini-Hochberg across pairs) and exceeds 200 molecules — 39 of 42
candidate pairs on the pancreas dataset, 13 on the breast crop, 27 on NSCLC,
with no manual curation.

A correction is then scored by how much of the excess it eliminated
(Figure 1, blue): a value of 1 means exposed T cells now look exactly like
unexposed ones. Because the "after" rate is compared to the corrected
zero-exposure baseline, indiscriminate deletion cannot masquerade as cleanup
of the exposure gradient. The bottom panel of Figure 1 shows the
resulting per-pair profile for a standard single-fit cleanup (ls-NMF,
membrane scoring): effectiveness is highly heterogeneous across pairs, and
the largest pair by leakage mass (exocrine → ductal, over 200,000 molecules)
is missed entirely — a coverage failure invisible to aggregate diagnostics.

![Figure 1](figures/benchmark_fig1.png)

**Figure 1. The neighbor benchmark.** *Top:* pooled strict-tier
source-marker rates in target cells, stratified by the number of source-type
neighbors, before (red, dashed) and after (blue) a standard cleanup (bare
ls-NMF fit, membrane scoring, pancreas dataset); the dotted line marks the
zero-exposure baseline. The rise with exposure is contamination made
visible; cleanup quality is the degree to which the blue curve flattens to
baseline — compare the near-complete flattening of endocrine → endothelial
with the untouched fibroblast → immune pair, where the curves coincide.
*Bottom:* the same score for every detected pair (bars), with each pair's
estimated leaked-molecule count overlaid (orange, log scale). Take-home:
cleanup effectiveness is measurable without molecule-level ground truth, and
a standard single-fit correction is highly uneven across cell-type pairs —
including a near-zero score on one of the largest leakage pairs
(exocrine → ductal).

## From excess removal to sensitivity and specificity

The benchmark's two axes are estimators of the standard classification
quantities, each computed on a stratum where the truth is nearly known:

- **Estimated sensitivity.** For strict-tier markers, the excess above the
  zero-exposure baseline estimates the number of truly leaked molecules, so
  the fraction of excess removed is an estimate of TP/(TP+FN) on that
  stratum. Its bias: it samples leakage carried by source-specific genes —
  the *easiest* leakage to attribute — so it is an optimistic stratum for any
  marker-driven method.
- **Estimated false-positive rate.** Molecules of a cell type's own marker
  genes, inside cells of that type, are almost surely genuine; their removal
  rate is a direct FPR estimate (specificity = 1 − FPR). Its bias runs the
  other way: it samples the molecules easiest to keep. A worst-case
  companion — the most-affected pair's retention — tracks the shared-gene
  erasures (e.g. CFTR, genuinely expressed by both exocrine and ductal cells)
  that pooled numbers hide.

The two strata bracket the unmeasurable middle (shared and non-distinctive
genes), and all headline results below are reported in these terms.

## Stochasticity of factorization and scoring

Rerunning the identical pipeline with different random seeds changes the
correction substantially (Figure 2). On the pancreas dataset, ten seeds of
the default membrane-scoring pipeline remove between 1.07 and 1.54 million
molecules (invsqrt KL-NMF; 0.42-0.80 million for ls-NMF), and the removal
*sets* of any two seeds share only ~57% of molecules (median Jaccard;
~45% for ls-NMF). At the pair level, sensitivity can swing from 0.9 to
near 0 across seeds (Figure 2, right).

The source of this variability is structural, not numerical. KL-type NMF is
a mixture (topic) model; with the rank set above the number of well-separated
expression programs, the surplus factors face many near-equivalent choices —
duplicate a large program, split one, or form diffuse blends — and
multiplicative updates lock each restart into a different discrete gene
partition (81% of invsqrt loadings are zero; matched factors across seeds
share only ~43% of their owned genes). Restart objectives span ~2% while the
factor structures differ qualitatively, and the objective does not predict
cleanup quality — so best-objective selection among restarts cannot resolve
the ambiguity, and neither can more restarts. Scoring inherits this
variability twice over: a pair's removal *decision* may fall below threshold
when its evidence is split across factors, and the removal *labels* may be
absent when no factor owns the leaked molecules. Dense ls-NMF factors are far
more reproducible (matched ownership correlation 0.90 vs 0.43) but commit
weakly to molecule labels, trading instability for insensitivity — the two
variants fail at different stages rather than one being uniformly better.

![Figure 2](figures/benchmark_fig2.png)

**Figure 2. Correction stochasticity across random seeds (pancreas,
membrane scoring).** *Left:* total molecules removed by ten reruns of the
identical pipeline. *Middle:* pairwise overlap (Jaccard) of the removal
sets, invsqrt KL-NMF. *Right:* per-pair estimated sensitivity across seeds.
Take-home: the single-fit correction is effectively a lottery — comparable
total removal, but only about half the individual molecules agree between
runs, and individual pairs flip between fully cleaned and untouched.

## Ensemble corrections by molecule voting

Since each seed's correction is a per-molecule decision, an ensemble is
natural: run the pipeline N times, count for each molecule how many runs
removed it, and remove molecules whose vote count reaches a threshold. The
threshold is a sensitivity/specificity dial, and sweeping it traces an
estimated ROC curve (Figure 3).

Two properties of the resulting curves are noteworthy. First, the operating
knee sits at a *minority* vote — around 3 of 10 — not at majority: because
factor ownership is a lottery, the factor structure that correctly cleans a
given pair arises in only a minority of restarts, so demanding majority
agreement discards genuine cleanup (sensitivity halves between thresholds
3 and 5 on pancreas membrane scoring). Second, the vote count separates
systematic removals from idiosyncratic ones: worst-case pair specificity
rises from ~0.2 at the union threshold to ~1.0 at unanimity, meaning the
shared-gene erasures are low-vote events contributed by one or two seeds,
while genuine leakage removal accumulates votes up to its lottery ceiling.

A lighter-weight alternative — pooling only the removal *decisions* across
seeds and applying them with each fit's own labels, with imported decisions
re-vetted by the importing fit's native-factor check — repairs
decision-stage variance (it rescued every catastrophic seed on pancreas
membrane scoring: 0.40 → 0.85) but cannot supply labels a fit never
produced. The molecule-level vote subsumes it whenever at least the
threshold number of runs label the molecules correctly.

![Figure 3](figures/benchmark_fig3.png)

**Figure 3. Vote-threshold ROC for ensemble corrections (10 seeds).**
Estimated sensitivity (strict tier) versus estimated native-stratum FPR as
the required vote count varies from 1 (union) to 10 (unanimity); selected
thresholds annotated. Take-home: a vote fraction near 30% retains nearly all
of the union's sensitivity at meaningfully lower false-positive cost, while
majority and stricter thresholds are miscalibrated because correct removals
are typically minority events across restarts.

## Benchmarks across datasets and scoring methods

FIGURE4_AND_TABLE_PLACEHOLDER

## Conclusions

The neighbor benchmark makes admixture cleanup measurable on real data, at
the resolution of cell-type pairs, in estimated sensitivity/specificity
terms, and with no manual curation — which in turn makes method comparison
and tuning an empirical exercise rather than a judgment call. Applied to
cellAdmix it yields three conclusions:

1. **Single-fit corrections are lotteries.** Identical settings produce
   removal sets sharing only about half their molecules across seeds, with
   pair-level effectiveness flipping between complete and absent. Neither
   factorization variant escapes this: invsqrt KL-NMF gambles on factor
   ownership, ls-NMF on scoring decisions.
2. **Ensembling the correction fixes what selection cannot.** The restart
   objective does not identify good corrections, but the restarts jointly
   contain them: molecule-level voting across seeds dominates every
   single-fit arm, and its threshold is a transparent, per-dataset-calibrable
   sensitivity/specificity dial.
3. **Recommended solution:** derive corrections by molecule-level voting
   across the fit's restarts (which the pipeline already computes for its
   stability diagnostic), with the vote threshold defaulting to ~30% of
   restarts and exposed as a user-facing dial; retain the per-fit
   native-factor check within each contributing run. RECOMMEND_DETAIL_PLACEHOLDER
