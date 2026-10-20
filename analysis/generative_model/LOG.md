# OODA loop log: generative admixture model

Running log of observe/orient/decide/act cycles. Summarized in REPORT.md at
the end; kept verbatim here.

## Loop 0 - observe (setup)

Read: audit REPORT.md, 01_split_eval.R (harness), 08_phase_a.R (regression
corrector v3 - the bar: 0.906 validation / 0.992 production weighted removal
power, zero own-marker removal, shuffled control 0.007, kNN purity 0.967),
14_t7_spikeins.py (four-class spike-in), metrics.R (bench_* functions),
exported tables under audit_guided/results/.

Key observations that shape the design:

- The count matrix is 377 genes x 140,335 cells, 6.19M molecules; per-type
  cell counts 2.2k (Endocrine) to 36k (Exocrine). Dense per-type EM is
  feasible in numpy.
- Factor->type alignment (factor_types.csv): 1 Fibroblast/CAF, 2 Exocrine,
  3 Immune, 4 Ductal/tumor, 6 Endothelial, 7 Endocrine; factors 5, 8, 9
  unaligned. Factor 5 is a second ductal/exocrine state (22% of ductal
  content). Mural/pericyte has no own factor (factor 1 dominates it).
- Real induced-gene candidates (T3 test): CXCL6, CFB, PPP1R1B in
  ductal-source pairs, residual z 6-45. These sit INSIDE the evaluation
  marker pools, so retaining them costs measured removal power - an honest
  tension between gate 2 (>=0.99 power) and gate 3 (retain induced).
- The regression corrector's shuffled control was run with the ambient tier
  OFF (08_phase_a.R line 182). Any exposure-independent fractional removal
  of pool genes shows up as "power" under shuffling (excess_a = (1-f) *
  excess_b), so the shuffled control must disable exposure-independent
  components; I mirror that convention.
- eval_split re-estimates the exposure-0 baseline after correction, so
  removing the flat baseline component from all cells (exposed and not)
  still reads as power only via the gradient; fractional gene-level removal
  applied equally to all cells reads as power f. This is why the production
  regression arm turns the ambient (strict-gene) tier on.

Decisions (Loop 0):

- Implement the model in Python (numpy) on a count matrix exported from R
  via the exact build_split_defs() harness, so "before" counts match the
  evaluation bit-for-bit; corrections go back to R as per-entry retained
  fractions and are evaluated with the audit's own eval_split().
- Per target type T: Poisson mixture with (a) own expression = mixture of
  programs initialized from factor-labeled molecule profiles of T's own
  factors (pseudobulk fallback for Mural), (b) one contamination component
  per detected source pair with per-cell fraction alpha_cS given a Gamma
  prior whose mean is the isotonic exposure dose-response measured on guide
  genes (A halves in validation, full pools in production) scaled to the
  whole source cytoplasmic profile, (c) an ambient component with a global
  profile zeroed on the target's own genes, (d) a sparse per-pair induced
  term: support screened by a per-gene proportionality test (excess vs
  source-profile share, z threshold), per-gene rate refined by EM, with a
  per-cell activity multiplier so cell-specific induction is retained,
  (e) a tiny uniform floor so removal fractions are well-behaved when own
  expression is zero.
- Validation arm: alpha/dose evidence restricted to non-B genes (scaled by
  the source-profile share of allowed genes); no ambient; no strict-gene
  constraints; induced support restricted to A-side. Production arm: full
  pools, ambient on, own-program support zeroed on strict genes (their
  native baseline is near zero by definition), induced term on.
- Source cytoplasmic profiles from overlaps_nucleus == 0 molecules of
  source cells; refined over 2-3 outer rounds by multiplying with the
  fitted own-fraction per gene of the source type's own model
  (decontamination as required by the spec).

## Loop 1 - act (first end-to-end run)

Exported harness definitions: 39 pairs, all pools disjoint from target
native markers, 6,186,500 molecules in the count matrix (matches the run's
molecule table exactly - no quality-value filtering in the store).

Design decisions locked before seeing gate results:

- alpha update is a Gamma-posterior step (prior strength 300
  pseudo-molecules) with the posterior clipped to at most 3x the prior mean,
  and exactly 0 where the prior mean is 0. Rationale: the shuffled-exposure
  control (gate 4) requires that removal be anchored to the measured
  exposure signal; an unbounded evidence-driven alpha would keep finding the
  true contamination under shuffled exposure (it is genuinely there), which
  would show up as power against the original exposure bins and fail the
  pre-registered control. The bounded update still lets each cell's
  expression move its own contamination estimate up to 3-fold.
- The ambient component is restricted to the strict-gene union (genes with
  near-zero native baseline in the target), profile and scale taken from
  target cells >250 um from any source cell; this mirrors the regression
  corrector's ambient tier (production only, off in validation and in the
  shuffled control, following 08_phase_a.R's own convention).
- Induced-expression support is selected by a per-gene proportionality
  z test (z > 6 and > 50 excess molecules) - this is the sparsity mechanism
  (a hard analog of an L1 threshold); the per-gene rates and a per-cell
  activity multiplier are then refined by EM and the content is retained.
- Validation arm honesty: the alpha evidence sums exclude all held-out
  B-half genes (rescaled by the source-profile share of the remaining
  genes), the dose-response uses A halves only, and the induced screen
  cannot select B genes. Own-program profiles are free on all genes, which
  makes the test conservative: if the A-derived dose underestimates
  contamination, the own component absorbs B-gene excess and power_B drops.

## Loop 2 - observe (first gate numbers, three arms)

First full evaluation (before the validation arm finished; it was killed and
relaunched after the loop-3 changes below):

- shuffled-exposure control: weighted power 0.016, removal 22,487 molecules
  (0.36% of production) - the control collapses as required (gate 4 PASS at
  this stage).
- production: weighted power_B 0.877 (median 0.953), 10 pairs < 0.8,
  own-marker false removal 1.4% pooled / 4.4% worst - below both gate-2
  thresholds.
- production without the induced term: 0.955 (median 1.000), 2 pairs < 0.8;
  the induced term as first configured cost 0.078 weighted power.

Diagnosis from the diagnostics tables:

1. The induced screen fired on 153 genes. A pure z threshold is scale-
   dependent: AMY2A passed with excess 148,714 vs a proportional expectation
   of 138,248 (1.08-fold) simply because z grows with sqrt(N). Retained
   "induced" content reached 33.5k molecules on Ductal -> Exocrine alone -
   more than the pair's removed contamination - dragging measured power.
   Many flagged genes (SEMA3C, TNC, ALDH1A3, VCAN...) look like activated-
   fibroblast markers: source cells NEAR tumor have a shifted composition
   relative to the global source profile, which a global-profile
   proportionality test misreads as target-side induction.
2. Own-marker false removal (1.4-4.4%) comes entirely from the source
   cytoplasmic profiles carrying weight on target-owned genes (source cells
   are themselves contaminated by the target type).
3. Ductal -> Exocrine sits at 0.34-0.39 in BOTH arms with power_strictB =
   1.0: the shortfall is on shared genes, where the exocrine own component -
   whose second program was initialized from factor-5-labeled molecules, a
   state shared with ductal - absorbs the exposure-linked ductal content.
4. Exocrine -> Ductal runs at the prior dose cap (mean lambda 0.37 vs cap
   0.4).
5. The far-field ambient estimate fell back to the half-baseline
   approximation for every target: no pancreas target type has enough cells
   >250 um from all of its sources. Expected in dense tissue; the fallback
   (far-field is 0.5-0.7x the pooled zero-exposure baseline per the audit)
   is the honest substitute.

## Loop 3 - decide/act (five changes, all pre-specified before re-running)

a. Source profiles are zeroed on target-owned genes (no renormalization):
   contamination there is unidentifiable from own expression and is left in
   place; own-marker false removal becomes structurally zero, matching the
   regression corrector's gset convention.
b. Contamination profiles are now PER PAIR: the cytoplasmic profile of the
   source cells within 30 um of the nearest target-type cell (fallback 100
   um, then all cells, if under 30k cytoplasmic molecules) - the actual
   donor population - so source-composition gradients (CAF activation near
   tumor) live in the contamination profile, not in the induced term.
c. Induced screen gains a disproportionality requirement: excess > 3x the
   proportional expectation, in addition to z > 6 and > 50 molecules.
d. Own-expression programs are learned predominantly from lightly dosed
   cells (cell weight 1/(1 + 20 * prior contamination fraction) in the
   program M-step), so admixed content in heavily exposed cells cannot be
   re-absorbed as a "cell state". This targets the Ductal -> Exocrine
   factor-5 shield.
e. Prior dose cap raised from 0.4 to 0.6 (Exocrine -> Ductal was clipped).

## Loop 4 - observe/act (dose top-up)

Loop-3 results: validation 0.965 (gate-1 bar of 0.90 passed; regression
corrector: 0.906), production 0.975 (median 1.000, one pair < 0.8),
production without induced term 0.977, shuffled 0.017, own-marker false
removal exactly 0 everywhere. Per-pair comparison against the regression
corrector's validation arm: 25 pairs better, 2 slightly worse, 12 ties.
The remaining production shortfall concentrates in shared-gene pairs
(Fibroblast -> Ductal/Immune/Endothelial, Ductal -> Exocrine). Added the
analog of the regression corrector's top-up: after EM convergence,
re-measure the residual isotonic dose-response on the corrected counts
(excluding induced-support genes, whose gradient is retained by design),
add it to the prior dose, continue EM; raised the per-cell adaptation cap
to 5x. Result: production 0.9826, validation 0.9755. Freezing the own
programs during top-up (loop 5) changed nothing (0.9819) - the model's
remaining disagreement with the harness is not program drift.

## Loop 5 - orient (the stuck pair is the induced term working as designed)

Endothelial -> Ductal (0.73-0.75) decomposed gene by gene: every B-half
gene at removal 0.92-1.00 except GNG11 (0.25), which the screen flags at
5.4-fold disproportionality and deliberately retains; without the induced
term the pair scores 0.84 and no pair is below 0.8. The pure-removal
machinery's honest ceiling on this harness is ~0.985 ("admixture-only"
yardstick, i.e. excluding flagged genes: 0.9854, zero pairs < 0.8): the
last ~1.5% is exposure-linked excess on shared genes that the model
attributes to the target's own expression.

## Loop 6 - observe/decide (spike-in exposes a count-level limit; the
##                          molecule level carries gate 3)

First spike-in run (count level): fragments removed 0.83/0.85, ambient
0.85 - but planted induction retention 0.000. The planted genes are AQP8,
GATM, AMY2A - the largest acinar channels: natural proportional
contamination on them in exposed endothelial cells is 3.0k/12.4k/26.4k
molecules, the planted total is 3.75k. A 3-fold disproportionality guard
can never fire there; and with any profile-uncertainty allowance (CV >=
0.10) the overdispersed z is ~1. Only a Poisson-exact test "sees" the
spikes (z=21) - the same assumption that flooded loop 2 with 153 false
induced genes. Conclusion: population count-level proportionality cannot
detect induction riding on a dominant contamination channel; the spatial
signature (planted induction is not in gene-diverse source-matching
patches) is the only discriminator. Actions:
- replaced the ratio guard with the overdispersed z screen (profile CV
  0.15), threshold calibrated to z > 8: keeps the strong deviants (CFTR,
  SPIB 22x, FXYD2 20x, TM4SF18 16x, TNC 10x, SEMA3C 7.7x, GNG11 5.4x) and
  drops a fold-~2 band that includes obvious composition artifacts (GCG
  "induced" in fibroblasts = islet-periphery alpha-cell enrichment).
  CXCL6 (z 6.5, fold 2.1 against the NEAR-source profile) now falls below
  threshold: most of its exposure-linked excess matches what bordering
  ductal cells express - admixture from locally activated sources rather
  than target induction. The audit's global-profile test could not make
  this distinction.
- implemented the molecule-level realization (the stretch goal, now
  necessary): per-entry spike-and-slab posterior - removable molecules
  are Poisson with the fitted contamination+ambient mean, retained
  molecules are the fitted own+induced mean plus a small-probability
  (2%) geometric burst of unmodeled induction - combined with a
  per-molecule spatial-coherence feature (number of distinct other
  source-owned genes within 1.5 um), whose class-conditional Poisson
  means are estimated from reference sets (real strict-gene molecules in
  exposed cells vs native-marker molecules). Exact subset posterior via
  elementary symmetric polynomials gives each molecule a removal
  probability: gene identity sets the prior, spatial coherence moves it.

## Loop 7 - final calibration and the spike-in boundary

- Screen recalibrated to overdispersed z > 8 (profile CV 0.15): 28
  pair-gene rows / 25 unique genes, all with 2.3-22x disproportionality.
  Final harness numbers: validation 0.9736, production 0.9820 (the single
  sub-0.8 pair is Endothelial -> Ductal at 0.735, entirely the GNG11
  retention), production without induced term 0.9850 with all pairs >=
  0.8, admixture-only yardstick 0.9855 (min pair 0.915), shuffled 0.0180,
  own-marker removal 0.000. Downstream kNN purity 0.963 (original 0.804
  reproduced exactly; ensemble 0.968; regression 0.967).
- Real-data induced retention: the flagged genes keep 79.5% of their
  exposure-linked excess (CFTR: 89k of 111k molecules) while non-flagged
  pool content is removed at 96-99%.
- Molecule-level scorer upgraded to two spatial features at 1.0 um
  (distinct other source-owned genes; native-marker neighbors) with
  class-conditional Poisson likelihoods estimated from reference sets.
  Measured contrasts are weak (0.29 vs 0.18; 0.17 vs 0.50): single-
  molecule geometry at this density cannot overturn a strict-gene
  identity prior - consistent with the audit's T7 conclusion that patch
  statistics corroborate but cannot separate induction from fragments.
- Spike-in outcome, primary design (induction planted on the pair's top-3
  pool genes AQP8/GATM/AMY2A): fragments removed 0.826/0.846 (count) and
  0.809/0.830 (molecule), ambient removed 0.851/0.782 - but induction
  retention 0.000/0.102: FAIL, with an impossibility argument (the
  planted excess is a ~14% bump on the largest contamination channels;
  any test robust to >=10% profile error reads it as noise).
- Spike-in, supplementary variant (identical procedure, induction planted
  on mid-rank pool genes CA4/DIRAS3/PROX1): induction retention 0.949
  (count) / 0.942 (molecule), fragments removed 0.817/0.836 and
  0.803/0.821, ambient 0.749/0.694 - the retention machinery passes all
  gate-3 targets in the regime where induction is statistically
  detectable at all.

## Loop 8 - the dose prior: joint estimation and the distance ladder

- Observe: two recorded weaknesses of the dose prior - marginal per-pair
  estimation (joint/marginal ratio median 0.88, mural -> ductal 0.58 in
  07_robustness) and a prior of exactly zero in unexposed cells, where
  the audit's ladder shows contamination persists.
- Orient: both are properties of the covariate and of the fitting, not
  of the model's components. An ordered stratum axis (exposure count for
  exposed cells; ladder depth of the nearest source cell for unexposed
  ones, deepest first) makes one isotonic fit yield gradient and decay,
  with its reference level set by the audit's rule. Backfitting of monotone
  components across all sources of a target gives each pair its own
  component only.
- Decide/act: implemented in C++ behind `dose_mode`, `dose_ladder`,
  `dose_min_z`, with a built-in covariate-permutation control
  (`permute_covariates_seed`); both bindings pass the options through
  and the pair summaries report `prior_unexposed_molecules` and
  `cross_attributed_molecules`. Default options reproduce the previous
  behavior except for the guard (below).
- Observe (tests): the confounded fixture failed first because the
  "pure A behind a B stripe" geometry left A beyond the 14-neighbor
  window (A exposure zero everywhere) - redesigned as a randomly mixed
  A/C zone so the exposures co-vary across columns yet vary independently
  within one. The joint+ladder permutation control removed 17% of the
  planted transfer: with ~12 strata and a bimodal cell population, the
  isotonic fit of a permuted covariate drifts upward (prior at a tenth
  of the real fit's) and the top-up passes ratchet it. Orient: a
  maximum-type bias, not a bug; the audit guards pair detection with a
  significance test, the model had none. Act: `dose_min_z` - the pooled
  above-reference versus reference comparison on cell-level rates with an
  empirical standard error, applied to the partial residual under the
  joint fit and in every top-up pass. The control then collapses to
  exactly zero; real-data default removal changes by 0.17%, and the exact
  parity with the Python reference is recovered with the guard off.
- Observe (synthetic truth, first design): uninformative - all planted
  genes were strict, and the ambient tier's bounded per-cell update
  already removed 98% of the hidden-zone content in every arm; and the
  confounded A -> B pair was never detected because its only A-marker
  content in B was the shared gene, which the audit's screen flagged as
  disproportionate. Act: give B a low native level of two A genes (pool
  but not strict - the ambient tier cannot reach them), plant a small
  genuine A transfer so the pair exists, and add the variant in which C
  shares a proportional slice of A's profile (nothing to flag).
- Observe (final synthetic): shared-gene hidden-zone contamination 0.4%
  removed without the ladder, 53% with it, 0.1% of native content
  removed in reference cells; the joint fit turns a 1.97x over-credit of
  the small pair (shared profile) into 0.75x, and 0.93x -> 0.88x where
  only one gene is shared. Collinearity shrinkage is the joint fit's
  cost; it never inflates.
- Next: the real-data 2x2 on pancreas through the gate harness and the
  audit (11_dose_modes.py, 02_eval_gates.R, 12_dose_purity.R).
- Observe (pancreas, three yardsticks): the ladder removes all of the
  audit-measured unexposed-bin excess (55.7k pool molecules; 1.5% left
  by the marginal fit), raises leakage removal 0.942 -> 0.983, retention
  0.795 -> 0.806, purity 0.963 -> 0.967, own-marker removal 0; the
  harness's held-out power falls 0.980 -> 0.973. Orient: the harness
  references the post-correction zero-neighbor bin; cleaning that bin
  exposes residuals it hid. Checked per pair: of six pairs losing > 0.02,
  five show zero residual against the ambient reference; fibroblast ->
  mural keeps 7-11%. The joint fit is neutral on every aggregate
  (removal -0.2%, power -0.3/+0.1 points, purity -0.0003) and moves the
  attribution of the robustness check's pairs (mural -> ductal 0.58);
  its cost concentrates on fibroblast -> mural (0.91 -> 0.84).
- Decide: both on by default, pending the NSCLC and breast replication
  (same arms through the harness; audit scoring where the matrix fits).
- Observe (NSCLC, same arms, before the reference cap): the ladder lifts the
  audit's leakage removal 0.785 -> 0.970 and the unexposed-bin excess
  removed 0.55 -> 1.00 (removal +24%), purity 0.893 -> 0.902, own-marker
  0; the harness's held-out power falls 0.964 -> 0.919. Per pair, most of
  the drop is the reference artifact again (endothelial -> macrophage,
  endothelial -> malignant, immune other -> malignant: harness drops of
  0.15-0.18 with zero residual against the ambient reference), but
  endothelial -> immune other is a genuine regression: harness 0.999 ->
  0.405, audit residual 0 -> 12% of the exposed excess.
- Orient: endothelium is everywhere, so "no endothelial cell among the
  240 nearest" selects a biologically distinct population (immune
  aggregates without vessels) whose marker rate exceeds the pooled
  unexposed rate; the monotone fit then pools it with the next rungs,
  the reference level rises, and the exposed strata's prior shrinks. The audit
  guards its reference against exactly this by never letting it exceed
  the pooled unexposed rate; the model's reference level lacked the cap.
- Act: the reference level is now the smaller of the fit at the reference stratum
  and the pooled rate of all unexposed strata, in both the marginal and
  the joint path (the joint fit takes the cap as the weighted mean of its
  component over the unexposed levels). Synthetic and unit results
  unchanged (the cap never binds there); all arms rerun on the three
  datasets, with every arm's options set explicitly - the breast arms of
  the first pass had silently run one configuration after the package
  defaults changed mid-run.
- Observe (rerun with the cap): NSCLC unchanged - the four pairs that lose
  the most held-out power under the ladder (endothelial -> immune other,
  macrophage, malignant; immune other -> malignant) have a prior of
  exactly zero in the ladder arms. Per-stratum rates for endothelial ->
  immune other: deep stratum 7.3 per 1,000, the near-unexposed rungs
  5.7-6.2, exposed cells 7.3 (one neighbor) to 13.7 (six). The cap
  bounds the reference level, but the trend guard compared the rows above the
  reference stratum (near-unexposed rungs and exposed cells together)
  against the deep stratum alone - whose distinct population sits at the
  exposed level - and rejected the whole dose-response.
- Act: the guard now compares the exposed cells against all unexposed
  cells pooled (the split is the last unexposed stratum), which is how
  the audit detects a pair; without a ladder the split is the unexposed
  stratum, as before. The reference level keeps the cap. Rerun from the start.
- Observe (rerun with the cap and the exposed-vs-unexposed guard): the
  NSCLC regression is gone (endothelial -> immune other: harness 1.00 ->
  0.97 under the default, prior in unexposed cells 8.6k, zero residual
  against the ambient reference), and on all three datasets the ladder
  removes the whole audit-measured unexposed-bin excess at zero
  own-marker removal, with purity up (pancreas 0.962 -> 0.966, NSCLC
  0.893 -> 0.902, breast 0.810 -> 0.822). But the covariate-permutation
  control is not clean: it keeps 1.6% of the production removal on
  pancreas (three pairs), 6% on NSCLC (two pairs, fibroblast -> immune
  other alone 112k molecules), 0 on breast.
- Orient: the trend guard's standard error was the per-cell variance
  divided by the number of cells, which is the error of an unweighted
  mean; the guard compares molecule-weighted mean rates, whose error is
  the sum of squared weights times squared residuals over the squared
  total weight. With cell totals of coefficient of variation 0.6-0.7
  (effective sample size 0.67-0.76 of n) the z statistic of a permuted
  covariate had standard deviation 1.1-1.6 instead of 1.0 (300 random
  permutations per pair, 8 pairs on NSCLC and pancreas), so a threshold
  of 3 let through 0.7-2% of null tests per pair and pass - enough, over
  27-39 pairs and the top-up passes, for the survivors seen. With the
  weighted-mean error the permuted z has standard deviation 0.95-1.1;
  the real pairs' z falls (e.g. 92 -> 64) but stays far above 3.
- Act: the guard uses the weighted-mean standard error. Unit, Python and
  R suites rerun; all three datasets and the example notebooks rerun from
  the start.
- Observe (rerun with the weighted-mean error): production arms unchanged
  to three decimals on all three datasets; the control keeps 1.1% of the
  default's removal on pancreas (two endothelial-target pairs; harness
  0.025, audit 0.011) and exactly nothing on NSCLC and breast. NSCLC:
  the ladder lifts the audit's leakage removal 0.785 -> 0.976 under the
  default, purity 0.893 -> 0.902, retention 0.887 -> 0.923, with the
  harness at 0.950 against 0.976 and eight of the nine pairs losing
  more than 0.02 showing zero residual against the ambient reference;
  endothelial -> immune other is ordinary (187k against 196k marginal).
  Breast: no audit scoring; the ladder raises purity 0.810 -> 0.828 and
  removal 25%; the joint fit costs most on this dataset (removal -4%,
  harness -0.5, two small fibroblast pairs dropped where co-located
  T/NK and myeloid components explain the fibroblast content); the
  default stands at 0.927 against 0.942 on the harness and 0.822 against
  0.810 on purity.
- Decide: joint + ladder + guard confirmed as the defaults; the costs are
  recorded per pair in REPORT.md. Open question for the user: whether the
  joint fit, whose benefit shows only on the synthetic truth, should stay
  the default or ship as an option.
