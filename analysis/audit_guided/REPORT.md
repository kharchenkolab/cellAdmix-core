# Audit-guided correction: pancreas evaluation

Evaluation of two audit-guided correction strategies against the current
defaults, under a gene-split validation harness. Scripts `00`-`07` in this
directory; per-pair CSVs under `results/`.

## Harness

Each of the 39 detected pairs' 20-gene marker pools is split into a guide
half (A) and a validation half (B), alternating along the contrast ranking.
Guidance (anchor gene sets, exposure budgets, vote-threshold selection) sees
A only; performance is reported on B with B's own re-estimated exposure-0
baseline. Metrics per pair: `power_B` (fraction of exposure-linked B-excess
removed), `power_strictB` (same on the strict, near-zero-baseline subset),
`depth_B` (how much of the exposure-0 B baseline was removed - erosion when
the baseline is native, desirable ambient cleanup when it is not). Per type:
own-marker false removal on native markers disjoint from all guide pools.

## Scoreboard (excess-weighted over 39 pairs)

| arm | power_B | strict-B | depth_B (med) | false removal |
|---|---|---|---|---|
| membrane ensemble (default, vote 0.3) | 0.961 | 0.998 | 0.97 | 7.5% |
| membrane single fit | 0.670 | 0.869 | 0.81 | 5.5% |
| bridge ensemble (ls) | 0.503 | 0.527 | 0.43 | 1.7% |
| bridge ensemble (invsqrt) | 0.666 | 0.482 | 0.70 | 6.8% |
| vote auto-calibrated (A-side, 5% budget -> 0.7) | 0.839 | - | 0.77 | 4.7% |
| + mural anchor, membrane (single) | 0.804 | 0.929 | 0.70 | 5.7% |
| + mural anchor + audit-attributed rules, bridge (single) | 0.425 | 0.453 | 0.39 | 5.3% |
| full audit-guided type anchors, bridge (single) | 0.376 | 0.519 | 0.23 | 5.3% |
| historical SPA anchors s1-s3, bridge (single) | 0.41-0.53 | - | - | 2.1-9.6% |
| exposure-regression corrector (E4 prototype) | 0.710 | 0.607 | **0.00** | **0.0%** |
| exposure corrector v3, validation arm (A-guided) | 0.906 | - | 0.00 | 0.0% |
| exposure corrector v3, production arm (Phase A) | **0.992** | - | - | **0.0%** |
| composite: regression + induced-gene retention (17) | 0.988 | 1.000 | 0.43 | **0.0%** |

Controls: junk anchors align to no type and produce zero rules (both
scoring methods); shuffled exposure collapses the E4 corrector to
power 0.001; E4 leaves zero-exposure cells untouched by construction.

## Findings

1. **The ensemble covers far more than its rules admit.** Six of the eight
   pairs flagged "no removal rule covers this pair" are in fact corrected by
   member votes - Immune -> Fibroblast/CAF loses 100% of its marker content
   (31,624 of 31,765 B-half molecules) with zero rows for the pair in
   `correction$rules`. The uncovered-pair warnings and the rules table
   describe the primary member, not realized behavior. Product fix needed:
   report realized per-pair removal (the audit measures it) and per-member
   rule provenance.
2. **The ensemble's real cost is baseline erosion, not coverage.** Median
   `depth_B` 0.97: on 36/39 pairs the default correction removes most of the
   unexposed-baseline marker signal along with the leakage. For strict-tier
   genes that baseline is ambient and its removal is a feature; for
   shared-gene pairs (e.g. Immune -> Fibroblast/CAF, whose strict tier is
   empty) it erases genuine expression that the audit says makes up ~half
   the pool content.
3. **The audit's per-pair dose-response supports calibrated removal.** The
   minimal exposure-regression corrector - per-cell budgets from A-half
   excess rates, spent across source-owned genes by likelihood ratio with
   waterfilling - reaches 0.71 power at exactly zero baseline erosion and
   zero own-marker false removal, with clean negative controls. No other arm
   touches that corner of the trade-off space. Its remaining gap is pairs
   where leakage sits on genes the budget weighting under-serves.
4. **Targeted anchor recovery works and its mechanism is verified.** The
   audit identified Mural/pericyte as the one type without an aligned factor.
   A soft anchor built from its one-vs-rest profile margin aligns cleanly
   (factor-source score 0.73), forms scoring rules, passes the native check,
   and fixes the whole mural cluster: Fibroblast -> Mural 0.02 -> 0.81
   (the missing native factor had been shielding fibro content in mural
   cells), Mural -> Endothelial 0.85, Mural -> Immune 0.95,
   Mural -> Fibroblast 0.49 - single-fit results at or above the 10-member
   ensemble on three of four pairs. Junk anchors are rejected by the
   existing gates, so anchoring does not make the pipeline credulous.
5. **Full-replacement type anchors are not the answer.** A perfect
   one-factor-per-type factorization (all 7 types aligned) underperforms
   the historical SPA anchors and the extension strategy - state-level
   substructure matters for bridge contrast. Extend the fitted basis;
   don't replace it.
6. **Vote-threshold auto-calibration works end to end.** Selecting on
   A-halves under a 5% false-removal budget picks vote=0.7 and delivers
   0.839 B-power at 4.7% - a legitimate, per-dataset operating-point
   chooser for `correct()`.
7. **Two rule-derivation gaps surfaced.** The anchor factor's automatic
   source call was wrong (Immune instead of Mural) - audit-side attribution
   fixed it; and Exocrine -> Immune has significant membrane evidence under
   an ambiguous factor that source-calling never credits to exocrine.
   Both argue for letting pair-conditional evidence (which the scoring
   already computes) drive rule sources when the audit detects the pair.
8. **A real pipeline bug found and fixed:** fixed-H fits sized cell factor
   fractions by the annotation-derived rank instead of H's row count,
   silently dropping extra factors from cell-level fractions (and thereby
   from the native check). Fixed in `pipeline_store.cpp` with a regression
   test; 245+4 tests pass.

## Phase A: exposure corrector v3

Three measured shortfalls of the E4 prototype were fixed: an isotonic
dose-response on raw exposure counts (the coarse "3+" bin had carried
65-95% of each pair's excess), exact stratum-level budget delivery with
gene waterfilling (per-cell budgets had died on 0.6-2.5 pool molecules per
cell), and an ambient tier that removes strict-gene content outright
(near-zero native baseline means that content is contamination at any
exposure). Results (script `08_phase_a.R`):

- validation arm (budgets from A halves only, no ambient): power_B 0.906 -
  the dose-response generalizes from guide genes to held-out genes;
- production arm (full pools, ambient on): power_B 0.992 (median 0.998),
  own-marker false removal 0.0%, budget delivery 98-99%, top-up residual
  <1% after one pass;
- coverage gate: zero pairs below 0.8 - every previously-stuck pair is
  corrected;
- shuffled-exposure control: 0.007;
- downstream cell-state kNN purity: original 0.804, membrane ensemble
  0.968, exposure v3 0.967 - identical downstream benefit while removing
  1.16M molecules versus the ensemble's 1.99M and touching no own-marker
  content.

The regression tier meets every gate on pancreas; the full guided-NMF/EM
tier (E5) remains unjustified here and is deferred to the breast phase.

## Composite corrector: regression removal with induced-gene retention

The generative-model evaluation (`analysis/generative_model/FINDINGS.md`)
predicted that its two transferable ideas - the interface-local
overdispersed proportionality screen and per-gene retention of the
disproportionate (induced) excess share - could be grafted onto the
regression corrector to dominate both parents. Script `17_composite.R`
implements exactly that: the screen (now the package's audit screen,
applied to every source-owned gene per pair) yields flagged (pair, gene)
combinations with an induced share f = (excess - expected) / excess; the
corrector measures its dose on unflagged guide genes only, caps each
flagged gene's removable content at the proportional share 1 - f
(waterfilling the remainder to other source-owned genes), and the ambient
tier retains the f share of flagged strict genes. On pancreas, 26
pair-gene combinations are flagged across 14 pairs.

Production configuration (pool-guided, ambient on): weighted power_B
0.988 (strict tier 1.000); on the admixture-only yardstick (B halves
excluding flagged genes) 0.992 weighted with every pair at or above 0.836
- matching the plain regression corrector on transferred content. The
flagged genes retain 79.6% of their exposure-linked excess in aggregate
(CFTR 0.81, PROX1 0.74, SEMA3C 0.86), versus 8.4% when the retention
filter is ablated with the same machinery - the filter, not the
corrector, is what preserves the biology. The generative model's own
retention figure is 79.5%, so the count-level cap reproduces its
entry-level split almost exactly. Own-marker false removal 0.0% (pooled
and worst type), shuffled-exposure control 0.008, downstream kNN purity
0.964 (ensemble 0.968, regression 0.967, generative model 0.963). The
two pairs below 0.8 on the standard yardstick (Endothelial -> Ductal
0.67, Fibroblast -> Mural 0.73) are precisely the pairs whose B halves
carry retained flagged genes - deliberate retention, not missed removal.

The same script run on the CosMx NSCLC dataset (`COMPOSITE_DATASET=nsclc`;
88 flagged pair-gene combinations across 18 of 27 detected pairs,
dominated by stress and chemokine programs - FOS, JUNB, HSPA1A/B,
CCL3/4, CXCL2/3) reproduces every gate: production power_B 0.982
(median 1.000), admixture-only yardstick 0.998 with a per-pair minimum
of 0.966, retention 80.1%, shuffled control 0.013, own-marker false
removal 0.0%, purity 0.811 -> 0.885. The retention ablation is milder
there (73.9% kept even without the filter): NSCLC's flagged genes are
channels the corrector's delivery weighting touches lightly anyway,
where on pancreas the filter is what saves CFTR's 111k-molecule excess.

The composite's remaining gap to the generative model is generalization
of the dose: the A-guided validation arm scores 0.903 versus the model's
0.974, because a fixed stratum budget generalizes to held-out genes worse
than the model's per-cell posterior delivery. For product purposes the
composite offers the model's retention behavior at the regression
corrector's complexity and runtime.

## Verdicts against the pre-registered gates

- E1 (anchors): PASS - mean B-power over the mural cluster 0.77 (gate 0.5),
  no false-removal increase, junk control clean.
- E2 (vs SPA anchors): full replacement 0.38 vs floor 0.52 - FAIL for the
  replacement strategy, superseded by the extension strategy which was not
  part of the historical comparison.
- E3 (calibration): PASS - stable selection, honest B-side reporting.
- E4/Phase A (exposure corrector): PASS - v3 reaches 0.992 aggregate
  power with full pair coverage, clean controls, zero own-marker damage,
  and ensemble-parity downstream purity. The full exposure-informed EM
  (E5) is not justified by any remaining pancreas gap; revisit on breast.

## Molecule-level diagnostics: patchiness, out-of-plane leakage (T1/T2)

Scripts `09`-`11`. Two section-geometry concerns were tested directly on
the run's molecule table (continuous z, CRF factor labels, nucleus fields):

- **The fragment (CRF) hypothesis holds at molecule level.** Source-pool
  molecules inside target cells co-locate with molecules of *different*
  pool genes within 1um at 1.7x the within-cell permutation null (median
  over pairs; up to 2.9x for fibroblast sources). A radius of 2.5um
  saturates the statistic in dense cells - the nearest-different-gene
  formulation is the usable one.
- **The zero-exposure baseline is not ambient.** In target cells with zero
  observed source neighbors, strict-gene content decays 5-10x with lateral
  distance to the nearest source cell, patch enrichment is *higher* than
  in exposed cells (1.8x median), and strict molecules sit mildly toward
  the slab surfaces - the signature of out-of-plane/near-field leakage
  inside the reference level. Exocrine -> Ductal is the instructive
  exception: its unexposed baseline shows no patch enrichment (dispersed,
  genuinely ambient-like), so ambient and out-of-plane components separate
  per pair.
- **Repricing the reference is material.** Replacing the pooled
  zero-exposure reference with the >250um far-field rate roughly halves
  baselines (e.g. 24.2 -> 11.5 per 1k for exocrine -> ductal) and inflates
  excess estimates by a median 1.41x on the 16 measurable pairs (+12%
  aggregate; up to 1.9x per pair). Pairs whose targets never sit far from
  sources need a decay-curve extrapolation instead of a far bin.

Method consequences: the audit and the exposure corrector should use a
distance-decay-extrapolated reference (they currently under-remove by the
structured-baseline component); patch membership (CRF labels, already per
molecule) is the right molecule-selection prior and extends ambient-tier
removal to shared genes via patch evidence rather than gene identity; and
the non-saturating patch statistic enables per-gene induction flagging
(gradient without patch enrichment) ahead of the formal proportionality
test.

## Discriminator characterization (T3/T5/T7)

- **T3 (profile proportionality - the primary admixture/induction
  discriminator).** Median correlation of per-gene exposure-linked excess
  with the source profile: 0.82-0.84 across 38 pairs (cytoplasmic vs
  whole-cell vs nuclear profiles do not separate on a 377-gene panel);
  median 70% of each pair's excess is profile-proportional. 99 outlier
  genes flagged, led by CXCL6, CFB, PPP1R1B in ductal-source pairs -
  inflammation/activation genes with exposure-linked excess far above
  their source-profile share: the induction signature, present in real
  data. Patch composition corroborates: source-factor-labeled molecules in
  target cells match the source cytoplasmic profile (cos 0.63) and not the
  target profile (0.12).
- **T5 (control-probe null).** Even pure noise is structured: control
  molecules are clustered (0.37x dispersion vs random), surface-biased in
  z (|dev| 0.27 vs 0.18 null), and nucleus-enriched - so the strict-gene
  baseline's z behavior matches noise/surface material, and small-scale
  clumping alone cannot certify admixture.
- **T7 (four-class spike-in ground truth).** Fragments (in-plane and
  out-of-plane), induction, and ambient singletons injected into real
  exposed cells: the T3 test recovers exactly the three planted induction
  genes as top outliers (z 21-76); the patch statistic separates
  fragments/induction (0.34-0.44) from ambient (0.05) but NOT induction
  from fragments - co-induced genes co-locate by cell density alone.
  Patch absence is a corroborating clue only; profile proportionality is
  the primary distinction.

## Recommended next steps

1. Breast 5K phase 2: anchor endothelial + fibroblast (the two missing
   signatures there), same harness - the case anchoring was designed for.
2. Composite correction: ensemble removal for strict-tier(ambient-baseline)
   content, exposure-calibrated budgets for shared-gene pairs - the two
   halves are complementary corners of the measured trade-off surface.
3. Productize: realized-removal reporting in `evaluate()`/`correction$rules`,
   `fit$recover_factors(audit)`, and `correct(vote = "auto", ...)`.
