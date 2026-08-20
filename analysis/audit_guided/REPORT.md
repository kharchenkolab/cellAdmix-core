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
| exposure-regression corrector (E4) | 0.710 | 0.607 | **0.00** | **0.0%** |

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

## Verdicts against the pre-registered gates

- E1 (anchors): PASS - mean B-power over the mural cluster 0.77 (gate 0.5),
  no false-removal increase, junk control clean.
- E2 (vs SPA anchors): full replacement 0.38 vs floor 0.52 - FAIL for the
  replacement strategy, superseded by the extension strategy which was not
  part of the historical comparison.
- E3 (calibration): PASS - stable selection, honest B-side reporting.
- E4 (exposure corrector): concept validated; 0.71 power at perfect
  specificity qualifies it as a fallback/composite tier. The full
  exposure-informed EM (E5) is not yet justified by the gap it would close
  on pancreas; revisit after breast.

## Recommended next steps

1. Breast 5K phase 2: anchor endothelial + fibroblast (the two missing
   signatures there), same harness - the case anchoring was designed for.
2. Composite correction: ensemble removal for strict-tier(ambient-baseline)
   content, exposure-calibrated budgets for shared-gene pairs - the two
   halves are complementary corners of the measured trade-off surface.
3. Productize: realized-removal reporting in `evaluate()`/`correction$rules`,
   `fit$recover_factors(audit)`, and `correct(vote = "auto", ...)`.
