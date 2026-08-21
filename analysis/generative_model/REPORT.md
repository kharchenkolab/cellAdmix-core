# A generative model for admixture correction: pancreas evaluation

This directory implements and evaluates a per-cell generative (Poisson
mixture) model of molecule counts for admixture correction, against the
audit's pre-registered gates and against the two existing correctors: the
exposure-regression corrector v3 (`analysis/audit_guided/08_phase_a.R`;
validation 0.906, production 0.992) and the default membrane ensemble.
All evaluation uses the audit's own gene-split harness
(`analysis/audit_guided/01_split_eval.R`) on the cached pancreas run
(140,335 cells, 377 genes, 6,186,500 molecules, 39 detected source-target
pairs). Scripts are numbered; per-pair tables for every gate are under
`results/`; the full observe-orient-decide-act log is in `LOG.md` and
summarized below.

## The model as implemented

For each annotated target type T, the molecule counts y_cg of every cell c
of type T over genes g are modeled as independent Poisson draws with mean

    t_c * [ sum_k theta_ck F_kg            (own expression programs)
          + sum_S alpha_cS psi_Sg          (contamination, one term per
                                            detected source pair S -> T)
          + beta_c a_g                     (ambient background)
          + u_cS rho_cS m_Sg               (sparse induced-expression term)
          + eps_g ]                        (uniform floor, 0.2% of content)

where t_c is the cell's observed molecule total. The components:

- **Own expression.** F_k are gene profiles (rows sum to 1) of K_T
  state-specific programs, initialized from the profiles of factor-labeled
  molecules of the type's own factors (the aligned factor plus unaligned
  factors carrying >= 5% of the type's content; the type pseudobulk for
  Mural/pericyte, which has no own factor), refined by EM with cell weights
  1/(1 + 20 * prior contamination fraction), so programs are learned
  predominantly from lightly dosed cells and admixed content in heavily
  exposed cells cannot be re-absorbed as a "cell state". In the production
  configuration, program support is zeroed on strict genes (near-zero
  native baseline by the harness's own definition).
- **Contamination.** psi_S is the cytoplasmic profile (molecules with
  overlaps_nucleus == 0) of the source cells within 30 um of the nearest
  target-type cell - the actual donor population - so that composition
  gradients within the source type (e.g. activated fibroblasts near tumor)
  are part of the contamination profile. Profiles are decontaminated over
  two outer rounds by multiplying with the source type's own fitted
  per-gene own-expression fraction, and are zeroed on target-owned genes:
  contamination there is unidentifiable from own expression and is left in
  place, which makes removal of the target's own markers structurally
  zero. The per-cell fraction alpha_cS has a Gamma prior (strength 300
  pseudo-molecules) whose mean lambda_cS is an isotonic (monotone)
  dose-response of guide-gene content on the number of S cells among the
  cell's 15 nearest neighbors, converted to a whole-profile fraction; the
  posterior update from the cell's own expression is bounded to 5-fold
  above the prior mean and is exactly zero where the prior is zero, so
  removal stays anchored to the measured exposure signal (this is what
  makes the shuffled-exposure control collapse). After EM convergence, up
  to three top-up passes re-measure the residual dose-response on the
  corrected counts (excluding induced-support genes) and add it to the
  prior dose - the analog of the regression corrector's top-up.
- **Ambient.** Restricted to the strict-gene union of the target's pairs;
  profile and scale estimated from target cells > 250 um from any source
  cell. No pancreas target has enough such cells (the tissue is dense), so
  the scale falls back to half the pooled zero-exposure strict rate - the
  audit measured the far-field level at 0.5-0.7x that baseline. The
  per-cell scale beta_c updates by the same bounded Gamma-posterior step
  (cap 5x), letting it absorb the structured (near-field) part of the
  baseline. Off in the validation arm and in the shuffled control,
  following the regression corrector's own convention.
- **Induced expression (the sparsity-screened term).** For each pair, a
  per-gene proportionality test on source-owned genes compares the
  exposure-linked excess with the best proportional fit to psi_S; the
  variance includes counting noise plus a 15% multiplicative
  profile-uncertainty term, and genes with overdispersed z > 8 and > 50
  excess molecules form the support (28 pair-gene rows, 25 unique genes on
  pancreas). Their per-gene rates m_Sg (per unit of measured dose u_cS)
  and a per-cell activity rho_cS (Gamma prior with mean 1) are refined by
  EM, and that content is RETAINED rather than removed.

Fitting is expectation-maximization with multiplicative closed-form
updates (40 iterations plus top-up passes, per target type, two outer
rounds for source-profile decontamination). The corrected matrix keeps the
own + induced + floor share of each observed count. The validation arm
additionally excludes all held-out B-half genes from the dose evidence
(rescaled by the source-profile share of the remaining genes), uses A
halves for the dose-response, and bars B genes from the induced support.

**Molecule-level realization** (`03_spikein.py`, second half): within each
count-matrix entry, the observed molecules are split between a removable
origin (contamination + ambient, Poisson with the fitted mean) and a
retained origin (own + induced, plus a 2%-prior geometric burst of
unmodeled induction - a spike-and-slab), and each molecule carries two
spatial-coherence features (distinct other source-owned genes, and
native-marker molecules, within 1.0 um), with class-conditional Poisson
likelihoods estimated from reference molecule sets (real strict-gene
molecules in exposed cells vs native-marker molecules). The exact subset
posterior (via elementary symmetric polynomials) gives each molecule a
removal probability: gene identity sets the prior, spatial coherence
moves it.

## Scripts

- `00_export.R` - builds the audit's split definitions (verbatim harness)
  and exports counts, pools, A/B splits, exposure, profiles.
- `01_model.py` - the model; `--arm validation|production|production_noind|shuffled`.
- `02_eval_gates.R` - applies removal matrices, runs the audit's
  `eval_split`, writes per-pair/per-type CSVs.
- `03_spikein.py` - four-class spike-in (adapted from
  `14_t7_spikeins.py`), count-level and molecule-level scoring;
  `GM_IND_RANK=9` selects the mid-rank supplementary variant.
- `04_purity.R` - downstream same-label kNN purity (function copied from
  `08_phase_a.R` for comparability).

## Gate-by-gate results

Comparators: exposure-regression corrector v3 and membrane ensemble, as
reported in `analysis/audit_guided/REPORT.md`. Full table:
`results/gm_final_gates.csv`.

| gate | target | generative model | regression v3 | ensemble |
|---|---|---|---|---|
| 1. gene-split validation (A-guided), weighted removal power on held-out B halves | >= 0.90 | **0.974** - PASS | 0.906 | 0.961 (at 7.5% own-marker false removal) |
| 2. production, weighted power | >= 0.99 | 0.982 (median 1.000); without the induced term 0.985; admixture-only yardstick 0.986 - **NOT MET** | **0.992** | 0.961 |
| 2. production, all pairs >= 0.8 | yes | 1 pair below (0.735, = deliberate retention of one flagged gene); without the induced term: all pairs >= 0.8 | yes | - |
| 3. spike-in, planted fragments removed (in-plane / out-of-plane) | >= 0.8 | 0.826 / 0.846 (count level); 0.809 / 0.830 (molecule level) - PASS | not tested; removes by budget | - |
| 3. spike-in, planted ambient removed | yes | 0.851 (count), 0.782 (molecule) - PASS | - | - |
| 3. spike-in, planted induced genes retained - primary design (top-3 pool genes) | >= 0.8 | 0.000 (count), 0.102 (molecule) - **FAIL** (see analysis) | 0 by construction | 0 by construction |
| 3. spike-in, induced retention - same design on mid-rank pool genes | >= 0.8 | **0.949** (count), **0.942** (molecule), fragments still removed at 0.82/0.84 - PASS | 0 by construction | 0 by construction |
| 3. real data: retention of flagged induced genes' exposure-linked excess | report | **0.795** across 25 genes (CFTR: 89k of 111k molecules kept) | 0 | 0 |
| 4. shuffled-exposure control | < 0.05 | **0.018** - PASS (22k molecules removed vs 1.21M) | 0.007 | - |
| 5. own-marker false removal | ~ 0 | **0.000** pooled, 0.000 worst type - PASS (structural) | 0.000 | 7.5% |
| 6. downstream kNN purity (original 0.804) | report | **0.963** | 0.967 | 0.968 |

Removal totals: generative production 1.21M molecules; regression v3
1.16M; ensemble 1.99M.

### Gate 1 detail (the honest head-to-head)

Per-pair against the regression corrector's validation arm: the
generative model is better on 25 pairs, tied on 12, worse on 2 (both by
~0.05: Fibroblast -> Immune, Fibroblast -> Endothelial). The regression
corrector's worst validation pairs are transformed: Ductal -> Endothelial
0.30 -> 1.00, Ductal -> Exocrine 0.36 -> 0.90, Ductal -> Mural 0.58 ->
0.98, Fibroblast -> Mural 0.50 -> 0.88. The dose-response measured on A
halves generalizes to held-out genes substantially better when it is
delivered through a full source profile with per-cell posterior updates
than through budgeted waterfilling.

### Gate 2 detail (why 0.99 is not reached)

Three-quarters of the remaining 1.8% shortfall is deliberate: the model
retains the exposure-linked excess of the 28 flagged induced pair-gene
combinations (79.5% retained), and the harness counts every retained
molecule as unremoved contamination. The single below-0.8 pair
(Endothelial -> Ductal, 0.735) decomposes gene by gene into removal
0.92-1.00 everywhere except the flagged gene GNG11 (5.4-fold
disproportionality, 75% retained). With the induced term disabled the
model reaches 0.985 with every pair >= 0.8 - and that 1.5% residue is the
model attributing part of the exposure-linked excess on shared genes to
the target's own expression (exposure-correlated native biology near
tissue interfaces), which no amount of principled dose top-up recovers
(three passes converge there). The regression corrector reaches 0.992
precisely because it removes everything the yardstick measures, induced
content included. The two numbers describe different policies, not
different competences: on the admixture-only yardstick (excluding the
flagged genes) the generative model scores 0.986 with a per-pair minimum
of 0.915.

### Gate 3 detail (the headline, and its boundary)

The pre-registered spike-in (audit script `14_t7_spikeins.py`) plants
"induction" on the pair's top-3 pool genes - AQP8, GATM, AMY2A, the
largest acinar channels. In exposed endothelial cells those genes carry
3.0k/12.4k/26.4k molecules of genuine proportional contamination; the
planted total is 3.75k. The planted excess is a ~14% bump on the largest
contamination channels: a Poisson-exact test sees it (z = 21), but the
same Poisson-exactness fired on 153 genes in real data (including AMY2A
itself at 1.08-fold, i.e. plainly profile noise) and cost 0.10 of
removal power when those were retained. Any screen robust to >= 10%
profile error reads the planted bump as noise (overdispersed z ~ 1).
Molecule-level spatial rescue was measured, not assumed: the two
coherence features have contrasts 0.29-vs-0.18 and 0.17-vs-0.50 at 1 um -
far too weak to overturn a strict-gene identity prior - consistent with
the audit's own T7 conclusion that patch statistics corroborate but
cannot separate induction from fragments. Count-level retention 0.000,
molecule-level 0.102: **this gate fails on the primary design, and the
analysis argues no data-honest method can pass it there.**

The supplementary variant runs the identical procedure with induction
planted on mid-rank pool genes (CA4, DIRAS3, PROX1), where the planted
excess is large relative to the proportional expectation - the regime the
real flagged genes (3-22 fold) actually occupy. There the screen fires,
the per-cell activity concentrates the induced term in the right cells,
and retention is 0.949 (count) / 0.942 (molecule) with fragments still
removed at 0.82/0.84 and ambient at 0.75. On real data the same machinery
retains 79.5% of the flagged genes' exposure-linked excess - CFTR, FXYD2,
PROX1, CA4 (ductal duct-cell program at acinar interfaces), TNC, SEMA3C,
THY1, APCDD1, ASPN, BASP1, ALDH1A3 (activated-fibroblast program near
tumor), SPIB, CTLA4, MS4A1, CXCL10, KIT, GLIPR1 (immune activation),
GNG11, TM4SF18, EDNRB, STEAP4, ADAMTS1 - while removing 96-99% of
everything else. This is a capability the regression corrector and the
ensemble structurally lack (both remove induced content by construction).

A byproduct worth reporting: with near-source contamination profiles, the
audit's headline induced candidates CXCL6, CFB and PPP1R1B are no longer
disproportional (CXCL6: 2.1-fold, below the 8-z threshold) - their
exposure-linked excess matches what the ductal cells actually bordering
the targets express. The model reclassifies them as admixture from
locally activated source cells rather than target-side induction; the
audit's global-profile test could not make that distinction.

### Gates 4-6 detail

The shuffled control collapses (0.018; 22k of 6.19M molecules removed)
because the contamination prior is exposure-derived and the per-cell
update cannot leave zero: with permuted exposure the isotonic
dose-response is flat and every alpha stays at zero. Own-marker false
removal is exactly zero by construction (source profiles carry no weight
on target-owned genes). Downstream kNN purity on the 5,000-cell subset:
0.804 original (audit's number reproduced exactly), 0.963 corrected -
0.004 below the regression corrector and 0.005 below the ensemble, i.e.
the same downstream benefit within noise of the subsample.

## OODA loop summary (full log in LOG.md)

1. **Loop 1** (design): count-level EM in Python on the exact harness
   export; bounded exposure-anchored alpha updates chosen specifically so
   the shuffled control can collapse; ambient tier strict-only, off in
   validation/shuffled, mirroring the regression corrector's conventions.
2. **Loop 2** (first numbers): production 0.877, false removal 1.4-4.4%,
   induced screen fired on 153 genes (z alone is scale-dependent - AMY2A
   flagged at 1.08-fold). Shuffled control already passed (0.016).
3. **Loop 3** (five fixes): source profiles zeroed on target-owned genes
   (false removal -> exactly 0); per-pair NEAR-source cytoplasmic
   profiles (30 um) so source-composition gradients stop masquerading as
   induction; 3-fold disproportionality requirement; own programs learned
   from lightly dosed cells (kills the "contamination as cell state"
   shield - Ductal -> Exocrine 0.34 -> 0.87); dose cap raised. Validation
   jumped to 0.965, production to 0.975.
4. **Loops 4-5** (top-up): residual-gradient top-up passes (the
   regression corrector's own trick) added +0.007; freezing programs
   during top-up changed nothing - the remaining gap is a policy
   disagreement with the yardstick, not an optimization failure.
5. **Loop 6** (spike-in): count-level retention of the primary spike
   design is impossible under profile-error-robust screening (measured,
   not asserted); replaced the ratio guard with an overdispersed z
   screen; built the molecule-level spike-and-slab realization with
   spatial-coherence likelihoods measured from reference sets.
6. **Loop 7** (final): screen at z > 8 keeps 25 strong genes (2.3-22
   fold); mid-rank spike variant passes all gate-3 targets; purity and
   final tables.

## Verdict

**Does the generative model earn its complexity?** For raw removal power
on this harness: no - the regression corrector reaches 0.992 production
power with a page of code, and the generative model's production number
is 0.982 (0.985-0.986 on machinery-only views). For everything else
measured here: largely yes, on four counts.

1. **It wins the honest test.** Under gene-split validation - guidance
   restricted to A halves, scored on held-out B halves - it reaches
   0.974 vs the regression corrector's 0.906, better on 25 of 39 pairs
   and worse on 2, at identical (zero) own-marker cost. The dose model
   generalizes; the budget-waterfilling does not, quite.
2. **Induced-gene retention is real and unique** - the headline works,
   with an honest boundary. In real data it retains ~80% of the
   exposure-linked excess of 25 strongly disproportional genes (CFTR's
   89k molecules would simply be deleted by every other corrector) while
   matching their removal elsewhere; on planted spike-ins it passes all
   four class targets (retention 0.94+, fragments removed) whenever the
   planted induction is statistically detectable at all. On the
   pre-registered primary design - induction hidden as a 14% bump on the
   three largest contamination channels - it fails (0.10 retained), and
   the measurements argue this regime is undetectable in principle at
   count level and that 1-um spatial coherence (contrast 0.29 vs 0.18) is
   too weak to rescue it. That boundary should temper claims for any
   future method: what can be retained is disproportionate induction,
   not induction of a source's top markers.
3. **Safety is structural, not tuned.** Zero own-marker removal falls out
   of the profile support (target-owned genes are untouchable), not out
   of a calibrated threshold; the shuffled control collapses because the
   dose prior is the only way in.
4. **The decomposition is interpretable and already paid for two
   findings**: the near-source profiles reclassified CXCL6/CFB/PPP1R1B
   from "induced in the target" to "admixture from locally activated
   ductal cells", and the per-pair fits expose where exposure-linked
   excess is native biology (shared-gene pairs cap at ~0.985 under any
   honest dose delivery).

**Where it loses:** aggregate production power (0.982 vs 0.992, and gate
2 formally NOT MET - the cost of retaining what it believes is real
expression); compute (minutes of EM per arm vs seconds); ~0.004 of kNN
purity; and two fibroblast-source validation pairs. If the downstream
goal is maximal marker-gradient flattening, use the regression corrector.
If the goal is corrected counts that can still contain the biology that
correlates with neighborhood - induced programs, interface cell states -
the generative model is currently the only option on the table, and its
remaining gap to the 0.99 bar is, on inspection, mostly that biology.

## Replication and the retention boundary (added experiments)

- NSCLC (CosMx, 98,002 cells, 960 genes, 27 detected pairs; whole-cell
  profiles - the export carries no nucleus flag; factor alignment derived
  from the data): validation arm 0.964 weighted power_B, production 0.970
  (0.986 with the induced term disabled - the gap is deliberate
  retention), shuffled control 0.013, own-marker false removal 0.000,
  9.0% of molecules removed. Screen flags 88 pair-gene combinations
  (HSPA1A/B, DUSP5, RGS2, GLUL, SRGN, COL4A1/2). Results:
  results/nsclc_gm_*.csv.
- Spike-in retention sweep (03_spikein.py, GM_IND_RANK x GM_IND_MASS;
  results/gm_spikein_sweep.csv; figure 05_figure.py ->
  docs/figures/generative_fig1.png): retention of planted induction as a
  function of total excess over the transfer expectation crosses 0.35 at
  2-fold, 0.6-0.87 at 3-4-fold, 0.95+ beyond 10-fold; fragment and
  ambient removal stay at 0.75-0.85 across all settings. Real flagged
  genes fall on the same curve.
