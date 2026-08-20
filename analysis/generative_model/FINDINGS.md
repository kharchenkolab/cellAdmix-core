# Findings: a generative model for admixture correction (pancreas, Xenium 377-gene panel)

Self-contained record of the model, the measurements, and the conclusions.
Companion files: `REPORT.md` (gate summary), `LOG.md` (chronological
observe-orient-decide-act log), `results/` (all tables cited below),
numbered scripts `00`-`04` (implementation and evaluation). Everything was
evaluated on the cached pancreas run
`examples/xenium_pancreas_membrane_377_full/out/runs/fit_manual_rank9_invsqrt_kl`
(140,335 segmented cells, 377 genes, 6,186,500 cell-assigned molecules)
with the audit's gene-split harness
(`analysis/audit_guided/01_split_eval.R`) and its 39 detected
source-to-target cell-type pairs. Throughout, "admixture" means molecules
physically present inside a segmented cell that originate from another
cell's expression (spilled fragments, out-of-plane material, or freely
diffusing ambient molecules); "induced expression" means transcription by
the target cell itself that happens to correlate with the neighborhood
(e.g. a gene switched on near another cell type).

---

## 1. The model as implemented

### 1.1 Generative form

For each annotated target cell type T, the observed molecule count
y_cg of every cell c of that type at gene g is modeled as an independent
Poisson draw with expected value

    t_c * [ sum_k theta_ck F_kg          own expression programs
          + sum_S alpha_cS psi_Sg        contamination, one term per
                                         detected source pair S -> T
          + beta_c a_g                   ambient background
          + u_cS rho_cS m_Sg             sparse induced-expression term
          + eps_g ]                      uniform floor

where t_c is the cell's total observed molecule count, so every quantity
in brackets is a fraction of the cell's content. The observed count at
each (cell, gene) entry is then split among the components in proportion
to their fitted rates (the standard Poisson-mixture posterior), and the
corrected count keeps the own + induced + floor share. Corrections were
delivered to the evaluation as a matrix of removed molecule mass and
applied to the identical "before" matrix the harness itself uses.

### 1.2 Components, priors, and where each prior comes from

**Own expression.** F_k are gene profiles (each summing to 1 over genes)
of K_T expression programs per type, with per-cell weights theta_ck. They
were initialized from the empirical profiles of factor-labeled molecules:
the fitted run assigns each molecule a factor label from the underlying
9-factor factorization, and for each type the programs are the molecule
profiles of the type's own factors - the factor aligned to that type in
the audit's factor-to-type map, plus any unaligned factor carrying at
least 5% of the type's molecule content (this admits factor 5, a second
ductal/exocrine state carrying 22% of ductal content). Mural/pericyte has
no aligned factor and was initialized from its type pseudobulk profile.
Programs are re-estimated inside the fit (see 1.3), with two protections
that mattered:

- *Dose-weighted learning.* In the program update, each cell's
  contribution is weighted by 1/(1 + 20 * prior contamination fraction),
  so programs are learned predominantly from lightly contaminated cells.
  Without this, the second exocrine program (initialized from factor-5
  molecules, a state shared with ductal cells) absorbed the ductal
  material inside exposed exocrine cells, and the Ductal -> Exocrine pair
  scored 0.34-0.39 removal power; with it, 0.87-0.91.
- *Strict-gene support.* In the production configuration, program values
  are fixed to zero on "strict" genes - the harness's own designation for
  pool genes with near-zero native expression in the target (baseline
  under 5% of the source level). A gene the target does not express
  cannot be part of an own program.

**Contamination.** psi_S is the source type's cytoplasmic expression
profile - estimated from source-cell molecules with overlaps_nucleus == 0
- and, critically, it is the profile of the source cells that actually
border the target type: cells within 30 um of the nearest target-type
cell (widened to 100 um, then to all cells, if fewer than 30,000
cytoplasmic molecules are available). Section 3 explains why this
near-interface restriction matters. Profiles are refined over two outer
rounds: after all seven types are fitted, each type's per-gene "own
fraction" (the fitted share of that gene's counts attributed to the
type's own expression) multiplies its profile, so a source's profile is
progressively cleaned of the contamination it itself contains. Finally,
each pair's profile is zeroed on target-owned genes (genes whose
highest-expressing type is the target): contamination there cannot be
distinguished from own expression, and the model deliberately leaves it
in place. This makes removal of the target's own marker genes
structurally impossible rather than merely calibrated to be small.

The per-cell contamination fraction alpha_cS is the model's central
random variable. Its prior mean lambda_cS is a monotone (isotonic)
dose-response measured from the data, exactly in the spirit of the
audit's exposure-regression corrector: target cells are stratified by
exposure (the number of source-type cells among the cell's 15 nearest
neighbors, the harness's own exposure definition), the pooled content of
guide genes per stratum is fitted by weighted isotonic regression, the
stratum-0 level is subtracted, and the excess rate is converted to a
whole-profile contamination fraction by dividing by the guide genes'
share of psi_S (capped at 0.6). Guide genes are the A halves of the
harness's marker pools in the validation arm and the full pools in
production. The prior has strength 300 pseudo-molecules (comparable to
the median cell's total, so evidence and prior carry similar weight for
a typical cell), and the posterior update from the cell's own expression
is clipped to at most 5-fold the prior mean - and to exactly zero where
the prior mean is zero. The clipping is not cosmetic: it is what makes
the shuffled-exposure control collapse (section 2.4). An unbounded
evidence-driven alpha keeps finding the true contamination even when the
exposure covariate is destroyed - the contamination is genuinely there -
and would therefore fail the pre-registered control; the bounded update
keeps removal anchored to the measured exposure signal while still
letting each cell's expression move its own estimate several-fold.

After the main fit converges, up to three "top-up" passes re-measure the
residual dose-response on the corrected counts (excluding
induced-support genes, whose gradient is retained by design) and add it
to the prior dose, then continue the fit. This is the direct analog of
the regression corrector's top-up and contributed about +0.007 weighted
power.

**Ambient background.** Restricted to the union of the target's strict
genes. The intended estimate - profile and level of target cells more
than 250 um from any source cell - turned out to be unavailable: no
pancreas target type has even 20,000 molecules in such cells (the tissue
is dense; every target sits near its sources somewhere). The fallback
encodes the audit's own measurement that the far-field level is 0.5-0.7
of the pooled zero-exposure baseline: the ambient level prior is half the
zero-exposure strict-gene rate, with the profile taken from the same
cells. The per-cell ambient scale beta_c is updated by the same bounded
Gamma-posterior step (cap 5-fold), which lets it absorb the structured,
near-interface part of the baseline that the audit showed is not truly
ambient. The ambient component is on only in the production
configuration; the validation arm and the shuffled control run without
it, following the regression corrector's own conventions (an
exposure-independent removal of a fixed gene set registers as "power"
under any exposure permutation, so controls must exclude it).

**Induced expression.** A sparse per-pair, per-gene term. Support is
selected by a proportionality screen: for each pair, the per-gene
exposure-linked excess (pooled exposed-cell counts minus the
zero-exposure rate scaled to exposed totals) is regressed on the source
profile, and the residual is standardized by a variance that includes
counting (Poisson) noise plus a 15% multiplicative profile-uncertainty
term. Genes with standardized residual above 8 and excess above 50
molecules enter the support. The overdispersion term is the sparsity
mechanism that survived testing: a pure Poisson z threshold flagged 153
genes on real data, including AMY2A with excess 148,714 against a
proportional expectation of 138,248 - a 1.08-fold deviation that is
plainly profile noise but reaches z = 27 because z grows with the square
root of the channel size. On the final screen, 28 pair-gene combinations
(25 unique genes) are flagged, all between 2.3-fold and 22-fold above
their proportional expectation. Flagged genes receive a per-gene rate
m_Sg (per unit of measured dose) and a per-cell activity multiplier
rho_cS with a Gamma prior of mean 1 (shape 2, capped at 30), so
induction concentrated in a subset of cells is retained where it occurs.
Flagged content is retained, not removed - this is the model's
distinguishing policy.

**Floor.** A uniform 0.2%-of-content rate spread over all genes keeps
posterior fractions well-behaved at entries where every structured
component is near zero.

### 1.3 Fitting

Expectation-maximization with multiplicative closed-form updates: 40
iterations per target type, then the top-up passes (15 iterations each),
all repeated for two outer rounds to refine the source profiles. Cell
parameters (theta, alpha, beta, rho) update per iteration under their
priors; gene-side parameters (F, m) update by responsibility-weighted
re-estimation. The validation arm additionally excludes all held-out
B-half genes from the alpha evidence sums (rescaled by the source
profile's share on the remaining genes), so the dose is estimated
without ever seeing the genes it is scored on. Runtime is a few minutes
per arm on the full dataset with 4 threads.

### 1.4 Implementation choices that mattered (with reasons)

1. **Near-interface source profiles** (30 um): moves source-composition
   gradients out of the induced term and into the contamination term
   where they belong; see section 3.
2. **Zeroing source profiles on target-owned genes**: converted
   own-marker false removal from 1.4-4.4% (by type) to exactly 0.000,
   with no tuning, at no measured cost to removal power.
3. **Dose-weighted program learning**: prevented "contamination as a
   cell state"; worth ~0.5 of removal power on the worst pair.
4. **Bounded, exposure-anchored alpha with a zero trap**: the only
   configuration in which the shuffled-exposure control collapses while
   per-cell adaptation is preserved. An early analysis (LOG loop 0)
   showed that any exposure-independent fractional removal f of pool
   genes registers as measured power f under the harness - so the
   control's integrity has to come from the dose prior, not from
   downstream thresholds.
5. **Overdispersed proportionality screen**: one interpretable knob
   (profile coefficient of variation, 0.15) replaced a fold-change guard
   and separated real induction (2.3-22 fold) from channel-size
   artifacts (1.08 fold at z=27).
6. **Evaluating through the audit's own harness verbatim** (counts
   exported from the identical build, corrections returned as removal
   matrices): every number here is commensurable with the audit's
   published numbers.

---

## 2. Observations

### 2.1 Gene-split validation (the honest generalization test)

Setup: each pair's 20-gene marker pool is split into an A half (visible
to guidance) and a B half (held out); the model's dose-response,
evidence sums, and induced screen see only A-side and non-pool genes;
power is the fraction of exposure-linked B-half excess removed, weighted
across pairs by excess size.

- Generative model: **0.974** weighted (median 1.000), own-marker false
  removal 0.000.
- Exposure-regression corrector v3: 0.906, false removal 0.000.
- Membrane ensemble (default): 0.961, but at 7.5% own-marker false
  removal and deep baseline erosion (audit numbers).

Per-pair, the generative model beats the regression corrector on 25 of
39 pairs, ties (within 0.02) on 12, and loses on 2. The wins are where
the regression corrector's budget mechanism under-delivers:
Ductal -> Endothelial 1.000 vs 0.303, Ductal -> Exocrine 0.901 vs 0.364,
Ductal -> Mural 0.978 vs 0.578, Fibroblast -> Mural 0.862 vs 0.500,
Fibroblast -> Exocrine 0.999 vs 0.687, Endocrine -> Immune 1.000 vs
0.731, Mural -> Exocrine 0.975 vs 0.718, Immune -> Endocrine 1.000 vs
0.784. The two losses are fibroblast-source pairs with heavy shared-gene
content: Fibroblast -> Endothelial 0.963 vs 1.000 and Fibroblast ->
Immune 0.935 vs 0.967. Interpretation: a dose measured on half the
markers generalizes to the other half much better when it is delivered
through a full source profile with per-cell posterior updates than when
it is spent as a fixed budget by waterfilling; the fibroblast losses are
pairs where forced budget delivery on shared genes happens to match the
yardstick better than the model's own-vs-foreign split.

### 2.2 Production configuration and the retention ledger

Full pools guide the dose; ambient tier on; induced term on. Weighted
power **0.982** (median 1.000); removed 1,206,425 of 6,186,500 molecules
(19.5%; the regression corrector removes 1.16M, the ensemble 1.99M). One
pair falls below 0.8: Endothelial -> Ductal at 0.735. Decomposed gene by
gene, every B-half gene of that pair is removed at 0.92-1.00 except
GNG11 - flagged by the screen at 5.4-fold disproportionality and 75%
retained - which accounts for 603 of the pair's 613 residual excess
molecules. With the induced term disabled the model scores 0.985 with
every pair at or above 0.8; on the "admixture-only" yardstick (the same
metric computed on B halves excluding the 28 flagged genes) it scores
0.986 with a per-pair minimum of 0.915. The pre-registered 0.99 bar is
therefore NOT met, and the gap has two named parts: (i) ~0.4 points of
deliberately retained flagged-gene excess, and (ii) ~1.4 points of
exposure-linked excess on shared genes that the model attributes to the
target's own expression near tissue interfaces and that three top-up
passes do not (and, by design, will not) remove.

Real-data retention behavior (`results/gm_production_induced_retention.csv`):
the 28 flagged pair-gene combinations retain 79.5% of their
exposure-linked excess in aggregate while non-flagged pool content is
removed at 96-99%. Representative per-gene retained fractions: CFTR
0.805 (89,229 of 110,867 excess molecules in Exocrine -> Ductal - the
single largest retention decision; CFTR is the canonical pancreatic duct
gene and duct cells at acinar interfaces are exactly where it peaks),
SPIB 0.979, FXYD2 1.05 (fully retained), TM4SF18 0.808, GNG11 0.754,
SEMA3C 0.891, TNC 0.895, ACTG2 0.899, APOLD1 0.901, EDNRB 0.94-0.96,
CTLA4 0.829, PMP22 0.820, GPC3 0.812; the least-retained flagged genes
(ASPN 0.658, GLIPR1 0.663, BASP1 0.670) are those whose entries overlap
heavy proportional contamination, where the entry-level split hands part
back to removal. The flagged set reads as three biological programs:
a ductal duct-cell/interface program (CFTR, CA4, PROX1, FXYD2, TM4SF18,
GNG11), an activated-fibroblast program (SEMA3C, TNC, THY1, APCDD1,
ASPN, BASP1, ALDH1A3, PMP22, C5orf46), and immune activation (SPIB,
CTLA4, MS4A1, CXCL10, KIT, GLIPR1), plus endothelial/mural interface
genes (ADAMTS1, STEAP4, EDNRB, APOLD1, CAVIN1, ACTG2, GPC3).

### 2.3 Spike-in test (four planted molecule classes)

Reusing the audit's T7 generator verbatim (same pair Exocrine ->
Endothelial, same random seed, same placement rules): 400 exposed
endothelial cells per class receive (a) in-plane fragments - 4-12
molecules copied from a 1.5-um ball of a nearby exocrine cell, placed at
the target edge facing the source; (b) out-of-plane fragments - the same
material placed over the footprint with z pushed to the section
surfaces; (c) "induction" - molecules of three exocrine marker genes
placed like the target's own molecules, scaled with exposure; (d)
ambient singletons - single pool-gene molecules placed uniformly.
Planted molecule counts: 1,928 / 1,898 / 3,750 / 800. The model was
refitted on the spiked matrix in production configuration and every
planted molecule scored by the removal probability of its entry (count
level) and by the per-molecule posterior (molecule level).

Results, primary design (induction on the pool's top-3 genes AQP8,
GATM, AMY2A):

| class | count-level removal | molecule-level removal | target |
|---|---|---|---|
| in-plane fragments | 0.826 | 0.809 | >= 0.8: pass |
| out-of-plane fragments | 0.846 | 0.830 | >= 0.8: pass |
| ambient singletons | 0.851 | 0.782 | removed: pass |
| planted induction | 1.000 (retention 0.000) | 0.898 (retention 0.102) | retention >= 0.8: **fail** |

Why the primary retention test fails, measured rather than asserted.
The three planted genes are the pancreas's largest acinar channels. In
the exposed endothelial population they already carry 3,010 / 12,401 /
26,358 molecules of exposure-linked excess that is proportional to the
exocrine profile (genuine admixture); the planted addition is 3,750
molecules in total - a ~14% increment on those channels. A Poisson-exact
proportionality test detects it (AMY2A residual z = 21.1), but the same
Poisson-exactness is what flagged 153 genes on real data including the
1.08-fold AMY2A artifact; allowing any plausible profile error destroys
the detection - the standardized residual falls to 1.4 / 0.9 / 0.7 at
profile coefficients of variation 0.10 / 0.15 / 0.20. The molecule-level
rescue was also measured: the two spatial-coherence features have
class-conditional means 0.29 vs 0.18 (distinct other source-owned genes
within 1 um, admixed vs native reference molecules) and 0.17 vs 0.50
(native-marker neighbors within 1 um) - likelihood ratios per molecule
of order 0.7-1.4, hopeless against a prior of order 1000:1 that a
strict-gene molecule in an exposed cell is contamination. This agrees
with the audit's own T7 conclusion that patch statistics separate
fragments and induction from ambient but not from each other.

Results, supplementary variant - the identical procedure with induction
planted on mid-rank pool genes (CA4, DIRAS3, PROX1), where the planted
excess is several-fold the proportional expectation, i.e. the regime the
real flagged genes (2.3-22 fold) actually occupy:

| class | count-level | molecule-level |
|---|---|---|
| planted induction retention | **0.949** | **0.942** |
| in-plane fragment removal | 0.817 | 0.803 |
| out-of-plane fragment removal | 0.836 | 0.821 |
| ambient removal | 0.749 | 0.694 |

The screen fires, the per-cell activity multiplier concentrates the
induced term in the spiked cells, and retention passes with fragments
still removed. The retention machinery works wherever retention is
statistically possible at all.

### 2.4 Shuffled-exposure control

Each pair's exposure values permuted across its target cells; ambient
off (matching the regression corrector's control convention, since any
exposure-independent removal reads as power under permutation). The
permuted dose-response is flat, every prior contamination mean is ~0,
and the zero-trap keeps every alpha at zero: 22,332 molecules removed
(0.36% of the production arm's 1.21M), weighted power **0.018** against
the required < 0.05. The regression corrector's control is 0.007.

### 2.5 Own-marker false removal and downstream purity

False removal of the target types' own top-30 marker genes: **0.000
pooled, 0.000 worst type** - structural, as explained in 1.4.2. (The
ensemble's figure is 7.5%.) Downstream check - same-label k-nearest-
neighbor purity of corrected counts in principal-component space on the
audit's 5,000-cell subset, same code and seed: original 0.804 (the
audit's number reproduced exactly), generative model 0.963, regression
corrector 0.967, ensemble 0.968. The three correctors are within half a
point; the generative model achieves this while removing 1.21M molecules
against the ensemble's 1.99M and while keeping the flagged biology.

---

## 3. The interface-activation finding

The audit's per-gene proportionality test (T3) flagged CXCL6, CFB and
PPP1R1B - inflammation/activation genes - as induced: their
exposure-linked excess in cells near ductal/tumor cells is far above
their share of the ductal expression profile (global-profile folds 4-20,
z 9-45 across the ductal-source pairs). The generative model's screen,
which measures proportionality against the profile of the source cells
that actually border each target type, does not flag them. The
difference is in the source cells themselves:

Measured cytoplasmic profile shares (percent of molecules) in ductal
cells, globally vs within 30 um of each target type:

| gene | global share | near exocrine | near immune | near endothelial | near fibroblast |
|---|---|---|---|---|---|
| CXCL6 | 1.57% | 2.57% (1.64x) | 2.31% (1.47x) | 2.46% (1.57x) | 2.07% (1.32x) |
| PPP1R1B | 0.96% | 1.57% (1.64x) | 1.08% (1.13x) | 1.50% (1.56x) | 0.99% (1.03x) |
| CFB | 3.70% | 4.25% (1.15x) | 4.07% (1.10x) | 4.29% (1.16x) | 3.84% (1.04x) |

Ductal cells at tissue interfaces express these genes 1.1-1.6 times more
than the ductal average (they are, plausibly, the inflamed tumor edge).
When the proportionality test uses the interface-local profile, the
excess in the neighboring target cells is largely explained as
transferred material: for Ductal -> Exocrine, CXCL6 excess is 5,252
against an interface-proportional expectation of 6,254 (0.84-fold,
standardized residual -1.1 - i.e. fully proportional), CFB 13,786 vs
10,344 (1.33-fold, z 2.2), PPP1R1B 6,378 vs 3,818 (1.67-fold, z 4.1);
for Ductal -> Immune, CXCL6 1,334 vs 646 (2.07-fold, z 6.5), CFB
1,382 vs 1,141 (1.21-fold, z 1.3), PPP1R1B 533 vs 303 (1.76-fold,
z 4.3). All fall below the 8-z screening threshold that the genuinely
retained genes exceed comfortably.

What distinguishes the two readings - "the target cell was induced to
express CXCL6" vs "the target cell received CXCL6 molecules made by an
activated ductal neighbor"? Three pieces of evidence favor the second:
(i) the interface enrichment is measured in the source cells' own
cytoplasm, independent of any target; (ii) the residual excess in
targets after accounting for it is proportional across the trio in the
way transfer predicts (the same 1.1-1.6x factors appear in every target
direction, tracking the source-side enrichment rather than any
target-specific response); (iii) the audit's own patch-composition
measurement showed that source-labeled material inside target cells
matches the source's cytoplasmic profile (cosine 0.63) - and the profile
it matches is necessarily that of the *adjacent* source cells. A genuine
target-side induction should instead produce excess that scales with
exposure but not with the source-side composition shift, which is
exactly what the retained genes (CFTR at 4.6-fold, SPIB at 22-fold, etc.)
do. A residual mild induction of CXCL6 in immune cells (2.07-fold, z
6.5) cannot be excluded - it sits just under threshold - but the bulk of
the audit's flagged excess for this trio is transferred material.

Implication for the audit's induced-gene screening: the T3 test should
be re-based on interface-local source profiles (source cells within
~30 um of the target type, cytoplasmic molecules, decontaminated). On
this panel that one change removes the entire headline trio from the
induced list and replaces it with the CFTR/SEMA3C/SPIB-class genes at
2.3-22 fold - and it should be combined with an overdispersed residual
(profile coefficient of variation ~0.15) so that channel-size artifacts
like AMY2A (1.08-fold at Poisson z 27) never enter. Both changes are
implemented in `01_model.py` (`near_source_profile`, `induced_screen`)
and can be lifted into the audit directly.

---

## 4. The identifiability limit (general statement)

Any correction method of this family - one that models a target cell's
counts as its own expression plus admixture proportional to a source
expression profile, using exposure or distance as the driving covariate
- can only distinguish induced expression from admixture through one of
three signals, and each has a hard floor:

1. **Disproportionality.** Induced expression of gene g is detectable if
   its exposure-linked excess exceeds what the source profile predicts by
   more than the profile's own uncertainty. If the source profile is
   known to relative precision c (on this dataset, interface-composition
   effects alone are 10-60%, so c >= 0.1 is realistic), induction on a
   channel that already carries proportional excess E is detectable only
   when the induced amount is comparable to c * E. Induction of a
   source's top marker genes - where E is largest - is therefore
   invisible until it is enormous: the spike-in's planted 14% bump on
   AQP8/GATM/AMY2A is undetectable in principle (overdispersed z ~ 1),
   while the same planted mass on mid-rank genes is recovered at 94%.
   What is recoverable is *disproportionate* induction: genes expressed
   by the source weakly or not at all, or induced several-fold beyond
   their source share. All 25 genes retained on real data are in this
   regime (2.3-22 fold).
2. **Spatial coherence.** Admixed fragments are gene-diverse patches
   matching the source's cytoplasm; induced molecules sit among the
   target's own molecules. This signal exists (the audit's 1.7x
   co-location enrichment) but is weak at single-molecule resolution in
   dense tissue: measured class-conditional feature means of 0.29 vs
   0.18 (source-gene neighbors within 1 um) and 0.17 vs 0.50
   (native-marker neighbors) give per-molecule likelihood ratios well
   under an order of magnitude - enough to move a balanced prior (it
   contributes the 0.10 retention on the primary spike and sharpens the
   mid-rank result to 0.94), never enough to overturn a strict-gene
   identity prior of order 1000:1. Deeper segmentation-free spatial
   modeling might improve the contrast but faces the audit's T5 finding
   that even control-probe noise is spatially structured.
3. **Cell-level concentration.** Induction concentrated in a cell
   subpopulation, on few genes, deviates from the many-gene
   proportional signature of fragments. The per-cell activity
   multiplier exploits this once a gene is flagged; it cannot by itself
   flag a gene hidden under signal floor (1).

Consequently: no method in this family - however elaborate - can
simultaneously (a) be robust to realistic source-profile error and (b)
retain induction planted on the source's dominant expression channels.
The pre-registered spike-in retention gate, as instantiated by the T7
script's choice of the top-3 pool genes, sits on the wrong side of that
boundary; the same gate on mid-rank genes is passed at 0.94-0.95. Claims
for this model - and for any successor - should be stated as: induced
expression is preserved when it is disproportionate to the local source
profile; proportionate induction of a source's own top markers is
indistinguishable from admixture at count level and will be removed.

---

## 5. Conclusions and recommendations

### 5.1 When to use which corrector

- **Exposure-regression corrector v3**: when the goal is maximal
  flattening of exposure-linked marker gradients at minimal complexity
  and runtime (0.992 production power, seconds, simple to audit). Its
  known cost: it removes induced/interface expression by construction
  (CFTR loses its 89k-molecule interface excess) and its dose
  generalizes noticeably worse (0.906 gene-split validation).
- **Generative model**: when corrected counts feed biology that may
  legitimately correlate with neighborhood - interface cell states,
  induced programs, ligand-receptor analyses - or when per-cell, per-
  component decompositions (contamination fraction per source, ambient
  share, induced activity) are wanted as outputs in their own right.
  Costs: 0.982 aggregate production power (the retention policy itself),
  minutes of runtime, more moving parts.
- The ensemble's remaining niche (aggressive strict-tier cleanup) is
  matched by both at far lower own-marker cost (0.000 vs 7.5%).

### 5.2 A composite

The natural composite keeps each mechanism where it is strongest and is
straightforward to assemble from existing parts:

1. regression-style ambient/strict tier (or this model's beta component)
   for strict-gene baseline content;
2. the generative dose model (near-interface profiles, bounded per-cell
   alpha, top-up) for the exposure-linked removal - it is the better
   generalizer;
3. the overdispersed interface-local proportionality screen as a
   *retention filter* applied to ANY corrector's removal matrix: molecules
   of flagged (pair, gene) combinations get their removal scaled back by
   the fitted induced share. Applied to the regression corrector, this
   would preserve most of its 0.992 on non-flagged content while adding
   the retention capability - likely the best product configuration.

### 5.3 Follow-up experiments, in priority order

1. **Port the two screen upgrades into the audit** (interface-local
   source profiles + overdispersed residual) and re-run the audit's T3
   across pancreas and breast. Small effort (both functions exist in
   `01_model.py`); directly changes the audit's induced-gene reporting
   and re-tests CXCL6/CFB/PPP1R1B on independent tissue.
2. **Retention-filter composite** (5.2.3) evaluated under the same six
   gates. Moderate effort (a day): reuses the regression corrector and
   this model's screen; expected to dominate both parents on the gate
   table (production ~0.99 on non-flagged content, retention ~0.8).
3. **Breast 5K replication.** The whole evaluation is scripted; the
   breast run exercises the anchoring case (two missing native factors)
   and tests whether the interface-activation finding (and the flagged
   gene programs) replicate on an independent panel and tissue. Moderate
   effort, mostly compute.
4. **Spike-in boundary map.** Sweep the planted-induction gene rank and
   mass (the `GM_IND_RANK` parameterization already exists) to chart
   retention as a function of planted-to-proportional ratio, giving the
   quantitative version of section 4's boundary for the docs. Small
   effort.
5. **Molecule-level realization for production output.** The per-molecule
   posterior currently serves the spike-in; exporting posterior origin
   probabilities for all molecules of flagged genes (and optionally
   using the factor label as a third feature - unavailable for synthetic
   spikes but available for every real molecule) would give downstream
   users molecule-resolved retention. Moderate effort; the exact
   subset-posterior machinery is written.
6. **Distance-decay reference for the dose baseline.** The audit showed
   the zero-exposure baseline is contaminated (median 1.41x excess
   inflation when re-priced); the model currently absorbs this through
   the ambient tier's evidence updates, but an explicit lateral-distance
   decay term in the dose prior would price it per pair. Larger effort;
   revisit after the breast replication shows whether it limits anything.

