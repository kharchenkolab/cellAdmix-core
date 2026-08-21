# Admixture and induced expression: a generative model

Admixture correction removes source-marker molecules whose content in
target cells rises with the number of source-type neighbors. But two
different realities produce such exposure-linked excess: molecules
physically transferred from neighboring cells (spilled fragments,
out-of-plane material, ambient background), and transcription by the
target cell itself that responds to its neighborhood — a cell switching
on an activation program at a tissue interface. A correction that
removes both flattens every gradient, at the cost of erasing exactly the
biology that spatial data is often collected to see: interface cell
states, induced programs, ligand–receptor signaling. This document
describes a generative model that decomposes each cell's counts into
origin components so that induced expression can be told apart from
transferred material and retained, the evidence for where that
separation is possible and where it is not, and the practical correction
that resulted. The complete experimental record lives in
`analysis/generative_model/` (model, evaluation scripts, findings) and
`analysis/audit_guided/` (the shared evaluation harness); the audit
vocabulary used throughout — exposure, marker pools, the ambient
reference — is defined in [benchmarks.md](benchmarks.md).

## The model

For each annotated target cell type, the observed count of cell $c$ at
gene $g$ is modeled as a Poisson draw with expected value

$$y_{cg} \sim \mathrm{Poisson}\Big(t_c \big[ \textstyle\sum_k \theta_{ck} F_{kg}
\;+\; \sum_S \alpha_{cS}\, \psi_{Sg} \;+\; \beta_c\, a_g \;+\;
\sum_S u_{cS}\, \rho_{cS}\, m_{Sg} \;+\; \epsilon_g \big]\Big)$$

where $t_c$ is the cell's total molecule count, so every bracketed
quantity is a fraction of the cell's content. The terms are, in order:
the cell's own expression programs ($F_k$, gene profiles with per-cell
weights $\theta_{ck}$); contamination from each detected source type $S$
(profile $\psi_S$, per-cell fraction $\alpha_{cS}$); ambient background
($a_g$, per-cell scale $\beta_c$); a sparse induced-expression term
active only on screened gene support ($u_{cS}$; per-gene rate $m_{Sg}$
with a per-cell activity multiplier $\rho_{cS}$); and a small uniform
floor $\epsilon_g$. After fitting by expectation–maximization, each
observed count is split among the components in proportion to their
fitted rates, and the corrected count keeps the own, induced, and floor
shares.

Four structural choices carry most of the model's behavior, each adopted
for a measured reason:

- **Interface-local source profiles.** $\psi_S$ is estimated from the
  cytoplasmic molecules of the source cells that actually border the
  target type (within 30 µm of the nearest target-type cell, widening
  when sparse), not from the source type as a whole. Cell states are not
  uniform across a tissue: source cells at an interface can express
  activation genes several-fold above the source average, and material
  transferred from them carries that elevated share. Against the global
  profile such genes would be misread as induced in the target (see the
  interface-activation finding below).
- **Source profiles are zeroed on target-owned genes.** Contamination on
  a gene whose top expresser is the target cannot be distinguished from
  own expression, and the model deliberately leaves it in place. This
  makes removal of a target's own marker genes structurally impossible —
  measured own-marker false removal is exactly zero — rather than merely
  calibrated to be small.
- **The contamination fraction is anchored to the measured exposure
  dose.** The prior mean of $\alpha_{cS}$ is a monotone (isotonic)
  dose–response of marker content on the number of source cells among
  the cell's nearest neighbors, and the posterior update from the cell's
  own expression is bounded (at most five-fold the prior, and exactly
  zero where the prior is zero). The bound is what keeps removal tied to
  the demonstrated spatial signal: with an unbounded update the model
  keeps finding the true contamination even after the exposure covariate
  is destroyed, and the permutation control fails.
- **Own programs are learned from lightly contaminated cells.** Program
  updates weight each cell down by its prior contamination fraction;
  without this, a program shared between source and target absorbs the
  transferred material as if it were a cell state, and the affected
  pair's removal collapses.

## The induced-expression screen

The induced term's gene support is chosen by a proportionality test.
Transferred material samples the source transcriptome, so a gene's
exposure-linked excess in target cells should be proportional to its
share of the (interface-local) source profile. For each pair, the
per-gene excess is regressed on the profile, and the residual is
standardized by counting noise plus a multiplicative profile-uncertainty
term of 15%. The second term is the sparsity mechanism: a pure
counting-noise test flags large transfer channels for deviations of a
few percent — on pancreas it flagged AMY2A at 1.08-fold its expectation,
plain profile noise reaching formal significance only because the
channel is huge — while the real induced genes deviate 2.3- to 22-fold.
Genes past the threshold receive the induced term; their disproportionate
excess is retained where it occurs (the per-cell activity multiplier
concentrates it in the responding cells) while their proportional share
is still removed as transfer.

This same screen, applied to the audit's marker panels, is how the
package excludes induced genes from admixture estimation (see
[benchmarks.md](benchmarks.md) and the binding documentation); flagged
genes and their measured disproportionality are reported by
`markers()`.

## The interface-activation finding

Whether the comparison profile is global or interface-local changes the
verdict for a whole class of genes. On pancreas, the inflammation genes
CXCL6, CFB and PPP1R1B show exposure-linked excess in cells near
ductal/tumor tissue at 4- to 20-fold their share of the *global* ductal
profile — reading as induction in the neighbors.
But the enrichment sits in the source cells themselves: ductal cells
within 30 µm of another cell type carry 1.1–1.6 times the ductal-average
cytoplasmic share of these genes (they are, plausibly, the inflamed
tumor edge), and against that interface-local profile the neighbors'
excess is largely proportional — transferred material made by activated
source cells, not a response by the receiving cells. The genes the
screen does flag as induced behave differently: their excess tracks
exposure but not the source-side composition shift, at folds of 2.3–22
(on pancreas, led by CFTR — the canonical duct gene, peaking in duct
cells at acinar interfaces, at 3.9-fold its transfer expectation —
with interface programs of fibroblasts, immune and endothelial/mural
cells behind it). The same screen applied across datasets flags
inflammatory and myoepithelial contractile programs on breast 5K (CCL2,
IL6, LIF; MYH11, MYLK, CNN1) and the classic stress and chemokine
programs on NSCLC (FOS, JUNB, HSPA1A/B, CCL3/4, CXCL2/3).

## Validation

All evaluation runs through the audit's gene-split harness: each pair's
20-gene marker pool is split in half, everything that guides the model
sees only the guide half, and performance is the fraction of
exposure-linked excess removed on the held-out half, weighted across the
39 detected pancreas pairs.

- **Held-out-gene generalization**: 0.974 (median 1.000), versus 0.906
  for an exposure-regression corrector guided by the same genes and
  0.961 for the default molecule-vote ensemble — the ensemble reaching
  its number at 7.5% own-marker false removal, the model at 0.000. A
  dose measured on half the markers generalizes to the other half best
  when delivered through a full source profile with per-cell posterior
  updates.
- **Production behavior**: 0.982 weighted removal, taking out 1.21M of
  the dataset's 6.19M cell-assigned molecules (the ensemble removes
  1.99M for the same downstream benefit). The flagged genes retain 79.5%
  of their exposure-linked excess in aggregate — for CFTR, 89k of its
  111k excess molecules stay in the ductal cells — while unflagged
  marker content is removed at 96–99%.
- **Spike-in test**: with four classes of molecules planted in one
  pair's exposed cells, the model removes 81–85% of in-plane fragments,
  out-of-plane fragments, and ambient singletons, at every setting of the
  planted induction. The planted induction itself is retained or removed
  depending on how far it exceeds the transfer expectation of its genes —
  the full sweep over gene ranks and planted amounts is shown in
  Figure 1a and discussed under the identifiability limit below.
- **Replication**: on the CosMx NSCLC dataset (98,002 cells, 960-gene
  panel, 27 detected pairs) the same configuration reaches 0.964 on
  held-out genes, removes 9.0% of the molecules, keeps own-marker
  removal at zero and the permutation control at 0.013 — and its screen
  flags the same kind of interface biology: heat-shock and
  immediate-early programs (HSPA1A/B, DUSP5, RGS2), macrophage
  activation (GLUL, SRGN), and collagen programs at
  endothelial–fibroblast borders.
- **Controls**: permuting each pair's exposure values collapses removal
  to 0.018 (the dose prior, not a downstream threshold, carries the
  control); own-marker false removal is 0.000 pooled and for the worst
  type; downstream cell-type k-nearest-neighbor purity rises from 0.804
  to 0.963, within half a point of the ensemble's 0.968.

## The identifiability limit

Any method of this family — target counts modeled as own expression
plus transfer proportional to a source profile, driven by a spatial
covariate — can separate induction from admixture through only three
signals, and each has a floor. Disproportionality is bounded by how well
the source profile is known: with interface-composition effects alone
shifting profiles 10–60%, induction on a source's dominant expression
channels is invisible until it rivals the transfer itself (the planted
14% increment on the pancreas's top acinar genes is undetectable in
principle, while the same mass on mid-rank genes is recovered at 94%).
Spatial coherence — fragments arrive as gene-diverse patches, induced
molecules sit among the cell's own — is real but weak at
single-molecule resolution in dense tissue: per-molecule likelihood
ratios stay under an order of magnitude, against identity priors of
order a thousand to one. And cell-level concentration of induction can
sharpen a gene already flagged but cannot flag one hidden below the
disproportionality floor. The honest claim for this model and any
successor is therefore: induced expression is preserved when it is
disproportionate to the local source profile; induction proportionate
to a source's own top markers is indistinguishable from admixture at
count level and will be removed.

Figure 1a shows this boundary measured directly. Synthetic induced
expression was planted across a range of gene ranks and amounts, and
its retained fraction charted against the gene's total exposure-linked
excess relative to the transfer expectation: retention turns on where
the excess exceeds the expectation about three-fold and is essentially
complete beyond ten-fold, while excess below two-fold is removed as
transfer. The real flagged genes, plotted at their measured
disproportionality and their realized retained fractions, fall on the
same curve — the model treats real and planted induction alike.
Figure 1b shows the model's largest single retention decision as a
worked example: CFTR content of ductal cells rises fourteen-fold with
exocrine exposure, and the model splits that rise into a retained share
(native expression plus the induced duct-cell program, most of the
signal) and a removed transferred share, rather than flattening the
gradient wholesale as a marker-based correction would.

![Figure 1](figures/generative_fig1.png)

**Figure 1. The separation boundary, and a worked example.** **(a)**
Fraction of planted induced expression retained, as a function of the
gene's exposure-linked excess relative to the transfer expectation
(blue: spike-in sweep over planted gene ranks and amounts, pancreas);
red squares: the real flagged genes at their measured
disproportionality. The dotted line marks 80% retention. **(b)** CFTR
in ductal cells stratified by the number of exocrine cells among the 15
nearest: observed rate (red), the model's retained share (blue), and
the removed transferred share (shaded); the dotted line is the rate in
unexposed ductal cells.

## A practical composite

The model's two transferable mechanisms — the interface-local
overdispersed screen and per-gene retention of the disproportionate
excess share — graft directly onto the much simpler exposure-regression
corrector: measure the removal dose on unflagged genes only, and cap
each flagged gene's removable content at its proportional (transferred)
share, leaving the induced share with the cell. Under the same harness
this composite removes transferred content at the regression corrector's
full power (0.992 on held-out genes excluding flagged ones, every pair
covered) while retaining 79.6% of the flagged-gene excess — against
8.4% when the retention cap is ablated, and matching the full model's
79.5% — at zero own-marker removal and a clean permutation control. The
same script reproduces every gate on the CosMx NSCLC dataset (0.998 on
non-flagged content, 80.1% retention). Scripts and per-pair tables:
`analysis/audit_guided/17_composite.R` and
`analysis/audit_guided/results/`.

The three correctors order naturally by need. The exposure-regression
corrector is the simplest and fastest way to flatten exposure-linked
marker gradients, at the known cost of removing induced and interface
expression along with transfer. The composite adds retention of
identified induced genes at essentially no cost in removal power or
complexity. The full generative model additionally offers the best
held-out-gene generalization and per-cell, per-component decompositions
(contamination fraction per source, ambient share, induced activity) as
outputs in their own right, at the cost of more machinery and minutes of
runtime.
