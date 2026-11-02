# Multi-sample processing: OODA log

Working notes for the `dev_multisample` branch. Each loop records what was
observed, the decision taken, and what was done. Results and scripts live
next to this file; the design is in
`docs/development/multisample_proposal.md`.

## Loop 1 (2026-11-01): orientation and terms

Observed. Per-sample processing yields its own factors, rules, and
molecule labels per sample, so replicates diverge after correction. The
core already has the pieces for a shared factor space: persisted restart
member pools, a fixed-profile fit (`nmf_fixed_h`, R only), Hungarian
matching on gene-ownership profiles (`nmf_stability.cpp`), pair-level
score tables, and `correct(rules = )` with an explicit rule table. The
fixed-profile fit carries no member pool, so corrections on it fall back
to single-fit.

Terms fixed with the user. *Joint factorization*: one NMF on neighborhood
rows pooled from all samples. *Consensus factorization*: per-sample fits,
factors matched across samples and aggregated into shared profiles, each
sample projected onto the fixed profiles. *Reference-sample projection*:
one sample's profiles projected onto the others (comparison arm only).
*Fixed-profile projection*: the labeling step common to all.

Decisions. Namespace cell ids `<sample>_<cell_id>`; per-sample rule
support is reported, never gated; `sample_allocation` option,
proportional by default; tabular input and per-sample crops in scope;
consensus factorization first, joint factorization only if the
comparison calls for it.

Done. Worktree `wt-dev_multisample` with isolated C++ build, Python 3.11
venv, and R 4.4.3 library; all three suites green at the branch point
(C++ 1/1, Python 35/35, R 290/290). Breast replicate bundles downloaded
to `/home/pkharchenko/cellAdmix/datasets/` with the published
supervised labels for both sections. Found that the 2022 bundles mark
unassigned molecules with cell id `-1`, which the loader did not treat
as unassigned (2.8M molecules in replicate 1 would have formed one
cell); fixed in the Xenium and tabular loaders.

## Loop 2: comparison protocol and the consensus layer

Plan. Three arms per dataset: separate (today's pipeline per sample),
consensus (new collection layer), joint stand-in (both samples in one
molecule table with offset coordinates and a `sample_id` column, through
the tabular path). A fourth arm injects the joint stand-in's profiles
into per-sample fixed-profile runs, isolating pooled training from the
rest. Criteria: per-sample loss gap, cross-sample match quality against
within-sample restart stability, replicate divergence of per-type
pseudobulks after correction, per-sample audit sensitivity and
own-marker false removal, removal fraction per type per sample, rule
support.

Done. Core: `nmf_fixed_h_members` (a member pool persisted with a
fixed-profile run, the fixed profiles leading as member 0; member
agreement with the fixed profiles reported in the stability slots),
`match_factor_profiles` (Hungarian on ownership profiles, rectangular
inputs padded), `align_profiles_to_genes`, `training_loss_under_profiles`
(a run's sampled training rows under given profiles in the run's own
objective; weighted LS usages by multiplicative updates, KL by the fixed-
profile projector in the run's transformed feature space), fixed-profile
runs keep the given factor order, manifest fields `fixed_profiles` and
`fixed_profile_members`; `-1` and integer cell ids accepted by the
loaders. Bindings: `nmf_fixed_h` in Python, member pools in both, factor
matcher, training loss, member-pool reader. Collection layer in R
(`r/R/collection.R`) and Python (`collection.py`) with the same
structure: `cellAdmixCollection(datasets)` → `fit()` → `fit_consensus()`
or `fit_profiles()` → `score_bridge()` (per-sample caps split by
`sample_allocation`) → pooled `summary()` (Stouffer over samples) →
`rules()` with `sample_support` → per-sample ensemble `correct()` →
per-sample `audit_admixture()`/`evaluate()`. Consensus profiles are
averaged in raw gene space (KL-variant profiles divided by the run's gene
weights), unit-sum normalized, rescaled to the weighted mean of the
original row sums; unmatched factors appended as sample-specific. Tests:
C++ 1/1 (new cases for matching, alignment, LS usages, the `-1`
sentinel), Python 38/38, R all green (collection end to end, fixed
profiles with members, reference projection).

Arms running. Pancreas halves: separate (done, 4 min each), consensus
(10 consensus factors from 9 + 9; 26 pooled rules, 25 kept, 12 supported
by both halves), joint stand-in through the tabular path (6.19M molecules,
50 mm gap between halves). Breast replicates: separate arms (rank 23,
fits 2–3 min, bridge scoring 3–4 min, 213 and 235 kept rules).

Found. The tabular parquet reader (`read_parquet_arrow` in `src/tabular.cpp`)
segfaults inside libparquet's string-column reader on every parquet file
written by this venv's pyarrow 26 (plain or large strings, with or without
dictionaries, writer versions 1.0 to 2.6); the same rows as gzip CSV load
and process normally, and the Xenium parquet readers are unaffected. The
C++ tests never covered parquet tabular input. The joint stand-in runs on
CSV; the parquet reader is a separate defect to return to.

## Loop 3 (2026-11-01): pancreas halves, five arms, bridge scoring

Observed (`analysis/multisample/results/pancreas_*.csv`; left and right
halves of one section, `ls_nmf`, 10 restarts, bridge scoring, default
ensemble). Separate processing of the two halves already diverges in the
pipeline's own choices: 14 versus 20 kept rules, rule-pair Jaccard 0.55,
and the left half's correction removes almost none of its measured
admixture (median pair sensitivity 0.03 against 0.38 on the right).

| arm | sample | leakage removed | median pair sensitivity | own-marker false removal | molecules removed |
|---|---|---|---|---|---|
| separate | left | 0.30 | 0.03 | 0.5% | 6.6% |
| separate | right | 0.20 | 0.38 | 0.3% | 6.4% |
| consensus | left | 0.53 | 0.54 | 1.8% | 13.2% |
| consensus | right | 0.24 | 0.46 | 0.3% | 7.0% |
| joint stand-in | left | 0.38 | 0.55 | 0.7% | 8.5% |
| joint stand-in | right | 0.17 | 0.15 | 0.2% | 5.7% |
| joint profiles, per-sample scoring | left | 0.44 | 0.20 | 1.2% | 9.5% |
| joint profiles, per-sample scoring | right | 0.18 | 0.08 | 0.3% | 6.2% |
| reference (left) projection | left | 0.38 | 0.51 | 2.7% | 11.4% |
| reference (left) projection | right | 0.32 | 0.46 | 2.1% | 11.5% |

Consensus diagnostics: 10 consensus factors (9 matched, one left-only
factor appended; matched ownership correlations 0.36–0.88, median 0.81,
against within-sample restart stability 0.86 and 0.79); loss gap 0.7%
(left) and 0.6% (right) against restart spreads of 0.4% and 1.4%; the
joint stand-in's profiles sit at 1.0% on both halves, and the left half's
own profiles at 2.8% on the right half. Member agreement with the shared
profiles 0.88–0.95. Factor usage divergence is large for the ductal/tumor
type (0.45) and small elsewhere (0.05–0.16), as expected for two regions
of one tumor section rather than replicates.

Divergence between halves (cell-weighted mean per-type pseudobulk
correlation; genes differentially expressed between the halves' 1 mm
block pseudobulks): uncorrected 0.924 / 389; separate 0.846 / 540;
consensus 0.861 / 501; joint stand-in 0.918 / 387; joint profiles
0.907 / 420; reference projection 0.890 / 566.

Oriented. The consensus route meets the factor-space criteria (loss gap
within or near the restart spread, cross-sample matches comparable to
restart stability) and gives the best cleanup in both halves at low false
removal; the joint stand-in cleans the left half as well but leaves the
right half largely uncorrected (sensitivity 0.15, below the separate
arm). The pooled scoring, not the shared profiles, is what rescued the
left half: the same joint profiles with per-sample scoring and Stouffer
pooling (joint-profiles arm) clean poorly, while the consensus profiles
with the same pooling clean well, so the profile set and the pooling
interact and the per-sample scoring caps deserve a look. The divergence
metric is confounded on these halves: they are different tissue regions,
and removing admixture, which blurs cell types toward their neighborhood,
can legitimately increase differences between regions; the arms that
remove more (consensus, reference) show more divergence. The breast
replicates are the proper test of that criterion.

Decided. Keep the consensus route as the primary; run the breast
replicates through all arms; run the thinning ladder (spreads 0.1, 0.2,
0.4, 0.8) for batch robustness; look at why the joint profiles pool
poorly under per-sample scoring (rule counts 20 vs 23 for the same
profiles).

Also learned: Rscript reads a script incrementally while executing, so
editing `run_separate.R` during the breast replicate 1 run broke that run
(rerun started).

Observed afterwards. No bridge-scored arm, separate or shared, produces a
rule for exocrine → ductal/tumor, the right half's largest admixture pair
(341,000 estimated molecules, 29% of its total), nor for the other
exocrine-source pairs; the per-pair sensitivity for exocrine → ductal is
zero in every arm. This is the known weakness of bridge scoring on the
pancreas section (audit.md: diffuse exocrine leakage, best handled by
membrane scoring with `invsqrt_kl`), and it bounds what any pooling can
achieve under bridge scoring here. The membrane configuration is being
run on the halves (separate first, consensus and reference next).

Looked into the joint-profiles arm's smaller rule set (20 kept against the
joint stand-in's 23 for the same profiles). The seven rules the joint
stand-in alone keeps are all borderline (its p-values 0.003–0.07 on 6–166
cell pairs); under per-sample scoring the two halves disagree on them
(e.g. immune → fibroblast: left p = 0.05 on 137 pairs, right no effect
on 267 pairs), so Stouffer pooling leaves them above the threshold even
though the per-sample scoring examined more pairs in total (65,000
against 38,000). One pooled test over mixed pairs can be carried by one
half; pooling per-sample p-values is conservative when the halves
disagree, which is the intended behavior of the collection layer, and
the four rules the joint-profiles arm alone keeps are of the same
borderline kind. The difference is rule-flip noise near the threshold,
not a defect of either pooling.

## Loop 4 (2026-11-01): breast replicates, first diagnostics

Observed while the breast arms run. Consensus of the two replicates:
19 of replicate 2's 23 factors match replicate 1's (ownership correlation
0.52–0.96, median 0.80, against within-sample restart stabilities 0.85
and 0.81); the 4 unmatched factors (appended as sample-specific) have
within-sample restart stability 0.62, 0.04, 0.46, 0.45, i.e. three of
them are among replicate 2's least reproducible factors. So on true
replicates the unmatched factors are mostly restart noise rather than
sample-specific biology, which argues for reporting unmatched factors
together with their own-sample stability, and possibly for an option
that appends only factors stable within their own sample. Loss gap
−1.3% and −1.2% (the 27 consensus factors fit each replicate's training
rows better than its own 23, partly because there are more of them);
member agreement with the consensus 0.92 and 0.91. Pooled rules: 298,
250 kept, 217 supported by both replicates; the separate arms kept 213
and 235.

Also observed: the tabular parquet reader segfaults on a 10x-written
parquet file too, so the failure is in the `ReadTable` code path on this
machine rather than in any particular writer; a C++ regression test for
parquet tabular input is being added to locate it.

Resolved. The tabular parquet crash is the use-after-free already fixed on
the unmerged remote branch `fix/parquet-lifetime-and-determinism`
(38ce649: the streaming store source returned a record-batch reader whose
parquet FileReader had been destroyed; Arrow IO threads then read freed
memory). A standalone C++ driver loaded the same files fine only by timing.
The fix commit, with its many-row-group regression test, is cherry-picked
onto `dev_multisample` (8b44a7a); the C++ suite passes with it.

Breast separate arms (bridge, `ls_nmf`, rank 23, 10 restarts). Replicate
1: 172 detected pairs, 12.1M estimated admixed molecules (30% of 39.7M),
leakage removed 0.62, median pair sensitivity 0.85, own-marker false
removal 1.6%, 22% of molecules removed, 213 kept rules. Replicate 2: 149
pairs, 10.3M (35% of 29.1M), 0.68, 0.89, 1.5%, 24% removed, 235 kept
rules. The replicate 1 script run that was corrupted by the mid-run edit
had in fact finished its correction and audit before failing on its last
statement; the rerun wrote identical outputs over it, and the correction
directory is complete.

Thinning ladder, spread 0.2 (right half thinned gene-wise, realized
log-efficiency SD 0.19). Consensus: 11 factors (7 of the left half's 9
matched at 0.66–0.87, 2 appended), loss gap −0.2% and −0.5%; cleanup
left 0.40 / 0.45, right 0.40 / 0.57 (separate 0.30 / 0.03 and
0.20 / 0.53); joint stand-in left 0.43 / 0.58, right 0.42 / 0.44; false
removal ≤ 0.8% everywhere. Divergence (correlation / DE genes):
uncorrected 0.920 / 400, separate 0.852 / 531, consensus 0.876 / 478,
joint 0.926 / 382. The thinning did not produce split factors in the
joint stand-in (rank fixed at 9) and did not change the consensus
matching picture beyond one more unmatched left factor.

## Loop 5 (2026-11-01): pancreas halves, membrane scoring with `invsqrt_kl`

Observed (`results/pancreas_membrane_*.csv`; separate, consensus, and
reference arms; the joint stand-in cannot run membrane scoring through the
tabular path). The recommended configuration removes far more admixture
than bridge scoring here (separate: leakage removed 0.73 and 0.85, median
pair sensitivity 0.97 and 0.99, own-marker false removal 2.3% and 3.6%,
25% and 33% of molecules removed), and separate processing diverges the
halves strongly: rule-pair Jaccard 0.63 (26 versus 31 pairs), pseudobulk
correlation between halves 0.92 → 0.49, DE genes 389 → 910.

| arm | sample | leakage removed | median pair sensitivity | false removal | removed |
|---|---|---|---|---|---|
| separate | left | 0.73 | 0.97 | 2.3% | 24.6% |
| separate | right | 0.85 | 0.99 | 3.6% | 33.1% |
| consensus | left | 0.75 | 0.96 | 2.2% | 24.4% |
| consensus | right | 0.77 | 0.96 | 2.8% | 28.9% |
| reference (left) | left | 0.73 | 0.98 | 2.1% | 24.9% |
| reference (left) | right | 0.45 | 0.92 | 1.1% | 18.0% |

Consensus: 41 kept rules, 30 supported by both halves; divergence
0.74 / 549 (separate 0.49 / 910). Diagnostics show the predicted weak
link of this variant: within-sample restart stability 0.47 and 0.37, only
4 of the left half's 9 factors match the right half's (0.35–0.63), five
are appended (14 consensus factors), and the loss gap is 1.8% (left) and
3.5% (right) against restart spreads of 1.5% and 1.4%, so the consensus
represents the right half worse than one of its own restarts would. The
cleanup still holds (0.75 / 0.77 at 2.2–2.8% false removal), while the
reference projection fails on the other half (0.45). 

Oriented. For the stable `ls_nmf` variant the consensus passes every
factor-space criterion; for `invsqrt_kl` the per-sample fits themselves
disagree, and matching the selected runs alone is too thin. The natural
next step is to match and average over every restart of every sample
(the consensus NMF recipe over the member pools) before fixing the shared
profiles, which the stored restart pools make cheap; the joint
factorization remains the alternative if that is not enough. Decision
deferred until the breast replicates and the ladder are in.

## Loop 6 (2026-11-01): thinning ladder, pancreas halves, bridge scoring

Observed (`results/ladder_<spread>/`; the right half thinned gene by gene
with log-efficiency deviations of spread 0.1, 0.2, 0.4, 0.8, realized SDs
0.09, 0.19, 0.37, 0.75; separate, consensus, and joint stand-in arms;
`ls_nmf`, rank 9, 10 restarts).

Batch-separating factors: none at any spread. In the joint stand-in the
largest cosine between two factors' profiles stays at 0.17–0.38, far from
a split, and the per-factor concentration of molecules in one half
(0.13–2.0 times proportional) is the same at spread 0 as at 0.8: it
reflects the halves' different composition (tumor region on one side),
not the thinning. Consensus matching likewise: 7 of the left half's 9
factors match the right half's at 0.60–0.88 at every spread, 2 are
appended, and the loss gap stays within the restart spread throughout
(−0.5% to +0.7% against spreads of 0.4–1.1%). Gene-wise efficiency
differences up to a log SD of 0.75 therefore do not split factors under
either route on this tissue at this rank.

Cleanup varies more between ladder steps than between routes: each step is
a different random thinning and a different right-half fit, and bridge
scoring on this tissue is seed-sensitive (audit.md). Consensus leakage
removed, left/right: 0.56/0.34 (0.1), 0.40/0.40 (0.2), 0.35/0.60 (0.4),
0.54/0.37 (0.8); separate right 0.21–0.38 with left fixed at 0.30; joint
right 0.15–0.42. Own-marker false removal of the consensus rose to 2.6–4.8%
at spread 0.4 and 3.4–3.8% at 0.1 (0.6–2.1% elsewhere), the one soft spot,
against ≤ 0.8% for the joint stand-in; it tracks the amount removed (up to
21% of the thinned half's molecules), not the spread. Divergence between
halves is confounded here as before (consensus ≈ separate, joint ≈
uncorrected, with the arms that remove more diverging more).

Oriented. The batch-robustness question is settled for the realistic
range: no split factors and stable matching. What varies is the bridge
scoring's rule set from fit to fit, which the ensemble vote and the pooled
evidence do not fully absorb on this tissue; that is the single-sample
limitation documented in audit.md rather than a multi-sample effect.

Realistic spread. The gene-wise log-efficiency difference between the two
breast replicate sections, measured within shared cell types as the
count-weighted SD of the centered log ratio of each gene's share (genes
with at least 200 counts per type pair), is 0.06–0.16 for most types
(tumor 0.06–0.07, T cells 0.07–0.09, stroma 0.16), 0.19–0.25 for the
sparse perivascular, mast, and second macrophage types, and 0.43 for the
two DCIS types, whose difference includes genuine between-section state
differences and the admixture itself; pooled 0.22. The ladder's spreads
0.1–0.2 therefore bracket a true replicate pair, and 0.4–0.8 are well
beyond it.

## Loop 7 (2026-11-01): breast replicates, four arms, bridge scoring

Observed (`results/breast_*.csv`; two serial FFPE sections, published
labels, `ls_nmf` rank 23, 10 restarts, bridge scoring; 19 types measured).

| arm | sample | leakage removed | median pair sensitivity | false removal | removed | kept rules |
|---|---|---|---|---|---|---|
| separate | rep1 | 0.62 | 0.86 | 1.6% | 22.2% | 213 |
| separate | rep2 | 0.68 | 0.89 | 1.5% | 24.1% | 235 |
| consensus | rep1 | 0.60 | 0.85 | 1.4% | 21.6% | 250 (217 in both) |
| consensus | rep2 | 0.54 | 0.71 | 1.3% | 18.0% | 250 |
| reference (rep1) | rep1 | 0.59 | 0.82 | 1.4% | 21.2% | 246 |
| reference (rep1) | rep2 | 0.53 | 0.71 | 1.7% | 17.4% | 246 |
| joint stand-in | rep1 | 0.62 | 0.88 | 1.2% | 22.7% | 260 |
| joint stand-in | rep2 | 0.68 | 0.88 | 1.5% | 24.0% | 260 |

Replicate concordance (cell-weighted mean per-type pseudobulk correlation
between replicates; genes differentially expressed between the replicates'
1 mm block pseudobulks):

| arm | correlation | DE genes |
|---|---|---|
| uncorrected | 0.988 | 242 |
| separate | 0.902 | 1066 |
| consensus | 0.862 | 1069 |
| reference (rep1) | 0.901 | 638 |
| joint stand-in | 0.991 | 117 |

The replicate test is decisive. Separate processing makes true replicates
diverge (rule-pair Jaccard 0.58; stromal correlation 0.990 → 0.918,
B cells 3 → 73 DE genes; myoepithelial KRT15+ loses 26% of its molecules
in one replicate and 6% in the other). The joint stand-in removes the
divergence entirely: every cell type's correlation is at or above its
uncorrected value, DE genes fall below the uncorrected count (B, T, and
endothelial cells to 0–2), removal fractions agree between replicates to
within a few percent per type (stromal 45% / 41%, macrophages 43% / 41%),
and cleanup per replicate equals separate processing. Cleaning admixture
makes replicates more alike, as it should, when every molecule of both
sections is treated by one factorization, one scoring, and one ensemble.

The consensus arm, as implemented, does not achieve this: its factor space
is fine (loss gap −1.3% and −1.2%, 19 of 23 factors matched at median
0.80, member agreement 0.91–0.92), but replicate 2 is under-corrected
relative to replicate 1 by a third in most types (stromal 43% / 28%,
macrophages 44% / 28%, endothelial 42% / 30%), so the replicates diverge
as much as under separate processing. The mechanism is downstream of the
profiles: each sample's ensemble votes with its own aligned restarts, and
each member's rules are derived from that sample alone (the pooled rules
restrict members to pooled pairs, but the member-level rule sets, and so
the votes, differ by sample); replicate 2 also carries four appended
sample-specific factors that are restart noise (own-sample stability
0.04–0.62) and compete for its molecules. The reference projection shares
one member pool across samples and halves the DE count (638) with its
member rules still per sample, which points the same way: uniformity has
to hold at the ensemble-member level, not only for the primary rules.

Decided. (1) The joint factorization is needed and is the route that
meets the aim; the stand-in is the joint route minus engineering (per-
sample neighbor indices, namespaced ids, per-sample stain images, one
store), so the combined-store design in the proposal is next. (2) The
consensus layer stays as the lightweight route over existing per-sample
runs, with two changes before it is compared again: a shared member pool
across samples (aligned restarts of all samples, capped) and member rules
derived from evidence pooled across samples, so every member's vote is
the same decision in every sample; and unmatched factors kept only when
they recur across restarts rather than appended blindly (to be designed
without gating on the stability metric). (3) The per-sample caps under
proportional allocation gave the smaller replicate fewer candidate pairs;
with shared members this should be revisited.

## Loop 8 (2026-11-01): shared ensemble across samples

Orient. The user's correction: scoring and rule derivation are the most
seed-dependent steps of the pipeline, so leaving the ensemble members'
scoring and rules per sample (Loop 7) left the least stable step
un-pooled; pooling the primary rules alone was cosmetic.

Act. The collection layer now shares the whole molecule-vote ensemble
(commit 55ded96, both bindings): one member pool per consensus (the
aligned restarts of all samples, interleaved sample by sample, consensus
as member 0) carried by every sample's fixed-profile run; every member
scored in every sample with the primary method and parameters, its
summaries pooled by the same Stouffer combination, one member rule set
derived once (pooled native check on the member's own cell factor
fractions, restricted to the primary rules' source → target pairs), the
same member rules voting in every sample; rules carry `support` (fraction
of members keeping a pair). The run id hashes profiles plus member pool,
so the earlier projections are not reused. The consensus and reference
arms of pancreas (bridge), pancreas (membrane, `invsqrt_kl`) and breast
are rerunning (23:17); separate and joint arms stand.

Expected if the diagnosis is right: breast consensus replicate
concordance moves from 0.862 / 1069 DE genes toward the joint stand-in's
0.991 / 117, with removal fractions agreeing between replicates per type.
If it does not, the residual divergence is in the factor space or the
per-sample candidate caps, not in the ensemble.

Observed, pancreas membrane (`invsqrt_kl`, shared ensemble;
`results/pancreas_membrane_*.csv`; earlier per-sample-ensemble consensus
in parentheses):

| arm | sample | leakage removed | median pair sensitivity | false removal | removed |
|---|---|---|---|---|---|
| separate | left | 0.73 | 0.97 | 2.3% | 25% |
| separate | right | 0.85 | 0.99 | 3.6% | 33% |
| consensus | left | 0.81 (0.75) | 0.97 (0.96) | 1.9% (2.2%) | 26% |
| consensus | right | 0.78 (0.77) | 0.98 (0.96) | 3.0% (2.8%) | 33% |
| reference (left) | left | 0.78 | 0.98 | 2.2% | 27% |
| reference (left) | right | 0.74 | 0.95 | 3.6% | 33% |

Concordance between halves (cell-weighted per-type pseudobulk correlation;
DE genes between 1 mm block pseudobulks): uncorrected 0.924 / 389,
separate 0.489 / 910, consensus 0.939 / 338 (was 0.74 / 549 with
per-sample ensembles), reference 0.947 / 446. Every type's correlation
except exocrine (0.971 → 0.952) is now above its uncorrected value, and
the DE count falls below uncorrected for six of seven types. Rule
support: 41 of 43 pooled rules kept, median member support 1.0 (lower
quartile 0.9). The halves are different regions, so removal fractions
still differ by type (ductal 28% left vs 87% right, the right half being
exocrine-dominated with ductal cells embedded in acinar tissue; the audit
keeps false removal at 3.0% and the ductal pseudobulk correlation rises
from 0.85 to 0.89, so the removed material is exocrine contamination).
The shared ensemble restored the aim on this variant even though its
factor matching is weak (4 of 9 left factors matched); the vote now
decides uniformly over the same voters in both halves.

Observed, pancreas bridge (`ls_nmf`, shared ensemble;
`results/pancreas_*.csv`; earlier per-sample-ensemble consensus in
parentheses):

| arm | sample | leakage removed | median pair sensitivity | false removal | removed |
|---|---|---|---|---|---|
| separate | left | 0.38 | 0.03 | 0.5% | 6.6% |
| separate | right | 0.51 | 0.38 | 0.3% | 6.4% |
| consensus | left | 0.43 (0.53) | 0.66 (0.54) | 3.2% (1.8%) | 11.9% |
| consensus | right | 0.35 (0.24) | 0.56 (0.46) | 3.4% (1.6%) | 13.2% |
| reference (left) | left | 0.41 | 0.55 | 3.0% | 12.5% |
| reference (left) | right | 0.31 | 0.48 | 1.6% | 11.2% |
| joint stand-in | left | 0.38 | 0.55 | 0.7% | 8.5% |
| joint stand-in | right | 0.17 | 0.15 | 0.2% | 5.7% |

Concordance between halves: uncorrected 0.924 / 389, separate 0.846 /
540, consensus 0.914 / 573 (was 0.861 / 501), reference 0.907 / 548,
joint 0.918 / 387. The correlation recovers to near the uncorrected level
for every type, but the DE count rises through one type: exocrine cells
lose 34% of their molecules in the left half against 17% in the right
(exocrine DE 84 → 238, correlation 0.971 → 0.952), while the joint
stand-in removes 7.5% and 5.3%. The five pooled rules that target
exocrine cells (ductal, fibroblast, endothelial, immune, and endocrine
source factors) all carry member support ≥ 0.8, and two of them rest on
right-half evidence alone (`samples_significant = right`); with the
shared ensemble they are applied in the left half, where the few
exocrine cells sit inside tumor tissue, whereas the earlier per-sample
voters effectively withheld them there. Own-marker false removal rises
to 3.2% and 3.4% with exocrine the worst type on the left. The halves are
different regions, so this is uniform rules acting on unequal contexts
rather than a replicate divergence; the breast replicates decide. Pooled
rules with low member support (0.1 to 0.3: fibroblast ← immune factor 5,
endocrine ← immune, mural ← endothelial) are kept but remove little,
because the vote threshold (3 of 10) is not reached from the primary
member alone; this is the single-sample behavior.

Observed, breast replicates (`ls_nmf`, bridge, shared ensemble;
`results/breast_*.csv`; Loop 7 per-sample-ensemble consensus in
parentheses):

| arm | sample | leakage removed | median pair sensitivity | false removal | removed |
|---|---|---|---|---|---|
| separate | rep1 | 0.62 | 0.86 | 1.6% | 22.2% |
| separate | rep2 | 0.68 | 0.89 | 1.5% | 24.1% |
| consensus | rep1 | 0.61 (0.60) | 0.86 (0.85) | 1.3% (1.4%) | 23.1% (21.6%) |
| consensus | rep2 | 0.67 (0.54) | 0.86 (0.71) | 1.9% (1.3%) | 24.4% (18.0%) |
| reference (rep1) | rep1 | 0.61 | 0.85 | 1.4% | 22.8% |
| reference (rep1) | rep2 | 0.68 | 0.86 | 1.9% | 24.4% |
| joint stand-in | rep1 | 0.62 | 0.88 | 1.2% | 22.7% |
| joint stand-in | rep2 | 0.68 | 0.88 | 1.5% | 24.0% |

Replicate concordance (cell-weighted per-type pseudobulk correlation; DE
genes between the replicates' 1 mm block pseudobulks):

| arm | correlation | DE genes |
|---|---|---|
| uncorrected | 0.988 | 242 |
| separate | 0.902 | 1066 |
| consensus, shared ensemble | 0.992 (0.862) | 141 (1069) |
| reference (rep1), shared ensemble | 0.991 (0.901) | 137 (638) |
| joint stand-in | 0.991 | 117 |

The shared ensemble closes the gap completely. Every type's correlation
is at or above its uncorrected value (stromal 0.990 → 0.995, B cells
0.992 → 0.995), the immune, endothelial, myoepithelial, and tumor types
have 0 or 1 DE genes between replicates, and the residual DE genes are
in the three types where the uncorrected replicates already differed
(DCIS 1 and 2, stromal: 138 of 141). Removal fractions now agree between
replicates to within about two points per type (stromal 45% / 41%,
macrophages 45% / 44%, endothelial 42% / 42%, myoepithelial KRT15+ 7% /
8% against 26% / 6% under separate processing), and cleanup per
replicate equals separate processing. 250 of 298 pooled rules kept,
median member support 1.0, lower quartile 0.8. The second replicate's
under-correction in Loop 7 was entirely the per-sample ensemble.

Decided. (1) The consensus route with the shared ensemble meets the aim
on true replicates as well as the joint stand-in does, so the combined-
store joint factorization is not required for replicate concordance; it
stays specified as a later option (one factorization, no matching step)
rather than scheduled. (2) The consensus layer as committed (55ded96) is
the deliverable of this branch. (3) Open for the user: whether unmatched
factors that do not recur across a sample's restarts should be dropped
rather than appended (breast appended four, three of them restart noise);
whether the per-voter scoring across samples should run samples in
parallel (it doubled the correction time, 93 min on breast); and the
pancreas bridge observation that pooled rules resting on one sample's
evidence are applied in the other sample by design.

## Loop 9 (2026-11-02): pooled-clustering consensus

Act. The shared factors are now built by clustering the pooled factors of
every restart of every sample (average linkage on gene-ownership
correlation, cut at `match_min`, implemented once in C++ as
`cluster_factor_profiles` and exposed in both bindings). A cluster becomes
a shared factor when its components come from at least `min_restarts`
(default 2) distinct restarts of any samples; the rest are reported as
unreproduced. The shared profile is the sample-weighted mean of each
sample's mean unit-sum component profile, rescaled to the weighted mean
row sum. Factors are ordered by sample coverage, then restarts. New
tables: `factors` (coverage, restarts, components, agreement, restarts
contributing more than one component), `recurrence()` (per factor and
sample), `unreproduced`, `components`; `matching` keeps describing each
sample's primary factors (consensus factor, ownership correlation with it,
within-cluster agreement, own-sample stability). The `reference` and
`unmatched` arguments are gone; `fit_profiles()` (reference projection)
fills the same tables from its matching. Documentation rewritten
accordingly (`docs/multisample.md`, binding pages, README).

Expected on the regression arms: the same or better cleanup and replicate
concordance as Loop 8, the breast consensus without restart-noise factors
(the four appended factors of Loop 7 should appear as unreproduced
components or as low-coverage factors with recurrence across restarts),
and the pancreas membrane configuration (weak pairwise matching, 4 of 9)
with higher coverage now that restarts enter the construction.

## Loop 10 (2026-11-02): per-sample diagnostics

Act. Two per-sample tables, both bindings. `consensus$diagnostics()`: loss
gap against restart spread with `representation_flag`, the number and
fraction of the sample's primary factors in the shared space with
`shared_factors_flag` (below one half), within-sample stability against
consensus agreement, and the largest factor-usage divergence of any cell
type with its name (`usage_divergence(by_sample = TRUE)` gives the full
per-type-and-sample table). `correction$diagnostics(report)`: the applied
rules lacking evidence from the sample and the fraction of its removed
molecules attributable to them (a removed molecule is attributed to the
rule of its cell's type and primary factor label, computed from each
cell's per-factor molecule mass before and after the correction, which
reproduces the removed total exactly and never goes negative on the
breast data), plus the audit's leakage removed, median pair sensitivity,
and own-marker false removal with `false_removal_flag` above twice the
median. All flags report; none gates. Documented with the remedies in
order (lower `match_min`, `fit_profiles()` for the sample, joint
factorization). Positive control (a rank-10 profile set projected onto
the rank-23 breast data) pending the end of the regression arms.

Observed, pancreas (pooled clustering, `min_restarts` 2, cut 0.3; Loop 8
pairwise construction in parentheses).

Bridge, `ls_nmf` (9 factors per sample, 20 restarts pooled): 12 shared
factors, 7 of them carried by every restart of both halves (agreement 0.68
to 0.92), four more carried by 3 to 13 restarts (alternative
decompositions found by a minority of restarts), one carried by one half
only; a single unreproduced component. Loss gaps −0.1% and −0.4%
(spreads 0.4% and 1.4%); member agreement 0.92. Cleanup: left 0.41 /
0.59 at 1.1% false removal, 9.5% removed (0.43 / 0.66 at 3.2%, 11.9%);
right 0.17 / 0.53 at 0.4%, 5.8% removed (0.35 / 0.56 at 3.4%, 13.2%).
Concordance 0.920 / 415 (0.914 / 573; uncorrected 0.924 / 389; joint
0.918 / 387). Less is removed, with far fewer own-marker losses, and the
halves end as concordant as under the joint stand-in; the right half's
leakage removal falls to the joint stand-in's level (0.17).

Membrane, `invsqrt_kl`: 15 shared factors, 7 covering both halves (3 of
them in 18 to 20 restarts, the rest in 2 to 10), 8 covering one half only
(2 to 10 restarts), 58 unreproduced components of which 7 are primary
factors with own-sample stability 0.05 to 0.27. Loss gaps +5.7% and
+9.1% against spreads 1.5% and 1.4%: the representation flag fires for
both halves (Loop 8: 1.8% and 3.5%). Cleanup: left 0.74 / 0.95 at 1.6%
(0.81 / 0.97 at 1.9%); right 0.45 / 0.93 at 0.4%, 17% removed (0.78 /
0.98 at 3.0%, 33%). Concordance 0.941 / 364 (0.939 / 338).

Interpretation. For the unstable variant every restart is a different,
equally good local optimum (within-sample stability 0.37 to 0.47, restart
spread 1.4%). Average-linkage clusters at the 0.3 cut have agreement 0.33
to 0.55, and the mean of such loosely related components is a blurred
profile that reconstructs the neighborhoods worse than any single
restart would; the right half, whose primary factors are the least
reproducible, loses half its cleanup. The pairwise construction averaged
fewer, better-matched components (one per sample) and blurred less. The
diagnostics table caught it (gap far above spread), which is what it is
for. Next: compare, by the same loss gap, the cluster mean against the
per-gene median, the medoid component, and an agreement-weighted mean,
at cuts 0.3, 0.5, and 0.7 (`consensus_variants.R`); the pancreas bridge
configuration serves as the control that should not change much.

Observed, profile construction experiment (`consensus_variants.R`; loss gap
per half under the shared profiles, restart spreads 0.4% / 1.4% bridge and
1.5% / 1.4% membrane):

| configuration | cut | factors | mean | median | medoid | agreement-weighted |
|---|---|---|---|---|---|---|
| bridge `ls_nmf` | 0.3 | 12 | −0.1% / −0.4% | −0.5% / +0.3% | 0.0% / +0.3% | −0.2% / −0.2% |
| bridge `ls_nmf` | 0.5 | 14 | −0.8% / −0.7% | −1.1% / −0.1% | −0.3% / −0.2% | −0.8% / −0.5% |
| bridge `ls_nmf` | 0.7 | 21 | −3.1% / −3.3% | −3.6% / −3.4% | −3.7% / −3.2% | −3.2% / −3.3% |
| membrane `invsqrt_kl` | 0.3 | 15 | +5.7% / +9.1% | +4.0% / +4.3% | +21.6% / +19.8% | +3.7% / +6.1% |
| membrane `invsqrt_kl` | 0.5 | 12 | +10.7% / +10.3% | +39% / +35% | +39% / +38% | +10.5% / +11.0% |
| membrane `invsqrt_kl` | 0.7 | 5 | +98% / +270% | | | |

Each half under the other half's primary profiles: bridge +2.2% / +2.8%;
membrane, left under right's +8.7%, right under left's −2.6% (the left
decomposition fits the right half's neighborhoods better than the right's
own selected restart does).

Interpretation. For the stable variant every construction is within the
restart spread and the choice does not matter; more factors (a higher
cut) lower the loss by over-completeness, not by better programs. For
the unstable variant no construction from pooled clusters reaches the
spread: the median and the agreement-weighted mean blur less than the
mean, the medoid is worst because factors drawn from different restarts
do not form one coherent decomposition, and a tighter cut discards most
components as unreproduced (114 of 180 at 0.5) and leaves too few
factors. What does represent both halves is one coherent decomposition:
the left half's own profiles (gap 0 and −2.6%), which is the reference
projection; its cleanup in Loop 8 (0.78 / 0.98 and 0.74 / 0.95) was also
close to the pairwise construction's and well above the pooled mean's
(0.74 / 0.95 and 0.45 / 0.93). Conclusion: averaging components is sound
when the restarts agree (stable variants, the breast and bridge cases)
and unsound when they do not; in the latter case a single decomposition
chosen for how well it represents every sample, or a joint refinement of
the profiles on the pooled training rows, is needed. The diagnostics
table identifies the case (gap above spread; within-sample stability low).

Observed, positive control for the diagnostics (`diagnostics_check.R`): the
rank-10 primary profiles of replicate 1 projected onto both rank-23 breast
replicates give loss gaps +8.4% and +9.0% against spreads 0.5% (the
representation flag fires for both), 10 of 23 primary factors in the
shared space (shared fraction 0.43, the shared-factor flag fires), usage
divergence at most 0.05. Both flags fire when they should.

Observed, breast under pooled clustering at cut 0.3 (projection stage of
the Loop 9 arm): 24 shared factors, 20 covering both replicates, 4 from
one replicate only with agreement 0.95 to 0.99; no unmatched primary
factor; but loss gaps +1.5% and +1.2% against spreads 0.5%, so the
representation flag fires, where the pairwise construction had −1.3% and
−1.2%. The cause is visible in the factor table: several clusters gather
up to 36 components from 20 restarts (15 restarts contributing two or
more), with agreement 0.52 to 0.64; distinct programs that a restart
finds side by side (correlated at 0.3 to 0.5) were merged, and their mean
is a blur. The pairwise Hungarian construction cannot do this because it
matches one to one.

Act. The clustering now carries the one-to-one structure as a constraint:
two clusters that both hold a factor of the same restart never merge
(C++ `cluster_factor_profiles`, tested), so a shared factor is the
one-to-one correspondence of a program across restarts and samples, and
each restart contributes at most one component to it. The purity column
is gone since it is always zero. The construction experiment and the
breast loss check are rerun with the constraint; the running breast arm
(unconstrained) is left to finish as the test of whether a fired
representation flag predicts worse cleanup or concordance on true
replicates.

Observed, construction experiment with the one-to-one constraint (loss gap
per half; spreads 0.4% / 1.4% bridge, 1.5% / 1.4% membrane):

| configuration | cut | factors | mean | median | medoid | agreement-weighted |
|---|---|---|---|---|---|---|
| bridge `ls_nmf` | 0.3 | 15 | −1.9% / −1.8% | −1.9% / −1.3% | −1.8% / −1.2% | −2.0% / −1.2% |
| membrane `invsqrt_kl` | 0.3 | 16 | +4.0% / +8.8% | −0.4% / +1.6% | +16% / +17% | +2.1% / +5.8% |

Cuts 0.5 and 0.7 are unchanged (the constraint does not bind there) and
remain unusable for the unstable variant. With one component per restart
in each cluster, the per-gene median of the unit-sum components
represents both halves of the membrane configuration within or at their
restart spread, where the mean does not; on the stable bridge
configuration the choice is immaterial. The medoid stays worst: factors
drawn from different restarts do not form one coherent decomposition.

Decided. The shared profile of a cluster is the weighted per-gene median
of its unit-sum components (each sample's weight divided among its
components in the cluster), renormalized and rescaled as before; both
bindings, documented. The three consensus arms are rerun with the
constrained clustering and the median construction.
