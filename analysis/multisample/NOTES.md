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

Observed, breast shared-factor structure with the one-to-one constraint
(scratch script on the cached fits): 37 shared factors from rank-23 fits,
11 carried by all 20 restarts, then 19, 19, 18, 18, 17, 17, 16, 15, 13,
12, 10, 8, 8, 5, 4, 3 covering both replicates, and 7, 6, 4, 4, 3, 2, 2,
2, 2, 2 from one replicate only; 4 unreproduced components. Loss gaps
−4.3% and −4.1%, that is an over-complete basis. Only 8 pairs of shared
factors have ownership correlation at or above 0.3, and 6 of those
co-occur in the same restarts (distinct programs found side by side), so
duplication is not the issue; the per-restart factor sets simply vary:
every restart finds the 11 core programs plus a varying dozen from a pool
of about 26 less recurrent ones. Removing the two true alternatives
changes little (35 factors, −3.8% and −3.6%).

Also observed: the first weighted-median implementation took the lower
middle value for even component counts, which for two-component clusters
is the per-gene minimum of the pair; the membrane median arm projected
with loss gaps of +31% and +26% where the experiment's ordinary median
gave −0.4% and +1.6%. Fixed (midpoint of the lower and upper weighted
medians, equal to the ordinary median for equal weights; tested in both
bindings); the three median arms were stopped and are relaunched.

Act. A per-sample recurrence threshold, `min_recurrence` (default 0,
off), requires a program in at least that fraction of some sample's
restarts; at 0.5 only programs that are the majority solution of at least
one sample enter, which on breast would give about 24 shared factors for
rank 23. The arms run both ways (`min_recurrence` 0 and 0.5) so that the
choice of default rests on cleanup and concordance, not on the loss gap
alone (which favors over-completeness).

Observed, pancreas with the corrected median (component weights balanced
per sample, the first weighting) at recurrence thresholds 0 (`median`)
and 0.5 (`majority`); `results/pancreas*_arms_per_sample.csv`:

| configuration | arm | factors | loss gap | left leak / sens / false | right leak / sens / false | concordance |
|---|---|---|---|---|---|---|
| membrane | separate | 9 + 9 | | 0.73 / 0.97 / 2.3% | 0.85 / 0.99 / 3.6% | 0.489 / 910 |
| membrane | pairwise (Loop 8) | 14 | +1.8% / +3.5% | 0.81 / 0.97 / 1.9% | 0.78 / 0.98 / 3.0% | 0.939 / 338 |
| membrane | pooled mean, unconstrained | 15 | +5.7% / +9.1% | 0.74 / 0.95 / 1.6% | 0.45 / 0.93 / 0.4% | 0.941 / 364 |
| membrane | median, recurrence 0 | 16 | +5.9% / +3.0% | 0.81 / 0.96 / 1.3% | 0.80 / 0.98 / 0.8% | 0.947 / 328 |
| membrane | median, recurrence 0.5 | 9 | +72% / +63% | 0.82 / 0.99 / 3.8% | 0.90 / 0.99 / 3.8% | 0.940 / 354 |
| bridge | separate | 9 + 9 | | 0.30 / 0.03 / 0.5% | 0.20 / 0.38 / 0.3% | 0.846 / 540 |
| bridge | pairwise (Loop 8) | 10 | +0.7% / +0.6% | 0.43 / 0.66 / 3.2% | 0.35 / 0.56 / 3.4% | 0.914 / 573 |
| bridge | pooled mean, unconstrained | 12 | −0.1% / −0.4% | 0.41 / 0.59 / 1.1% | 0.17 / 0.53 / 0.4% | 0.920 / 415 |
| bridge | median, recurrence 0 | 15 | −1.4% / −2.3% | 0.56 / 0.60 / 4.4% | 0.36 / 0.55 / 1.8% | 0.915 / 470 |
| bridge | median, recurrence 0.5 | 9 | +1.9% / +0.9% | 0.38 / 0.48 / 0.9% | 0.18 / 0.40 / 0.2% | 0.918 / 428 |
| bridge | joint stand-in | 9 | | 0.38 / 0.55 / 0.7% | 0.17 / 0.15 / 0.2% | 0.918 / 387 |

Interpretation. (1) The median construction repairs the membrane case:
both halves clean as well as or better than separately, at lower false
removal, with the best concordance of any arm (0.947 / 328). (2) The
training-loss gap does not predict cleanup. The nine majority programs
reconstruct the membrane halves' neighborhoods 60 to 70% worse than
their own fits, yet clean them best of all arms (0.82 and 0.90 leakage
removed at sensitivity 0.99), at the price of 3.8% own-marker false
removal, the level of separate processing for this configuration. The
gap measures how well fixed profiles reconstruct neighborhood
compositions, which an unstable variant's idiosyncratic per-restart
programs do and reproducible programs need not; cleanup depends on
whether the admixture programs are present and sharp. The representation
flag therefore reports reconstruction, and the per-sample audit remains
the arbiter of cleanup; the documentation will say so. (3) The recurrence
threshold trades thoroughness for specificity in the stable bridge
configuration: all programs (15 factors) remove the most leakage with the
most own-marker loss (4.4% on the left), the majority programs (9) the
least of both. Concordance between the halves is nearly the same for
every shared construction (0.915 to 0.920 against 0.924 uncorrected).

Also observed: the balanced weighting (each sample's weight divided among
its components in a cluster) lets a lone component of one sample carry
half the median when the other sample contributes many; the package now
weights each component by its sample's weight alone, which reduces to the
ordinary median under equal allocation and matches the construction the
loss experiment validated (−0.4% / +1.6% on membrane). The pancreas median
and majority arms are rerun with it (the balanced results are kept as
`*_balanced`); the running breast arms still use the balanced weighting.

Correction to the loss-gap record. The training loss under fixed profiles
for the KL variants depended on the overall scale of the profiles: the
row projection returns usage fractions, and the loss used them as usages,
so halving every profile lowered the loss by a quarter. Every membrane
(`invsqrt_kl`) loss gap above, including the construction experiments
and the Loop 8 and Loop 9 tables, is affected; the `ls_nmf` gaps are not
(the least-squares projection refits the usages). Fixed in commit
f238e55: the usages are rescaled to the observed row totals, the
projection starts at the observed scale, and a unit test holds the loss
invariant to the profile scale. Recomputed membrane gaps (spreads 1.5% /
1.4%): pairwise Loop 8 consensus (14 factors) −0.1% / 0.0%; pooled mean,
unconstrained (15) +4.4% / +5.5%; pooled median (16) +18% / +17%;
majority (9) +128% / +122%; reference projection of the left half 0 /
+13.8%; each half under the other's primary profiles +20.6% / +13.8%.
Bridge gaps are unchanged: pooled median (15) −1.5% / −2.1%, majority (9)
+1.8% / +0.9%, pairwise +0.7% / +0.6%.

Interpretation, revised. Averaging across a sample's own restarts blurs
when the restarts disagree; the per-gene median does not fix this (it
sparsifies instead, reconstructing worse than the mean). The pairwise
construction averaged one component per sample, the primary one, and
reconstructs both membrane halves as well as their own fits. The
cleanup rankings are nearly flat across constructions apart from the
pooled mean's failure on the right half, so representation is the
criterion left to decide between them, together with scalability. Next
test: the cluster structure of all restarts for correspondence,
recurrence, and coverage, with the shared profile of a cluster formed
from each sample's primary component only.

Observed, primary-component construction (scratch): with the pooled
clusters as correspondence but each cluster's profile averaged from the
samples' primary components only, membrane gives 8 factors with gaps
+65% / +76% (the unstable variant's primary factors mostly fall into
unreproduced clusters and are lost), bridge 10 factors at +0.1% / +0.2%,
breast 27 factors at −1.3% / −1.2% (the pairwise numbers). Not pursued:
it reproduces the pairwise construction on stable variants and fails on
unstable ones. For unstable variants no construction of reproducible
programs reconstructs the samples' neighborhoods as their own idiosyncratic
solutions do, and the cleanup evidence says it does not need to.

Decided (pending the breast arms). Keep the pooled clustering with the
one-to-one constraint and the per-gene median; make the majority rule
(`min_recurrence` 0.5, programs that are the majority solution of at least
one sample) the default if the breast arms confirm its concordance and
cleanup, since it gives rank-sized shared spaces with the plainest
meaning; keep the representation flag as information.

Observed, breast under the unconstrained pooled mean (the Loop 9 arm left
to finish as the test of the representation flag; 24 factors, gaps +1.5%
and +1.2% against spreads 0.5%, flag fired for both replicates; `ls_nmf`,
so the gap is trustworthy): leakage removed 0.64 / 0.70 (separate 0.62 /
0.68; pairwise Loop 8 0.61 / 0.67), median pair sensitivity 0.86 / 0.87,
own-marker false removal 1.9% / 2.6% (Loop 8 1.3% / 1.9%), removal 25% /
26%, 209 rules kept (177 supported by both replicates); concordance
0.992 / 145 DE genes (Loop 8 0.992 / 141; joint 0.991 / 117; uncorrected
0.988 / 242). A fired representation flag therefore did not predict
worse cleanup or concordance on true replicates; its only trace is the
somewhat higher own-marker false removal, consistent with blurred merged
profiles labeling a few native molecules as admixture. The flag stays
informational, as documented.

Observed, breast replicates under the corrected median construction
(`results/breast_*.csv`; bridge, `ls_nmf`):

| arm | factors | rep1 leak / sens / false | rep2 leak / sens / false | rules (both-supported) | concordance |
|---|---|---|---|---|---|
| separate | 23 + 23 | 0.62 / 0.86 / 1.6% | 0.68 / 0.89 / 1.5% | 213, 235 | 0.902 / 1066 |
| pairwise (Loop 8) | 27 | 0.61 / 0.86 / 1.3% | 0.67 / 0.86 / 1.9% | 250 (217) | 0.992 / 141 |
| pooled mean, unconstrained | 24 | 0.64 / 0.86 / 1.9% | 0.70 / 0.87 / 2.6% | 209 (177) | 0.992 / 145 |
| median, recurrence 0 | 37 | 0.63 / 0.89 / 1.4% | 0.70 / 0.89 / 2.0% | 326 (268) | 0.992 / 138 |
| median, recurrence 0.5 | 26 | 0.63 / 0.88 / 1.1% | 0.69 / 0.89 / 1.7% | 270 (233) | 0.991 / 128 |
| joint stand-in | 23 | 0.62 / 0.88 / 1.2% | 0.68 / 0.88 / 1.5% | 260 | 0.991 / 117 |

Every shared construction holds replicate concordance at or above
uncorrected (0.988 / 242) with cleanup equal to separate processing; the
differences between them are at the level of run-to-run variation. The
majority rule (recurrence 0.5) has the lowest own-marker false removal
of all arms, the highest sensitivity among the consensus arms, a
rank-sized shared space (26 factors for rank 23; loss gaps −0.8% and
−0.4%), and the fewest DE genes between replicates of any consensus arm
(128; the joint stand-in 117). The residual DE genes sit in the three
types where the uncorrected replicates already differ (DCIS 1 and 2,
stromal).

Decided. `min_recurrence = 0.5` is the default in both bindings: the
shared factors are the programs that are the majority solution of at
least one sample. Loops 9 and 10 are closed. Loop 11 (eight-tile cohort)
is running; Loop 12's headline run (breast replicates, membrane scoring,
`invsqrt_kl`) is launched behind it.

## Loop 11 (2026-11-02): eight-tile cohort

Setup. Each breast replicate cut into four quantile tiles (41,945 cells
per replicate-1 tile, 29,688 per replicate-2 tile; 6.2 to 11.0 million
molecules), replicate-1 tile 3 and replicate-2 tile 2 thinned gene by
gene at a log-efficiency SD of 0.4; separate per-tile fits (`ls_nmf`, 10
restarts each, 80 restarts pooled), bridge scoring.

Observed, factor space (`cohort_factor_space.R`, tile restart spreads
0.2% to 0.8%):

| construction | shared factors | duplicate pairs (ownership r ≥ 0.3) | max pair r | in all 8 tiles | in one tile only | tiles with loss gap above spread |
|---|---|---|---|---|---|---|
| star (pairwise to the largest tile) | 57 | 60 | 0.82 | 11 | 34 | 1 of 8 (+11%) |
| pooled clustering, recurrence 0 | 61 | 30 | 0.72 | 19 | 6 | 0 of 8 (gaps −10%) |
| pooled clustering, majority rule (the arm) | 37 | 12 | 0.66 | 19 | 2 | 0 of 8 (gaps −3.2% to −4.1%) |

The star construction behaves as predicted for many samples: programs
the reference tile lacks are appended once per tile that has them (34
single-tile factors, 60 duplicate pairs), and the reference tile's own
space fits one tile 11% worse than its own restarts. The pooled
clustering finds 19 programs present in every tile; with the majority
rule the shared space has 37 factors for rank-23 fits, 26 of them in at
least 7 of the 8 tiles, 127 unreproduced components, and reconstructs
every tile within 3 to 4% below its own fit. Its pooled rules: 246 kept,
53 supported by all eight tiles' own evidence, 170 by at least four, 23 by
one tile only (the per-tile candidate budget is one eighth of a
single-sample budget under proportional allocation, so single-tile
support is expected to be sparse). The cohort consensus arm (scoring,
shared ensemble over 8 tiles, per-tile audits) and the grouped replicate
comparison are running.

## Loop 12 (2026-11-02): breast replicates with membrane scoring

Observed (`invsqrt_kl`, rank 23, 10 restarts per replicate, membrane
scoring; `results/breast_membrane_*.csv`):

| arm | factors | rep1 leak / sens / false | rep2 leak / sens / false | removed | rules | concordance |
|---|---|---|---|---|---|---|
| separate | 23 + 23 | 0.60 / 0.88 / 3.0% | 0.66 / 0.85 / 2.1% | 32% / 29% | 188, 191 (Jaccard 0.48) | 0.642 / 1242 |
| consensus, majority rule | 12 | 0.51 / 0.59 / 0.7% | 0.57 / 0.66 / 0.8% | 21% / 22% | 115 (92 both) | 0.992 / 129 |

Separate membrane processing diverges the replicates far more than
separate bridge processing did (correlation 0.988 → 0.642; 1242 DE
genes); the consensus restores concordance to the bridge level (0.992 /
129) but cleans less than separate processing (leakage removed 0.51 /
0.57 against 0.60 / 0.66; sensitivity 0.59 / 0.66 against 0.88 / 0.85),
with very low own-marker false removal (0.7% / 0.8%).

Diagnosis from the tables. The variant is unstable on this tissue:
within-sample restart stability 0.38 for both replicates. The majority
rule keeps 12 shared programs (6 in every restart, 6 in 50 to 70%); 277
components are unreproduced, among them 14 of replicate 1's and 12 of
replicate 2's 23 primary factors, all with own-sample stability 0.13 to
0.26. The shared-factor flag of `diagnostics()` fires for both replicates
(9 and 11 of 23 primary factors in the shared space), and the loss gap is
+339% / +330%. Admixture programs among the unstable primaries are lost
to the shared space, hence the lower sensitivity. The pancreas membrane
case (rank 9, stability 0.37 to 0.47) kept 9 of its programs under the
same rule and cleaned best; at rank 23 the majority rule bites harder.

Act. The same configuration runs with `min_recurrence = 0` (every
program recurring in two restarts) to see whether the lost programs come
back with their cleanup, which would make "lower `min_recurrence`" the
documented first remedy when the shared-factor flag fires.

Observed, cohort cleanup and concordance (`results/breast_tiles_*.csv`;
bridge, `ls_nmf`, majority-rule consensus with 37 shared factors):

| arm | leakage removed, mean (range) | median pair sensitivity | own-marker false removal, mean (range) | removed | rules per tile | concordance rep1 vs rep2 (pooled tiles) |
|---|---|---|---|---|---|---|
| uncorrected | | | | | | 0.987 / 272 |
| separate per tile | 0.64 (0.56 to 0.71) | 0.83 | 1.5% (1.0 to 2.4%) | 21% | 95 to 135 (pairwise Jaccard 0.45) | 0.904 / 688 |
| consensus, majority rule | 0.69 (0.55 to 0.77) | 0.91 | 4.1% (2.3 to 7.5%) | 26% | 246 (53 supported by all tiles, 23 by one) | 0.991 / 153 |

Concordance between the replicates is restored (0.991 / 153, against
0.987 / 272 uncorrected and 0.904 / 688 separate), immune, endothelial,
and tumor types have 0 to 2 DE genes left, and sensitivity rises, but
own-marker false removal is 2.7 times that of separate processing, above
the 2% acceptance level in every tile, and the two-replicate consensus
(1.1% / 1.7%) did not show this. The thinned tiles are not the cause
(rep1_t3 3.1%, rep2_t2 6.0%, but the unthinned rep2_t4 7.5%).

Attribution (removed molecules per tile credited to the rule of their
cell's type and primary factor): rules lacking the tile's own evidence
account for 1.5 to 11% of a tile's removal, so pooled power is not the
driver; rules on factors carried by fewer than 7 of the 8 tiles account
for 12 to 30%. Eleven of the 37 shared factors are present in 1 to 6
tiles only (regional programs: DCIS present in some tiles, a program of
one replicate's tiles). A rule for such a program applies in tiles that
do not carry it, where the molecules labeled with it are nearest-factor
artifacts rather than that program, and they are removed. Per type the
cohort removes more than both the two-replicate consensus and separate
processing from B cells (59% against 48% and 40%), stromal (51 against
44 and 41), CD8 T cells (44 against 33 and 31), and KRT15+ myoepithelial
cells (22% against 6% and 10%).

Also observed: the eight-tile correction took 68 minutes, most of it the
eight per-tile votes applied one tile after another; the voter scoring
itself (72 jobs) ran concurrently.

Interpretation. With two replicate sections nearly every shared program
is present in both, so uniform rules act on equal contexts; with eight
regional tiles a third of the shared programs are regional, and the
uniform application of their rules in tiles lacking the program is where
the extra own-marker removal comes from. Two levers, both decisions for
the user: apply a rule in a sample only where the sample's own restarts
carry the rule's program (a presence condition on the data, not an
evidence threshold; it leaves true replicates untouched because their
programs are shared), or require a minimum sample coverage for a program
to enter the shared space (which drops regional programs altogether).
The per-sample `diagnostics()` tables already show coverage per factor
and the removed fraction attributable to rules without local evidence;
a column for the removed fraction attributable to programs the sample
does not carry would make this visible directly.

Observed, with the new correction diagnostic computed on the cohort: each
tile's own restarts carry 28 to 32 of the 37 shared programs; the rules on
the programs a tile does not carry number 5 to 46 per tile and account for
0.2% to 18% of its removed molecules (rep2_t3 18%, rep1_t4 14%, rep2_t4
11%); across tiles this fraction correlates only weakly with own-marker
false removal (r = 0.24). So rules on absent programs explain part of the
extra removal in some tiles but not the general rise; the rest comes with
the larger rule set itself (246 pooled rules against 95 to 135 per
separate tile; B cells, stromal, and T cells lose 10 to 20 points more of
their molecules than under separate processing). Whether a cohort of
regional tiles should be cleaned this much more thoroughly, at 4% own-
marker loss, is a question for the user; the diagnostics now make both
components visible per sample.

Observed, breast membrane with `min_recurrence = 0` (`consensus_all`):
71 shared factors from rank-23 fits (134 unreproduced components, 12 of
them primary), loss gaps −8.9% / −8.7%, 421 rules kept (276 supported by
both replicates), leakage removed 0.59 / 0.64, sensitivity 0.77 / 0.86,
own-marker false removal 4.7% / 5.5%, removal 31% / 32%, concordance
0.991 / 142. Against the majority rule (12 factors: 0.51 / 0.57, 0.59 /
0.66, 0.7% / 0.8%, 0.992 / 129) and separate processing (0.60 / 0.66,
0.88 / 0.85, 3.0% / 2.1%, 0.642 / 1242).

Interpretation. Lowering the recurrence threshold recovers the lost
sensitivity on the unstable variant at the price of a tripled shared
space and own-marker removal above separate processing; no threshold
reproduces separate processing's balance because the variant's programs
are not reproducible across restarts. The same replicates under `ls_nmf`
with bridge scoring and the default consensus reach 0.63 / 0.69, 0.88 /
0.89, 1.1% / 1.7%, 0.991 / 128: better than every membrane arm on every
measure except separate processing's removal fraction, which it pays for
with lost concordance. For multi-sample work the stable variant is the
one to use; documented as the remedy when the shared-factor flag fires.

## Loop 13 (2026-11-02): concurrent per-sample correction steps

Observe. Stage timing of the eight-tile collection correction on the
cached consensus run (8 threads per dataset, thread budget 64, votes one
tile after another): primary scoring 10.4 min; member scoring of the 72
(sample, member) jobs with 8 concurrent workers 18.1 min; pooled member
rules with the native-factor check 4.2 min; the eight votes 0.6 to 1.1
min each, 6.6 min in all. The run removed exactly the molecules the Loop
11 run removed in every tile, so member scoring is deterministic.
Correction to Loop 11: the votes are not where the 68 minutes went. That
correction ran at a load near 180 on 64 cores alongside the eight
separate-tile arms, and member scoring is the dominant stage at any load.

Correction to Loop 11 and to the open decisions: the per-sample candidate
budget is not one eighth of a single-sample budget. The allocation shares
sum to the number of samples, so every sample receives a full
single-sample budget scaled by its relative size (rep1_t1 received 1122
candidate pairs per type pair against the single-sample 1000). Single-tile
support is sparse because a tile holds fewer cells, and fewer cells of
each type pair, than a section, not because of the cap.

Act. The three per-sample steps of the correction run as independent jobs
in forked worker processes within one thread budget: labeling each
sample's molecules with every member (the dataset's threads per job),
scoring every member in every sample (as before), and the vote in each
sample (one thread per job). The number of concurrent jobs is the budget
divided by the threads one job uses, capped at the job count; the budget
is `num_threads` on `correct()`, by default `CELLADMIX_NUM_THREADS` or the
available cores (the Python dataset default now honors the same variable,
as R did). One job runner per binding (`.celladmix_run_jobs`, `_run_jobs`)
raises a worker's error in the parent; tests cover the worker count, the
environment default, and error propagation.

Result. With the member summaries cached, the eight-tile correction takes
5.7 min (pooled member rules 4.2, the eight votes concurrently about 1.1)
against 11 min with sequential votes, and removes exactly the same
molecules in every tile (rep1_t1 1,690,200 through rep2_t4 2,619,253). A
cold correction is dominated by member scoring, which was already
concurrent; its workers each run one scoring job with the dataset's
threads, so a dataset created with fewer threads admits more concurrent
jobs.

## Loop 14 (2026-11-02): rules confined to the samples carrying the program

Decision. The user chose to make the presence condition an option
(`apply_rules = "where_present"` on `correct()`) and to keep uniform
application (`"everywhere"`) as the default.

Act. A rule is applied in a sample only when at least one of the sample's
own restarts carries the rule's program (recurrence above zero in the
consensus components); primary and member rules are filtered alike, so the
vote threshold is unchanged. The returned rules carry
`samples_without_program`; the correction records `params$apply_rules`,
and `diagnostics()` counts the rules applied in each sample and the rules
on programs it does not carry (applied under the default, skipped under
the option). An empty per-sample rule set writes the sample's run
unchanged. Both bindings, tests on the row filter and the end-to-end path,
docs, skill. The analysis driver names the correction after the mode so
the two arms do not overwrite one another's corrected runs (the first run
did, and the comparison tables were regenerated after the uniform arm's
corrected runs were restored; they came back with the identical removal
counts).

Observe (eight tiles, bridge, `ls_nmf`, majority consensus; the audit per
tile and the concordance between the pooled replicates):

| arm | leakage removed | pair sensitivity | own-marker false removal | removed fraction | concordance (r / DE genes) |
|---|---|---|---|---|---|
| separate | 0.64 | 0.83 | 1.5% | 0.21 | 0.904 / 688 |
| consensus, rules everywhere | 0.69 | 0.91 | 4.1% | 0.26 | 0.991 / 153 |
| consensus, rules where present | 0.65 | 0.83 | 2.7% | 0.22 | 0.938 / 362 |
| uncorrected | | | | | 0.987 / 272 |

Per tile, confining the rules lowered false removal everywhere (largest
drops rep2_t2 6.0% to 2.0%, rep1_t4 4.1% to 1.4%, rep2_t3 3.2% to 1.1%)
but also the leakage removed (rep2_t3 0.70 to 0.56, rep1_t4 0.65 to 0.57)
and the pair sensitivity (rep2_t3 0.91 to 0.65, rep1_t4 0.91 to 0.75):
programs that no restart of those tiles found still carried real
admixture there. The concordance between the replicates fell below the
uncorrected value because corresponding regions of the two replicates were
no longer cleaned alike; which programs a tile's restarts find depends on
the tile's composition and on the restarts themselves.

Interpretation. The presence condition does what it says, removing less
native expression, but on this cohort it costs sensitivity and undoes the
uniformity that the shared rules provide, which is the purpose of the
multi-sample analysis. Uniform application stays the default; the option
is documented, with these numbers, for cohorts whose samples are not
replicates and whose composition differences are real. My earlier
recommendation of the presence condition as the default was wrong on the
evidence.

## Closing summary, Loops 9 to 14 (2026-11-02)

Delivered on `dev_multisample` (both bindings, tests green, docs
current): consensus by pooled clustering of all restarts with the
one-to-one constraint, the per-gene weighted median, the majority
recurrence rule (default 0.5), the shared molecule-vote ensemble with
concurrent per-sample labeling, member scoring, and votes, per-sample
diagnostics on the fit and the correction, a scale-invariant training loss for the KL variants, the
agent skills, and the cohort analysis scripts.

Evidence. Two breast replicates: concordance 0.991 / 128 DE genes (0.988 /
242 uncorrected, 0.902 / 1066 separate), cleanup equal to separate
processing with the lowest false removal of any arm. Eight regional
tiles: concordance between replicates restored (0.991 / 153 against
0.904 / 688 separate), sensitivity up, own-marker false removal 4.1%
against 1.5% separate, part of it from rules on programs a tile does not
carry (0.2 to 18% of a tile's removal, now reported per sample). The star
construction duplicates programs on cohorts (60 duplicate pairs, 34
single-tile factors) where the pooled clustering gives 12 and 2.
Unstable variant (`invsqrt_kl`): the consensus restores concordance but
cannot match separate processing's cleanup balance at any recurrence
threshold; the stable variant should be used.

Resolved: (1) the presence condition is an option
(`apply_rules = "where_present"`), with uniform application the default
after the Loop 14 cohort test showed the option lowers false removal but
costs sensitivity and replicate concordance; a minimum sample coverage
for shared programs was not pursued; (2) the per-sample steps of the
correction run concurrently (Loop 13); (3) the candidate budget already
scales with the sample count, one single-sample budget per sample scaled
by relative size (Loop 13). The pull request to `dev` waits for the user.

## Loop 15 (2026-11-11): the stability gate, and what restart instability is

Trigger. A colleague's six-zone run produced two shared factors; the
zones' own fits had few stable factors after a coarser annotation (under
investigation on their side, rank fixed to the fine annotation's value).
The user set the terms: recurrence is across samples, stability is across
restarts within a sample; the consensus code had used "recurrence" for the
restart share and "coverage" for the sample share.

Observed, what the majority rule (restart share 0.5) drops, from the
pooled clustering at threshold 0 (`pooled_factors.csv` of Loop 11 and the
replicates' pools):

| data | factors at 0 | dropped at 0.5 | of those in 2+ samples | in 3+ | in 4 to 6 |
|---|---|---|---|---|---|
| breast replicates | 37 | 11 | 3 | | |
| eight tiles | 61 | 24 | 20 | 12 | 6 |

One dropped program occurs in 6 of 8 tiles (10 restarts in all, at most
2 per tile); the rule keeps 2 single-tile programs. The dropped programs
have lower within-cluster agreement (median 0.58 against 0.76). The
cleanup differences between the thresholds on stable fits (Loop 10) were
within run-to-run variation; the sensitivity cost showed on the pancreas
halves (0.60 / 0.55 → 0.48 / 0.40) and the unstable breast fits (Loop
12). More restarts do not raise a restart share; they only help the
two-restart minimum.

Observed, fixed-rank weighted NMF of the pooled unit-sum restart profiles
(rows weighted by sample weight over the sample's restarts; fitted on the
odd restarts, even restarts projected by nonnegative least squares):

| variant | construction | profiles | held-out unexplained |
|---|---|---|---|
| `ls_nmf` | pooled NMF rank 23 / 37 | 23 / 37 | 0.17 / 0.07 |
| `ls_nmf` | cluster medians, 0.5 / 0 | 26 / 37 | 0.13 / 0.07 |
| `invsqrt_kl` | pooled NMF rank 23 / 45 | 23 / 45 | 0.59 / 0.45 |
| `invsqrt_kl` | cluster medians, 0.5 / 0 | 12 / 71 | 0.72 / 0.36 |

On the stable variant the factorization rediscovers the cluster programs
(each profile matches one at ownership r 0.81 to 0.97; primaries load on
one profile, median share 0.91) and reconstructs held-out restarts as the
medians do at equal size; the held-out curve has no elbow at any rank. On
the unstable variant no common low-rank cone exists, and the fixed-rank
profiles are mixtures (matches 0.25 to 0.91, half the primaries spread
over several profiles), the blurring of the pooled mean (Loop 9).

Observed, restart structure (`restart_structure.R`: for each pair of
restarts, the share of a factor's unit-sum profile that nonnegative least
squares on the other restart's factors cannot represent; "noise factors"
are primary factors above 0.3 in the median over the other restarts):

| fit | objective spread | one-to-one stability | unrepresentable, mass-weighted | noise factors (mass) |
|---|---|---|---|---|
| breast `ls_nmf` rank 23 | 0.5% / 0.5% | 0.82 / 0.79 | 2% / 2% | 3 / 4 of 23 (1% / 1%) |
| breast `invsqrt_kl` rank 23 | 1.7% / 1.1% | 0.31 / 0.31 | 24% / 23% | 15 / 17 of 23 (24% / 27%) |
| breast `invsqrt_kl` rank 23, 600 iterations | 2.0% | 0.32 | 25% | 16 of 23 (29%) |
| breast `invsqrt_kl` rank 8 | 1.6% / 2.1% | 0.47 / 0.52 | 8% / 6% | 4 / 2 of 8 (14% / 5%) |
| pancreas halves `invsqrt_kl` rank 9 | 1.5% / 1.4% | 0.41 / 0.29 | 8% / 9% | 4 / 3 of 9 (10% / 6%) |

On the pancreas right half 8 of 9 primary factors have stability below
0.5 while 92% of the molecules sit in factors representable in every
other restart at under 10% residual (INS/GCG at 1%): re-partitionings of
the same programs. The noise factors are small (0.5 to 3% of molecules
each), nearly orthogonal across restarts (85 to 95% unrepresentable), and
their top genes are not programs (THAP2/EIF4EBP1/LPL/DMKN,
MAP3K8/FSTL3/EDNRB/RAPGEF3, LARS/JUP/CLDN4/KLF5).

Observed, rank ladder (`rank_ladder.R`; `results/breast_rank_ladder.csv`;
10 restarts, 150 iterations):

| variant | rank | stability mean (rep1 / rep2) | stable factors | unrepresentable, mass-weighted | noise mass |
|---|---|---|---|---|---|
| `ls_nmf` | 8 | 0.94 / 0.89 | 8 / 7 | 0% / 2% | 0% / 2% |
| `ls_nmf` | 12 | 0.80 / 0.85 | 11 / 11 | 2% / 2% | 1% / 2% |
| `ls_nmf` | 16 | 0.80 / 0.80 | 15 / 14 | 3% / 6% | 0% / 8% |
| `ls_nmf` | 23 | 0.82 / 0.79 | 21 / 20 | 2% / 2% | 1% / 1% |
| `invsqrt_kl` | 5 | 0.53 / 0.58 | 2 / 3 | 4% / 4% | 4% / 0% |
| `invsqrt_kl` | 8 | 0.47 / 0.52 | 4 / 5 | 8% / 6% | 14% / 5% |
| `invsqrt_kl` | 12 | 0.42 / 0.40 | 5 / 4 | 11% / 11% | 10% / 9% |
| `invsqrt_kl` | 16 | 0.36 / 0.34 | 6 / 4 | 16% / 17% | 18% / 16% |
| `invsqrt_kl` | 23 | 0.31 / 0.31 | 6 / 6 | 24% / 23% | 24% / 27% |

Interpretation. (1) The majority rule was a stability gate, against the
standing decision that stability is report-only, and it discards programs
that recur across samples while keeping single-sample ones. (2) A
fixed-rank factorization of the pooled factor space is not a robustness
step: it reproduces the clustering where the fits are stable and blurs
where they are not, and no held-out criterion picks its rank. (3) Two
things hide under low stability. Re-partitioning of real programs across
restarts (the same cone, a different basis) happens on the pancreas
halves and is what the one-to-one metric under-reports. Noise factors,
small, nearly orthogonal across restarts, with gene-soup profiles, are the
other, and `invsqrt_kl` produces them in proportion to the rank: it
resolves 4 to 7 programs on the breast replicates at any rank and fills
the rest of the rank with interchangeable fits of the residual (objective
spread 1 to 2%, unchanged by four times the iterations), which is why its
stability falls monotonically with the rank while `ls_nmf` keeps the
unrepresentable mass at 2% up to rank 23. The ensemble vote tolerates the
noise factors because their labels differ from member to member and
cancel in the vote. The likely cause is the inverse-square-root gene
weighting handing the objective to the low-count genes, whose Poisson
noise the extra factors fit; untested.

Decided. `min_stability = 0` is the default in both bindings: every
program that two restarts agree on enters, with recurrence (share of
samples) and stability (largest restart share in any sample) reported per
factor and ordered by recurrence; the terms are renamed throughout
(`fit_consensus(min_stability)`, `factors$recurrence`,
`factors$stability`, `factor_stability()`), and the docs no longer
prescribe a more stable variant or more restarts as the remedy for a
fired shared-factor flag. The example pages are re-rendered under the
new default (37 shared factors on the replicates). Open: a KL fit without
the gene weighting to test the cause of the noise tail; across-sample
recurrence as an optional admission gate.

## Loop 16 (2026-11-11): the KL noise tail is a sparse-row effect, not the gene weighting

Trigger. The user asked for the KL fit without the inverse-square-root
gene weighting and for the comparison with `ls_nmf`; the weighting had
been chosen from examples, to avoid factors dominated by single genes,
not by a robust evaluation. The core already has `kl` (raw neighborhood
counts) and `sqrt_kl` (square-root counts) next to `invsqrt_kl`.

Observed, rank ladder (`rank_ladder.R`, 10 restarts, 150 iterations,
20-molecule neighborhoods, 10,000 training rows;
`results/variant_ladder.csv`; mean one-to-one stability / share of the
molecules in factors more than 30% unrepresentable in other restarts):

| data | rank | `ls_nmf` | `kl` | `sqrt_kl` | `invsqrt_kl` |
|---|---|---|---|---|---|
| breast rep1 | 5 | | 0.77 / 0% | 0.73 / 17% | 0.53 / 4% |
| breast rep1 | 8 | 0.94 / 0% | 0.69 / 10% | 0.66 / 15% | 0.47 / 14% |
| breast rep1 | 12 | 0.80 / 1% | 0.55 / 42% | 0.53 / 24% | 0.42 / 10% |
| breast rep1 | 16 | 0.80 / 0% | 0.53 / 52% | 0.49 / 59% | 0.36 / 18% |
| breast rep1 | 23 | 0.82 / 1% | 0.43 / 56% | 0.40 / 47% | 0.31 / 24% |
| breast rep2 | 23 | 0.79 / 1% | 0.44 / 43% | 0.41 / 73% | 0.31 / 27% |
| pancreas left | 9 | 0.85 / 0% | 0.69 / 30% | 0.56 / 32% | 0.41 / 10% |
| pancreas right | 9 | 0.78 / 5% | 0.58 / 12% | 0.48 / 18% | 0.30 / 6% |
| pancreas left | 14 | 0.80 / 2% | 0.56 / 20% | 0.50 / 24% | 0.31 / 9% |
| pancreas right | 14 | 0.85 / 0% | 0.47 / 33% | 0.39 / 28% | 0.25 / 9% |

Observed, factor character of the selected restart (`factor_character.R`;
`results/variant_factor_character.csv`, per factor in
`results/variant_factors.csv`): single-gene factors are factors with more
than half of their unit-sum loading on one gene; a cell type has a home
factor when one factor holds at least half of its molecules.

| fit | single-gene factors (molecule share) | median effective genes | mass-weighted purity | cell types with a home factor |
|---|---|---|---|---|
| breast rep1 `ls_nmf` 23 | 0 | 36 | 0.53 | 13 of 19 |
| breast rep1 `invsqrt_kl` 23 | 0 | 17 | 0.48 | 5 of 19 |
| breast rep1 `sqrt_kl` 23 | 0 | 11 | 0.45 | 1 of 19 |
| breast rep1 `kl` 23 | 2 (9%) | 10 | 0.45 | 0 of 19 |
| pancreas right `ls_nmf` 9 | 0 | 50 | 0.73 | 5 of 7 |
| pancreas right `invsqrt_kl` 9 | 0 | 41 | 0.65 | 4 of 7 |
| pancreas right `kl` 9 | 3 (58%) | 26 | 0.65 | 4 of 7 |

On the pancreas right half the raw-count fit puts AMY2A (97% of its
factor's loading), GATM (76%) and CFTR (58%) into single-gene factors
holding 58% of the molecules, the acinar program in pieces; `ls_nmf` does
the same for AMY2A (79%) and INS (61%) at rank 14. On the breast the
raw-count fit splits the tumor program into GATA3, CCND1, TACSTD2,
ERBB2 (79%) and ANKRD30A pieces with 2 to 13 effective genes, where
`invsqrt_kl` keeps one tumor factor of 28 effective genes and `ls_nmf`
one of 34.

Observed, gene counts of the noise factors (`gene_weight_check.R`):
`invsqrt_kl` noise factors load on genes with mean log10 total 4.8
against 5.2 for its program factors (3 to 4 times the loading on the
bottom-quartile genes); `kl` noise factors load on high-count genes (5.3
against 5.35). The weighting moves which genes the spare factors fit; it
does not create them.

Observed, mechanism tests on `invsqrt_kl` at rank 23 on the breast
replicates (rep1 / rep2):

| change | stability | stable factors | unrepresentable mass | noise mass |
|---|---|---|---|---|
| default (20 molecules per row, 10,000 rows) | 0.31 / 0.31 | 6 / 6 | 24% / 23% | 24% / 27% |
| 600 iterations | 0.32 | 7 | 25% | 29% |
| 50,000 training rows | 0.33 / 0.33 | 5 / 5 | 25% / 20% | 28% / 23% |
| 50 molecules per row | 0.45 / 0.42 | 10 / 7 | 10% / 10% | 6% / 7% |
| 100 molecules per row | 0.51 / 0.54 | 11 / 13 | 6% / 5% | 7% / 4% |
| `ls_nmf`, 50 molecules per row | 0.83 / 0.80 | 22 / 20 | 4% / 2% | 2% / 1% |
| `ls_nmf`, 50,000 rows | 0.83 / 0.81 | 22 / 20 | 7% / 5% | 1% / 7% |

With 50-molecule rows the rank-23 `invsqrt_kl` fit recovers 13 programs
of 28 to 45 effective genes (tumor, two stromal, macrophages, T cells,
myoepithelial, endothelial, DCIS, B cells, dendritic cells,
proliferation, CD8 T cells, a second tumor program) and 7 to 9 of 19 cell
types get a home factor; the remaining noise factors hold 6% of the
molecules. On the pancreas halves 50-molecule rows cut the `invsqrt_kl`
noise mass from 10% / 6% to 3% / 4% at rank 9 and raise stability from
0.41 / 0.30 to 0.55 / 0.38; the rest of its instability there is the
re-partitioning of Loop 15.

Interpretation. (1) Removing the gene weighting does not remove the
noise tail; raw counts make it larger at every rank and bring back the
single-gene factors the weighting was introduced against. (2) The tail
is the KL objective fitting the sampling noise of sparse neighborhood
rows (20 molecules over 321 genes) with the spare rank: denser rows
remove most of it, more rows and more iterations do nothing, and every
count transform shows it. The least-squares variant with
inverse-frequency gene weights does not reward fitting the noise of
single rows and is unaffected by row density, row count and rank. (3)
`ls_nmf` is the most stable and representable variant on both datasets at
every rank, with the broadest programs and the most cell types resolved;
`invsqrt_kl` with 50- to 100-molecule rows comes close in noise mass but
not in one-to-one stability. The automatic neighborhood size
(20 × sqrt(genes / 400), capped at half the median cell size, never below
20) trades locality for counts and was set for the least-squares variant.

Decided. No default changes: the variant and neighborhood defaults stay,
because a larger neighborhood costs sub-cellular locality and only the
cleanup benchmark (docs/audit.md) can price that. The multi-sample page
names the lever for an unstable KL fit (a larger `ncv_k`), next to more
restarts. Open: the cleanup benchmark of `invsqrt_kl` at `ncv_k = 50`
against `ls_nmf`; the colleague's six-zone case, where `rank_ladder.R`
on one zone separates noise factors from re-partitioning.

## Loop 17 (2026-11-11): choosing the neighborhood size

Trigger. The user asked for strategies to pick the neighborhood size
`ncv_k`: it has to stay below the cell size to capture sub-cellular
admixture patterns, and it depends on the cell size, the panel size and
more.

Observed, the automatic rule and what it gives (`neighborhood_reach.py`;
k = 20 × sqrt(genes / 400), capped at half the median cell's molecules,
never below 20; reach is the median distance to the k-th within-cell
neighbor as a share of the cell's equivalent radius, 2D from the molecule
spread, over 1,500 cells of 30 or more molecules):

| dataset | genes | molecules per cell (p25 / p50 / p75) | automatic k | k as share of the median cell | reach at k (share of the radius) |
|---|---|---|---|---|---|
| breast rep1 | 321 | 121 / 206 / 323 | 20 | 9% | 0.29 |
| pancreas left | 377 | 15 / 30 / 52 | 20 | 38% (cap active) | 0.66 |
| breast 5k crop | 5,101 | 73 / 143 / 325 | 71 | 33% | about 0.62 |

Reach at 10 / 20 / 50 / 100 molecules: 0.20 / 0.29 / 0.48 / 0.71 on the
breast replicate, 0.42 / 0.66 / 1.06 / 1.26 on the pancreas half, 0.23 /
0.34 / 0.56 / 0.69 on the 5k crop.

Observed, the spatial structure of admixture (`marker_patch_scale.py`;
breast rep1; the molecules of genes whose share in the target type is 30
times below the source type's, inside target cells within 20 µm of a
source cell, five pairs): they sit at the periphery (median radial position
0.72 to 0.84 of the radius against 0.63 to 0.66 for all molecules), on the
half of the cell facing the source (16 to 22 points more than the cell's
molecules), and are weakly clustered: among a foreign molecule's k nearest
within-cell neighbors the foreign share is 2 to 3 times the cell's share at
k = 10, 2 times at k = 20, 1.1 to 1.7 times at k = 50. Admixture is a
source-facing rim whose foreign molecules are mixed with native ones at the
10- to 20-molecule scale. The molecules the shared correction removes show
the same decay (1.5 times at k = 5 to 10, none at k = 100), but that set is
shaped by the 20-molecule neighborhoods and the label smoothing.

Observed, the information side (`information_floor.R`; Bhattacharyya
distance per molecule between the fitted programs, molecules for an error
bound of exp(-3)): 2 to 3 molecules separate the median pair of programs on
the 321- and 377-gene panels and 4 to 8 on the 5,101-gene panel (hardest
tenth of pairs). Labeling against known programs needs far fewer molecules
than a row holds; the counts matter for finding the programs, and only
under the KL objective.

Observed, fit structure against k (`results/variant_ladder.csv`; stability
mean / noise mass):

| data | k | `ls_nmf` | `invsqrt_kl` |
|---|---|---|---|
| breast rep1, rank 23 | 10 | 0.78 / 1% | 0.23 / 53% |
| | 20 | 0.82 / 1% | 0.31 / 24% |
| | 50 | 0.83 / 2% | 0.45 / 6% |
| | 100 | | 0.51 / 7% |
| breast 5k crop, rank 8 | 20 | 0.60 / 2% | 0.09 / 23% |
| | 35 | 0.78 / 0% | 0.18 / 8% |
| | 71 (automatic) | 0.84 / 1% | 0.32 / 4% |
| | 140 | 0.72 / 0% | 0.39 / 1% |
| breast 5k crop, rank 16 | 20 | 0.67 / 0% | 0.05 / 45% |
| | 71 | 0.75 / 2% | 0.20 / 7% |
| | 140 | 0.68 / 2% | 0.25 / 3% |

Observed, cleanup against k (`cleanup_vs_k.R`;
`results/breast_rep1_cleanup_vs_ncv_k.csv`; `ls_nmf` rank 23, bridge
scoring, the 10-member vote at 0.3, audited over 172 detected pairs):

| k | stability | rules kept | leakage removed overall | median pair sensitivity | pairs below 0.2 | own-marker false removal |
|---|---|---|---|---|---|---|
| 10 | 0.78 | 207 | 0.61 | 0.89 | 30 | 0.9% |
| 20 | 0.82 | 213 | 0.63 | 0.86 | 22 | 1.6% |
| 50 | 0.83 | 248 | 0.64 | 0.84 | 24 | 2.5% |
| 100 | 0.79 | 250 | 0.61 | 0.79 | 23 | 2.8% |

The pairs move in both directions: DCIS_2 → Myoepi_ACTA2+ (393,000 admixed
molecules, the largest pair) is cleaned at 50 and 100 (0.82, 0.86) and not
at 10 and 20 (0.05, 0.28); DCIS_1 → Myoepi_ACTA2+ (28,000) the reverse
(0.86 at 20, 0.25 at 100); from 20 to 100, 29 pairs lose more than 0.1 and
7 gain.

Interpretation. (1) Two bounds. Locality bounds k from above: the
admixture is a facing rim mixed at the 10- to 20-molecule scale, so the
neighborhood belongs near a tenth of the cell's molecules (0.3 of its
radius); the breast cleanup agrees, with false removal falling
monotonically as k shrinks and the median pair sensitivity rising, while
overall leakage removal stays within three points from 10 to 100 and the
thin cell types lose their pairs at 100. Counts bound k from below, and the
bound belongs to the fit, not to the labeling: 2 to 8 molecules separate
the programs once they exist, but the KL objective needs 50 to 140 per row
to find them without a noise tail (a third to all of the cell on the 5k
panel and all of it on the pancreas), where the least-squares objective
finds them at 10 to 20 on every panel tested. (2) The current rule serves
the KL variant's count need (k grows with the square root of the panel)
and meets the locality bound only through the cap at half the median
cell: a third of the median cell on the 5k panel, two thirds on the
pancreas. (3) No single k serves every pair: the largest breast pair needs
50 or more, the thin-cell pairs 20 or less, and 10 leaves more pairs
without a rule.

Strategies, in order of preference:

1. Locality rule with a variant floor: k = f × median molecules per cell,
   f about 0.1, clamped to a floor that is the fit's count need (10 for
   `ls_nmf`; for the KL variants 50 at 300 genes and 70 to 140 at 5,000)
   and reported when the floor exceeds the locality bound instead of
   silently capping. On the datasets here this gives 20 (breast), 14 (5k
   crop) and 3 (pancreas, floored to 10) for `ls_nmf`.
2. Multi-scale ensemble: pool the molecule vote over members labeled at
   two or three neighborhood sizes (10, 20, 50) the way it pools over
   restarts, so the pairs that need different scales are both covered and
   the vote threshold keeps the specificity; the pair table shows both
   directions of movement, so this is the one that can beat every single
   k. The cost is one labeling pass per extra scale.
3. Decoupled training and labeling neighborhoods for the KL variants: fit
   the programs on dense rows (50 or more molecules, or whole-cell rows)
   and label molecules with small rows projected onto the fixed profiles;
   the fixed-profile path makes this testable without new code (a fit at
   50, then a fixed-profile run at 20), and a `ncv_k_train` option would
   make it a setting.
4. Per-cell neighborhoods, k_cell = max(floor, f × N_cell), so that the
   pancreas's 30-molecule cells and the breast's 200-molecule cells get
   the same relative reach; the within-cell search supports it directly.
5. A physical radius instead of a count equalizes the scale in µm across
   cell sizes but ties the counts to the local density, which reintroduces
   the count floor as a density floor; less attractive than 4.

Decided. Nothing changed in the code or the defaults. Open: the cleanup
ladder on the 5k crop, where the rule picks 71; strategies 2 and 3 cost one
audited cleanup run each on the breast replicate (about 1.5 hours per arm).

## Loop 18 (2026-11-11): the neighborhood ladder for the KL variant

Trigger. Loop 17 ran the cleanup ladder for `ls_nmf` only; the user pointed
out that the instability is the KL variant's. Same protocol
(`cleanup_vs_k.R`, now taking the variant; breast rep1, rank 23, bridge
scoring, the 10-member vote at 0.3, audited over 172 detected pairs), plus
the decoupled arms of `decoupled_kl.R`: the programs of the 50- or
100-molecule fit (its selected restart as the fixed profiles, all its
restarts as the member pool) with every molecule labeled from 20-molecule
rows.

Observed (`results/breast_rep1_cleanup_vs_ncv_k.csv`; strict-tier
sensitivity is 0.92 to 1.0 in every arm, so the arms separate on
specificity):

| fit | stability | rules kept | molecules removed (share of all) | leakage removed overall | median pair sensitivity | pairs below 0.2 | own-marker false removal |
|---|---|---|---|---|---|---|---|
| `ls_nmf`, 10 | 0.78 | 207 | 8.9M (23%) | 0.61 | 0.89 | 30 | 0.9% |
| `ls_nmf`, 20 | 0.82 | 213 | 8.9M (22%) | 0.63 | 0.86 | 22 | 1.6% |
| `ls_nmf`, 50 | 0.83 | 248 | 9.1M (23%) | 0.64 | 0.84 | 24 | 2.5% |
| `ls_nmf`, 100 | 0.79 | 250 | 8.7M (22%) | 0.61 | 0.79 | 23 | 2.8% |
| `invsqrt_kl`, 10 | 0.23 | 301 | 19.0M (48%) | 0.77 | 1.00 | 8 | 11.6% |
| `invsqrt_kl`, 20 | 0.31 | 271 | 14.4M (36%) | 0.70 | 1.00 | 16 | 7.1% |
| `invsqrt_kl`, 50 | 0.45 | 242 | 10.8M (27%) | 0.64 | 0.95 | 24 | 3.2% |
| `invsqrt_kl`, 100 | 0.51 | 253 | 10.3M (26%) | 0.65 | 0.88 | 25 | 3.4% |
| `invsqrt_kl`, programs at 50, labels at 20 | 0.45 | 242 | 10.7M (27%) | 0.64 | 0.95 | 27 | 2.6% |
| `invsqrt_kl`, programs at 100, labels at 20 | 0.51 | 253 | 10.2M (26%) | 0.65 | 0.89 | 25 | 3.3% |

Interpretation. (1) For the KL variant the ladder is the fit's noise, not
locality: the noise factors of the sparse-row fits carry removal rules, and
the 10- and 20-molecule fits remove half and a third of all molecules at
12% and 7% own-marker false removal; denser training rows bring the
removal to the least-squares level (26 to 27% of molecules) with false
removal of about 3%. (2) The labeling neighborhood hardly matters: the
decoupled arms reproduce their training fit's cleanup almost exactly, with
false removal 0.6 points lower at 50 and the same at 100. The neighborhood
size acts through the programs, so for the KL variant the training rows
can be dense while the labels stay local, and strategy 3 of Loop 17 costs
nothing in sensitivity. (3) Even so, the best KL arm (programs at 50, 2.6%)
removes own markers at 1.6 times the rate of `ls_nmf` at 20 (1.6%) and
three times `ls_nmf` at 10 (0.9%), for the same overall leakage removal
(0.64 against 0.63); the least-squares fit remains the better default on
this dataset, and its own ladder is governed by locality alone.

Decided. Nothing changed in the code or defaults. For a KL fit, the
recommendation is the decoupled form: 50 or more molecules per training
row (70 to 140 on a 5,000-gene panel), labels at the locality rule's size.
Open: the same two ladders on the 5k crop and on the pancreas halves, where
the count floor and the cell size collide.
