---
name: celladmix-generative
description: The generative admixture model of cellAdmix - how a target cell's counts are decomposed into own expression, per-source contamination, ambient background, and neighborhood-induced expression; how to fit it from the audit, derive a correction that retains induced genes, and read its outputs. Use when a correction must preserve interface biology or when working with audit$fit_generative().
---

Use this skill when a cleanup should remove transferred molecules but keep expression that the target cell itself induces at a tissue interface, when you call `audit$fit_generative()` / `audit.fit_generative()` or `correct_generative()`, or when you need to interpret the model's per-cell decomposition and its list of induced genes. The authoritative page is `docs/generative.md`; the implementation is `src/generative.cpp` with `include/celladmix/generative.hpp`, wrapped by `r/R/audit.R` (`CellAdmixGenerativeModel`, `CellAdmixGenerativeCorrection`) and `python/src/celladmix/generative.py` (`GenerativeModel`, `GenerativeCorrection`). The audit vocabulary (exposure, marker pools, reference) is in `celladmix-audit-and-stability`.

## Why a model

Exposure-linked excess of source genes in target cells has two causes: material physically transferred from neighbors (spilled fragments, out-of-plane material, ambient background) and the target cell's own transcription responding to its neighborhood (an activation program at an interface). A rule-based correction removes both. The generative model separates them so that induced expression can be retained while transferred material is removed, and it does so without any factorization: the NMF enters only as an optional initializer of the expression programs.

## The model

For each annotated target type, the count of cell c at gene g is Poisson with mean equal to the cell's total molecule count times the sum of:

- **Own expression**: the type's expression programs (gene profiles with per-cell weights), re-estimated during fitting mostly from lightly contaminated cells.
- **Contamination**, one term per detected source type: the cytoplasmic profile of the source cells bordering the target type (within 30 um, widening when sparse), zeroed on target-owned genes (a gene is owned by the type where it makes up the largest share of the transcriptome), so removal of the target's own markers is structurally impossible. The per-cell fraction has a prior from a monotone dose-response of marker content on source-neighbor count and a bounded posterior update (at most five-fold, zero where the prior is zero), which keeps removal anchored to the spatial signal.
- **Ambient background** on the strict genes (essentially absent from the target natively), with level measured on target cells far from any source.
- **Induced expression**: a sparse term on genes selected by the proportionality screen (below), with a per-cell activity that concentrates the induced content in the responding cells. This content is retained; the genes' proportional transferred share is still removed.
- A uniform floor of 0.2% of cell content.

Each observed count is split among the components in proportion to their fitted rates; the removed mass is the contamination plus ambient share.

## The induced-expression screen

Transferred material samples the source transcriptome, so a gene's exposure-linked excess in target cells should be proportional to its share of the interface-local source profile. Per pair, the per-gene excess is regressed on that profile and the residual standardized by counting noise plus a 15% multiplicative profile-uncertainty allowance; genes whose excess exceeds the expectation beyond the threshold are flagged as induced (real flagged genes deviate 2.3- to 22-fold; a 1.08-fold deviation on a large transfer channel is not flagged). The screen examines only source-owned genes, the ones whose excess could be mistaken for transfer. The comparison profile comes from the bordering source cells because interface source cells can genuinely express activation genes above the source average; against a global profile such genes would be misread as induced. The same screen cleans the audit's marker panels (`markers()$induced`).

Identifiability limit: induction is preserved when it is disproportionate to the local source profile (retention turns on at about three-fold the transfer expectation and is essentially complete beyond ten-fold); induction proportionate to a source's own top markers is indistinguishable from admixture at count level and is removed.

## Usage

```r
ds <- cellAdmix(source, output_dir = "out", annotation = annotation)
audit <- ds$audit_admixture()            # measurement: counts, positions, annotation
nmf_fit <- ds$fit()                      # optional: initializes the expression programs
model <- audit$fit_generative(init = nmf_fit, n_programs = 4L, seed = 1L, num_threads = 8)
model$pairs                              # per pair: removed and induced molecule totals
model$induced                            # retained induced genes with their disproportionality
model$composition("Exocrine epithelial", "Ductal/tumor epithelial")
                                         # per target cell: dose, contamination fraction, induced activity
correction <- model$correct()            # corrected counts; induced expression kept
alt <- model$correct(retain_induced = FALSE)   # induced share removed along with contamination
audit$evaluate(correction)$summary()     # verified like any other correction
correction$counts()                      # sparse gene x cell counts
```

```python
audit = ds.audit_admixture()
model = audit.fit_generative(init=nmf_fit, num_threads=8)
model.pairs; model.induced; model.composition("Exocrine epithelial", "Ductal/tumor epithelial")
correction = model.correct()             # or model.correct(retain_induced=False)
audit.evaluate(correction).summary()
```

- `init` selects the expression-program initialization: an NMF fit (its factor-labeled molecules; the recommended configuration and the default for an audit made from a fit via `fit$audit_admixture()`), `"clusters"` (clustering each type's cells weighted toward lightly dosed ones; the default for an audit made from a dataset), `"pseudobulk"` (one profile per type), or a named list (R) / dict (Python) of programs-by-genes matrices per cell type for explicit programs, for example from an external reference. The initialization matters little for removal (agreement 0.91-1.00 across datasets) and more for retention on large panels (on breast 5K the factor initialization flags 180 induced genes and keeps 78% of their excess versus 129 genes and 62% otherwise).
- `n_programs` (4) is the number of expression programs per type; further keyword arguments override the constants of `GenerativeOptions` in `include/celladmix/generative.hpp` (contamination prior strength and update cap, ambient settings, screen threshold and profile-uncertainty allowance, interface radius, dose weighting, EM schedule, exposure neighbor count, component switches).
- `audit$correct_generative(name = "generative", ...)` / `audit.correct_generative(**kwargs)` is the one-call form; it discards the fitted model's decomposition.
- The fit runs in the C++ core: about 20 s on the pancreas dataset, 50 s on CosMx NSCLC, about four minutes on breast 5K. R and Python marshal the same inputs to one implementation and give identical results.

The returned correction (`CellAdmixGenerativeCorrection` / `GenerativeCorrection`) exposes `counts()`, `rules` (the pairs it acted on), `retained_induced`, and `cell_state_umap()`; it is accepted by `audit$evaluate()`.

## Reading the outputs

- `model$pairs`: per source -> target pair, the molecules removed as contamination plus ambient and the molecules retained as induced.
- `model$induced`: the flagged genes per pair with excess, transfer expectation, fold disproportionality, and retained share; inspect these before trusting a retention decision, since material from a source with no detected pair (an unannotated or rare type) has no term to load onto and could be misread as induced if its arrival tracks a modeled source.
- `model$composition(source, target)`: per target cell, the dose (source neighbor count), the prior and fitted contamination fraction, and the induced activity; stratifying by dose shows how the attributed composition changes with exposure.
- `audit$evaluate(correction)`: sensitivity per pair and own-marker false removal as for any correction; the model's design makes own-marker removal zero by construction, so compare it to the rule-based ensemble on sensitivity and on how much induced excess survives (`docs/generative.md`, Validation).

## Validation record

Held-out-gene generalization 0.974 on pancreas (ensemble 0.961 at 7.5% own-marker false removal; model 0.000), replication on NSCLC (0.964) and breast 5K (0.863, with strict-tier leakage removed completely and the shortfall on genes shared between source and target), a spike-in test, and an exposure-permutation control (removal collapses to 0.018) are described in `docs/generative.md`; scripts and tables are in `analysis/generative_model/` and `analysis/audit_guided/`. Worked examples: `examples/xenium_pancreas_membrane_377_full/pancreas_generative_minimal.md`, `pancreas_generative.md`, `examples/cosmx_nsclc_giotto/nsclc_generative.md`, `examples/xenium_breast_membrane_5k_full/breast_5k_generative.ipynb`.
