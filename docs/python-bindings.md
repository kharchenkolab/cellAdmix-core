# Python Bindings

The Python package exposes the same C++ analysis store and persisted run
artifacts used by the R package. The current implementation is focused on
Xenium bundles and the scverse ecosystem; tabular and other input adapters can
reuse the same store abstraction later.

## Install

From the repository root:

```bash
python -m pip install -e python
```

SpatialData integration requires a Python 3.11 or newer environment because
current `spatialdata`, `spatialdata-io`, and `spatialdata-plot` releases
require Python `>=3.11`. The core Python bindings themselves currently support
Python 3.8+.

See [install.md](install.md) for C++ system dependencies, optional extras, and
verification commands.

## Standalone Xenium API

```python
import celladmix as ca

ds = ca.CellAdmix("data", output_dir="out", annotation=cell_annotation)
fit = ds.fit()
score = fit.score_membrane()
rules = score.rules()
correction = score.correct(rules=rules)
counts, genes, cells = correction.counts()
```

`CellAdmix()` validates configuration and owns an on-disk output/cache
directory. The full molecule table is not loaded into Python memory. Fitting,
scoring, and correction call the C++ core through pybind11.

Fitted runs are cached on disk by run id and reused only when the requested
parameters match the cached run: `fit()` compares the requested parameters
(and the annotation content) against the run manifest, reuses on a match with
a `Reusing cached run ...` message, and refits automatically when anything
differs, listing the changes. Omitted auto-resolved values (`ncv_k`,
`nmf_n_runs`) match whatever the cached run recorded; `overwrite=True` forces
a refit. Scores and corrections recompute on every call and never go stale.

`score.rules()` applies a native-factor false-positive check by default:
rules whose factor persists in target cells that have no source-type cells
among their nearest neighbors are flagged `keep=False` with the reason in the
`native_check` column, and `correct()` skips them. Pass `native_check=False`
to disable, or tune `native_median_thresh` and related thresholds. See the
Native-Factor Check section in
[docs/scoring_methods.md](scoring_methods.md).

`score.correct()` applies the correction as a molecule-vote ensemble by
default: up to `ensemble` (default 10) of the fit's NMF restarts are scored
and vetted independently (each with its own native-factor check), and a
molecule is removed when at least `vote` (default 0.3) of the members remove
it. The
vote stabilizes the seed-dependence of single-fit corrections and acts as a
sensitivity/specificity dial (see [benchmarks.md](benchmarks.md)). Pass
`ensemble=1` for a single-fit correction from the selected restart;
`correction.ensemble()` reports the member count, minimum votes, and vote
histogram, and the returned `correction.rules` carry a `support` column with
the fraction of members keeping each source→target pair. The first ensemble
correction computes and caches per-member molecule labelings in the run
directory; runs fitted with `nmf_n_runs=1` (or cached runs whose stores lack
a restart member pool) fall back to the single-fit correction with a message.

`correct()` warns when a rule set removes more than 60% of a cell type's
molecules — usually a sign that the factorization has no native factor for
that type; the minimum-molecule gate for the warning is
`celladmix.fit.OVERREMOVAL_MIN_MOLECULES`.

Default behavior mirrors the R API:

- NMF method: `ls_nmf`.
- Molecule node potentials: gene-loading based scoring.
- Rank: `ceil(1.2 * number_of_annotation_labels)`, capped at 30.
- Xenium control/codeword/non-gene features are excluded during input-store
  construction. Pass `keep_non_gene=True` only for control-feature diagnostics.
- NMF restarts: at least 10 by default (more on higher thread counts), with
  every restart's loadings kept as the ensemble member pool.

### Loading Existing Runs

Persisted runs can be reattached without refitting or manually globbing the
private run directory layout:

```python
fit = ca.CellAdmixFit.load("out/runs/fit_rank8_ls_nmf",
                          source="data",
                          annotation=cell_annotation)
```

or, when working from an output directory:

```python
ds = ca.CellAdmix.attach_existing("out", source="data", annotation=cell_annotation)
fit = ds.load_fit("fit_rank8_ls_nmf")
```

The `source` argument is optional when the input-store manifest records a usable
source path. It is required for methods that need original bundle-side files,
such as membrane image discovery or cell-boundary plotting, if the path cannot
be inferred.

### Factor Indexing

The Python public API uses one-based factors throughout:

- integer factor ids are in `1..K` and are stored in a `factor` column;
- display labels are strings such as `F1`, `F2`, ... in `factor_label`;
- cell-level columns remain `factor_1_fraction`, `factor_2_fraction`, ...
  and `dominant_factor` is one-based.

The on-disk `molecules.parquet` file stores the implementation label as a
zero-based integer. This is an internal detail. Use `fit.molecules(raw=True)`
only when debugging the persisted parquet representation directly.

### Scoring Additional Molecules

Small same-fit or synthetic molecule tables can be scored against the learned
factor loadings without refitting:

```python
scored = fit.score_molecules(
    new_molecules,              # columns include at least "gene"
    return_scores=True,
)
```

The returned table includes `factor`, `factor_label`, and `factor_margin`.
With `return_scores=True`, it also includes pre-smoothing soft scores in
columns `F1..FK`. These scores are the gene-loading molecule potentials used by
the default pipeline.

This helper only scores the molecule rows supplied by the caller. Molecules
removed during input-store construction, for example by a strict QV threshold
or by excluding unassigned molecules, are not recoverable from a completed run
unless they are supplied separately or the store was built with permissive
filters.

To also apply the same within-cell hard smoothing used by the fit pipeline:

```python
scored = fit.score_molecules(new_molecules, smooth=True)
```

`smooth=True` requires `cell_id`, `x`, and `y` columns. The smoothing step is
iterative hard-label ICM, not marginal-probability CRF inference, so the
`F1..FK` columns, when requested, remain the pre-smoothing soft potentials.

Cell-level `factor_K_fraction` values are fractions of molecules assigned to
each hard factor label after smoothing. They are not soft cell probabilities.

### Auditing Admixture and Verifying Cleanup

Independently of factorization and scoring, admixture can be estimated
from the dataset's spatial structure: source-marker content in target
cells rises with source-type neighbor exposure, while unexposed target
cells provide an internal negative control (see
[benchmarks.md](benchmarks.md) for the methodology). Because a section
shows only a slab of the tissue, cells with zero observed source-type
neighbors can still carry material from source cells above or below the
section plane; the audit therefore takes its reference level from target
cells with zero source-type cells among a progressively larger set of
nearest neighbors (up to 240, spanning only a ~100 um radius, so the
comparison stays within the local tissue neighborhood), using the largest
neighborhood that retains enough reference molecules. The `reference_kind`
and `reference_inflation` columns of `pairs()` record the neighborhood used
and how much contamination this removed from the comparison group. Marker panels are screened for likely
induced genes — exposure-linked excess far above the gene's share of the
source expression profile indicates a transcriptional response to
proximity rather than transferred material; such genes are excluded,
replaced by the next-ranked source-specific genes, and listed in
`markers()["induced"]`. The audit measures every ordered cell-type pair
and verifies corrections against the same measurements:

```python
audit = fit.audit_admixture()
audit.pairs()                       # per-pair admixture rates and molecule estimates
audit.plot_map()                    # admixture-rate map (% of target-type molecules)
audit.plot_exposure()               # pooled excess-exposure profile, 95% intervals
audit.plot_remaining({"membrane": correction})

report = audit.evaluate(correction) # per-pair sensitivity, false removal,
report.summary()                    # warnings for pairs left largely
report.plot_cleanup()               # uncorrected and for severe own-marker
                                    # over-removal (>25%)
```

`evaluate()` warns from the corrected counts themselves — a pair is
flagged when its molecules were measurably not removed, regardless of what
the correction's rule list claims. The reported rates and molecule counts
include the contamination present in zero-neighbor cells above the ambient
reference, and extrapolate the marker-panel measurement to the full
transcriptome by the markers' share of the source expression profile, kept
as the `coverage` column of `audit.pairs()`. Pair detection itself remains
based on the rise of exposed cells over unexposed ones, which no reference
choice can inflate.

## SpatialData API

SpatialData integration keeps SpatialData as the Python-side source of
cell-level state while cellAdmix continues to stream molecules and stains from
the original Xenium bundle:

```python
import spatialdata_io as sd_io
import celladmix as ca

sdata = sd_io.xenium("data")
ds = ca.from_spatialdata(
    sdata,
    xenium_dir="data",
    output_dir="out_spatialdata",
    annotation="cell_type",
)
fit = ds.fit()
ca.add_fit_to_spatialdata(sdata, fit)

score = fit.score_membrane()
correction = score.correct()
ca.add_corrected_counts_to_spatialdata(
    sdata,
    correction,
    layer="celladmix_corrected_counts",
)
```

The explicit `xenium_dir` argument is intentional. SpatialData objects do not
reliably preserve enough original-bundle path information to support optimized
streaming reads and membrane-image discovery.

## Examples

- [Standalone Xenium pancreas workflow](../examples/xenium_pancreas_membrane_377_full/pancreas_python_standalone.ipynb) -
  loads a Xenium bundle through the Python API, fits factors, scores membrane
  admixture, applies correction, and shows standard diagnostics.
- [SpatialData Xenium pancreas integration](../examples/xenium_pancreas_membrane_377_full/pancreas_spatialdata_integration.ipynb) -
  reads the same dataset as a `SpatialData` object, adds cellAdmix factors and
  corrected counts back into the SpatialData table, and illustrates Scanpy /
  SpatialData plotting on the outputs.
- [Standalone Xenium breast 5K workflow](../examples/xenium_breast_membrane_5k_full/breast_5k_python_standalone.ipynb) -
  the same core workflow on a full 5K-panel dataset (82.1M molecules,
  688K cells after QC).

These notebooks are intended as starting templates. The SpatialData notebook
needs a Python 3.11+ environment with `spatialdata`, `spatialdata-io`, and
`spatialdata-plot`.

## Batch Processing

The Python package also includes a CLI wrapper for reproducible batch runs:

```bash
python scripts/celladmix_batch.py --input data --output out \
  --annotation annotations/annotation.csv.gz \
  --annotation-col merged_annotation \
  --score auto --report --threads 10
```

The script writes factor top genes, score summaries, correction rules,
correction summaries, a JSON/CSV batch summary, and optionally a compact
self-contained HTML report. See [python-batch.md](python-batch.md) for the
option list and current limitations.
