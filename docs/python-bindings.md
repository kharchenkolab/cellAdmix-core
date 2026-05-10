# Python Bindings

The Python package exposes the same C++ analysis store and persisted run
artifacts used by the R package. The current implementation is focused on
Xenium bundles and the scverse ecosystem; tabular and other input adapters can
reuse the same store abstraction later.

## Install

From the repository root:

```bash
python -m pip install -e python --no-build-isolation
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

Default behavior mirrors the R API:

- NMF method: `invsqrt_kl`.
- Molecule node potentials: gene-loading based scoring.
- Rank: `ceil(1.2 * number_of_annotation_labels)`, capped at 30.
- Xenium control/codeword/non-gene features are excluded during input-store
  construction. Pass `keep_non_gene=True` only for control-feature diagnostics.
- NMF restarts: default to the dataset thread count.

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
