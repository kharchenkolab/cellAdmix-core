# celladmix Python bindings

This package exposes the C++ cellAdmix core to Python. The initial API mirrors
the R workflow while returning Python-native objects such as pandas DataFrames,
SciPy sparse matrices, and AnnData objects.

The package is currently experimental. The primary target is scverse workflows,
including SpatialData integration.

```python
import pandas as pd
import celladmix as ca

annotation = pd.read_csv("annotations/annotation.csv.gz")
labels = annotation.set_index("cell_id")["merged_annotation"]

ds = ca.CellAdmix("data", output_dir="out", annotation=labels, num_threads=10)
fit = ds.fit()
score = fit.score_membrane()
rules = score.rules()
correction = score.correct(rules=rules)
counts, genes, cells = correction.counts()
```

Existing output directories can be reattached without refitting:

```python
fit = ca.CellAdmixFit.load("out/runs/fit_rank8_invsqrt_kl",
                          source="data",
                          annotation=labels)
```

Public factor ids are one-based. Molecule tables returned by
`fit.molecules()` and `fit.score_molecules()` use `factor` (`1..K`) and
`factor_label` (`F1..FK`); pass `raw=True` only to inspect internal parquet
labels.

SpatialData integration is available through `ca.from_spatialdata()`.
Rendering SpatialData-native figures also uses `spatialdata-plot`; current
SpatialData packages require Python 3.11 or newer.

For command-line Xenium processing, use:

```bash
python scripts/celladmix_batch.py --input data --output out \
  --annotation annotations/annotation.csv.gz \
  --annotation-col merged_annotation \
  --score auto --report --threads 10
```

See `docs/python-batch.md` for details.
