# Python Batch Processing

`scripts/celladmix_batch.py` runs the core cellAdmix workflow from the Python
bindings without opening a notebook. It is intended for command-line Xenium
processing and mirrors the high-level behavior of `scripts/celladmix_batch.R`
where the Python bindings currently expose the needed backend features.

## Minimal Run

```bash
python scripts/celladmix_batch.py --input data --output out \
  --annotation annotations/annotation.csv.gz \
  --annotation-col merged_annotation \
  --score auto --report --threads 10
```

This performs:

- dataset construction and on-disk input-store reuse/building;
- NMF factor fitting with the default `invsqrt_kl` variant;
- membrane scoring when a membrane image is discoverable, otherwise bridge
  scoring;
- score-rule generation at `--p-thresh`;
- molecule correction unless `--no-correct` is supplied;
- summary table output and optional HTML report generation.

If no annotation is available, pass `--auto-annotate` to cluster cell-level
expression profiles and use the resulting clusters as provisional labels:

```bash
python scripts/celladmix_batch.py --input data --output out \
  --auto-annotate --score auto --report --threads 10
```

With `--report`, auto-annotation also writes `annotation_report.html`.
`--annotation-report` can be used to request that report explicitly.

## Outputs

The output directory receives:

- `batch_summary.json` and `batch_summary.csv`: run-level metadata;
- `fit_summary.csv`: compact fit metadata;
- `factor_top_genes.csv`: top loading genes per factor;
- `score_summary.csv`: factor/source/target score summaries;
- `score_rules.csv`: factor/target correction rules;
- `correction_summary.csv`: compact per-cell-type correction summary when
  correction is run;
- `report_minimal.html`: compact report when `--report` is supplied.

For `--auto-annotate`, the output directory also receives:

- `annotation_clusters.csv`: cell ids and generated cluster labels;
- `annotation_cluster_markers.csv`: quick one-vs-rest marker summaries for
  the generated clusters when an annotation report is requested;
- `annotation_report.html`: UMAP/spatial cluster overview when `--report` or
  `--annotation-report` is supplied.

The fit, molecule labels, corrected run, and input-store caches are written
under the same `output_dir` run-store layout used by the notebook API.

## Common Options

```bash
python scripts/celladmix_batch.py --help
```

Important options:

- `--threads N`: default worker count for fit, scoring, and NMF restarts.
- `--rank auto|N`: factor rank; `auto` uses the annotation-based default.
- `--nmf-variant invsqrt_kl|kl|sqrt_kl|ls_nmf`: NMF formulation.
- `--nmf-runs auto|N`: multiseed NMF restarts; `auto` uses `--threads`.
- `--score auto|membrane|bridge`: scoring method.
- `--targets A,B,C`: restrict correction rules to selected target cell types.
- `--auto-annotate`: cluster cells when no annotation is supplied.
- `--annotation-report`: write an annotation report even without the main
  batch report.
- `--no-correct`: fit and score only.
- `--report`: write a compact self-contained HTML report.

## Current Python Limitations

Python v1 focuses on Xenium bundle processing. The CLI exposes tabular and
related schema options for parity with the R script, but generic tabular input
currently fails with an explicit message:

- Use `scripts/celladmix_batch.R` for generic tabular inputs until the Python
  tabular-store constructor is exposed.
- Coherence scoring is not exposed in the Python batch script yet.

The Python standalone and SpatialData notebooks cover interactive workflows:

- `examples/xenium_pancreas_membrane_377_full/pancreas_python_standalone.ipynb`
- `examples/xenium_pancreas_membrane_377_full/pancreas_spatialdata_integration.ipynb`
