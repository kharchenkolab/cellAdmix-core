---
name: celladmix-python-api
description: The Python user API of celladmix - CellAdmix, CellAdmixFit, CellAdmixScore, CellAdmixCorrection, the audit objects, run reattachment, factor indexing, SpatialData helpers, and the command-line batch script. Use when writing or reading Python code that drives the pipeline.
---

Use this skill when writing Python against the `celladmix` package, reattaching to runs produced by either binding, integrating with SpatialData or AnnData, or running `scripts/celladmix_batch.py`. Sources: `python/src/celladmix/*.py` (docstrings), `python/README.md`, `docs/python-bindings.md`, `docs/python-batch.md`. Scoring semantics are in `celladmix-scoring-and-correction`; audit details in `celladmix-audit-and-stability`.

## Minimal workflow

```python
import pandas as pd, celladmix as ca
labels = pd.read_csv("annotation.csv.gz").set_index("cell_id")["merged_annotation"]
ds = ca.CellAdmix("xenium_bundle", output_dir="out", annotation=labels, num_threads=10)
fit = ds.fit()
audit = fit.audit_admixture()
score = fit.score_membrane()               # or fit.score_bridge()
rules = score.rules()
correction = score.correct(rules=rules)    # molecule-vote ensemble, vote=0.3
counts, genes, cells = correction.counts() # scipy CSC, genes x cells
report = audit.evaluate(correction); report.summary()
```

Public names (`celladmix/__init__.py`): `CellAdmix`, `CellAdmixFit`, `CellAdmixScore`, `CellAdmixCorrection`, `CellAdmixFactorSourceScore`, `CellAdmixAudit`, `CellAdmixCleanupReport`, `CellAdmixCollection`, `CellAdmixConsensusFit`, `CellAdmixCollectionScore`, `CellAdmixCollectionCorrection`, `CellAdmixCollectionAudit`, `CellAdmixCollectionCleanupReport`, `read_annotation`, `plot_spatial`, `annotation_knn_purity`, `from_spatialdata`, `add_fit_to_spatialdata`, `add_corrected_counts_to_spatialdata`. Python input support is Xenium bundles only; tabular tables go through R.

## CellAdmix (dataset)

`CellAdmix(source=None, *, output_dir="out", format="xenium", annotation=None, annotation_col=None, cell_id_col="cell_id", num_threads=None, cell_filter=None, gene_filter=None, min_qv=-1.0, keep_unassigned=False, keep_non_gene=False)`.

- `attach_existing(output_dir, *, source=None, format=None, annotation=None, ...)` (classmethod): reattach to an output directory; `source` defaults to the path recorded in `input_store/store.json`.
- `ensure_store(*, force=False, materialize_molecules=True, verbose=True)`, `counts()` -> `(matrix, genes, cells)`, `cells()` (cell id, x, y), `set_annotation(annotation, *, annotation_col, cell_id_col)`.
- `fit(*, rank=None, rank_multiplier=1.2, rank_cap=30, nmf_variant="ls_nmf", nmf_init="auto", nmf_n_runs=None, molecule_scoring="gene_loadings", run_id=None, overwrite=False, num_threads=None, verbose=True, **kwargs)`. `kwargs` reach the core: `ncv_k` (0 = auto), `graph_k=10`, `same_label_ratio=5.0`, `nmf_iterations=150`, `nmf_train_max_rows=10000`, `nmf_min_molecules=10`, `seed=1`, `training_scope_cell_types`, `nmf_fixed_h` (factor-by-gene DataFrame with gene columns), `nmf_fixed_h_members` (list of such frames). Rank defaults to `ceil(1.2 * n_labels)` capped at 30; `nmf_n_runs` to `max(num_threads, 10)`; run id to `fit_rank<K>_<variant>`.
- `load_fit(run_id_or_path)`, `audit_admixture(**kwargs)` (no fit needed), `cell_state_umap(...)`.

Caching follows the R rules: a cached `run.json` is reused when `rank`, `nmf_variant`, `nmf_init`, `molecule_scoring`, explicitly passed core parameters, and the annotation hash match (`Reusing cached run ...`); otherwise the differences are printed and the run is refitted; `overwrite=True` forces it.

Reattaching a run directly:

```python
fit = ca.CellAdmixFit.load("out/runs/fit_rank8_ls_nmf", source="xenium_bundle", annotation=labels)
```

## CellAdmixFit

Attributes: `dataset`, `run_path`, `manifest` (the parsed `run.json`, including `nmf_diagnostics`).

- Data: `cells()`, `cell_factors()`, `molecules(columns=None, *, raw=False)`, `region(*, bbox, cache=True)`, `factor_loadings()` (gene x factor DataFrame, columns `F1..FK`), `counts()`.
- Scoring: `score_membrane(*, image_path=None, pixel_size=None, num_threads=None, verbose=True, **kwargs)` (image discovered from the bundle when omitted), `score_bridge(*, num_threads=None, verbose=True, **kwargs)`, `score_factor_sources(*, annotation=None, counts=None, **kwargs)`. Keyword defaults are those of the core bindings: bridge `candidate_mode="cell_center"`, `candidate_pairs_per_type_pair=1000`, `crossing_k=20`, `min_factor_molecules=5`, `min_pairs=5`, `max_cells_per_type_pair=400`, `null_iterations=3`, `compute_null=True`; membrane `candidate_pairs_per_type_pair=400`, `max_cells_per_type_pair=400`, `control_distance_fraction=0.025`, `line_samples=16`, `x_offset=0`, `y_offset=0`. Coherence scoring is not exposed in Python.
- `score_molecules(molecules, *, gene_col="gene", cell_col="cell_id", x_col="x", y_col="y", return_scores=False, smooth=False, ...)`: label additional molecule rows against the learned loadings without refitting; `smooth=True` applies the same within-cell hard smoothing.
- `correct(rules, *, name="clean", rule_members=None, min_votes=1)`: the low-level single rule-table correction; use `CellAdmixScore.correct()` for the ensemble.
- `audit_admixture(**kwargs)`, `plot_loadings(*, n_genes=8)`, `plot_stability(**kwargs)`, `stain(stain="membrane")`, `stain_image()`, `stains()`, `stain_crop(image, *, bbox, max_pixels=512)`, `cell_boundaries(*, boundary_path, cells, bbox)`, `prepare_cell_example()`, `plot_cell_example()`.

`score_neighbor_enrichment()` and `score_neighbor_frequency()` import a module `celladmix.neighbor` that does not exist in the package; do not call them.

### Factor indexing

Public tables use one-based factors: `factor` in `1..K`, `factor_label` as `F1..FK`, cell columns `factor_<k>_fraction`, `dominant_factor` one-based. `molecules.parquet` stores zero-based labels; `fit.molecules(raw=True)` exposes them. Cell-level fractions are fractions of hard labels after smoothing, not soft probabilities.

## CellAdmixScore

Attributes: `fit`, `method`, `summary` (DataFrame), `pairs` (DataFrame), `params`.

- `annotation(*, p_thresh=0.1, adjust_p=False)`; `rules(*, p_thresh=0.1, adjust_p=False, targets=None, max_rules=None, native_check=True, native_median_thresh=0.1, native_expr_thresh=0.05, native_outlier_min_frac=0.1, neighbor_k=15)` -> columns `factor`, `source_cell_type`, `target_cell_type`, `p_value`, `neg_log10_p`, `rule_id`, plus `keep`, `native_check`, and the supporting `native_*` columns.
- `correct(rules=None, *, p_thresh=0.1, name=None, ensemble=None, vote=0.3)`: ensemble of up to 10 members by default; `ensemble=1` is a single-fit correction; explicit `rules` restrict every member to their source/target pairs; the returned `correction.rules` carries `support`.
- `plot_heatmap(*, p_thresh=0.1, adjust_p=False, ...)` (vetoed calls marked `(*)`), `plot_pairs(*, factors=None, factor=None, **kwargs)`, `examples(**kwargs)`, `plot_example(example, *, p_thresh, adjust_p, **kwargs)`, `plot_examples(examples=None, *, cells=None, **kwargs)`.

## CellAdmixCorrection

Attributes: `run_path`, `manifest`, `rules`. Methods: `counts()` -> `(csc_matrix, genes, cells)`, `to_anndata()`, `summary()` (per cell type; first row `all`), `cell_summary()`, `ensemble()` -> `members`, `min_votes`, `n_removed`, `vote_histogram`; `plot_removed_molecules()`, `cell_state_umap(...)`. `fit.correct()` warns when more than 60% of a cell type's molecules are removed (gate `celladmix.fit.OVERREMOVAL_MIN_MOLECULES = 10000`).

## Audit objects

`CellAdmixAudit` (built by `fit.audit_admixture()` or `ds.audit_admixture()`; constructor keywords `neighbor_k=15, n_pool=20, q_thresh=0.01, min_excess=200, min_target_cells=200, min_reference_cells=100`): `pairs(detected_only=False)`, `markers(source, target)`, `plot_map(value="rate", detected_only=True, ax=None)`, `plot_exposure(source=None, target=None, correction=None, ...)`, `plot_remaining(corrections=None, ax=None)`, `exposure(source=None, target=None, correction=None, strict=True)` and `remaining(corrections=None)` (the tables behind the two plots), `plot_reference(source, target, ax=None)`, `evaluate(correction, warn_uncovered=True)`, `fit_generative(**kwargs)`, `correct_generative(**kwargs)`. `CellAdmixCleanupReport`: `pairs()`, `false_removal()`, `summary()`, `plot_cleanup(ax=None)`, `plot_map(value="rate", axes=None)` (before and after the correction, one color scale). See `celladmix-audit-and-stability` and `celladmix-generative`.

## SpatialData integration

```python
import spatialdata_io as sd_io, celladmix as ca
sdata = sd_io.xenium("xenium_bundle")
ds = ca.from_spatialdata(sdata, xenium_dir="xenium_bundle", output_dir="out_sdata",
                         annotation="cell_type")       # table=, cell_id_key=, num_threads=
fit = ds.fit(); ca.add_fit_to_spatialdata(sdata, fit)   # obs columns celladmix_factor_<k>_fraction, celladmix_dominant_factor
correction = fit.score_membrane().correct()
ca.add_corrected_counts_to_spatialdata(sdata, correction, layer="celladmix_corrected_counts")
```

`xenium_dir` is required because SpatialData objects do not retain enough bundle path information for streaming reads and image discovery. With one table in `sdata.tables`, `table=` may be omitted; `annotation` defaults to the first of `cell_type`, `annotation`, `cluster`, `cluster_id` found in `obs`. SpatialData extras need Python 3.11+ (`pip install -e "python[plot,spatialdata]"`).

## Batch script

```bash
python scripts/celladmix_batch.py --input data --output out \
  --annotation annotations/annotation.csv.gz --annotation-col merged_annotation \
  --score auto --report --threads 10
```

Stages: dataset and store, fit (`--rank auto|N`, `--nmf-variant ls_nmf|invsqrt_kl|kl|sqrt_kl`, `--nmf-runs auto|N`, `--nmf-init`, `--seed`), scoring (`--score auto|membrane|bridge`; `auto` uses membrane when a stain is discoverable; `coherence` is accepted by the parser but raises `NotImplementedError`), rules (`--p-thresh 0.1`, `--adjust-p`, `--targets A,B`), correction (`--ensemble 10`, `--vote 0.3`, `--no-correct`), reports (`--report`, `--report-file`, `--top-genes 10`), clustering when no annotation is given (`--auto-annotate`, `--cluster-*`, `--annotation-report`). Outputs in `--output`: `batch_summary.json/.csv`, `fit_summary.csv`, `factor_top_genes.csv`, `score_summary.csv`, `score_rules.csv`, `correction_summary.csv`, `report_minimal.html`, and with auto-annotation `annotation_clusters.csv`, `annotation_cluster_markers.csv`, `annotation_report.html`. Tabular options exist for parity with the R script but fail with an explicit message; use `scripts/celladmix_batch.R` for tables. Full option list in `docs/python-batch.md`.

## Tests

`python -m pytest python/tests` (unittest-style modules covering caching, native check, ensemble, audit, generative model, collections, stains).
