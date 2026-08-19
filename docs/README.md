# Documentation

This folder contains user-facing documentation for `cellAdmix-core`.

## User Documentation

- [Installation](install.md): system dependencies and R/Python binding
  installation.
- [R bindings](r-bindings.md): high-level R workflow, common operations, batch
  processing, and example notebooks.
- [Inputs](inputs.md): supported input formats and input-specific options.
- [Scoring methods](scoring_methods.md): bridge, membrane, and coherence
  scoring definitions.
- [NMF restart stability](nmf_stability.md): the per-factor stability
  diagnostic, restart recipe, and how to use it for factor trust and rank
  selection.
- [Cleanup benchmarks](benchmarks.md): the neighbor-benchmark methodology
  and the factorization/scoring comparisons behind the recommended
  workflow settings.
- [Python bindings](python-bindings.md): Python API and SpatialData integration.
- [Python batch processing](python-batch.md): CLI workflow for Xenium batch
  processing through the Python bindings.

## Examples

- [Pancreas quickstart](../examples/xenium_pancreas_membrane_377_full/pancreas_quickstart.md)
- [Minimal CosMx NSCLC tutorial](../examples/cosmx_nsclc_giotto/celladmix_cosmx_minimal.md)
- [Detailed Xenium pancreas tutorial](../examples/xenium_pancreas_membrane_377_full/pancreas_membrane_scoring_clean.md)
- [Seurat Xenium integration tutorial](../examples/xenium_pancreas_membrane_377_full/pancreas_seurat_integration.md)
- [Xenium breast 5K membrane scoring tutorial](../examples/xenium_breast_membrane_5k_full/breast_5k_membrane_scoring.md)
- [Xenium breast 5K medium-crop tutorial](../examples/xenium_breast_membrane_5k_full/breast_5k_membrane_scoring_medium_crop.md)
- [Standalone Xenium pancreas tutorial](../examples/xenium_pancreas_membrane_377_full/pancreas_python_standalone.ipynb)
- [SpatialData Xenium integration tutorial](../examples/xenium_pancreas_membrane_377_full/pancreas_spatialdata_integration.ipynb)
- [Standalone Xenium breast 5K tutorial](../examples/xenium_breast_membrane_5k_full/breast_5k_python_standalone.ipynb)
