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

- [Pancreas quickstart](../examples/xenium_pancreas_membrane_377_full/pancreas_quickstart.ipynb)
- [Minimal CosMx NSCLC tutorial](../examples/cosmx_nsclc_giotto/celladmix_cosmx_minimal.ipynb)
- [Xenium pancreas membrane/bridge scoring tutorial](../examples/xenium_pancreas_membrane_377_full/pancreas_membrane_scoring_clean.ipynb)
- [Seurat Xenium integration tutorial](../examples/xenium_pancreas_membrane_377_full/pancreas_seurat_integration.ipynb)
- [Xenium breast 5K membrane scoring tutorial](../examples/xenium_breast_membrane_5k_full/breast_5k_membrane_scoring.ipynb)
- [Python Xenium standalone tutorial](../examples/xenium_pancreas_membrane_377_full/pancreas_python_standalone.ipynb)
- [Python SpatialData integration tutorial](../examples/xenium_pancreas_membrane_377_full/pancreas_spatialdata_integration.ipynb)
- [Python Xenium breast 5K standalone tutorial](../examples/xenium_breast_membrane_5k_full/breast_5k_python_standalone.ipynb)
