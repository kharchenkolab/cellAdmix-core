# Figure helpers for the generative model

`generative_plots.R` and `generative_plots.py` hold the analysis and figure
functions used by the generative example notebooks: the designs of
`docs/generative.md` (the composition of a cell type's content by exposure,
the transfer-proportionality screen with both expectations, the per-gene
rate contrast with exact intervals, the neighborhood-change volcano) and
the supporting views (burden by type, contamination maps, example-cell
composition and gene-origin rendering, induced-gene retention, state-space
UMAPs). They take the package's public objects and return tables or plots,
and they live here rather than in the package until the model's place in
it is settled.

Inputs, by name, in both languages:

- `before`, `after`: gene x cell count matrices with cell ids, from
  `ds$counts()` / `ds.counts()` and `correction$counts()` /
  `correction.counts()` (Python functions also take the gene and cell lists
  those calls return).
- `annotation`: cell types by cell id (named vector in R, Series in Python).
- `cells_xy`: the cell table from `ds$cells()` / `ds.cells()`.
- `exposure`: source-type neighbor counts per cell id, from
  `generative_exposure(cells_xy, annotation, source)`.
- `model`: a fitted `audit$fit_generative()` / `audit.fit_generative()`
  object; `nmf_fit`: the NMF fit, needed only for the molecule-level
  cell rendering.

`generative_induced_changes()` is the general screen of neighborhood-
associated expression changes on the corrected counts (every gene of a
target type, all sources' exposures fitted jointly, quasi-Poisson errors,
with the removed-fraction and transfer-profile guards); it returns one row
per gene and source, and `generative_plot_volcano()` draws one source.

Usage from a notebook in `examples/<dataset>/`:

```r
source("../../analysis/generative_model/tools/generative_plots.R")
```

```python
import sys; sys.path.insert(0, "../../analysis/generative_model/tools")
from generative_plots import *
```
