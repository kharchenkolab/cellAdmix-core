# NMF Restart Stability

When a fit uses multiple NMF restarts (`nmf_n_runs > 1`), cellAdmix-core
reports a per-factor stability diagnostic alongside the selected
factorization. Stability answers two questions:

- **Which factors can be trusted?** A factor that is re-found by independent
  random restarts reflects structure in the data; a factor that appears in
  only one restart is a seed artifact and should not drive scoring or
  correction decisions.
- **Is the rank appropriate?** The number of stable factors levels off at the
  number of factors the data can resolve. Raising the rank past that point
  adds unstable factors rather than new structure.

## Definition

Stability is measured on **gene-ownership profiles**: the selected run's
loading matrix H is renormalized so that each gene's loadings sum to one
across factors. A factor is thereby described by which genes it owns rather
than by raw loading magnitudes. This cancels gene abundance, so two unrelated
factors score near zero even when all factors share a strong
abundance backbone, and the values are comparable across `nmf_variant`
settings and across ranks.

Factors of the selected run are paired one-to-one with the factors of each
other restart by a maximum-similarity (Hungarian) assignment on the Pearson
correlations of ownership profiles. The per-factor stability is the mean
matched correlation across the comparison restarts:

- values near 1 mean the factor is re-found with the same gene ownership in
  every restart;
- values near 0 mean nothing specific to this factor reproduces;
- **0.3** is the default threshold for counting a factor as stable, with
  values above roughly 0.6 indicating solidly re-found factors. Two genuinely
  similar factors (for example, two states of one cell type) do not lower each
  other's score: each is compared against its own matched partner.

Note that stability certifies reproducibility across restarts, not biological
validity, and restarts share the same training subsample, so it does not
capture sensitivity to training-data resampling.

## Restart recipe

Restarts differ only in their random seed. When cluster initialization is
active (`nmf_init = "auto"` with an annotation available, or
`nmf_init = "cluster"`), only the first restart is cluster-initialized; the
remaining restarts always use fully random initialization. The best final
objective across all restarts is still selected — the cluster-initialized
candidate frequently wins — but stability is averaged over the independent
random restarts only. This keeps the diagnostic meaningful: it reads as "do
independent random starts re-find the selected solution?", not "do
identically-seeded runs stay where they were put?".

Stability requires at least two comparison restarts to be informative;
`nmf_n_runs` of 5–10 is recommended when the diagnostic matters. With
`nmf_n_runs = 1` the reported stability is a placeholder of 1.0 and
`stability_comparison_runs` is 0.

## Reported fields

The run manifest's `nmf_diagnostics` block records:

| field | meaning |
|---|---|
| `selected_factor_stability` | per-factor matched ownership correlation, ordered like the factors |
| `stable_factor_count` | number of factors at or above `stability_threshold` |
| `stability_threshold` | threshold used for the count (default 0.3) |
| `stability_comparison_runs` | how many independent restarts the average uses |
| `candidate_best_match_correlations` | per-restart mean matched correlation against the selected run |
| `stability_metric` | `"ownership_matched"`; manifests written by older versions load as `"best_match_legacy"` |

Values recorded under the legacy metric (best-match Pearson on raw loadings,
without one-to-one matching) are not comparable to ownership-matched values:
the legacy kernel scores unrelated factors well above zero for
abundance-dominated variants, and inflates comparisons between `nmf_variant`
settings. Refit with the current version when stability matters.

## Usage

R:

```r
fit <- ds$fit(nmf_n_runs = 10)
fit$plot_stability()          # per-factor stability vs. factor importance
fit$stability()               # the underlying table
```

Python:

```python
fit = ds.fit(nmf_n_runs=10)
fit.plot_stability()
fit.manifest["nmf_diagnostics"]["selected_factor_stability"]
```

The verbose fit log summarizes the same information, for example
`stable_factors=8/24` together with the restart objective spread. A large gap
between the rank and the stable-factor count is the primary signal that the
rank is set higher than the data support; factors below the threshold should
be treated as noise when interpreting factors or building correction rules.

Low stability across the board on a large panel usually indicates
undersized NCV neighborhoods rather than an intrinsic property of the data:
spreading few neighborhood draws over thousands of genes leaves almost no
gene co-occurrence signal per neighborhood. The automatic `ncv_k` default
scales the neighborhood with the panel size to avoid this; see the fitting
documentation in [r-bindings.md](r-bindings.md).
