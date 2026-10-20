Generative Admixture Correction on CosMx NSCLC
================

-   <a href="#the-admixture-pattern-between-cell-types"
    id="toc-the-admixture-pattern-between-cell-types">The admixture pattern
    between cell types</a>
-   <a href="#correcting-with-the-generative-model"
    id="toc-correcting-with-the-generative-model">Correcting with the
    generative model</a>
-   <a href="#the-induction-pattern" id="toc-the-induction-pattern">The
    induction pattern</a>
-   <a href="#verification" id="toc-verification">Verification</a>

The generative model corrects admixture by decomposing each cell’s
counts into the cell’s own expression, contamination from each
neighboring cell type, ambient background, and expression the cell
switched on in response to its neighborhood; the contamination and
ambient shares are removed, the induced expression retained
([docs/generative.md](../../docs/generative.md)). This page runs it on
the CosMx NSCLC dataset and characterizes both sides of the separation:
the pattern of material moving between cell types, and the induced
programs the correction preserves. One platform note applies throughout:
this CosMx export carries no nucleus-overlap flag, so the model uses
whole-cell rather than cytoplasmic source profiles.

``` r
cell_meta <- read.csv(file.path("prepared", "cell_metadata_all.csv.gz"),
  stringsAsFactors = FALSE)
cell_annotation <- setNames(cell_meta$cell_type_coarse, cell_meta$cell)

ds <- cellAdmix(file.path("prepared", "molecules_all.csv.gz"),
  output_dir = "out", annotation = cell_annotation)
```

## The admixture pattern between cell types

The audit is computed straight from the dataset — counts, cell
positions, and the annotation; no factorization is involved:

``` r
audit <- ds$audit_admixture()
```

    ## Excluded 72 likely induced genes from marker panels (exposure-linked excess far above the source-profile expectation): ADIRF, B2M, C5AR2, CAV1, CCL21, CCL3, CCL3L3, CCL4

``` r
audit$plot_map()
```

<img src="nsclc_generative_files/figure-gfm/audit-map-1.png" alt="" width="624" style="display: block; margin: auto;" />

``` r
pairs <- audit$pairs(detected_only = TRUE)
head(pairs[order(-pairs$rate),
  c("source", "target", "admixed_molecules", "rate")], 8)
```

    ##          source       target admixed_molecules       rate
    ## 13   fibroblast immune other           1325145 0.12202791
    ## 17 immune other   fibroblast            474902 0.12076415
    ## 16 immune other  endothelial            214239 0.09316072
    ## 22   macrophage immune other            747404 0.06882583
    ## 11   fibroblast  endothelial            107226 0.04662678
    ## 12   fibroblast   epithelial             58450 0.03815173
    ## 25    malignant   epithelial             51489 0.03360857
    ## 21   macrophage   fibroblast            130021 0.03306331

Aggregated over targets, the shedding is dominated by the tumor and
myeloid compartments:

``` r
by_source <- aggregate(admixed_molecules ~ source, pairs, sum)
ggplot(by_source, aes(reorder(source, admixed_molecules),
    admixed_molecules / 1e3)) +
  geom_col(fill = "#c0392b", alpha = 0.85) +
  coord_flip() +
  labs(x = NULL, y = "estimated admixed molecules shed (thousands)") +
  theme_classic(base_size = 11)
```

<img src="nsclc_generative_files/figure-gfm/audit-sources-1.png" alt="" width="576" style="display: block; margin: auto;" />

## Correcting with the generative model

The model is fitted on the audit’s detected pairs. The NMF fit enters
only to initialize the model’s expression programs — the recommended
configuration — and the correction derives from the fitted model:

``` r
nmf_fit <- ds$fit()
```

    ## Reusing cached run fit_manual_rank8_ls_nmf (parameters match)

``` r
model <- audit$fit_generative(init = nmf_fit, num_threads = 8)
model
```

    ## cellAdmix generative model
    ##   pairs: 27 
    ##   induced genes retained: 88 
    ##   removed molecules (expected): 3,390,467 
    ##   initialization: nmf_factors

``` r
correction <- model$correct()
```

The per-cell decomposition of the largest pair’s target type, averaged
by exposure to the source, and the pair’s contamination fraction in
space (the figure helpers are those of the model documentation, in
[analysis/generative_model/tools/generative_plots.R](../../analysis/generative_model/tools/generative_plots.R)):

``` r
top_pair <- pairs[order(-pairs$admixed_molecules), ][1, ]
cells_xy <- ds$cells()
expo <- generative_exposure(cells_xy, cell_annotation, top_pair$source)
generative_plot_composition(model, top_pair$target, expo,
  sprintf("%s cells among 15 nearest neighbors", top_pair$source))
```

<img src="nsclc_generative_files/figure-gfm/composition-1.png" alt="" width="576" style="display: block; margin: auto;" />

``` r
generative_plot_composition_map(model, top_pair$source, top_pair$target,
  cells_xy, cap = 0.5)
```

<img src="nsclc_generative_files/figure-gfm/composition-map-1.png" alt="" width="720" style="display: block; margin: auto;" />

## The induction pattern

On this tissue the retained induced set reads as the classic stress,
immediate-early and chemokine programs of cells at tissue interfaces:

``` r
induced <- model$induced
head(induced[order(-induced$excess),
  c("source", "target", "gene", "excess", "fold", "z")], 10)
```

    ##         source       target   gene    excess     fold         z
    ## 36 endothelial immune other HSPA1A 10822.071 4.115143 17.490986
    ## 37 endothelial immune other HSPA1B  9197.786 4.069398 16.710038
    ## 43  macrophage immune other   GLUL  7484.991 2.531707  9.664492
    ## 46  macrophage immune other   SRGN  7152.531 5.427592 20.806496
    ## 41  macrophage immune other   NPPC  6637.024 2.434147  8.265181
    ## 84  macrophage    malignant  DUSP5  5885.255 7.805129 24.148237
    ## 20 endothelial   fibroblast COL4A2  5441.595 4.089366 15.346591
    ## 5   macrophage  endothelial  DUSP5  5257.256 6.922909 16.729099
    ## 28  macrophage   fibroblast  DUSP5  5209.542 3.296008 10.064316
    ## 16 endothelial   fibroblast COL4A1  4278.497 3.739442 14.421808

The contrast between a transferred and an induced gene, for the pair
with the strongest induced signal — the observed rise of the transferred
gene is flattened to the unexposed level, while the induced gene keeps
its gradient:

``` r
ind_top <- induced[order(-induced$excess), ][1, ]
mk <- audit$markers(ind_top$source, ind_top$target)
transfer_gene <- setdiff(mk$pool, induced$gene)[1]
before <- ds$counts()
after <- correction$counts()
t_cells <- model$composition(ind_top$source, ind_top$target)$cell_id
expo_t <- generative_exposure(cells_xy, cell_annotation, ind_top$source)
generative_plot_gene_contrast(before, after, t_cells, expo_t,
  c(transfer_gene, ind_top$gene),
  sprintf("%s cells among 15 nearest neighbors", ind_top$source),
  sprintf("rate in %s cells (per 1,000 molecules)", ind_top$target))
```

<img src="nsclc_generative_files/figure-gfm/gene-contrast-1.png" alt="" width="816" style="display: block; margin: auto;" />

Beyond the genes the model retains by its own screen, the corrected
counts support the general question of which genes — of any ownership,
in either direction — change in the target cells with the source’s
proximity, tested jointly across all sources
([docs/generative.md](../../docs/generative.md), Figure 4). Filled
points are genes with essentially none of their content removed as
transfer and no meaningful share of the source’s transfer profile; open
grey points are entangled with the transfer split:

``` r
changes <- generative_induced_changes(before, after, model, cells_xy,
  cell_annotation, ind_top$target)
generative_plot_volcano(changes, ind_top$source)
```

<img src="nsclc_generative_files/figure-gfm/volcano-1.png" alt="" width="652.8" style="display: block; margin: auto;" />

## Verification

``` r
report <- audit$evaluate(correction)
report$summary()
```

    ## $detected_pairs
    ## [1] 27
    ## 
    ## $estimated_admixed_molecules
    ## [1] 4381180
    ## 
    ## $leakage_removed_overall
    ## [1] 0.9688949
    ## 
    ## $median_pair_sensitivity
    ## [1] 0.997736
    ## 
    ## $own_marker_false_removal
    ## [1] 0
    ## 
    ## $worst_false_removal
    ## [1] "endothelial"

``` r
report$plot_cleanup()
```

<img src="nsclc_generative_files/figure-gfm/verify-1.png" alt="" width="768" style="display: block; margin: auto;" />

The model never removes a target type’s own marker genes, and induced
expression is only recognizable when disproportionate to the source’s
transfer profile — see the identifiability limit in
[docs/generative.md](../../docs/generative.md).
