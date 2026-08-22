Generative Admixture Correction in Brief: Xenium Pancreas
================

The generative model corrects admixture by decomposing each cell’s
counts into the cell’s own expression, contamination from each
neighboring cell type, ambient background, and expression induced by the
neighborhood — removing the contamination and ambient shares while
keeping the induced ones. This page shows the shortest path from data to
corrected counts; the [full example](pancreas_generative.md)
characterizes what the model measures along the way, and
[docs/generative.md](../../docs/generative.md) describes the model
itself.

``` r
annotation <- read.csv(file.path("annotations", "annotation.csv.gz"))
cell_annotation <- setNames(annotation$merged_annotation, annotation$cell_id)

ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)
fit <- ds$fit(nmf_variant = "invsqrt_kl")
```

    ## Reusing cached run fit_manual_rank9_invsqrt_kl (parameters match)

The audit measures the admixture; the generative correction removes it:

``` r
audit <- fit$audit_admixture()
```

    ## Excluded 23 likely induced genes from marker panels (exposure-linked excess far above the source-profile expectation): ACTG2, ADAMTS1, APCDD1, APOLD1, BASP1, C5orf46, CA4, CAVIN1

``` r
correction <- audit$correct_generative(num_threads = 8)
correction
```

    ## cellAdmix generative correction
    ##   pairs: 39 
    ##   removed molecules (expected): 1,082,659 
    ##   induced genes retained: 28

The corrected gene-by-cell matrix drops into any downstream analysis:

``` r
corrected_counts <- correction$counts()
dim(corrected_counts)
```

    ## [1]    377 140335

The correction keeps the genes whose exposure-linked excess reflects the
cells’ own response to their neighbors rather than transferred material:

``` r
head(correction$induced[order(-correction$induced$excess), ], 5)
```

    ##                 source                  target   gene     excess   expected
    ## 4  Exocrine epithelial Ductal/tumor epithelial   CFTR 110866.738 23889.0283
    ## 6  Exocrine epithelial Ductal/tumor epithelial  PROX1  11709.747  3526.6096
    ## 8     Fibroblast / CAF Ductal/tumor epithelial SEMA3C   4181.181   540.1331
    ## 25    Fibroblast / CAF                  Immune  PMP22   3183.480  1029.5018
    ## 9     Fibroblast / CAF Ductal/tumor epithelial  BASP1   2750.508  1179.0843
    ##            z     fold
    ## 4  24.154536 4.640906
    ## 6  15.121410 3.320398
    ## 8  22.477063 7.741020
    ## 25 10.279826 3.092253
    ## 9   8.124244 2.332749

And it is verified like any other correction, against the audit’s own
measurements:

``` r
report <- audit$evaluate(correction)
report$summary()
```

    ## $detected_pairs
    ## [1] 39
    ## 
    ## $estimated_admixed_molecules
    ## [1] 1574854
    ## 
    ## $leakage_removed_overall
    ## [1] 0.9405623
    ## 
    ## $median_pair_sensitivity
    ## [1] 0.9587951
    ## 
    ## $own_marker_false_removal
    ## [1] 0
    ## 
    ## $worst_false_removal
    ## [1] "Ductal/tumor epithelial"
