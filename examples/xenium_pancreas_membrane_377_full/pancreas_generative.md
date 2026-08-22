Generative Admixture Correction on Xenium Pancreas
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

Molecules leak between neighboring cells, so a cell’s counts mix its own
expression with material from the cells around it. But not every count
that rises near a neighbor is leakage: cells also change their own
transcription in response to their neighborhood, and a correction that
flattens every neighbor-linked gradient erases exactly that biology. The
generative model separates the two: it decomposes each cell’s counts
into the cell’s own expression programs, contamination from each
neighboring cell type, ambient background, and neighborhood-induced
expression, then removes the contamination and ambient shares while
keeping the induced ones ([docs/generative.md](../../docs/generative.md)
describes the model and its validation; [the minimal
version](pancreas_generative_minimal.md) of this page is four calls
long).

This example characterizes both sides on the pancreas dataset: how much
material moves between which cell types, and which genes the model
recognizes as the cells’ own induced expression.

``` r
annotation <- read.csv(file.path("annotations", "annotation.csv.gz"))
cell_annotation <- setNames(annotation$merged_annotation, annotation$cell_id)

ds <- cellAdmix("data", output_dir = "out", annotation = cell_annotation)
fit <- ds$fit(nmf_variant = "invsqrt_kl")
```

    ## Reusing cached run fit_manual_rank9_invsqrt_kl (parameters match)

## The admixture pattern between cell types

The audit measures, for every ordered pair of cell types, how many
molecules leaked from the source type into cells of the target type,
using the rise of source-marker content with the number of source-type
neighbors:

``` r
audit <- fit$audit_admixture()
```

    ## Excluded 23 likely induced genes from marker panels (exposure-linked excess far above the source-profile expectation): ACTG2, ADAMTS1, APCDD1, APOLD1, BASP1, C5orf46, CA4, CAVIN1

``` r
audit$plot_map()
```

<img src="pancreas_generative_files/figure-gfm/audit-map-1.png" alt="" width="624" style="display: block; margin: auto;" />

The magnitudes are far from uniform. Sorting the detected pairs by their
estimated admixture rate — the fraction of the target type’s molecules
that originated in the source:

``` r
pairs <- audit$pairs(detected_only = TRUE)
head(pairs[order(-pairs$rate),
  c("source", "target", "admixed_molecules", "rate")], 8)
```

    ##                 source                  target admixed_molecules      rate
    ## 21 Exocrine epithelial             Endothelial             67068 0.2107178
    ## 19 Exocrine epithelial Ductal/tumor epithelial            382869 0.1995796
    ## 29    Fibroblast / CAF                  Immune             97795 0.1928331
    ## 35              Immune        Fibroblast / CAF            206179 0.1886504
    ## 27    Fibroblast / CAF             Endothelial             44954 0.1412410
    ## 23 Exocrine epithelial                  Immune             69990 0.1380071
    ## 36              Immune        Mural / pericyte              6666 0.1173330
    ## 30    Fibroblast / CAF        Mural / pericyte              5979 0.1052369

Aggregating over targets shows which cell types shed the most material —
on this tissue, the exocrine compartment dominates as a source:

``` r
by_source <- aggregate(admixed_molecules ~ source, pairs, sum)
ggplot(by_source, aes(reorder(source, admixed_molecules),
    admixed_molecules / 1e3)) +
  geom_col(fill = "#c0392b", alpha = 0.85) +
  coord_flip() +
  labs(x = NULL, y = "estimated admixed molecules shed (thousands)") +
  theme_classic(base_size = 11)
```

<img src="pancreas_generative_files/figure-gfm/audit-sources-1.png" alt="" width="576" style="display: block; margin: auto;" />

## Correcting with the generative model

One call fits the model on all detected pairs (about twenty seconds on
this dataset):

``` r
correction <- audit$correct_generative(num_threads = 8)
correction
```

    ## cellAdmix generative correction
    ##   pairs: 39 
    ##   removed molecules (expected): 1,082,659 
    ##   induced genes retained: 28

The model’s central quantity is the per-cell contamination fraction: how
much of each cell’s content arrived from each source. Its dependence on
the number of source neighbors, for the largest pair:

``` r
comp <- correction$composition("Exocrine epithelial",
  "Ductal/tumor epithelial")
cells_xy <- fit$cell_factors()
expo <- setNames(
  cellAdmixCore:::.celladmix_source_exposure_counts(
    cells_xy, cell_annotation, 15L)$counts[, "Exocrine epithelial"],
  as.character(cells_xy$cell_id))
comp$exposure <- pmin(expo[comp$cell_id], 4)
ggplot(comp, aes(factor(exposure), contamination)) +
  geom_boxplot(outlier.size = 0.3, fill = "#f4a582") +
  labs(x = "exocrine cells among 15 nearest neighbors",
    y = "fraction of the ductal cell's molecules\nattributed to exocrine contamination") +
  theme_classic(base_size = 11)
```

<img src="pancreas_generative_files/figure-gfm/composition-1.png" alt="" width="576" style="display: block; margin: auto;" />

The same fractions drawn in space show the admixture concentrating at
the tissue interfaces:

``` r
xy <- cells_xy[match(comp$cell_id, as.character(cells_xy$cell_id)), ]
ggplot(cbind(comp, x = xy$x, y = xy$y), aes(x, y,
    color = pmin(contamination, 0.6))) +
  geom_point(size = 0.25) +
  scale_color_viridis_c(option = "inferno",
    name = "exocrine share\nof ductal-cell\ncontent") +
  coord_equal() +
  labs(x = NULL, y = NULL,
    title = "Exocrine contamination of ductal cells") +
  theme_void(base_size = 11)
```

<img src="pancreas_generative_files/figure-gfm/composition-map-1.png" alt="" width="672" style="display: block; margin: auto;" />

## The induction pattern

The model retains genes whose exposure-linked excess exceeds what
transfer from the source can deliver — expression the target cells
switched on themselves. On this tissue the retained set is led by the
duct-cell program that ductal cells activate at acinar interfaces:

``` r
induced <- correction$induced
head(induced[order(-induced$excess),
  c("source", "target", "gene", "excess", "fold", "z")], 8)
```

    ##                 source                  target    gene     excess     fold
    ## 4  Exocrine epithelial Ductal/tumor epithelial    CFTR 110866.738 4.640906
    ## 6  Exocrine epithelial Ductal/tumor epithelial   PROX1  11709.747 3.320398
    ## 8     Fibroblast / CAF Ductal/tumor epithelial  SEMA3C   4181.181 7.741020
    ## 25    Fibroblast / CAF                  Immune   PMP22   3183.480 3.092253
    ## 9     Fibroblast / CAF Ductal/tumor epithelial   BASP1   2750.508 2.332749
    ## 19    Mural / pericyte             Endothelial ADAMTS1   2045.614 3.725193
    ## 7     Fibroblast / CAF Ductal/tumor epithelial  APCDD1   2035.355 3.465643
    ## 5  Exocrine epithelial Ductal/tumor epithelial     CA4   1681.018 4.159392
    ##            z
    ## 4  24.154536
    ## 6  15.121410
    ## 8  22.477063
    ## 25 10.279826
    ## 9   8.124244
    ## 19 15.015948
    ## 7  11.999946
    ## 5  16.955261

The distinction is visible gene by gene. AMY2A is the largest exocrine
transfer channel: its rise with exocrine exposure in ductal cells is
transferred material, and the correction flattens it to the unexposed
level. CFTR carries an induced duct-cell program: its rise is the ductal
cells’ own expression, and the correction preserves it while removing
only the transferred share:

``` r
before <- fit$counts()
after <- correction$counts()
t_cells <- comp$cell_id
tot <- Matrix::colSums(before[, t_cells])
bins <- pmin(expo[t_cells], 3)
rate_by_bin <- function(m, gene) {
  tapply(m[gene, t_cells], bins, sum) / tapply(tot, bins, sum) * 1e3
}
df <- do.call(rbind, lapply(c("AMY2A", "CFTR"), function(g) {
  rbind(
    data.frame(gene = g, counts = "observed", bin = 0:3,
      rate = as.numeric(rate_by_bin(before, g))),
    data.frame(gene = g, counts = "corrected", bin = 0:3,
      rate = as.numeric(rate_by_bin(after, g))))
}))
ggplot(df, aes(bin, rate, color = counts)) +
  geom_line() + geom_point() +
  facet_wrap(~gene, scales = "free_y") +
  scale_color_manual(values = c(observed = "#c0392b",
    corrected = "#2980b9")) +
  labs(x = "exocrine cells among 15 nearest neighbors",
    y = "rate in ductal cells (per 1,000 molecules)") +
  theme_classic(base_size = 11)
```

<img src="pancreas_generative_files/figure-gfm/gene-contrast-1.png" alt="" width="816" style="display: block; margin: auto;" />

## Verification

The correction is verified against the audit’s own measurements, like
any other correction in the package:

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

``` r
report$plot_cleanup()
```

<img src="pancreas_generative_files/figure-gfm/verify-1.png" alt="" width="768" style="display: block; margin: auto;" />

Two limits are worth keeping in mind when reading the results. The model
never removes a target type’s own marker genes (their contamination is
indistinguishable from native expression and is deliberately left in
place), so shared-gene admixture is under-removed by construction. And
induced expression is only recognizable when it is disproportionate to
the source’s transfer profile: induction of a source’s own top marker
genes cannot be told apart from admixture at count level and is removed
with it ([docs/generative.md](../../docs/generative.md), the
identifiability limit).
