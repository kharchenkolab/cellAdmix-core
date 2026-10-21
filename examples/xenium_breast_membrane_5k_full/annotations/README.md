# Cell-type annotation of the Xenium breast 5K example

Three tables over the 470,549 cells of the bundle, keyed by the Xenium
`cell_id`.

- `annotation.csv.gz`: per cell, the cluster id (`cluster`), a
  descriptive cluster label (`cluster_label`), the cell-type label used
  by the notebooks (`merged_annotation`) and the markers behind it
  (`marker_rationale`).
- `cluster_annotation.csv`: the same label and rationale once per
  cluster.
- `domain_annotation.csv.gz`: per cell, the tissue domain and region from
  `identify_domains()` (scales 25, 75 and 200 µm, ten domains), with the
  domain named after the lineage-labelled cells it contains
  (tumor_epithelial, immune_stroma, immune_rich, mixed_region).

The labels come from a cellAdmixCore clustering of the bundle at
resolution 3 (cells with at least 30 molecules and 15 genes,
15-nearest-neighbor graph), with the lymphoid clusters sub-clustered
separately to separate T cells, B cells and plasma cells; clusters were
labelled from their top detection-AUC markers and the per-cluster rates
of canonical lineage markers. Cluster ids of the form `L<n>` are the
lymphoid sub-clusters. Labels and counts: Epithelial 144,214; Ambiguous
/ low-quality 76,811; T cell 57,491; Inflammatory fibroblast 54,230;
Matrix fibroblast 31,849; B cell 24,723; Macrophage 24,603; Monocyte /
inflammatory myeloid 23,979; Endothelial 15,539; Plasma cell 12,886;
Mature dendritic cell 4,224.
