.celladmix_probe_xenium <- function(path) {
  .Call("_cellAdmixCore_celladmix_probe_xenium", path, PACKAGE = "cellAdmixCore")
}

.celladmix_probe_tabular <- function(molecules_path, segmentation_mask_path = NULL) {
  .Call(
    "_cellAdmixCore_celladmix_probe_tabular",
    as.character(molecules_path[[1]]),
    if (is.null(segmentation_mask_path)) NULL else as.character(segmentation_mask_path[[1]]),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_read_input_store <- function(store_dir) {
  .Call("_cellAdmixCore_celladmix_read_input_store", as.character(store_dir[[1]]), PACKAGE = "cellAdmixCore")
}

.celladmix_collect_input_store_counts <- function(store_dir) {
  .Call("_cellAdmixCore_celladmix_collect_input_store_counts", as.character(store_dir[[1]]), PACKAGE = "cellAdmixCore")
}

.celladmix_identify_domains_store <- function(
    store_dir,
    cell_types = NULL,
    scales = c(15L, 50L, 150L),
    transform = c("clr", "hellinger", "none"),
    include_self = TRUE,
    self_weight = 1,
    include_center_type = FALSE,
    center_type_weight = 0.25,
    n_domains = 8L,
    smooth_k = 20L,
    smooth_lambda = 0.5,
    max_edge_distance = NA_real_,
    min_region_size = 25L,
    seed = 1L,
    num_threads = 1L
  ) {
  transform <- match.arg(transform)
  .Call(
    "_cellAdmixCore_celladmix_identify_domains_store",
    as.character(store_dir[[1]]),
    if (is.null(cell_types)) NULL else stats::setNames(as.character(cell_types), names(cell_types)),
    as.integer(scales),
    as.character(transform[[1]]),
    isTRUE(include_self),
    as.numeric(self_weight),
    isTRUE(include_center_type),
    as.numeric(center_type_weight),
    as.integer(n_domains),
    as.integer(smooth_k),
    as.numeric(smooth_lambda),
    as.numeric(max_edge_distance),
    as.integer(min_region_size),
    as.integer(seed),
    as.integer(num_threads),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_build_xenium_store <- function(
    path,
    crops = NULL,
    cell_filter = NULL,
    gene_filter = NULL,
    min_qv = NA_real_,
    keep_unassigned = FALSE,
    keep_non_gene = FALSE,
    read_cells = TRUE,
    prefer_parquet = TRUE,
    store_dir,
    materialize_molecules = TRUE,
    force = FALSE,
    num_threads = 1L,
    parquet_row_group_size = 65536L,
    verbose = FALSE
  ) {
  .Call(
    "_cellAdmixCore_celladmix_build_xenium_store",
    path,
    crops,
    if (is.null(cell_filter)) NULL else as.character(cell_filter),
    if (is.null(gene_filter)) NULL else as.character(gene_filter),
    as.numeric(min_qv),
    isTRUE(keep_unassigned),
    isTRUE(keep_non_gene),
    isTRUE(read_cells),
    isTRUE(prefer_parquet),
    as.character(store_dir[[1]]),
    isTRUE(materialize_molecules),
    isTRUE(force),
    as.integer(num_threads),
    as.integer(parquet_row_group_size),
    isTRUE(verbose),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_build_tabular_store <- function(
    molecules_path,
    crops = NULL,
    min_qv = NA_real_,
    keep_unassigned = FALSE,
    x_col = "x",
    y_col = "y",
    z_col = NULL,
    gene_col = "gene",
    qv_col = NULL,
    cell_id_col = NULL,
    segmentation_mask_path = NULL,
    cell_type_col = NULL,
    cell_metadata_path = NULL,
    cell_metadata_cell_id_col = "cell",
    cell_metadata_cell_type_col = NULL,
    sample_id_col = NULL,
    fov_id_col = NULL,
    sample_id = NULL,
    fov_id = NULL,
    store_dir,
    materialize_molecules = TRUE,
    force = FALSE,
    num_threads = 1L,
    parquet_row_group_size = 65536L,
    verbose = FALSE
  ) {
  .Call(
    "_cellAdmixCore_celladmix_build_tabular_store",
    as.character(molecules_path[[1]]),
    crops,
    as.numeric(min_qv),
    isTRUE(keep_unassigned),
    as.character(x_col[[1]]),
    as.character(y_col[[1]]),
    if (is.null(z_col)) NULL else as.character(z_col[[1]]),
    as.character(gene_col[[1]]),
    if (is.null(qv_col)) NULL else as.character(qv_col[[1]]),
    if (is.null(cell_id_col)) NULL else as.character(cell_id_col[[1]]),
    if (is.null(segmentation_mask_path)) NULL else as.character(segmentation_mask_path[[1]]),
    if (is.null(cell_type_col)) NULL else as.character(cell_type_col[[1]]),
    if (is.null(cell_metadata_path)) NULL else as.character(cell_metadata_path[[1]]),
    as.character(cell_metadata_cell_id_col[[1]]),
    if (is.null(cell_metadata_cell_type_col)) NULL else as.character(cell_metadata_cell_type_col[[1]]),
    if (is.null(sample_id_col)) NULL else as.character(sample_id_col[[1]]),
    if (is.null(fov_id_col)) NULL else as.character(fov_id_col[[1]]),
    if (is.null(sample_id)) NULL else as.character(sample_id[[1]]),
    if (is.null(fov_id)) NULL else as.character(fov_id[[1]]),
    as.character(store_dir[[1]]),
    isTRUE(materialize_molecules),
    isTRUE(force),
    as.integer(num_threads),
    as.integer(parquet_row_group_size),
    isTRUE(verbose),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_cluster_store <- function(
    store_dir,
    min_molecules = 10L,
    min_genes = 5L,
    cells_max = NA_integer_,
    n_variable_genes = 1000L,
    pca_dims = 30L,
    graph_k = 15L,
    cluster_resolution = 1,
    compute_umap = TRUE,
    umap_neighbors = 15L,
    umap_epochs = 200L,
    num_threads = 1L,
    umap_parallel_optimization = TRUE,
    normalization_scale = 5000,
    seed = 1L,
    clusters_out,
    embedding_out,
    verbose = FALSE
  ) {
  .Call(
    "_cellAdmixCore_celladmix_cluster_store",
    as.character(store_dir[[1]]),
    as.integer(min_molecules),
    as.integer(min_genes),
    as.integer(cells_max),
    as.integer(n_variable_genes),
    as.integer(pca_dims),
    as.integer(graph_k),
    as.numeric(cluster_resolution),
    isTRUE(compute_umap),
    as.integer(umap_neighbors),
    as.integer(umap_epochs),
    as.integer(num_threads),
    isTRUE(umap_parallel_optimization),
    as.numeric(normalization_scale),
    as.integer(seed),
    as.character(clusters_out[[1]]),
    as.character(embedding_out[[1]]),
    isTRUE(verbose),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_cluster_run_counts <- function(
    run_dir,
    min_molecules = 10L,
    min_genes = 5L,
    cells_max = NA_integer_,
    n_variable_genes = 1000L,
    pca_dims = 30L,
    graph_k = 15L,
    cluster_resolution = 1,
    compute_umap = TRUE,
    umap_neighbors = 15L,
    umap_epochs = 200L,
    num_threads = 1L,
    umap_parallel_optimization = TRUE,
    normalization_scale = 5000,
    seed = 1L,
    verbose = FALSE
  ) {
  .Call(
    "_cellAdmixCore_celladmix_cluster_run_counts",
    as.character(run_dir[[1]]),
    as.integer(min_molecules),
    as.integer(min_genes),
    as.integer(cells_max),
    as.integer(n_variable_genes),
    as.integer(pca_dims),
    as.integer(graph_k),
    as.numeric(cluster_resolution),
    isTRUE(compute_umap),
    as.integer(umap_neighbors),
    as.integer(umap_epochs),
    as.integer(num_threads),
    isTRUE(umap_parallel_optimization),
    as.numeric(normalization_scale),
    as.integer(seed),
    isTRUE(verbose),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_fit_store_run <- function(
    store_dir,
    analysis_crop = NULL,
    ncv_k = 0L,
    rank = 2L,
    graph_k = 10L,
    same_label_ratio = 5,
    nmf_iterations = 150L,
    nmf_init = "auto",
    nmf_variant = "invsqrt_kl",
    molecule_scoring = "gene_loadings",
    nmf_n_runs = 1L,
    nmf_train_max_rows = NA_integer_,
    nmf_min_molecules = 10L,
    num_threads = 1L,
    training_labels_path = NULL,
    use_cell_type_training = TRUE,
    training_scope_cell_types = NULL,
    seed = 1L,
    out_dir,
    tile_size = 100,
    parquet_row_group_size = 65536L,
    report_ncv_umap = FALSE,
    verbose = FALSE
  ) {
  .Call(
    "_cellAdmixCore_celladmix_fit_store_run",
    as.character(store_dir[[1]]),
    if (is.null(analysis_crop)) NULL else as.character(analysis_crop[[1]]),
    as.integer(ncv_k),
    as.integer(rank),
    as.integer(graph_k),
    as.numeric(same_label_ratio),
    as.integer(nmf_iterations),
    as.character(nmf_init[[1]]),
    as.character(nmf_variant[[1]]),
    as.character(molecule_scoring[[1]]),
    as.integer(nmf_n_runs),
    as.integer(nmf_train_max_rows),
    as.integer(nmf_min_molecules),
    as.integer(num_threads),
    if (is.null(training_labels_path)) NULL else as.character(training_labels_path[[1]]),
    isTRUE(use_cell_type_training),
    if (is.null(training_scope_cell_types)) NULL else as.character(training_scope_cell_types),
    as.integer(seed),
    as.character(out_dir[[1]]),
    as.numeric(tile_size),
    as.integer(parquet_row_group_size),
    isTRUE(report_ncv_umap),
    isTRUE(verbose),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_build_report_data <- function(
    path,
    build_ncv_umap = TRUE,
    force = FALSE,
    umap_neighbors = 15L,
    umap_epochs = 200L,
    normalization_scale = 5000,
    pca_dims = 30L,
    num_threads = 1L,
    seed = 1L,
    verbose = FALSE
  ) {
  .Call(
    "_cellAdmixCore_celladmix_build_report_data",
    as.character(path[[1]]),
    isTRUE(build_ncv_umap),
    isTRUE(force),
    as.integer(umap_neighbors),
    as.integer(umap_epochs),
    as.numeric(normalization_scale),
    as.integer(pca_dims),
    as.integer(num_threads),
    as.integer(seed),
    isTRUE(verbose),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_read_run_manifest <- function(path) {
  .Call("_cellAdmixCore_celladmix_read_run", path, PACKAGE = "cellAdmixCore")
}

.celladmix_collect_run_transcripts <- function(
    path,
    analysis_crop = NULL,
    bbox = NULL,
    sample_n = NA_integer_,
    seed = 1L
  ) {
  .Call(
    "_cellAdmixCore_celladmix_collect_run_transcripts",
    path,
    if (is.null(analysis_crop)) NULL else as.character(analysis_crop[[1]]),
    bbox,
    as.integer(sample_n),
    as.integer(seed),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_collect_run_cells <- function(path) {
  .Call("_cellAdmixCore_celladmix_collect_run_cells", path, PACKAGE = "cellAdmixCore")
}

.celladmix_collect_training_molecules <- function(path) {
  .Call("_cellAdmixCore_celladmix_collect_training_molecules", path, PACKAGE = "cellAdmixCore")
}

.celladmix_collect_run_counts <- function(path, analysis_crop = NULL, bbox = NULL) {
  .Call(
    "_cellAdmixCore_celladmix_collect_run_counts",
    path,
    if (is.null(analysis_crop)) NULL else as.character(analysis_crop[[1]]),
    bbox,
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_collect_run_counts_sparse <- function(path, analysis_crop = NULL, bbox = NULL) {
  .Call(
    "_cellAdmixCore_celladmix_collect_run_counts_sparse",
    path,
    if (is.null(analysis_crop)) NULL else as.character(analysis_crop[[1]]),
    bbox,
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_bridge_scores_run <- function(
    path,
    analysis_crop = NULL,
    bbox = NULL,
    cell_types = NULL,
    candidate_mode = "molecule_global",
    candidate_k = 10L,
    cell_candidate_k = 50L,
    candidate_pairs_per_type_pair = 1000L,
    crossing_k = 20L,
    min_type_pair_contacts = 5L,
    min_factor_molecules = 5L,
    min_pairs = 5L,
    max_cells_per_type_pair = 200L,
    null_iterations = 3L,
    null_max_iterations = 50L,
    cell_candidate_halo = -1,
    null_step = 0.02,
    null_pool_fraction = 0.1,
    num_threads = 1L,
    seed = 1L,
    compute_null = TRUE,
    verbose = FALSE
  ) {
  .Call(
    "_cellAdmixCore_celladmix_bridge_scores_run",
    path,
    if (is.null(analysis_crop)) NULL else as.character(analysis_crop[[1]]),
    bbox,
    if (is.null(cell_types)) NULL else stats::setNames(as.character(cell_types), names(cell_types)),
    as.character(candidate_mode[[1]]),
    as.integer(candidate_k),
    as.integer(cell_candidate_k),
    as.integer(candidate_pairs_per_type_pair),
    as.integer(crossing_k),
    as.integer(min_type_pair_contacts),
    as.integer(min_factor_molecules),
    as.integer(min_pairs),
    as.integer(max_cells_per_type_pair),
    as.integer(null_iterations),
    as.integer(null_max_iterations),
    as.numeric(cell_candidate_halo),
    as.numeric(null_step),
    as.numeric(null_pool_fraction),
    as.integer(num_threads),
    as.integer(seed),
    isTRUE(compute_null),
    isTRUE(verbose),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_membrane_scores_run <- function(
    path,
    analysis_crop = NULL,
    bbox = NULL,
    image_path,
    pixel_size = 1,
    x_offset = 0,
    y_offset = 0,
    epsilon = 1,
    cell_types = NULL,
    cell_candidate_k = 50L,
    candidate_pairs_per_type_pair = 400L,
    cell_candidate_halo = -1,
    min_factor_molecules = 5L,
    min_pairs = 5L,
    max_cells_per_type_pair = 200L,
    control_distance_fraction = 0.025,
    line_samples = 16L,
    num_threads = 1L,
    seed = 1L,
    verbose = FALSE
  ) {
  .Call(
    "_cellAdmixCore_celladmix_membrane_scores_run",
    path,
    if (is.null(analysis_crop)) NULL else as.character(analysis_crop[[1]]),
    bbox,
    as.character(image_path[[1]]),
    as.numeric(pixel_size),
    as.numeric(x_offset),
    as.numeric(y_offset),
    as.numeric(epsilon),
    if (is.null(cell_types)) NULL else stats::setNames(as.character(cell_types), names(cell_types)),
    as.integer(cell_candidate_k),
    as.integer(candidate_pairs_per_type_pair),
    as.numeric(cell_candidate_halo),
    as.integer(min_factor_molecules),
    as.integer(min_pairs),
    as.integer(max_cells_per_type_pair),
    as.numeric(control_distance_fraction),
    as.integer(line_samples),
    as.integer(num_threads),
    as.integer(seed),
    isTRUE(verbose),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_coherence_scores_run <- function(
    path,
    analysis_crop = NULL,
    bbox = NULL,
    cell_types = NULL,
    stain_mode = "none",
    image_path = "",
    pixel_size = 1,
    x_offset = 0,
    y_offset = 0,
    stain_max_pixels = 2048L,
    k_neighbors = 10L,
    min_factor_molecules = 5L,
    min_cells = 5L,
    max_neighbor_distance = -1,
    distance_sigma = -1,
    lambda_coherence = 1,
    beta_margin = 1,
    score_threshold = 0,
    source_pseudocount = 0.5,
    min_source_log_enrichment = 0,
    include_self_source = FALSE,
    compute_null = TRUE,
    null_method = "label_permutation",
    null_iterations = 20L,
    null_exclude_factor = TRUE,
    null_match_nucleus = TRUE,
    null_match_density = TRUE,
    null_nucleus_distance_bins = 3L,
    null_density_bins = 3L,
    seed = 1L,
    normalize_membrane = TRUE,
    membrane_low_quantile = 0.05,
    membrane_high_quantile = 0.995,
    membrane_alpha = 3,
    line_samples = 12L,
    patch_edge_weight_min = 0.1,
    num_threads = 1L,
    verbose = FALSE
  ) {
  .Call(
    "_cellAdmixCore_celladmix_coherence_scores_run",
    path,
    if (is.null(analysis_crop)) NULL else as.character(analysis_crop[[1]]),
    bbox,
    if (is.null(cell_types)) NULL else stats::setNames(as.character(cell_types), names(cell_types)),
    as.character(stain_mode[[1]]),
    as.character(image_path[[1]]),
    as.numeric(pixel_size),
    as.numeric(x_offset),
    as.numeric(y_offset),
    as.integer(stain_max_pixels),
    as.integer(k_neighbors),
    as.integer(min_factor_molecules),
    as.integer(min_cells),
    as.numeric(max_neighbor_distance),
    as.numeric(distance_sigma),
    as.numeric(lambda_coherence),
    as.numeric(beta_margin),
    as.numeric(score_threshold),
    as.numeric(source_pseudocount),
    as.numeric(min_source_log_enrichment),
    isTRUE(include_self_source),
    isTRUE(compute_null),
    as.character(null_method[[1]]),
    as.integer(null_iterations),
    isTRUE(null_exclude_factor),
    isTRUE(null_match_nucleus),
    isTRUE(null_match_density),
    as.integer(null_nucleus_distance_bins),
    as.integer(null_density_bins),
    as.integer(seed),
    isTRUE(normalize_membrane),
    as.numeric(membrane_low_quantile),
    as.numeric(membrane_high_quantile),
    as.numeric(membrane_alpha),
    as.integer(line_samples),
    as.numeric(patch_edge_weight_min),
    as.integer(num_threads),
    isTRUE(verbose),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_read_stain_image_crop <- function(
    image_path,
    pixel_size,
    x_offset = 0,
    y_offset = 0,
    bbox,
    max_pixels = 512L
  ) {
  .Call(
    "_cellAdmixCore_celladmix_read_stain_image_crop",
    as.character(image_path[[1]]),
    as.numeric(pixel_size),
    as.numeric(x_offset),
    as.numeric(y_offset),
    bbox,
    as.integer(max_pixels),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_write_cell_labels <- function(path, labels, row_group_size = 65536L) {
  .Call(
    "_cellAdmixCore_celladmix_write_cell_labels",
    as.character(path[[1]]),
    stats::setNames(as.character(labels), names(labels)),
    as.integer(row_group_size),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_correct_run <- function(path, rules, out_dir, cell_types = NULL) {
  .Call(
    "_cellAdmixCore_celladmix_correct_run",
    path,
    rules,
    as.character(out_dir[[1]]),
    cell_types,
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_collect_correction_summary <- function(path) {
  .Call(
    "_cellAdmixCore_celladmix_collect_correction_summary",
    path,
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_weighted_nmf_matrix <- function(
    x,
    column_weights = NULL,
    rank = 2L,
    max_iterations = 150L,
    init_mode = c("random", "cluster"),
    row_groups = NULL,
    n_runs = 1L,
    num_threads = 1L,
    tolerance = 1e-6,
    update_epsilon = 1e-10,
    renormalize_each_iteration = FALSE,
    lee_style_epsilon = FALSE,
    random_init = c("legacy", "rnmf"),
    seed = 1L
  ) {
  init_mode <- match.arg(init_mode)
  random_init <- match.arg(random_init)
  if (!is.matrix(x)) {
    stop("x must be a numeric matrix")
  }
  if (!is.null(column_weights) && length(column_weights) != ncol(x)) {
    stop("column_weights must match ncol(x)")
  }
  if (!is.null(row_groups) && length(row_groups) != nrow(x)) {
    stop("row_groups must match nrow(x)")
  }

  out <- .Call(
    "_cellAdmixCore_celladmix_weighted_nmf_matrix",
    x,
    if (is.null(column_weights)) NULL else as.numeric(column_weights),
    as.integer(rank),
    as.integer(max_iterations),
    as.character(init_mode),
    if (is.null(row_groups)) NULL else as.integer(row_groups),
    as.integer(n_runs),
    as.integer(num_threads),
    as.numeric(tolerance),
    as.numeric(update_epsilon),
    as.logical(renormalize_each_iteration),
    as.logical(lee_style_epsilon),
    as.character(random_init),
    as.integer(seed),
    PACKAGE = "cellAdmixCore"
  )

  factor_names <- paste0("factor_", seq_len(ncol(out$w)))
  dimnames(out$w) <- list(rownames(x), factor_names)
  dimnames(out$h) <- list(factor_names, colnames(x))
  out
}

.celladmix_weighted_nmf_sparse_matrix <- function(
    x,
    column_weights = NULL,
    rank = 2L,
    max_iterations = 200L,
    init_mode = c("random", "cluster"),
    row_groups = NULL,
    n_runs = 1L,
    num_threads = 1L,
    tolerance = 1e-6,
    update_epsilon = 1e-10,
    renormalize_each_iteration = FALSE,
    lee_style_epsilon = FALSE,
    random_init = c("legacy", "rnmf"),
    seed = 1L
  ) {
  init_mode <- match.arg(init_mode)
  random_init <- match.arg(random_init)
  out <- .Call(
    "_cellAdmixCore_celladmix_weighted_nmf_sparse_matrix",
    x,
    if (is.null(column_weights)) NULL else as.numeric(column_weights),
    as.integer(rank),
    as.integer(max_iterations),
    as.character(init_mode),
    if (is.null(row_groups)) NULL else as.integer(row_groups),
    as.integer(n_runs),
    as.integer(num_threads),
    as.numeric(tolerance),
    as.numeric(update_epsilon),
    as.logical(renormalize_each_iteration),
    as.logical(lee_style_epsilon),
    as.character(random_init),
    as.integer(seed),
    PACKAGE = "cellAdmixCore"
  )

  factor_names <- paste0("factor_", seq_len(ncol(out$w)))
  dimnames(out$w) <- list(rownames(x), factor_names)
  dimnames(out$h) <- list(factor_names, colnames(x))
  out
}

.celladmix_sparse_nmf_matrix <- function(
    x,
    column_weights = NULL,
    rank = 2L,
    max_iterations = 200L,
    init_mode = c("random", "cluster"),
    row_groups = NULL,
    n_runs = 1L,
    num_threads = 1L,
    tolerance = 1e-6,
    update_epsilon = 1e-10,
    random_init = c("legacy", "rnmf"),
    loss_mode = c("euclidean", "kl"),
    h_l1_penalty = 0,
    h_diversity_penalty = 0,
    seed = 1L
  ) {
  init_mode <- match.arg(init_mode)
  random_init <- match.arg(random_init)
  loss_mode <- match.arg(loss_mode)
  out <- .Call(
    "_cellAdmixCore_celladmix_sparse_nmf_matrix",
    x,
    if (is.null(column_weights)) NULL else as.numeric(column_weights),
    as.integer(rank),
    as.integer(max_iterations),
    as.character(init_mode),
    if (is.null(row_groups)) NULL else as.integer(row_groups),
    as.integer(n_runs),
    as.integer(num_threads),
    as.numeric(tolerance),
    as.numeric(update_epsilon),
    as.character(random_init),
    as.character(loss_mode),
    as.numeric(h_l1_penalty),
    as.numeric(h_diversity_penalty),
    as.integer(seed),
    PACKAGE = "cellAdmixCore"
  )

  factor_names <- paste0("factor_", seq_len(ncol(out$w)))
  dimnames(out$w) <- list(rownames(x), factor_names)
  dimnames(out$h) <- list(factor_names, colnames(x))
  out
}

.celladmix_build_reference_ncv_matrix <- function(
    transcripts,
    k,
    query_ids = NULL,
    row_sample_n = NA_integer_,
    seed = 1L
  ) {
  .Call(
    "_cellAdmixCore_celladmix_build_reference_ncv_matrix",
    transcripts,
    as.integer(k),
    if (is.null(query_ids)) NULL else as.character(query_ids),
    as.integer(row_sample_n),
    as.integer(seed),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_build_pipeline_ncv_matrix <- function(transcripts, k, query_ids = NULL) {
  .Call(
    "_cellAdmixCore_celladmix_build_pipeline_ncv_matrix",
    transcripts,
    as.integer(k),
    if (is.null(query_ids)) NULL else as.character(query_ids),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_load_run_training_transcript_ids <- function(path) {
  .Call(
    "_cellAdmixCore_celladmix_load_run_training_transcript_ids",
    as.character(path),
    PACKAGE = "cellAdmixCore"
  )
}

.celladmix_export_training_ncv_matrix <- function(path) {
  .Call(
    "_cellAdmixCore_celladmix_export_training_ncv_matrix",
    as.character(path),
    PACKAGE = "cellAdmixCore"
  )
}
