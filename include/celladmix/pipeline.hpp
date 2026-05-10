// End-to-end molecule-level factorization, projection, and smoothing pipeline.

#pragma once

#include <string>
#include <vector>

#include "celladmix/graph.hpp"
#include "celladmix/ncv.hpp"
#include "celladmix/nmf_kl.hpp"
#include "celladmix/types.hpp"

namespace celladmix {

struct BasicPipelineOptions {
  int ncv_k = 20;
  int rank = 2;
  int graph_k = 10;
  double same_label_ratio = 5.0;
  int nmf_iterations = 150;
  int nmf_n_runs = 1;
  std::string nmf_init = "auto";
  std::string nmf_variant = "kl";
  std::string molecule_scoring = "gene_loadings";
  int nmf_train_max_rows = -1;
  int nmf_min_molecules = 10;
  int num_threads = 1;
  bool return_ncv = true;
  unsigned int seed = 1;
  std::vector<std::string> training_cell_strata;
  // Optional cell types used only to select NMF training rows. Projection and
  // downstream labeling still run over the full analysis scope.
  std::vector<std::string> training_scope_cell_types;
};

// Stage-level timings reported from the core fit pipeline.
struct BasicPipelineTiming {
  double training_query_sampling_sec = 0.0;
  double training_ncv_sec = 0.0;
  double nmf_fit_sec = 0.0;
  double projection_sec = 0.0;
  double label_smoothing_sec = 0.0;
  double total_sec = 0.0;
};

// Outputs of one full pipeline run on an in-memory transcript table.
struct BasicPipelineResult {
  DenseMatrix ncv;
  NcvFeatureTransform nmf_transform;
  SparseNmfResult nmf;
  DenseMatrix factor_scores;
  std::vector<int> labels;
  std::vector<int> training_query_indices;
  BasicPipelineTiming timing;
};

// Run the full molecule-level factorization and labeling pipeline.
BasicPipelineResult run_basic_pipeline(
    const TranscriptTable& table,
    const BasicPipelineOptions& options = {},
    bool verbose = false);

// Build a keep-mask by removing one factor from one target cell type.
std::vector<bool> apply_removal_rule(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    int factor,
    const std::string& target_cell_type);

// Collapse transcript labels to a dense cell-by-gene count matrix.
DenseMatrix build_cell_gene_counts(
    const TranscriptTable& table,
    const std::vector<bool>& keep_mask = {});

// Pick the factor with the largest loading sum across marker genes.
int select_factor_by_marker_sum(
    const DenseMatrix& h,
    const std::vector<std::string>& genes,
    const std::vector<std::string>& marker_genes);

}  // namespace celladmix
