// Local molecular coherence scoring for factor-assigned transcripts.

#pragma once

#include <functional>
#include <string>
#include <vector>

#include "celladmix/membrane.hpp"
#include "celladmix/types.hpp"

namespace celladmix {

// Controls within-cell molecular coherence scoring. The score asks whether
// factor-assigned molecules form locally coherent same-factor patches, with
// optional membrane-stain barriers reducing support across image boundaries.
struct CoherenceTestOptions {
  int k_neighbors = 10;
  int min_factor_molecules = 5;
  int min_cells = 5;
  double max_neighbor_distance = -1.0;
  double distance_sigma = -1.0;
  double lambda_coherence = 1.0;
  double beta_margin = 1.0;
  double score_threshold = 0.0;
  double source_pseudocount = 0.5;
  double min_source_log_enrichment = 0.0;
  bool include_self_source = false;
  bool compute_null = true;
  std::string null_method = "label_permutation";
  int null_iterations = 20;
  bool null_exclude_factor = true;
  bool null_match_nucleus = true;
  bool null_match_density = true;
  int null_nucleus_distance_bins = 3;
  int null_density_bins = 3;
  unsigned int seed = 1;
  bool use_membrane_barrier = false;
  bool normalize_membrane = true;
  double membrane_low_quantile = 0.05;
  double membrane_high_quantile = 0.995;
  double membrane_alpha = 3.0;
  int line_samples = 12;
  double patch_edge_weight_min = 0.1;
  int num_threads = 1;
  std::function<void(const std::string&, double)> progress;
};

// Per-target-cell score for one source type and factor.
struct CoherenceCellScore {
  int target_cell = -1;
  int target_type = -1;
  int source_type = -1;
  int factor = -1;
  int factor_count = 0;
  int active_count = 0;
  double mean_raw_score = 0.0;
  double mean_coherence_score = 0.0;
  double mean_score = 0.0;
  double mean_null_score = 0.0;
  double score_delta = 0.0;
  double q75_score = 0.0;
  double active_fraction = 0.0;
  double mean_edge_weight = 0.0;
  int largest_patch_count = 0;
  double largest_patch_fraction = 0.0;
  double patch_score = 0.0;
  double mean_null_patch_score = 0.0;
  double patch_score_delta = 0.0;
  double source_log_enrichment = 0.0;
  double source_probability = 0.0;
  bool used_in_summary = false;
};

// Statistical summary for one target/source/factor group.
struct CoherenceSummary {
  int target_type = -1;
  int source_type = -1;
  int factor = -1;
  int n_cells = 0;
  int n_molecules = 0;
  int n_active_molecules = 0;
  double active_fraction = 0.0;
  double mean_score = 0.0;
  double mean_null_score = 0.0;
  double mean_delta_score = 0.0;
  double q75_score = 0.0;
  double mean_patch_score = 0.0;
  double mean_null_patch_score = 0.0;
  double mean_delta_patch_score = 0.0;
  double mean_largest_patch_fraction = 0.0;
  double patch_p_value = 1.0;
  double patch_neg_log10_p = 0.0;
  double p_value = 1.0;
  double neg_log10_p = 0.0;
  double source_log_enrichment = 0.0;
  double source_probability = 0.0;
};

struct CoherenceTestResult {
  std::vector<std::string> cell_types;
  std::vector<double> source_log_enrichment;
  std::vector<double> source_probability;
  bool membrane_normalized = false;
  double membrane_scale_low = 0.0;
  double membrane_scale_high = 1.0;
  std::string null_method = "none";
  bool null_match_nucleus_overlap_used = false;
  bool null_match_nucleus_distance_used = false;
  bool null_match_density_used = false;
  int null_nucleus_distance_bins = 0;
  int null_density_bins = 0;
  int n_factors = 0;
  std::vector<CoherenceCellScore> cell_scores;
  std::vector<CoherenceSummary> summaries;
};

// Run local coherence scoring over fitted transcript factor labels.
CoherenceTestResult run_coherence_test(
    const TranscriptTable& table,
    const CellTable& cells,
    const std::vector<int>& labels,
    const std::vector<double>& factor_margin,
    const Image2D* membrane_image,
    const CoherenceTestOptions& options = {});

}  // namespace celladmix
