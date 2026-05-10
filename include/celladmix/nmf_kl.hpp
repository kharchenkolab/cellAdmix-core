// Sparse NMF implementation supporting KL and Euclidean objectives, with the
// KL path used by the production pipeline.

#pragma once

#include <string>
#include <vector>

#include "celladmix/matrix.hpp"

namespace celladmix {

struct SparseNmfOptions {
  int rank = 2;
  int max_iterations = 200;
  double tolerance = 1e-6;
  double update_epsilon = 1e-10;
  unsigned int seed = 1;
  int num_threads = 1;
  int n_runs = 1;
  std::string init_mode = "random";
  std::string random_init = "legacy";
  std::string loss_mode = "euclidean";
  double h_l1_penalty = 0.0;
  double h_diversity_penalty = 0.0;
  std::vector<int> init_groups;
};

// Result of a sparse NMF solve, including multirun selection metadata.
struct SparseNmfResult {
  DenseMatrix w;
  DenseMatrix h;
  std::vector<double> losses;
  double final_objective = 0.0;
  unsigned int selected_seed = 1;
  int selected_run = 0;
  std::vector<double> candidate_final_objectives;
  std::vector<double> candidate_best_match_correlations;
  std::vector<double> selected_factor_stability;
  double candidate_final_objective_mean = 0.0;
  double candidate_final_objective_sd = 0.0;
  double candidate_best_match_correlation_mean = 0.0;
};

// Fit sparse NMF under the requested objective and initialization strategy.
SparseNmfResult sparse_nmf(
    const SparseRowMatrix& x,
    const std::vector<double>& column_weights,
    const SparseNmfOptions& options = {});

// Canonicalize factor numbering by descending selected-run training W mass.
void order_nmf_factors_by_training_importance(SparseNmfResult& result);

}  // namespace celladmix
