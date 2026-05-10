// Weighted Euclidean/least-squares NMF used by the ls_nmf production option
// and by comparison utilities.

#pragma once

#include <string>
#include <vector>

#include "celladmix/matrix.hpp"

namespace celladmix {

struct WeightedNmfOptions {
  int rank = 2;
  int max_iterations = 200;
  double tolerance = 1e-6;
  double update_epsilon = 1e-10;
  unsigned int seed = 1;
  int num_threads = 1;
  int n_runs = 1;
  bool renormalize_each_iteration = false;
  bool lee_style_epsilon = false;
  std::string random_init = "legacy";
  std::vector<int> init_groups;
};

// Result of one weighted-Euclidean NMF solve, possibly selected from multirun.
struct WeightedNmfResult {
  DenseMatrix w;
  DenseMatrix h;
  std::vector<double> losses;
  unsigned int selected_seed = 1;
  int selected_run = 0;
  std::vector<double> candidate_final_losses;
  std::vector<double> selected_factor_stability;
};

// Compute default inverse column-frequency weights for a dense matrix.
std::vector<double> default_column_weights(const DenseMatrix& x);
// Compute default inverse column-frequency weights for a sparse matrix.
std::vector<double> default_column_weights(const SparseRowMatrix& x);

// Evaluate weighted Euclidean reconstruction loss on a dense matrix.
double weighted_reconstruction_loss(
    const DenseMatrix& x,
    const DenseMatrix& reconstruction,
    const std::vector<double>& column_weights);

// Evaluate weighted Euclidean reconstruction loss on a sparse matrix.
double weighted_reconstruction_loss(
    const SparseRowMatrix& x,
    const DenseMatrix& reconstruction,
    const std::vector<double>& column_weights);

// Fit weighted Euclidean NMF on a dense input matrix.
WeightedNmfResult weighted_nmf(
    const DenseMatrix& x,
    const std::vector<double>& column_weights,
    const WeightedNmfOptions& options = {});

// Fit weighted Euclidean NMF on a sparse input matrix.
WeightedNmfResult weighted_nmf(
    const SparseRowMatrix& x,
    const std::vector<double>& column_weights,
    const WeightedNmfOptions& options = {});

// Project dense rows into Euclidean factor space with fixed H.
DenseMatrix project_to_factors(const DenseMatrix& x, const DenseMatrix& h);
// Project sparse rows into Euclidean factor space with fixed H.
DenseMatrix project_to_factors(const SparseRowMatrix& x, const DenseMatrix& h);

}  // namespace celladmix
