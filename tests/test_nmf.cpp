#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "celladmix/nmf_kl.hpp"
#include "celladmix/nmf_euclidean.hpp"
#include "test_framework.hpp"

using namespace celladmix;

namespace {

double row_correlation(const DenseMatrix& x, int row_a, const DenseMatrix& y, int row_b) {
  double mean_a = 0.0;
  double mean_b = 0.0;
  for (int j = 0; j < x.cols(); ++j) {
    mean_a += x(row_a, j);
    mean_b += y(row_b, j);
  }
  mean_a /= static_cast<double>(x.cols());
  mean_b /= static_cast<double>(y.cols());

  double num = 0.0;
  double den_a = 0.0;
  double den_b = 0.0;
  for (int j = 0; j < x.cols(); ++j) {
    const double da = x(row_a, j) - mean_a;
    const double db = y(row_b, j) - mean_b;
    num += da * db;
    den_a += da * da;
    den_b += db * db;
  }
  return num / std::sqrt(std::max(den_a * den_b, 1e-12));
}

}  // namespace

TEST_CASE("Weighted NMF reconstructs a simple two-factor matrix") {
  DenseMatrix x(6, 4, 0.0);
  const double raw[6][4] = {
      {4.8, 4.5, 0.1, 0.2},
      {4.3, 3.9, 0.2, 0.1},
      {0.2, 0.3, 5.1, 4.7},
      {0.3, 0.2, 4.6, 4.2},
      {2.2, 2.0, 2.0, 1.8},
      {2.6, 2.2, 1.6, 1.4},
  };
  for (int i = 0; i < x.rows(); ++i) {
    for (int j = 0; j < x.cols(); ++j) {
      x(i, j) = raw[i][j];
    }
  }

  WeightedNmfOptions options;
  options.rank = 2;
  options.max_iterations = 300;
  options.seed = 7;

  const auto result = weighted_nmf(x, default_column_weights(x), options);
  const auto reconstruction = multiply(result.w, result.h);

  REQUIRE_LT(weighted_reconstruction_loss(x, reconstruction, default_column_weights(x)), 0.2);
  REQUIRE_GE(result.losses.size(), 2U);
  REQUIRE_LE(result.losses.back(), result.losses.front());

  const double factor0_marker_sum = result.h(0, 0) + result.h(0, 1);
  const double factor1_marker_sum = result.h(1, 0) + result.h(1, 1);
  const double factor0_other_sum = result.h(0, 2) + result.h(0, 3);
  const double factor1_other_sum = result.h(1, 2) + result.h(1, 3);

  const bool factor0_is_left = factor0_marker_sum > factor0_other_sum;
  const bool factor1_is_left = factor1_marker_sum > factor1_other_sum;
  REQUIRE(factor0_is_left != factor1_is_left);
}

TEST_CASE("Sparse weighted NMF reconstructs the same two-factor structure") {
  SparseRowMatrix x(
      6,
      4,
      {0, 2, 4, 6, 8, 12, 16},
      {0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3},
      {4.8, 4.5, 4.3, 3.9, 5.1, 4.7, 4.6, 4.2, 2.2, 2.0, 2.0, 1.8, 2.6, 2.2, 1.6, 1.4});

  WeightedNmfOptions options;
  options.rank = 2;
  options.max_iterations = 300;
  options.seed = 7;

  const auto weights = default_column_weights(x);
  const auto result = weighted_nmf(x, weights, options);
  const auto reconstruction = multiply(result.w, result.h);

  REQUIRE_LT(weighted_reconstruction_loss(x, reconstruction, weights), 0.2);
  REQUIRE_GE(result.losses.size(), 2U);
  REQUIRE_LE(result.losses.back(), result.losses.front());

  const auto projected = project_to_factors(x, result.h);
  REQUIRE_EQ(projected.rows(), 6);
  REQUIRE_EQ(projected.cols(), 2);
}

TEST_CASE("Sparse weighted NMF supports multirun selection") {
  SparseRowMatrix x(
      6,
      4,
      {0, 2, 4, 6, 8, 12, 16},
      {0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3},
      {4.8, 4.5, 4.3, 3.9, 5.1, 4.7, 4.6, 4.2, 2.2, 2.0, 2.0, 1.8, 2.6, 2.2, 1.6, 1.4});

  WeightedNmfOptions options;
  options.rank = 2;
  options.max_iterations = 300;
  options.seed = 7;
  options.num_threads = 2;
  options.n_runs = 4;

  const auto weights = default_column_weights(x);
  const auto result = weighted_nmf(x, weights, options);
  const auto reconstruction = multiply(result.w, result.h);

  REQUIRE_LT(weighted_reconstruction_loss(x, reconstruction, weights), 0.2);
  REQUIRE_GE(result.losses.size(), 2U);
  REQUIRE_LE(result.losses.back(), result.losses.front());
  REQUIRE_EQ(result.candidate_final_losses.size(), 4U);
  REQUIRE_GE(result.selected_run, 0);
  REQUIRE_LT(result.selected_run, 4);
  REQUIRE_EQ(
      result.losses.back(),
      *std::min_element(result.candidate_final_losses.begin(), result.candidate_final_losses.end()));
}

TEST_CASE("Group-guided initialization stabilizes factor loadings across seeds") {
  DenseMatrix x(6, 4, 0.0);
  const double raw[6][4] = {
      {6.0, 5.5, 0.2, 0.1},
      {5.8, 5.1, 0.3, 0.2},
      {5.6, 5.3, 0.2, 0.3},
      {0.2, 0.3, 6.1, 5.7},
      {0.1, 0.2, 5.8, 5.3},
      {0.3, 0.1, 5.6, 5.5},
  };
  for (int i = 0; i < x.rows(); ++i) {
    for (int j = 0; j < x.cols(); ++j) {
      x(i, j) = raw[i][j];
    }
  }

  WeightedNmfOptions options_a;
  options_a.rank = 2;
  options_a.max_iterations = 200;
  options_a.seed = 3;
  options_a.init_groups = {0, 0, 0, 1, 1, 1};

  auto options_b = options_a;
  options_b.seed = 17;

  const auto weights = default_column_weights(x);
  const auto fit_a = weighted_nmf(x, weights, options_a);
  const auto fit_b = weighted_nmf(x, weights, options_b);

  const double corr_00 = std::abs(row_correlation(fit_a.h, 0, fit_b.h, 0));
  const double corr_01 = std::abs(row_correlation(fit_a.h, 0, fit_b.h, 1));
  const double corr_10 = std::abs(row_correlation(fit_a.h, 1, fit_b.h, 0));
  const double corr_11 = std::abs(row_correlation(fit_a.h, 1, fit_b.h, 1));
  const double best_alignment = std::max(corr_00 + corr_11, corr_01 + corr_10);

  REQUIRE_GT(best_alignment, 1.9);
}

TEST_CASE("Sparse NMF supports generalized KL loss") {
  SparseRowMatrix x(
      6,
      4,
      {0, 2, 4, 6, 8, 12, 16},
      {0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3},
      {4.8, 4.5, 4.3, 3.9, 5.1, 4.7, 4.6, 4.2, 2.2, 2.0, 2.0, 1.8, 2.6, 2.2, 1.6, 1.4});

  SparseNmfOptions options;
  options.rank = 2;
  options.max_iterations = 150;
  options.loss_mode = "kl";
  options.seed = 11;

  const auto result = sparse_nmf(x, {}, options);
  REQUIRE_EQ(result.w.rows(), 6);
  REQUIRE_EQ(result.w.cols(), 2);
  REQUIRE_EQ(result.h.rows(), 2);
  REQUIRE_EQ(result.h.cols(), 4);
  REQUIRE_GE(result.losses.size(), 2U);
  REQUIRE_LE(result.losses.back(), result.losses.front());
  REQUIRE(std::isfinite(result.final_objective));
}

TEST_CASE("NMF factor ordering follows selected training W mass") {
  SparseNmfResult result;
  result.w = DenseMatrix(3, 3, 0.0);
  result.h = DenseMatrix(3, 2, 0.0);
  result.selected_factor_stability = {0.1, 0.2, 0.3};

  result.w(0, 0) = 1.0;
  result.w(0, 1) = 2.0;
  result.w(1, 1) = 3.0;
  result.w(2, 2) = 3.0;

  result.h(0, 0) = 10.0;
  result.h(0, 1) = 11.0;
  result.h(1, 0) = 20.0;
  result.h(1, 1) = 21.0;
  result.h(2, 0) = 30.0;
  result.h(2, 1) = 31.0;

  order_nmf_factors_by_training_importance(result);

  const auto masses = result.w.col_sums();
  REQUIRE_EQ(masses[0], 5.0);
  REQUIRE_EQ(masses[1], 3.0);
  REQUIRE_EQ(masses[2], 1.0);
  REQUIRE_EQ(result.h(0, 0), 20.0);
  REQUIRE_EQ(result.h(1, 0), 30.0);
  REQUIRE_EQ(result.h(2, 0), 10.0);
  REQUIRE_EQ(result.selected_factor_stability[0], 0.2);
  REQUIRE_EQ(result.selected_factor_stability[1], 0.3);
  REQUIRE_EQ(result.selected_factor_stability[2], 0.1);
}

TEST_CASE("Sparse NMF supports H regularization") {
  SparseRowMatrix x(
      6,
      4,
      {0, 2, 4, 6, 8, 12, 16},
      {0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3},
      {4.8, 4.5, 4.3, 3.9, 5.1, 4.7, 4.6, 4.2, 2.2, 2.0, 2.0, 1.8, 2.6, 2.2, 1.6, 1.4});

  SparseNmfOptions base_options;
  base_options.rank = 2;
  base_options.max_iterations = 150;
  base_options.seed = 11;
  const auto weights = default_column_weights(x);
  const auto baseline = sparse_nmf(x, weights, base_options);

  auto penalized_options = base_options;
  penalized_options.h_l1_penalty = 0.05;
  penalized_options.h_diversity_penalty = 0.05;
  const auto penalized = sparse_nmf(x, weights, penalized_options);

  const double h_sum_baseline = std::accumulate(
      baseline.h.data().begin(),
      baseline.h.data().end(),
      0.0);
  const double h_sum_penalized = std::accumulate(
      penalized.h.data().begin(),
      penalized.h.data().end(),
      0.0);
  REQUIRE_LT(h_sum_penalized, h_sum_baseline);
  REQUIRE(std::isfinite(penalized.final_objective));
}
