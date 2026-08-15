#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "celladmix/nmf_kl.hpp"
#include "celladmix/nmf_stability.hpp"
#include "test_framework.hpp"

using namespace celladmix;

namespace {

DenseMatrix from_rows(const std::vector<std::vector<double>>& rows) {
  DenseMatrix out(static_cast<int>(rows.size()), static_cast<int>(rows.front().size()), 0.0);
  for (int i = 0; i < out.rows(); ++i) {
    for (int j = 0; j < out.cols(); ++j) {
      out(i, j) = rows[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    }
  }
  return out;
}

double brute_force_best_assignment(const DenseMatrix& similarity) {
  const int n = similarity.rows();
  std::vector<int> perm(static_cast<std::size_t>(n));
  std::iota(perm.begin(), perm.end(), 0);
  double best = -1e300;
  do {
    double total = 0.0;
    for (int i = 0; i < n; ++i) {
      total += similarity(i, perm[static_cast<std::size_t>(i)]);
    }
    best = std::max(best, total);
  } while (std::next_permutation(perm.begin(), perm.end()));
  return best;
}

}  // namespace

TEST_CASE("Hungarian match maximizes total similarity") {
  const DenseMatrix a = from_rows({
      {0.9, 0.1, 0.4},
      {0.8, 0.7, 0.2},
      {0.5, 0.6, 0.3},
  });
  const DenseMatrix b = from_rows({
      {-0.2, 0.9, 0.3, 0.1},
      {0.7, -0.5, 0.6, 0.2},
      {0.1, 0.4, -0.1, 0.8},
      {0.9, 0.2, 0.5, -0.3},
  });
  for (const DenseMatrix* m : {&a, &b}) {
    const auto match = hungarian_match(*m);
    std::vector<char> used(static_cast<std::size_t>(m->rows()), 0);
    double total = 0.0;
    for (int i = 0; i < m->rows(); ++i) {
      const int j = match[static_cast<std::size_t>(i)];
      REQUIRE_GE(j, 0);
      REQUIRE_LT(j, m->cols());
      REQUIRE(!used[static_cast<std::size_t>(j)]);
      used[static_cast<std::size_t>(j)] = 1;
      total += (*m)(i, j);
    }
    REQUIRE_NEAR(total, brute_force_best_assignment(*m), 1e-9);
  }
}

TEST_CASE("Ownership profiles normalize genes across factors") {
  const DenseMatrix h = from_rows({
      {3.0, 0.0, 1.0, 5.0},
      {1.0, 0.0, 3.0, 5.0},
  });
  const auto profiles = ownership_profiles(h);
  REQUIRE_NEAR(profiles(0, 0), 0.75, 1e-12);
  REQUIRE_NEAR(profiles(1, 0), 0.25, 1e-12);
  // Zero-loading gene stays zero rather than becoming NaN.
  REQUIRE_EQ(profiles(0, 1), 0.0);
  REQUIRE_EQ(profiles(1, 1), 0.0);
  for (int j = 0; j < h.cols(); ++j) {
    if (j == 1) {
      continue;
    }
    REQUIRE_NEAR(profiles(0, j) + profiles(1, j), 1.0, 1e-12);
  }
}

TEST_CASE("Matched ownership stability is invariant to factor permutation") {
  const DenseMatrix h1 = from_rows({
      {5.0, 4.0, 0.1, 0.2, 0.1, 0.3},
      {0.2, 0.1, 6.0, 5.0, 0.2, 0.1},
      {0.1, 0.3, 0.2, 0.1, 4.0, 5.0},
  });
  // Same factors in a different order, lightly perturbed.
  const DenseMatrix h2 = from_rows({
      {0.1, 0.3, 0.2, 0.1, 4.2, 4.8},
      {5.2, 3.9, 0.1, 0.2, 0.1, 0.3},
      {0.2, 0.1, 5.9, 5.1, 0.2, 0.1},
  });
  const std::vector<const DenseMatrix*> runs = {&h1, &h2};
  const auto result = matched_ownership_stability(runs, 0, {1});
  REQUIRE_EQ(result.comparison_runs, 1);
  REQUIRE_EQ(static_cast<int>(result.factor_stability.size()), 3);
  for (const double value : result.factor_stability) {
    REQUIRE_GT(value, 0.95);
  }
  REQUIRE_EQ(result.stable_factor_count, 3);
  REQUIRE_NEAR(result.run_matched_means[0], 1.0, 1e-12);
  REQUIRE_GT(result.run_matched_means[1], 0.95);
}

TEST_CASE("Unrelated factors sharing an abundance backbone score near zero") {
  // Every factor carries the same strong backbone on genes 0-3; the specific
  // genes differ between the two runs, so nothing is genuinely re-found.
  const double b = 50.0;
  const DenseMatrix h1 = from_rows({
      {b, b, b, b, 9.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
      {b, b, b, b, 0.1, 9.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
      {b, b, b, b, 0.1, 0.1, 9.0, 0.1, 0.1, 0.1, 0.1, 0.1},
      {b, b, b, b, 0.1, 0.1, 0.1, 9.0, 0.1, 0.1, 0.1, 0.1},
  });
  const DenseMatrix h2 = from_rows({
      {b, b, b, b, 0.1, 0.1, 0.1, 0.1, 9.0, 0.1, 0.1, 0.1},
      {b, b, b, b, 0.1, 0.1, 0.1, 0.1, 0.1, 9.0, 0.1, 0.1},
      {b, b, b, b, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 9.0, 0.1},
      {b, b, b, b, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 9.0},
  });
  const std::vector<const DenseMatrix*> runs = {&h1, &h2};
  const auto result = matched_ownership_stability(runs, 0, {1});
  for (const double value : result.factor_stability) {
    REQUIRE_LT(std::abs(value), 0.3);
  }
  REQUIRE_EQ(result.stable_factor_count, 0);
}

TEST_CASE("Stability averages only the requested comparison runs") {
  const DenseMatrix h1 = from_rows({
      {9.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
      {0.1, 9.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
      {0.1, 0.1, 9.0, 0.1, 0.1, 0.1, 0.1, 0.1},
  });
  const DenseMatrix same = from_rows({
      {8.8, 0.1, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1},
      {0.2, 9.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1},
      {0.1, 0.2, 8.9, 0.1, 0.1, 0.1, 0.1, 0.1},
  });
  const DenseMatrix different = from_rows({
      {0.1, 0.1, 0.1, 9.0, 0.1, 0.1, 0.1, 0.1},
      {0.1, 0.1, 0.1, 0.1, 9.0, 0.1, 0.1, 0.1},
      {0.1, 0.1, 0.1, 0.1, 0.1, 9.0, 0.1, 0.1},
  });
  const std::vector<const DenseMatrix*> runs = {&h1, &different, &same};
  // Only run 2 is in the comparison set; run 1 still gets a matched mean.
  const auto result = matched_ownership_stability(runs, 0, {2});
  REQUIRE_EQ(result.comparison_runs, 1);
  for (const double value : result.factor_stability) {
    REQUIRE_GT(value, 0.9);
  }
  REQUIRE_LT(result.run_matched_means[1], 0.3);
  REQUIRE_GT(result.run_matched_means[2], 0.9);
}

TEST_CASE("Sparse NMF multirun reports ownership stability diagnostics") {
  // Two clean blocks so restarts agree on the factors.
  std::vector<int> indptr = {0};
  std::vector<int> indices;
  std::vector<double> values;
  for (int i = 0; i < 12; ++i) {
    const bool left = i < 6;
    indices.push_back(left ? 0 : 2);
    values.push_back(5.0);
    indices.push_back(left ? 1 : 3);
    values.push_back(4.0);
    indptr.push_back(static_cast<int>(indices.size()));
  }
  const SparseRowMatrix x(12, 4, std::move(indptr), std::move(indices), std::move(values));

  SparseNmfOptions options;
  options.rank = 2;
  options.n_runs = 4;
  options.loss_mode = "kl";
  options.max_iterations = 200;
  options.seed = 3;

  const auto result = sparse_nmf(x, {}, options);
  REQUIRE_EQ(static_cast<int>(result.candidate_best_match_correlations.size()), 4);
  REQUIRE_EQ(static_cast<int>(result.selected_factor_stability.size()), 2);
  REQUIRE_EQ(result.stability_comparison_runs, 3);
  REQUIRE_NEAR(result.stability_threshold, kStableFactorThreshold, 1e-12);
  for (const double value : result.selected_factor_stability) {
    REQUIRE_GT(value, 0.8);
  }
  REQUIRE_EQ(result.stable_factor_count, 2);
}

TEST_CASE("Cluster-initialized multirun keeps run zero out of the comparison set") {
  std::vector<int> indptr = {0};
  std::vector<int> indices;
  std::vector<double> values;
  std::vector<int> groups;
  for (int i = 0; i < 12; ++i) {
    const bool left = i < 6;
    groups.push_back(left ? 0 : 1);
    indices.push_back(left ? 0 : 2);
    values.push_back(5.0);
    indices.push_back(left ? 1 : 3);
    values.push_back(4.0);
    indptr.push_back(static_cast<int>(indices.size()));
  }
  const SparseRowMatrix x(12, 4, std::move(indptr), std::move(indices), std::move(values));

  SparseNmfOptions options;
  options.rank = 2;
  options.n_runs = 4;
  options.loss_mode = "kl";
  options.init_mode = "cluster";
  options.init_groups = groups;
  options.max_iterations = 200;
  options.seed = 3;

  const auto result = sparse_nmf(x, {}, options);
  // Run 0 (cluster init) is excluded from the comparison set; the selected
  // run is excluded as well, so at most n_runs - 1 comparisons remain.
  REQUIRE_GE(result.stability_comparison_runs, 2);
  REQUIRE_LE(result.stability_comparison_runs, 3);
  REQUIRE_EQ(static_cast<int>(result.candidate_best_match_correlations.size()), 4);
  for (const double value : result.selected_factor_stability) {
    REQUIRE_GT(value, 0.8);
  }
}
