#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "celladmix/nmf_kl.hpp"
#include "celladmix/nmf_stability.hpp"
#include "celladmix/pipeline.hpp"
#include "celladmix/pipeline_common.hpp"
#include "celladmix/simulation.hpp"
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

TEST_CASE("Factor matching recovers a permutation and pads rectangular inputs") {
  const DenseMatrix reference = from_rows({
      {5.0, 4.0, 0.1, 0.0, 0.0, 0.2},
      {0.1, 0.0, 6.0, 5.0, 0.0, 0.1},
      {0.0, 0.2, 0.0, 0.1, 7.0, 6.0}});
  const DenseMatrix permuted = from_rows({
      {0.0, 0.3, 0.1, 0.0, 6.5, 6.5},
      {5.5, 3.5, 0.0, 0.1, 0.0, 0.2},
      {0.2, 0.0, 5.5, 5.5, 0.1, 0.0}});
  const auto match = match_factor_profiles(permuted, reference);
  REQUIRE_EQ(match.match.size(), 3U);
  REQUIRE_EQ(match.match[0], 2);
  REQUIRE_EQ(match.match[1], 0);
  REQUIRE_EQ(match.match[2], 1);
  for (const double value : match.correlation) {
    REQUIRE(value > 0.9);
  }

  // Fewer reference factors than profiles: one profile stays unmatched.
  const DenseMatrix two = from_rows({
      {5.0, 4.0, 0.1, 0.0, 0.0, 0.2},
      {0.1, 0.0, 6.0, 5.0, 0.0, 0.1}});
  const auto partial = match_factor_profiles(permuted, two);
  REQUIRE_EQ(partial.match[0], -1);
  REQUIRE_EQ(partial.correlation[0], 0.0);
  REQUIRE_EQ(partial.match[1], 0);
  REQUIRE_EQ(partial.match[2], 1);
}

TEST_CASE("Profile alignment places genes by name, drops unknown ones, and zero-fills") {
  const std::vector<double> profiles = {1.0, 2.0, 3.0,   // factor 0 over a, b, c
                                        4.0, 5.0, 6.0};  // factor 1
  int matched = 0;
  const auto aligned = align_profiles_to_genes(
      profiles, 2, {"a", "b", "c"}, {"c", "x", "a"}, &matched);
  REQUIRE_EQ(matched, 2);
  REQUIRE_EQ(aligned.size(), 6U);
  REQUIRE_EQ(aligned[0], 3.0);
  REQUIRE_EQ(aligned[1], 0.0);
  REQUIRE_EQ(aligned[2], 1.0);
  REQUIRE_EQ(aligned[3], 6.0);
  REQUIRE_EQ(aligned[4], 0.0);
  REQUIRE_EQ(aligned[5], 4.0);
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

TEST_CASE("Automatic ncv_k grows with panel size and respects cell caps") {
  // Anchor: a ~400-gene panel keeps the historical default.
  REQUIRE_EQ(resolve_auto_ncv_k(400, 1000.0), 20);
  // Small panels never drop below the default.
  REQUIRE_EQ(resolve_auto_ncv_k(100, 1000.0), 20);
  // sqrt growth: 960 genes -> ~31, 5000 genes -> ~71.
  REQUIRE_EQ(resolve_auto_ncv_k(960, 1000.0), 31);
  REQUIRE_EQ(resolve_auto_ncv_k(5000, 1000.0), 71);
  // Median cell size caps the growth at half a median cell.
  REQUIRE_EQ(resolve_auto_ncv_k(5000, 100.0), 50);
  // Tiny cells fall back to the historical default rather than below it.
  REQUIRE_EQ(resolve_auto_ncv_k(5000, 30.0), 20);
  REQUIRE_EQ(resolve_auto_ncv_k(400, 30.0), 20);
}

TEST_CASE("Pipeline resolves the automatic ncv_k sentinel into options") {
  SimulationParams params;
  params.transcripts_per_cell = 60;
  params.admixture_per_target_cell = 4;
  params.seed = 5U;
  const auto sim = simulate_nsclc_admixture(params);

  BasicPipelineOptions options;
  options.ncv_k = 0;
  options.rank = 2;
  options.graph_k = 4;
  options.nmf_iterations = 30;
  options.nmf_train_max_rows = 40;
  options.return_ncv = false;
  options.seed = 5U;
  const auto fit = run_basic_pipeline(sim.transcripts, options);
  // The simulation panel is tiny, so auto resolution lands on the floor of 20
  // and the resolved value is recorded on the result.
  REQUIRE_EQ(fit.resolved_options.ncv_k, 20);

  options.ncv_k = 6;
  const auto explicit_fit = run_basic_pipeline(sim.transcripts, options);
  REQUIRE_EQ(explicit_fit.resolved_options.ncv_k, 6);
}

TEST_CASE("Profile clustering groups recurring programs across factorizations") {
  // Three factorizations of the same three programs in different orders, with
  // small perturbations, plus one program that appears in a single
  // factorization only (gene 6 belongs to it alone).
  const DenseMatrix pooled = from_rows({
      {5.0, 4.0, 0.1, 0.0, 0.0, 0.2, 0.0},
      {0.1, 0.0, 6.0, 5.0, 0.0, 0.1, 0.0},
      {0.0, 0.2, 0.0, 0.1, 7.0, 6.0, 0.0},
      {0.0, 0.3, 0.1, 0.0, 6.5, 6.5, 0.0},
      {5.5, 3.5, 0.0, 0.1, 0.0, 0.2, 0.0},
      {0.2, 0.0, 5.5, 5.5, 0.1, 0.0, 0.0},
      {0.0, 0.1, 6.2, 4.8, 0.2, 0.0, 0.0},
      {0.1, 0.0, 0.0, 0.2, 7.2, 5.8, 0.0},
      {4.8, 4.2, 0.2, 0.0, 0.1, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0}});
  const std::vector<int> group = {0, 0, 0, 1, 1, 1, 2, 2, 2, 2};
  const auto result = cluster_factor_profiles(pooled, group, 0.3);
  REQUIRE_EQ(result.n_clusters, 4);
  REQUIRE_EQ(result.cluster.size(), 10U);
  // Program 1 in rows 0, 4, 8; program 2 in rows 1, 5, 6; program 3 in rows 2, 3, 7.
  REQUIRE_EQ(result.cluster[0], result.cluster[4]);
  REQUIRE_EQ(result.cluster[0], result.cluster[8]);
  REQUIRE_EQ(result.cluster[1], result.cluster[5]);
  REQUIRE_EQ(result.cluster[1], result.cluster[6]);
  REQUIRE_EQ(result.cluster[2], result.cluster[3]);
  REQUIRE_EQ(result.cluster[2], result.cluster[7]);
  REQUIRE(result.cluster[0] != result.cluster[1]);
  REQUIRE(result.cluster[0] != result.cluster[2]);
  REQUIRE(result.cluster[1] != result.cluster[2]);
  // Clusters are numbered by first appearance.
  REQUIRE_EQ(result.cluster[0], 0);
  REQUIRE_EQ(result.cluster[1], 1);
  REQUIRE_EQ(result.cluster[2], 2);
  REQUIRE_EQ(result.cluster[9], 3);
  for (int i = 0; i < 9; ++i) {
    REQUIRE(result.within_correlation[static_cast<std::size_t>(i)] > 0.8);
  }
  REQUIRE(std::isnan(result.within_correlation[9]));

  // Ownership is normalized within each factorization, so rescaling one
  // factorization leaves the clustering unchanged.
  DenseMatrix scaled = pooled;
  for (int j = 0; j < scaled.cols(); ++j) {
    for (int i = 3; i < 6; ++i) {
      scaled(i, j) *= 10.0;
    }
  }
  const auto rescaled = cluster_factor_profiles(scaled, group, 0.3);
  REQUIRE(rescaled.cluster == result.cluster);

  // A permissive threshold merges everything; a strict one merges nothing.
  REQUIRE_EQ(cluster_factor_profiles(pooled, group, -1.5).n_clusters, 1);
  REQUIRE_EQ(cluster_factor_profiles(pooled, group, 1.5).n_clusters, 10);
}
