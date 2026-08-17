#include "celladmix/nmf_euclidean.hpp"

#include "celladmix/nmf_stability.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_map>

#include "subpar/range.hpp"

namespace celladmix {

// Weighted Euclidean/least-squares NMF implementation backing the ls_nmf path.

// Derive inverse column-frequency weights from a dense matrix.
std::vector<double> default_column_weights(const DenseMatrix& x) {
  auto sums = x.col_sums();
  for (double& value : sums) {
    value = 1.0 / std::max(value, 1e-8);
  }
  return sums;
}

// Derive inverse column-frequency weights from a sparse matrix.
std::vector<double> default_column_weights(const SparseRowMatrix& x) {
  auto sums = x.col_sums();
  for (double& value : sums) {
    value = 1.0 / std::max(value, 1e-8);
  }
  return sums;
}

// Evaluate weighted squared reconstruction error for dense inputs.
double weighted_reconstruction_loss(
    const DenseMatrix& x,
    const DenseMatrix& reconstruction,
    const std::vector<double>& column_weights) {
  if (x.rows() != reconstruction.rows() || x.cols() != reconstruction.cols()) {
    throw std::runtime_error("Matrices must match for weighted_reconstruction_loss");
  }
  if (static_cast<int>(column_weights.size()) != x.cols()) {
    throw std::runtime_error("column_weights length must match matrix columns");
  }

  double loss = 0.0;
  for (int i = 0; i < x.rows(); ++i) {
    for (int j = 0; j < x.cols(); ++j) {
      const double diff = x(i, j) - reconstruction(i, j);
      loss += column_weights[static_cast<std::size_t>(j)] * diff * diff;
    }
  }
  return loss;
}

// Evaluate weighted squared reconstruction error for sparse inputs.
double weighted_reconstruction_loss(
    const SparseRowMatrix& x,
    const DenseMatrix& reconstruction,
    const std::vector<double>& column_weights) {
  if (x.rows() != reconstruction.rows() || x.cols() != reconstruction.cols()) {
    throw std::runtime_error("Matrices must match for weighted_reconstruction_loss");
  }
  if (static_cast<int>(column_weights.size()) != x.cols()) {
    throw std::runtime_error("column_weights length must match matrix columns");
  }

  double loss = 0.0;
  for (int i = 0; i < reconstruction.rows(); ++i) {
    for (int j = 0; j < reconstruction.cols(); ++j) {
      const double value = reconstruction(i, j);
      loss += column_weights[static_cast<std::size_t>(j)] * value * value;
    }
  }
  for (int i = 0; i < x.rows(); ++i) {
    for (int p = x.indptr()[static_cast<std::size_t>(i)];
         p < x.indptr()[static_cast<std::size_t>(i + 1)];
         ++p) {
      const int j = x.indices()[static_cast<std::size_t>(p)];
      const double value = x.values()[static_cast<std::size_t>(p)];
      loss += column_weights[static_cast<std::size_t>(j)] *
          (value * value - 2.0 * value * reconstruction(i, j));
    }
  }
  return loss;
}

namespace {

// Estimate a dense matrix scale for legacy initialization.
double dense_mean(const DenseMatrix& x) {
  return std::accumulate(x.data().begin(), x.data().end(), 0.0) /
      static_cast<double>(std::max(1, x.rows() * x.cols()));
}

// Estimate a dense matrix maximum for RNMF-style initialization.
double dense_max(const DenseMatrix& x) {
  return x.data().empty() ? 0.0 : *std::max_element(x.data().begin(), x.data().end());
}

// Estimate a sparse matrix scale for legacy initialization.
double sparse_mean(const SparseRowMatrix& x) {
  return std::accumulate(x.values().begin(), x.values().end(), 0.0) /
      static_cast<double>(std::max(1, x.rows() * x.cols()));
}

// Estimate a sparse matrix maximum for RNMF-style initialization.
double sparse_max(const SparseRowMatrix& x) {
  return x.values().empty() ? 0.0 : *std::max_element(x.values().begin(), x.values().end());
}

// Clamp the requested worker count to the available amount of work.
int effective_workers(int requested, int num_tasks) {
  if (num_tasks <= 0) {
    return 1;
  }
  return std::max(1, subpar::sanitize_num_workers(requested, num_tasks));
}

// Randomly initialize W under the requested initialization style.
DenseMatrix initialize_w(
    int rows,
    int rank,
    double scale_x,
    const WeightedNmfOptions& options,
    std::mt19937& rng) {
  DenseMatrix w(rows, rank, 0.0);
  if (options.random_init == "rnmf") {
    std::uniform_real_distribution<double> random_uniform(0.0, std::max(scale_x, 1.0));
    for (int i = 0; i < w.rows(); ++i) {
      for (int a = 0; a < w.cols(); ++a) {
        w(i, a) = std::max(random_uniform(rng), 1e-12);
      }
    }
  } else {
    std::uniform_real_distribution<double> random_unit(0.2, 1.0);
    for (int i = 0; i < w.rows(); ++i) {
      for (int a = 0; a < w.cols(); ++a) {
        w(i, a) = std::max(1e-3, scale_x * random_unit(rng));
      }
    }
  }
  return w;
}

// Randomly initialize H under the requested initialization style.
DenseMatrix initialize_h(
    int rank,
    int cols,
    double scale_x,
    const WeightedNmfOptions& options,
    std::mt19937& rng) {
  DenseMatrix h(rank, cols, 0.0);
  if (options.random_init == "rnmf") {
    std::uniform_real_distribution<double> random_uniform(0.0, std::max(scale_x, 1.0));
    for (int a = 0; a < h.rows(); ++a) {
      for (int j = 0; j < h.cols(); ++j) {
        h(a, j) = std::max(random_uniform(rng), 1e-12);
      }
    }
  } else {
    std::uniform_real_distribution<double> random_unit(0.2, 1.0);
    for (int a = 0; a < h.rows(); ++a) {
      for (int j = 0; j < h.cols(); ++j) {
        h(a, j) = std::max(1e-3, scale_x * random_unit(rng));
      }
    }
  }
  return h;
}

// Compute dense row sums used by group-aware initialization.
std::vector<double> dense_row_sums(const DenseMatrix& x) {
  std::vector<double> sums(static_cast<std::size_t>(x.rows()), 0.0);
  for (int i = 0; i < x.rows(); ++i) {
    double total = 0.0;
    for (int j = 0; j < x.cols(); ++j) {
      total += x(i, j);
    }
    sums[static_cast<std::size_t>(i)] = total;
  }
  return sums;
}

// Compute sparse row sums used by group-aware initialization.
std::vector<double> sparse_row_sums(const SparseRowMatrix& x) {
  std::vector<double> sums(static_cast<std::size_t>(x.rows()), 0.0);
  for (int i = 0; i < x.rows(); ++i) {
    double total = 0.0;
    for (int p = x.indptr()[static_cast<std::size_t>(i)];
         p < x.indptr()[static_cast<std::size_t>(i + 1)];
         ++p) {
      total += x.values()[static_cast<std::size_t>(p)];
    }
    sums[static_cast<std::size_t>(i)] = total;
  }
  return sums;
}

// Choose the largest initialization groups when there are more groups than factors.
std::vector<int> select_top_init_groups(const std::vector<int>& init_groups, int rank) {
  std::unordered_map<int, int> counts;
  for (const int group : init_groups) {
    if (group >= 0) {
      ++counts[group];
    }
  }
  std::vector<std::pair<int, int>> ordered(counts.begin(), counts.end());
  std::sort(
      ordered.begin(),
      ordered.end(),
      [](const std::pair<int, int>& lhs, const std::pair<int, int>& rhs) {
        if (lhs.second != rhs.second) {
          return lhs.second > rhs.second;
        }
        return lhs.first < rhs.first;
      });
  if (static_cast<int>(ordered.size()) > rank) {
    ordered.resize(static_cast<std::size_t>(rank));
  }
  std::vector<int> out;
  out.reserve(ordered.size());
  for (const auto& item : ordered) {
    out.push_back(item.first);
  }
  return out;
}

// Decide whether group-informed initialization is informative enough to use.
bool should_use_group_init(const std::vector<int>& init_groups, int rank) {
  return select_top_init_groups(init_groups, rank).size() >= 2U;
}

// Seed H from dense group averages when cluster-aware initialization is requested.
DenseMatrix group_init_h_dense(
    const DenseMatrix& x,
    const std::vector<int>& init_groups,
    DenseMatrix h) {
  const auto top_groups = select_top_init_groups(init_groups, h.rows());
  if (top_groups.size() < 2U) {
    return h;
  }

  std::unordered_map<int, int> factor_by_group;
  for (std::size_t factor = 0; factor < top_groups.size(); ++factor) {
    factor_by_group.emplace(top_groups[factor], static_cast<int>(factor));
  }

  std::vector<int> counts(static_cast<std::size_t>(h.rows()), 0);
  for (int i = 0; i < x.rows(); ++i) {
    const auto it = factor_by_group.find(init_groups[static_cast<std::size_t>(i)]);
    if (it == factor_by_group.end()) {
      continue;
    }
    const int factor = it->second;
    ++counts[static_cast<std::size_t>(factor)];
    for (int j = 0; j < x.cols(); ++j) {
      h(factor, j) += x(i, j);
    }
  }

  for (int factor = 0; factor < static_cast<int>(top_groups.size()); ++factor) {
    const double denom = static_cast<double>(std::max(counts[static_cast<std::size_t>(factor)], 1));
    for (int j = 0; j < h.cols(); ++j) {
      h(factor, j) = std::max(h(factor, j) / denom, 1e-10);
    }
  }
  return h;
}

// Seed H from sparse group averages when cluster-aware initialization is requested.
DenseMatrix group_init_h_sparse(
    const SparseRowMatrix& x,
    const std::vector<int>& init_groups,
    DenseMatrix h) {
  const auto top_groups = select_top_init_groups(init_groups, h.rows());
  if (top_groups.size() < 2U) {
    return h;
  }

  std::unordered_map<int, int> factor_by_group;
  for (std::size_t factor = 0; factor < top_groups.size(); ++factor) {
    factor_by_group.emplace(top_groups[factor], static_cast<int>(factor));
  }

  std::vector<int> counts(static_cast<std::size_t>(h.rows()), 0);
  for (int i = 0; i < x.rows(); ++i) {
    const auto it = factor_by_group.find(init_groups[static_cast<std::size_t>(i)]);
    if (it == factor_by_group.end()) {
      continue;
    }
    const int factor = it->second;
    ++counts[static_cast<std::size_t>(factor)];
    for (int p = x.indptr()[static_cast<std::size_t>(i)];
         p < x.indptr()[static_cast<std::size_t>(i + 1)];
         ++p) {
      h(factor, x.indices()[static_cast<std::size_t>(p)]) += x.values()[static_cast<std::size_t>(p)];
    }
  }

  for (int factor = 0; factor < static_cast<int>(top_groups.size()); ++factor) {
    const double denom = static_cast<double>(std::max(counts[static_cast<std::size_t>(factor)], 1));
    for (int j = 0; j < h.cols(); ++j) {
      h(factor, j) = std::max(h(factor, j) / denom, 1e-10);
    }
  }
  return h;
}

// Initialize W by projecting dense rows onto a seeded H matrix.
DenseMatrix initialized_w_from_h(
    const DenseMatrix& x,
    const DenseMatrix& h,
    const WeightedNmfOptions& options,
    std::mt19937& rng,
    double scale_x) {
  DenseMatrix w = initialize_w(x.rows(), h.rows(), scale_x, options, rng);
  const DenseMatrix scores = project_to_factors(x, h);
  const auto totals = dense_row_sums(x);
  for (int i = 0; i < w.rows(); ++i) {
    const double scale = std::max(totals[static_cast<std::size_t>(i)], 1.0);
    for (int a = 0; a < w.cols(); ++a) {
      w(i, a) = std::max(scores(i, a) * scale, 1e-10);
    }
  }
  return w;
}

// Initialize W by projecting sparse rows onto a seeded H matrix.
DenseMatrix initialized_w_from_h(
    const SparseRowMatrix& x,
    const DenseMatrix& h,
    const WeightedNmfOptions& options,
    std::mt19937& rng,
    double scale_x) {
  DenseMatrix w = initialize_w(x.rows(), h.rows(), scale_x, options, rng);
  const DenseMatrix scores = project_to_factors(x, h);
  const auto totals = sparse_row_sums(x);
  for (int i = 0; i < w.rows(); ++i) {
    const double scale = std::max(totals[static_cast<std::size_t>(i)], 1.0);
    for (int a = 0; a < w.cols(); ++a) {
      w(i, a) = std::max(scores(i, a) * scale, 1e-10);
    }
  }
  return w;
}

// Renormalize factor scales without changing the W*H reconstruction.
void renormalize_factors(DenseMatrix& w, DenseMatrix& h) {
  for (int a = 0; a < h.rows(); ++a) {
    double row_sum = 0.0;
    for (int j = 0; j < h.cols(); ++j) {
      row_sum += h(a, j);
    }
    if (row_sum > 0.0) {
      for (int j = 0; j < h.cols(); ++j) {
        h(a, j) /= row_sum;
      }
      for (int i = 0; i < w.rows(); ++i) {
        w(i, a) *= row_sum;
      }
    }
  }
}

// Apply one guarded multiplicative-update step.
double multiplicative_update(
    double current,
    double numerator,
    double denominator,
    const WeightedNmfOptions& options) {
  if (options.lee_style_epsilon) {
    return std::max(current * numerator, options.update_epsilon) /
        (denominator + options.update_epsilon);
  }
  const double updated = current * numerator / std::max(denominator, options.update_epsilon);
  return std::max(updated, options.update_epsilon);
}

// Form the dense Gram matrix X^T X over rows.
DenseMatrix gram_rows(const DenseMatrix& x) {
  DenseMatrix out(x.cols(), x.cols(), 0.0);
  for (int i = 0; i < x.rows(); ++i) {
    for (int a = 0; a < x.cols(); ++a) {
      const double left = x(i, a);
      for (int b = 0; b < x.cols(); ++b) {
        out(a, b) += left * x(i, b);
      }
    }
  }
  return out;
}

// Form a column-weighted Gram matrix for H updates.
DenseMatrix weighted_gram_rows(const DenseMatrix& h, const std::vector<double>& column_weights) {
  DenseMatrix out(h.rows(), h.rows(), 0.0);
  for (int a = 0; a < h.rows(); ++a) {
    for (int b = 0; b < h.rows(); ++b) {
      double sum = 0.0;
      for (int j = 0; j < h.cols(); ++j) {
        sum += column_weights[static_cast<std::size_t>(j)] * h(a, j) * h(b, j);
      }
      out(a, b) = sum;
    }
  }
  return out;
}

// Update H for dense weighted Euclidean NMF.
void update_h_dense(
    const DenseMatrix& x,
    const std::vector<double>& column_weights,
    const DenseMatrix& w,
    const DenseMatrix& h,
    DenseMatrix& out_h,
    const WeightedNmfOptions& options) {
  DenseMatrix wh = multiply(w, h);
  for (int a = 0; a < h.rows(); ++a) {
    for (int j = 0; j < h.cols(); ++j) {
      double numerator = 0.0;
      double denominator = 0.0;
      for (int i = 0; i < x.rows(); ++i) {
        const double weight = column_weights[static_cast<std::size_t>(j)];
        numerator += w(i, a) * weight * x(i, j);
        denominator += w(i, a) * weight * wh(i, j);
      }
      out_h(a, j) = multiplicative_update(out_h(a, j), numerator, denominator, options);
    }
  }
}

// Update W for dense weighted Euclidean NMF.
void update_w_dense(
    const DenseMatrix& x,
    const std::vector<double>& column_weights,
    const DenseMatrix& w,
    const DenseMatrix& h,
    DenseMatrix& out_w,
    const WeightedNmfOptions& options) {
  DenseMatrix wh = multiply(w, h);
  for (int i = 0; i < w.rows(); ++i) {
    for (int a = 0; a < w.cols(); ++a) {
      double numerator = 0.0;
      double denominator = 0.0;
      for (int j = 0; j < x.cols(); ++j) {
        const double weight = column_weights[static_cast<std::size_t>(j)];
        numerator += h(a, j) * weight * x(i, j);
        denominator += h(a, j) * weight * wh(i, j);
      }
      out_w(i, a) = multiplicative_update(out_w(i, a), numerator, denominator, options);
    }
  }
}

// Compute the sparse numerator needed for Euclidean H updates.
DenseMatrix sparse_numerator_h(
    const SparseRowMatrix& x,
    const DenseMatrix& w,
    const std::vector<double>& column_weights) {
  DenseMatrix numerator(w.cols(), x.cols(), 0.0);
  for (int i = 0; i < x.rows(); ++i) {
    for (int p = x.indptr()[static_cast<std::size_t>(i)];
         p < x.indptr()[static_cast<std::size_t>(i + 1)];
         ++p) {
      const int j = x.indices()[static_cast<std::size_t>(p)];
      const double value = x.values()[static_cast<std::size_t>(p)];
      const double weight = column_weights[static_cast<std::size_t>(j)];
      for (int a = 0; a < w.cols(); ++a) {
        numerator(a, j) += w(i, a) * weight * value;
      }
    }
  }
  return numerator;
}

// Compute the sparse numerator needed for Euclidean W updates.
DenseMatrix sparse_numerator_w(
    const SparseRowMatrix& x,
    const DenseMatrix& h,
    const std::vector<double>& column_weights) {
  DenseMatrix numerator(x.rows(), h.rows(), 0.0);
  for (int i = 0; i < x.rows(); ++i) {
    for (int p = x.indptr()[static_cast<std::size_t>(i)];
         p < x.indptr()[static_cast<std::size_t>(i + 1)];
         ++p) {
      const int j = x.indices()[static_cast<std::size_t>(p)];
      const double value = x.values()[static_cast<std::size_t>(p)];
      const double weight = column_weights[static_cast<std::size_t>(j)];
      for (int a = 0; a < h.rows(); ++a) {
        numerator(i, a) += h(a, j) * weight * value;
      }
    }
  }
  return numerator;
}

// Run one dense Euclidean NMF fit from a single seed.
WeightedNmfResult weighted_nmf_single(
    const DenseMatrix& x,
    const std::vector<double>& column_weights,
    const WeightedNmfOptions& options) {
  std::mt19937 rng(options.seed);
  const double init_scale = options.random_init == "rnmf" ? dense_max(x) : dense_mean(x);
  DenseMatrix h = initialize_h(options.rank, x.cols(), init_scale, options, rng);
  const bool use_group_init = !options.init_groups.empty() &&
      should_use_group_init(options.init_groups, options.rank);
  if (!options.init_groups.empty()) {
    if (static_cast<int>(options.init_groups.size()) != x.rows()) {
      throw std::runtime_error("init_groups length must match NMF matrix rows");
    }
    if (use_group_init) {
      h = group_init_h_dense(x, options.init_groups, std::move(h));
    }
  }
  DenseMatrix w = !use_group_init
      ? initialize_w(x.rows(), options.rank, init_scale, options, rng)
      : initialized_w_from_h(x, h, options, rng, init_scale);

  WeightedNmfResult result;
  double previous_loss = std::numeric_limits<double>::infinity();

  for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
    DenseMatrix wh = multiply(w, h);
    const double loss = weighted_reconstruction_loss(x, wh, column_weights);
    result.losses.push_back(loss);
    if (std::abs(previous_loss - loss) / std::max(loss, 1e-8) < options.tolerance) {
      break;
    }
    previous_loss = loss;

    update_h_dense(x, column_weights, w, h, h, options);
    if (options.renormalize_each_iteration) {
      renormalize_factors(w, h);
    }
    update_w_dense(x, column_weights, w, h, w, options);
  }

  result.w = std::move(w);
  result.h = std::move(h);
  result.selected_seed = options.seed;
  result.selected_run = 0;
  result.candidate_final_losses = {result.losses.empty() ? std::numeric_limits<double>::infinity() : result.losses.back()};
  result.candidate_matched_correlations = {1.0};
  result.selected_factor_stability = std::vector<double>(static_cast<std::size_t>(result.h.rows()), 1.0);
  result.stability_comparison_runs = 0;
  result.stable_factor_count = 0;
  result.stability_threshold = kStableFactorThreshold;
  return result;
}

// Run one sparse Euclidean NMF fit from a single seed.
WeightedNmfResult weighted_nmf_single(
    const SparseRowMatrix& x,
    const std::vector<double>& column_weights,
    const WeightedNmfOptions& options) {
  std::mt19937 rng(options.seed);
  const double init_scale = options.random_init == "rnmf" ? sparse_max(x) : sparse_mean(x);
  DenseMatrix h = initialize_h(options.rank, x.cols(), init_scale, options, rng);
  const bool use_group_init = !options.init_groups.empty() &&
      should_use_group_init(options.init_groups, options.rank);
  if (!options.init_groups.empty()) {
    if (static_cast<int>(options.init_groups.size()) != x.rows()) {
      throw std::runtime_error("init_groups length must match NMF matrix rows");
    }
    if (use_group_init) {
      h = group_init_h_sparse(x, options.init_groups, std::move(h));
    }
  }
  DenseMatrix w = !use_group_init
      ? initialize_w(x.rows(), options.rank, init_scale, options, rng)
      : initialized_w_from_h(x, h, options, rng, init_scale);

  WeightedNmfResult result;
  double previous_loss = std::numeric_limits<double>::infinity();

  // Constant data term of the weighted squared loss.
  double x_squared_weighted = 0.0;
  for (int i = 0; i < x.rows(); ++i) {
    for (int p = x.indptr()[static_cast<std::size_t>(i)];
         p < x.indptr()[static_cast<std::size_t>(i + 1)];
         ++p) {
      const double value = x.values()[static_cast<std::size_t>(p)];
      x_squared_weighted +=
          column_weights[static_cast<std::size_t>(x.indices()[static_cast<std::size_t>(p)])] *
          value * value;
    }
  }

  for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
    // Gram-form loss: sum_ij w_j (x - WH)^2 =
    // tr((W^T W)(H diag(w) H^T)) - 2 sum_nnz w_j x_ij (WH)_ij + sum_nnz w_j x_ij^2,
    // which avoids materializing the dense reconstruction each iteration.
    const DenseMatrix gram_w_loss = gram_rows(w);
    const DenseMatrix gram_h_loss = weighted_gram_rows(h, column_weights);
    double fit_energy = 0.0;
    for (int a = 0; a < h.rows(); ++a) {
      for (int b = 0; b < h.rows(); ++b) {
        fit_energy += gram_w_loss(a, b) * gram_h_loss(a, b);
      }
    }
    double cross = 0.0;
    for (int i = 0; i < x.rows(); ++i) {
      for (int p = x.indptr()[static_cast<std::size_t>(i)];
           p < x.indptr()[static_cast<std::size_t>(i + 1)];
           ++p) {
        const int j = x.indices()[static_cast<std::size_t>(p)];
        double fitted = 0.0;
        for (int a = 0; a < h.rows(); ++a) {
          fitted += w(i, a) * h(a, j);
        }
        cross += column_weights[static_cast<std::size_t>(j)] *
            x.values()[static_cast<std::size_t>(p)] * fitted;
      }
    }
    const double loss = fit_energy - 2.0 * cross + x_squared_weighted;
    result.losses.push_back(loss);
    if (std::abs(previous_loss - loss) / std::max(loss, 1e-8) < options.tolerance) {
      break;
    }
    previous_loss = loss;

    const DenseMatrix numerator_h = sparse_numerator_h(x, w, column_weights);
    const DenseMatrix gram_w = gram_w_loss;
    const DenseMatrix old_h = h;
    for (int a = 0; a < h.rows(); ++a) {
      for (int j = 0; j < h.cols(); ++j) {
        double denominator = 0.0;
        for (int b = 0; b < h.rows(); ++b) {
          denominator += gram_w(a, b) * old_h(b, j);
        }
        denominator *= column_weights[static_cast<std::size_t>(j)];
        h(a, j) = multiplicative_update(old_h(a, j), numerator_h(a, j), denominator, options);
      }
    }

    if (options.renormalize_each_iteration) {
      renormalize_factors(w, h);
    }

    const DenseMatrix numerator_w = sparse_numerator_w(x, h, column_weights);
    const DenseMatrix gram_h = weighted_gram_rows(h, column_weights);
    const DenseMatrix old_w = w;
    for (int i = 0; i < w.rows(); ++i) {
      for (int a = 0; a < w.cols(); ++a) {
        double denominator = 0.0;
        for (int b = 0; b < w.cols(); ++b) {
          denominator += old_w(i, b) * gram_h(a, b);
        }
        w(i, a) = multiplicative_update(old_w(i, a), numerator_w(i, a), denominator, options);
      }
    }
  }

  result.w = std::move(w);
  result.h = std::move(h);
  result.selected_seed = options.seed;
  result.selected_run = 0;
  result.candidate_final_losses = {result.losses.empty() ? std::numeric_limits<double>::infinity() : result.losses.back()};
  result.candidate_matched_correlations = {1.0};
  result.selected_factor_stability = std::vector<double>(static_cast<std::size_t>(result.h.rows()), 1.0);
  result.stability_comparison_runs = 0;
  result.stable_factor_count = 0;
  result.stability_threshold = kStableFactorThreshold;
  return result;
}

// Execute multiple Euclidean NMF restarts and keep the best final loss.
// Cluster initialization (init_groups) is applied to run 0 only; the other
// restarts are fully random and form the stability comparison set.
template<class Runner>
WeightedNmfResult weighted_nmf_multirun(
    const WeightedNmfOptions& options,
    Runner run_once) {
  const int n_runs = std::max(options.n_runs, 1);
  if (n_runs == 1) {
    return run_once(options);
  }

  const bool cluster_run0 = !options.init_groups.empty() &&
      should_use_group_init(options.init_groups, options.rank);

  std::vector<WeightedNmfResult> runs(static_cast<std::size_t>(n_runs));
  std::vector<double> final_losses(static_cast<std::size_t>(n_runs), std::numeric_limits<double>::infinity());
  const int workers = effective_workers(options.num_threads, n_runs);

  auto execute_range = [&](int start, int length) {
    for (int run = start; run < start + length; ++run) {
      auto run_options = options;
      run_options.n_runs = 1;
      run_options.num_threads = 1;
      run_options.seed = options.seed + static_cast<unsigned int>(run);
      if (run > 0) {
        run_options.init_groups.clear();
      }
      auto fit = run_once(run_options);
      fit.selected_seed = run_options.seed;
      fit.selected_run = run;
      final_losses[static_cast<std::size_t>(run)] =
          fit.losses.empty() ? std::numeric_limits<double>::infinity() : fit.losses.back();
      runs[static_cast<std::size_t>(run)] = std::move(fit);
    }
  };

  if (workers <= 1) {
    execute_range(0, n_runs);
  } else {
    subpar::parallelize_range<true>(workers, n_runs, [&](int, int start, int length) {
      execute_range(start, length);
    });
  }

  int best_run = 0;
  double best_loss = final_losses.front();
  for (int run = 1; run < n_runs; ++run) {
    const double candidate = final_losses[static_cast<std::size_t>(run)];
    if (candidate < best_loss) {
      best_loss = candidate;
      best_run = run;
    }
  }

  std::vector<const DenseMatrix*> candidate_h;
  candidate_h.reserve(static_cast<std::size_t>(n_runs));
  std::vector<int> comparison_runs;
  comparison_runs.reserve(static_cast<std::size_t>(n_runs));
  for (int run = 0; run < n_runs; ++run) {
    candidate_h.push_back(&runs[static_cast<std::size_t>(run)].h);
    if (run != best_run && !(cluster_run0 && run == 0)) {
      comparison_runs.push_back(run);
    }
  }
  const auto stability =
      matched_ownership_stability(candidate_h, best_run, comparison_runs);

  std::vector<DenseMatrix> retained_h;
  std::vector<unsigned int> retained_seeds;
  retained_h.reserve(static_cast<std::size_t>(n_runs));
  retained_seeds.reserve(static_cast<std::size_t>(n_runs));
  for (int run = 0; run < n_runs; ++run) {
    retained_h.push_back(runs[static_cast<std::size_t>(run)].h);
    retained_seeds.push_back(options.seed + static_cast<unsigned int>(run));
  }

  auto best = std::move(runs[static_cast<std::size_t>(best_run)]);
  best.candidate_h = std::move(retained_h);
  best.candidate_seeds = std::move(retained_seeds);
  best.selected_run = best_run;
  best.selected_seed = options.seed + static_cast<unsigned int>(best_run);
  best.candidate_final_losses = std::move(final_losses);
  best.candidate_matched_correlations = stability.run_matched_means;
  best.selected_factor_stability = stability.factor_stability;
  best.stability_comparison_runs = stability.comparison_runs;
  best.stable_factor_count = stability.stable_factor_count;
  best.stability_threshold = stability.threshold;
  return best;
}

}  // namespace

// Validate arguments and fit weighted Euclidean NMF on a dense matrix.
WeightedNmfResult weighted_nmf(
    const DenseMatrix& x,
    const std::vector<double>& column_weights,
    const WeightedNmfOptions& options) {
  if (x.rows() == 0 || x.cols() == 0) {
    throw std::runtime_error("weighted_nmf requires a non-empty matrix");
  }
  if (options.rank <= 0) {
    throw std::runtime_error("weighted_nmf rank must be positive");
  }
  if (static_cast<int>(column_weights.size()) != x.cols()) {
    throw std::runtime_error("column_weights length must match matrix columns");
  }
  return weighted_nmf_multirun(options, [&](const WeightedNmfOptions& run_options) {
    return weighted_nmf_single(x, column_weights, run_options);
  });
}

// Validate arguments and fit weighted Euclidean NMF on a sparse matrix.
WeightedNmfResult weighted_nmf(
    const SparseRowMatrix& x,
    const std::vector<double>& column_weights,
    const WeightedNmfOptions& options) {
  if (x.rows() == 0 || x.cols() == 0) {
    throw std::runtime_error("weighted_nmf requires a non-empty matrix");
  }
  if (options.rank <= 0) {
    throw std::runtime_error("weighted_nmf rank must be positive");
  }
  if (static_cast<int>(column_weights.size()) != x.cols()) {
    throw std::runtime_error("column_weights length must match matrix columns");
  }
  return weighted_nmf_multirun(options, [&](const WeightedNmfOptions& run_options) {
    return weighted_nmf_single(x, column_weights, run_options);
  });
}

// Project dense rows to factor scores with row normalization.
DenseMatrix project_to_factors(const DenseMatrix& x, const DenseMatrix& h) {
  DenseMatrix scores = transpose_multiply_right(x, h);
  scores.normalize_rows();
  return scores;
}

// Project sparse rows to factor scores with row normalization.
DenseMatrix project_to_factors(const SparseRowMatrix& x, const DenseMatrix& h) {
  if (x.cols() != h.cols()) {
    throw std::runtime_error("Matrix dimensions do not match for project_to_factors");
  }
  DenseMatrix scores(x.rows(), h.rows(), 0.0);
  for (int i = 0; i < x.rows(); ++i) {
    for (int p = x.indptr()[static_cast<std::size_t>(i)];
         p < x.indptr()[static_cast<std::size_t>(i + 1)];
         ++p) {
      const int j = x.indices()[static_cast<std::size_t>(p)];
      const double value = x.values()[static_cast<std::size_t>(p)];
      for (int a = 0; a < h.rows(); ++a) {
        scores(i, a) += value * h(a, j);
      }
    }
  }
  scores.normalize_rows();
  return scores;
}

}  // namespace celladmix
