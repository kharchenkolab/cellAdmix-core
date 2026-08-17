#include "celladmix/nmf_kl.hpp"

#include "celladmix/nmf_stability.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "subpar/range.hpp"

namespace celladmix {
namespace {

// Sparse KL-NMF implementation used by the current production path.

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
    const SparseNmfOptions& options,
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
    const SparseNmfOptions& options,
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

// Decide whether cluster-aware initialization has enough groups to be useful.
bool should_use_group_init(const std::vector<int>& init_groups, int rank) {
  return select_top_init_groups(init_groups, rank).size() >= 2U;
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

// Project sparse rows onto H without row normalization to seed W.
DenseMatrix project_to_factors_unscaled(const SparseRowMatrix& x, const DenseMatrix& h) {
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
  return scores;
}

// Evaluate one W*H entry with numerical guarding.
double fitted_entry(
    const DenseMatrix& w,
    const DenseMatrix& h,
    int row,
    int col,
    double epsilon) {
  double fitted = 0.0;
  for (int a = 0; a < h.rows(); ++a) {
    fitted += w(row, a) * h(a, col);
  }
  return std::max(fitted, epsilon);
}

// Initialize W by projecting sparse rows onto a seeded H matrix.
DenseMatrix initialized_w_from_h(
    const SparseRowMatrix& x,
    const DenseMatrix& h,
    const SparseNmfOptions& options,
    std::mt19937& rng,
    double scale_x) {
  DenseMatrix w = initialize_w(x.rows(), h.rows(), scale_x, options, rng);
  const DenseMatrix scores = project_to_factors_unscaled(x, h);
  const auto totals = sparse_row_sums(x);
  for (int i = 0; i < w.rows(); ++i) {
    const double scale = std::max(totals[static_cast<std::size_t>(i)], 1.0);
    for (int a = 0; a < w.cols(); ++a) {
      w(i, a) = std::max(scores(i, a) * scale, 1e-10);
    }
  }
  return w;
}

// Apply one guarded multiplicative-update step.
double multiplicative_update(
    double current,
    double numerator,
    double denominator,
    const SparseNmfOptions& options) {
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

// Form a column-weighted Gram matrix for optional Euclidean updates.
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

// Precompute row sums of H for KL W updates.
std::vector<double> h_column_sums(const DenseMatrix& h) {
  std::vector<double> sums(static_cast<std::size_t>(h.cols()), 0.0);
  for (int j = 0; j < h.cols(); ++j) {
    double total = 0.0;
    for (int a = 0; a < h.rows(); ++a) {
      total += h(a, j);
    }
    sums[static_cast<std::size_t>(j)] = total;
  }
  return sums;
}

// Evaluate optional penalties that discourage broad or redundant factors.
double h_penalty_value(const DenseMatrix& h, const SparseNmfOptions& options) {
  double penalty = 0.0;
  if (options.h_l1_penalty > 0.0) {
    for (int a = 0; a < h.rows(); ++a) {
      for (int j = 0; j < h.cols(); ++j) {
        penalty += options.h_l1_penalty * h(a, j);
      }
    }
  }
  if (options.h_diversity_penalty > 0.0) {
    for (int j = 0; j < h.cols(); ++j) {
      for (int a = 0; a < h.rows(); ++a) {
        for (int b = a + 1; b < h.rows(); ++b) {
          penalty += options.h_diversity_penalty * h(a, j) * h(b, j);
        }
      }
    }
  }
  return penalty;
}

// Evaluate the weighted Euclidean fallback objective for the sparse solver.
double weighted_euclidean_objective(
    const SparseRowMatrix& x,
    const DenseMatrix& reconstruction,
    const std::vector<double>& column_weights,
    const DenseMatrix& h,
    const SparseNmfOptions& options) {
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
  return loss + h_penalty_value(h, options);
}

// Evaluate generalized KL divergence without materializing a dense reconstruction.
double generalized_kl_objective(
    const SparseRowMatrix& x,
    const DenseMatrix& w,
    const DenseMatrix& h,
    const SparseNmfOptions& options) {
  const auto w_sums = w.col_sums();
  const auto h_sums = h.row_sums();
  double loss = 0.0;
  for (int a = 0; a < h.rows(); ++a) {
    loss += std::max(w_sums[static_cast<std::size_t>(a)], 0.0) *
        std::max(h_sums[static_cast<std::size_t>(a)], 0.0);
  }

  const int workers = effective_workers(options.num_threads, x.rows());
  if (workers <= 1) {
    for (int i = 0; i < x.rows(); ++i) {
      for (int p = x.indptr()[static_cast<std::size_t>(i)];
           p < x.indptr()[static_cast<std::size_t>(i + 1)];
           ++p) {
        const int j = x.indices()[static_cast<std::size_t>(p)];
        const double value = x.values()[static_cast<std::size_t>(p)];
        const double fitted = fitted_entry(w, h, i, j, options.update_epsilon);
        loss += value * (std::log(value / fitted) - 1.0);
      }
    }
  } else {
    std::vector<double> partial_losses(static_cast<std::size_t>(workers), 0.0);
    subpar::parallelize_range<true>(workers, x.rows(), [&](int worker, int start, int length) {
      double partial = 0.0;
      for (int i = start; i < start + length; ++i) {
        for (int p = x.indptr()[static_cast<std::size_t>(i)];
             p < x.indptr()[static_cast<std::size_t>(i + 1)];
             ++p) {
          const int j = x.indices()[static_cast<std::size_t>(p)];
          const double value = x.values()[static_cast<std::size_t>(p)];
          const double fitted = fitted_entry(w, h, i, j, options.update_epsilon);
          partial += value * (std::log(value / fitted) - 1.0);
        }
      }
      partial_losses[static_cast<std::size_t>(worker)] = partial;
    });
    loss += std::accumulate(partial_losses.begin(), partial_losses.end(), 0.0);
  }
  return loss + h_penalty_value(h, options);
}

std::pair<double, double> mean_and_sd(const std::vector<double>& values) {
  if (values.empty()) {
    return {0.0, 0.0};
  }
  const double mean = std::accumulate(values.begin(), values.end(), 0.0) /
      static_cast<double>(values.size());
  double ss = 0.0;
  for (const double value : values) {
    const double delta = value - mean;
    ss += delta * delta;
  }
  return {mean, std::sqrt(ss / static_cast<double>(values.size()))};
}

// Build the KL numerator for H updates from sparse rows.
DenseMatrix generalized_kl_numerator_h(
    const SparseRowMatrix& x,
    const DenseMatrix& w,
    const DenseMatrix& h,
    const SparseNmfOptions& options) {
  DenseMatrix numerator_h(h.rows(), h.cols(), 0.0);
  const int workers = effective_workers(options.num_threads, x.rows());
  if (workers <= 1) {
    for (int i = 0; i < x.rows(); ++i) {
      for (int p = x.indptr()[static_cast<std::size_t>(i)];
           p < x.indptr()[static_cast<std::size_t>(i + 1)];
           ++p) {
        const int j = x.indices()[static_cast<std::size_t>(p)];
        const double value = x.values()[static_cast<std::size_t>(p)];
        const double ratio = value / fitted_entry(w, h, i, j, options.update_epsilon);
        for (int a = 0; a < w.cols(); ++a) {
          numerator_h(a, j) += w(i, a) * ratio;
        }
      }
    }
    return numerator_h;
  }

  std::vector<DenseMatrix> partials;
  partials.reserve(static_cast<std::size_t>(workers));
  for (int worker = 0; worker < workers; ++worker) {
    partials.emplace_back(h.rows(), h.cols(), 0.0);
  }
  subpar::parallelize_range<true>(workers, x.rows(), [&](int worker, int start, int length) {
    DenseMatrix& partial = partials[static_cast<std::size_t>(worker)];
    for (int i = start; i < start + length; ++i) {
      for (int p = x.indptr()[static_cast<std::size_t>(i)];
           p < x.indptr()[static_cast<std::size_t>(i + 1)];
           ++p) {
        const int j = x.indices()[static_cast<std::size_t>(p)];
        const double value = x.values()[static_cast<std::size_t>(p)];
        const double ratio = value / fitted_entry(w, h, i, j, options.update_epsilon);
        for (int a = 0; a < w.cols(); ++a) {
          partial(a, j) += w(i, a) * ratio;
        }
      }
    }
  });
  for (const auto& partial : partials) {
    for (int a = 0; a < numerator_h.rows(); ++a) {
      for (int j = 0; j < numerator_h.cols(); ++j) {
        numerator_h(a, j) += partial(a, j);
      }
    }
  }
  return numerator_h;
}

// Build the KL numerator for W updates from sparse rows.
DenseMatrix generalized_kl_numerator_w(
    const SparseRowMatrix& x,
    const DenseMatrix& w,
    const DenseMatrix& h,
    const SparseNmfOptions& options) {
  DenseMatrix numerator_w(x.rows(), h.rows(), 0.0);
  const int workers = effective_workers(options.num_threads, x.rows());
  if (workers <= 1) {
    for (int i = 0; i < x.rows(); ++i) {
      for (int p = x.indptr()[static_cast<std::size_t>(i)];
           p < x.indptr()[static_cast<std::size_t>(i + 1)];
           ++p) {
        const int j = x.indices()[static_cast<std::size_t>(p)];
        const double value = x.values()[static_cast<std::size_t>(p)];
        const double ratio = value / fitted_entry(w, h, i, j, options.update_epsilon);
        for (int a = 0; a < h.rows(); ++a) {
          numerator_w(i, a) += h(a, j) * ratio;
        }
      }
    }
    return numerator_w;
  }

  subpar::parallelize_range<true>(workers, x.rows(), [&](int, int start, int length) {
    for (int i = start; i < start + length; ++i) {
      for (int p = x.indptr()[static_cast<std::size_t>(i)];
           p < x.indptr()[static_cast<std::size_t>(i + 1)];
           ++p) {
        const int j = x.indices()[static_cast<std::size_t>(p)];
        const double value = x.values()[static_cast<std::size_t>(p)];
        const double ratio = value / fitted_entry(w, h, i, j, options.update_epsilon);
        for (int a = 0; a < h.rows(); ++a) {
          numerator_w(i, a) += h(a, j) * ratio;
        }
      }
    }
  });
  return numerator_w;
}

// Run one sparse NMF fit from a single seed under KL or Euclidean loss.
SparseNmfResult sparse_nmf_single(
    const SparseRowMatrix& x,
    const std::vector<double>& column_weights,
    const SparseNmfOptions& options) {
  std::mt19937 rng(options.seed);
  const double init_scale = options.random_init == "rnmf" ? sparse_max(x) : sparse_mean(x);
  DenseMatrix h = initialize_h(options.rank, x.cols(), init_scale, options, rng);
  const bool use_group_init = options.init_mode == "cluster" &&
      !options.init_groups.empty() && should_use_group_init(options.init_groups, options.rank);
  if (options.init_mode == "cluster") {
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

  SparseNmfResult result;
  double previous_loss = std::numeric_limits<double>::infinity();

  for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
    DenseMatrix wh;
    const double loss = options.loss_mode == "kl"
        ? generalized_kl_objective(x, w, h, options)
        : (wh = multiply(w, h), weighted_euclidean_objective(x, wh, column_weights, h, options));
    result.losses.push_back(loss);
    if (std::abs(previous_loss - loss) / std::max(loss, 1e-8) < options.tolerance) {
      break;
    }
    previous_loss = loss;

    if (options.loss_mode == "kl") {
      const std::vector<double> w_sums = w.col_sums();
      DenseMatrix numerator_h = generalized_kl_numerator_h(x, w, h, options);
      const auto col_sums = h_column_sums(h);
      for (int a = 0; a < h.rows(); ++a) {
        for (int j = 0; j < h.cols(); ++j) {
          const double old_value = h(a, j);
          double denominator = std::max(w_sums[static_cast<std::size_t>(a)], options.update_epsilon);
          denominator += options.h_l1_penalty;
          if (options.h_diversity_penalty > 0.0) {
            denominator += options.h_diversity_penalty *
                std::max(col_sums[static_cast<std::size_t>(j)] - old_value, 0.0);
          }
          h(a, j) = multiplicative_update(old_value, numerator_h(a, j), denominator, options);
        }
      }

      DenseMatrix numerator_w = generalized_kl_numerator_w(x, w, h, options);
      const auto h_sums = h.row_sums();
      for (int i = 0; i < w.rows(); ++i) {
        for (int a = 0; a < w.cols(); ++a) {
          const double old_value = w(i, a);
          w(i, a) = multiplicative_update(
              old_value,
              numerator_w(i, a),
              std::max(h_sums[static_cast<std::size_t>(a)], options.update_epsilon),
              options);
        }
      }
    } else {
      DenseMatrix numerator_h(h.rows(), x.cols(), 0.0);
      for (int i = 0; i < x.rows(); ++i) {
        for (int p = x.indptr()[static_cast<std::size_t>(i)];
             p < x.indptr()[static_cast<std::size_t>(i + 1)];
             ++p) {
          const int j = x.indices()[static_cast<std::size_t>(p)];
          const double value = x.values()[static_cast<std::size_t>(p)];
          const double weight = column_weights[static_cast<std::size_t>(j)];
          for (int a = 0; a < w.cols(); ++a) {
            numerator_h(a, j) += w(i, a) * weight * value;
          }
        }
      }

      const DenseMatrix gram_w = gram_rows(w);
      const DenseMatrix old_h = h;
      const auto col_sums = h_column_sums(old_h);
      for (int a = 0; a < h.rows(); ++a) {
        for (int j = 0; j < h.cols(); ++j) {
          double denominator = 0.0;
          for (int b = 0; b < h.rows(); ++b) {
            denominator += gram_w(a, b) * old_h(b, j);
          }
          denominator *= column_weights[static_cast<std::size_t>(j)];
          denominator += options.h_l1_penalty;
          if (options.h_diversity_penalty > 0.0) {
            denominator += options.h_diversity_penalty *
                std::max(col_sums[static_cast<std::size_t>(j)] - old_h(a, j), 0.0);
          }
          h(a, j) = multiplicative_update(old_h(a, j), numerator_h(a, j), denominator, options);
        }
      }

      DenseMatrix numerator_w(x.rows(), h.rows(), 0.0);
      for (int i = 0; i < x.rows(); ++i) {
        for (int p = x.indptr()[static_cast<std::size_t>(i)];
             p < x.indptr()[static_cast<std::size_t>(i + 1)];
             ++p) {
          const int j = x.indices()[static_cast<std::size_t>(p)];
          const double value = x.values()[static_cast<std::size_t>(p)];
          const double weight = column_weights[static_cast<std::size_t>(j)];
          for (int a = 0; a < h.rows(); ++a) {
            numerator_w(i, a) += h(a, j) * weight * value;
          }
        }
      }
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
  }

  result.w = std::move(w);
  result.h = std::move(h);
  result.selected_seed = options.seed;
  result.selected_run = 0;
  result.final_objective =
      result.losses.empty() ? std::numeric_limits<double>::infinity() : result.losses.back();
  result.candidate_final_objectives = {result.final_objective};
  result.candidate_best_match_correlations = {1.0};
  result.selected_factor_stability = std::vector<double>(static_cast<std::size_t>(result.h.rows()), 1.0);
  result.candidate_final_objective_mean = result.final_objective;
  result.candidate_final_objective_sd = 0.0;
  result.candidate_best_match_correlation_mean = 1.0;
  result.stability_comparison_runs = 0;
  result.stable_factor_count = 0;
  result.stability_threshold = kStableFactorThreshold;
  return result;
}

// Execute multiple sparse NMF restarts and keep the best final objective.
// When cluster initialization is active it is applied to run 0 only; the
// remaining restarts are fully random so that the stability diagnostic
// averages over independent starts (the cluster run stays a selectable
// candidate but is excluded from the stability comparison set).
template <class Runner>
SparseNmfResult sparse_nmf_multirun(
    const SparseNmfOptions& options,
    Runner run_once) {
  const int n_runs = std::max(options.n_runs, 1);
  if (n_runs == 1) {
    return run_once(options);
  }

  const bool cluster_run0 = options.init_mode == "cluster" &&
      !options.init_groups.empty() &&
      should_use_group_init(options.init_groups, options.rank);

  std::vector<SparseNmfResult> runs(static_cast<std::size_t>(n_runs));
  std::vector<double> final_objectives(
      static_cast<std::size_t>(n_runs),
      std::numeric_limits<double>::infinity());
  const int workers = effective_workers(options.num_threads, n_runs);

  auto execute_range = [&](int start, int length) {
    for (int run = start; run < start + length; ++run) {
      auto run_options = options;
      run_options.n_runs = 1;
      run_options.num_threads = 1;
      run_options.seed = options.seed + static_cast<unsigned int>(run);
      if (run > 0) {
        run_options.init_mode = "random";
        run_options.init_groups.clear();
      }
      auto fit = run_once(run_options);
      fit.selected_seed = run_options.seed;
      fit.selected_run = run;
      final_objectives[static_cast<std::size_t>(run)] = fit.final_objective;
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
  double best_objective = final_objectives.front();
  for (int run = 1; run < n_runs; ++run) {
    if (final_objectives[static_cast<std::size_t>(run)] < best_objective) {
      best_objective = final_objectives[static_cast<std::size_t>(run)];
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
  const auto objective_stats = mean_and_sd(final_objectives);
  const auto correlation_stats = mean_and_sd(stability.run_matched_means);

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
  best.final_objective = best_objective;
  best.candidate_final_objectives = std::move(final_objectives);
  best.candidate_best_match_correlations = stability.run_matched_means;
  best.selected_factor_stability = stability.factor_stability;
  best.candidate_final_objective_mean = objective_stats.first;
  best.candidate_final_objective_sd = objective_stats.second;
  best.candidate_best_match_correlation_mean = correlation_stats.first;
  best.stability_comparison_runs = stability.comparison_runs;
  best.stable_factor_count = stability.stable_factor_count;
  best.stability_threshold = stability.threshold;
  return best;
}

}  // namespace

// Validate arguments and fit the sparse production NMF solver.
SparseNmfResult sparse_nmf(
    const SparseRowMatrix& x,
    const std::vector<double>& column_weights,
    const SparseNmfOptions& options) {
  if (x.rows() == 0 || x.cols() == 0) {
    throw std::runtime_error("sparse_nmf requires a non-empty matrix");
  }
  if (options.rank <= 0) {
    throw std::runtime_error("sparse_nmf rank must be positive");
  }
  if (options.init_mode != "random" && options.init_mode != "cluster") {
    throw std::runtime_error("sparse_nmf init_mode must be 'random' or 'cluster'");
  }
  if (options.loss_mode != "euclidean" && options.loss_mode != "kl") {
    throw std::runtime_error("sparse_nmf loss_mode must be 'euclidean' or 'kl'");
  }
  if (options.random_init != "legacy" && options.random_init != "rnmf") {
    throw std::runtime_error("sparse_nmf random_init must be 'legacy' or 'rnmf'");
  }
  if (options.loss_mode == "euclidean" &&
      static_cast<int>(column_weights.size()) != x.cols()) {
    throw std::runtime_error("column_weights length must match matrix columns");
  }
  return sparse_nmf_multirun(options, [&](const SparseNmfOptions& run_options) {
    return sparse_nmf_single(x, column_weights, run_options);
  });
}

// Reorder factors so F1 is the largest component in the selected training W.
void order_nmf_factors_by_training_importance(SparseNmfResult& result) {
  const int rank = result.h.rows();
  if (rank <= 1) {
    return;
  }
  if (result.w.cols() != rank) {
    throw std::runtime_error("Cannot order NMF factors: W columns and H rows do not match");
  }

  const auto masses = result.w.col_sums();
  std::vector<int> order(static_cast<std::size_t>(rank));
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](int left, int right) {
    const double left_mass = masses[static_cast<std::size_t>(left)];
    const double right_mass = masses[static_cast<std::size_t>(right)];
    if (left_mass != right_mass) {
      return left_mass > right_mass;
    }
    return left < right;
  });

  DenseMatrix ordered_w(result.w.rows(), result.w.cols(), 0.0);
  DenseMatrix ordered_h(result.h.rows(), result.h.cols(), 0.0);
  for (int new_factor = 0; new_factor < rank; ++new_factor) {
    const int old_factor = order[static_cast<std::size_t>(new_factor)];
    for (int row = 0; row < result.w.rows(); ++row) {
      ordered_w(row, new_factor) = result.w(row, old_factor);
    }
    for (int col = 0; col < result.h.cols(); ++col) {
      ordered_h(new_factor, col) = result.h(old_factor, col);
    }
  }
  result.w = std::move(ordered_w);
  result.h = std::move(ordered_h);

  if (!result.selected_factor_stability.empty()) {
    if (static_cast<int>(result.selected_factor_stability.size()) != rank) {
      throw std::runtime_error("Cannot order NMF factors: stability vector length does not match rank");
    }
    std::vector<double> ordered_stability(static_cast<std::size_t>(rank), 0.0);
    for (int new_factor = 0; new_factor < rank; ++new_factor) {
      const int old_factor = order[static_cast<std::size_t>(new_factor)];
      ordered_stability[static_cast<std::size_t>(new_factor)] =
          result.selected_factor_stability[static_cast<std::size_t>(old_factor)];
    }
    result.selected_factor_stability = std::move(ordered_stability);
  }

  // Keep the member pool's selected entry identical to the reordered h so
  // that ensemble member selected_run shares the run's factor numbering.
  if (result.selected_run >= 0 &&
      static_cast<std::size_t>(result.selected_run) < result.candidate_h.size()) {
    result.candidate_h[static_cast<std::size_t>(result.selected_run)] = result.h;
  }
}

}  // namespace celladmix
