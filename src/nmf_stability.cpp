#include "celladmix/nmf_stability.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace celladmix {

DenseMatrix ownership_profiles(const DenseMatrix& h) {
  DenseMatrix profiles(h.rows(), h.cols(), 0.0);
  for (int j = 0; j < h.cols(); ++j) {
    double total = 0.0;
    for (int a = 0; a < h.rows(); ++a) {
      total += h(a, j);
    }
    if (total <= 0.0) {
      continue;
    }
    for (int a = 0; a < h.rows(); ++a) {
      profiles(a, j) = h(a, j) / total;
    }
  }
  return profiles;
}

namespace {

// Center and unit-normalize matrix rows; zero-variance rows become all-zero.
DenseMatrix standardize_rows(const DenseMatrix& x) {
  DenseMatrix out(x.rows(), x.cols(), 0.0);
  for (int i = 0; i < x.rows(); ++i) {
    double mean = 0.0;
    for (int j = 0; j < x.cols(); ++j) {
      mean += x(i, j);
    }
    mean /= static_cast<double>(std::max(x.cols(), 1));
    double ss = 0.0;
    for (int j = 0; j < x.cols(); ++j) {
      const double delta = x(i, j) - mean;
      out(i, j) = delta;
      ss += delta * delta;
    }
    if (ss <= 1e-30) {
      for (int j = 0; j < x.cols(); ++j) {
        out(i, j) = 0.0;
      }
      continue;
    }
    const double inv_norm = 1.0 / std::sqrt(ss);
    for (int j = 0; j < x.cols(); ++j) {
      out(i, j) *= inv_norm;
    }
  }
  return out;
}

}  // namespace

DenseMatrix row_correlation_matrix(const DenseMatrix& lhs, const DenseMatrix& rhs) {
  if (lhs.cols() != rhs.cols()) {
    throw std::runtime_error("row_correlation_matrix requires equal row widths");
  }
  const DenseMatrix a = standardize_rows(lhs);
  const DenseMatrix b = standardize_rows(rhs);
  DenseMatrix out(lhs.rows(), rhs.rows(), 0.0);
  for (int i = 0; i < a.rows(); ++i) {
    for (int k = 0; k < b.rows(); ++k) {
      double dot = 0.0;
      for (int j = 0; j < a.cols(); ++j) {
        dot += a(i, j) * b(k, j);
      }
      out(i, k) = dot;
    }
  }
  return out;
}

std::vector<int> hungarian_match(const DenseMatrix& similarity) {
  const int n = similarity.rows();
  if (similarity.cols() != n) {
    throw std::runtime_error("hungarian_match requires a square matrix");
  }
  if (n == 0) {
    return {};
  }

  // Potentials formulation of the Hungarian algorithm on cost = -similarity,
  // 1-based with a virtual row/column 0.
  const double kInf = std::numeric_limits<double>::infinity();
  std::vector<double> u(static_cast<std::size_t>(n) + 1, 0.0);
  std::vector<double> v(static_cast<std::size_t>(n) + 1, 0.0);
  std::vector<int> p(static_cast<std::size_t>(n) + 1, 0);
  std::vector<int> way(static_cast<std::size_t>(n) + 1, 0);
  for (int i = 1; i <= n; ++i) {
    p[0] = i;
    int j0 = 0;
    std::vector<double> min_v(static_cast<std::size_t>(n) + 1, kInf);
    std::vector<char> used(static_cast<std::size_t>(n) + 1, 0);
    do {
      used[static_cast<std::size_t>(j0)] = 1;
      const int i0 = p[static_cast<std::size_t>(j0)];
      int j1 = -1;
      double delta = kInf;
      for (int j = 1; j <= n; ++j) {
        if (used[static_cast<std::size_t>(j)]) {
          continue;
        }
        const double cur = -similarity(i0 - 1, j - 1) -
            u[static_cast<std::size_t>(i0)] - v[static_cast<std::size_t>(j)];
        if (cur < min_v[static_cast<std::size_t>(j)]) {
          min_v[static_cast<std::size_t>(j)] = cur;
          way[static_cast<std::size_t>(j)] = j0;
        }
        if (min_v[static_cast<std::size_t>(j)] < delta) {
          delta = min_v[static_cast<std::size_t>(j)];
          j1 = j;
        }
      }
      for (int j = 0; j <= n; ++j) {
        if (used[static_cast<std::size_t>(j)]) {
          u[static_cast<std::size_t>(p[static_cast<std::size_t>(j)])] += delta;
          v[static_cast<std::size_t>(j)] -= delta;
        } else {
          min_v[static_cast<std::size_t>(j)] -= delta;
        }
      }
      j0 = j1;
    } while (p[static_cast<std::size_t>(j0)] != 0);
    do {
      const int j1 = way[static_cast<std::size_t>(j0)];
      p[static_cast<std::size_t>(j0)] = p[static_cast<std::size_t>(j1)];
      j0 = j1;
    } while (j0 != 0);
  }

  std::vector<int> match(static_cast<std::size_t>(n), -1);
  for (int j = 1; j <= n; ++j) {
    if (p[static_cast<std::size_t>(j)] > 0) {
      match[static_cast<std::size_t>(p[static_cast<std::size_t>(j)]) - 1] = j - 1;
    }
  }
  return match;
}

FactorStabilityResult matched_ownership_stability(
    const std::vector<const DenseMatrix*>& candidate_h,
    int selected_run,
    const std::vector<int>& comparison_runs,
    double stable_threshold) {
  const int n_runs = static_cast<int>(candidate_h.size());
  if (selected_run < 0 || selected_run >= n_runs) {
    throw std::runtime_error("matched_ownership_stability: selected_run out of range");
  }
  const DenseMatrix reference = ownership_profiles(*candidate_h[selected_run]);
  const int rank = reference.rows();

  FactorStabilityResult result;
  result.threshold = stable_threshold;
  result.run_matched_means.assign(static_cast<std::size_t>(n_runs), 0.0);
  result.run_matched_means[static_cast<std::size_t>(selected_run)] = 1.0;

  std::vector<char> in_comparison(static_cast<std::size_t>(n_runs), 0);
  for (const int run : comparison_runs) {
    if (run < 0 || run >= n_runs || run == selected_run) {
      throw std::runtime_error("matched_ownership_stability: invalid comparison run");
    }
    in_comparison[static_cast<std::size_t>(run)] = 1;
  }

  std::vector<double> totals(static_cast<std::size_t>(rank), 0.0);
  int n_compared = 0;
  for (int run = 0; run < n_runs; ++run) {
    if (run == selected_run) {
      continue;
    }
    const DenseMatrix profiles = ownership_profiles(*candidate_h[run]);
    if (profiles.rows() != rank || profiles.cols() != reference.cols()) {
      throw std::runtime_error("matched_ownership_stability: candidate shape mismatch");
    }
    const DenseMatrix correlations = row_correlation_matrix(reference, profiles);
    const auto match = hungarian_match(correlations);
    double run_total = 0.0;
    for (int factor = 0; factor < rank; ++factor) {
      const double value = correlations(factor, match[static_cast<std::size_t>(factor)]);
      run_total += value;
      if (in_comparison[static_cast<std::size_t>(run)]) {
        totals[static_cast<std::size_t>(factor)] += value;
      }
    }
    result.run_matched_means[static_cast<std::size_t>(run)] =
        run_total / static_cast<double>(std::max(rank, 1));
    if (in_comparison[static_cast<std::size_t>(run)]) {
      ++n_compared;
    }
  }

  result.comparison_runs = n_compared;
  if (n_compared == 0) {
    result.factor_stability.assign(static_cast<std::size_t>(rank), 1.0);
  } else {
    result.factor_stability.resize(static_cast<std::size_t>(rank));
    for (int factor = 0; factor < rank; ++factor) {
      result.factor_stability[static_cast<std::size_t>(factor)] =
          totals[static_cast<std::size_t>(factor)] / static_cast<double>(n_compared);
    }
  }
  result.stable_factor_count = 0;
  for (const double value : result.factor_stability) {
    if (n_compared > 0 && value >= stable_threshold) {
      ++result.stable_factor_count;
    }
  }
  return result;
}

}  // namespace celladmix
