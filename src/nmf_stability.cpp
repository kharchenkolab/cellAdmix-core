#include "celladmix/nmf_stability.hpp"

#include <algorithm>
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

FactorMatchResult match_factor_profiles(
    const DenseMatrix& profiles,
    const DenseMatrix& reference) {
  if (profiles.cols() != reference.cols()) {
    throw std::runtime_error("match_factor_profiles: gene dimensions differ");
  }
  const DenseMatrix lhs = ownership_profiles(profiles);
  const DenseMatrix rhs = ownership_profiles(reference);
  const DenseMatrix correlations = row_correlation_matrix(lhs, rhs);
  const int n = std::max(profiles.rows(), reference.rows());
  DenseMatrix padded(n, n, 0.0);
  for (int i = 0; i < profiles.rows(); ++i) {
    for (int j = 0; j < reference.rows(); ++j) {
      padded(i, j) = correlations(i, j);
    }
  }
  const auto assignment = hungarian_match(padded);
  FactorMatchResult out;
  out.match.assign(static_cast<std::size_t>(profiles.rows()), -1);
  out.correlation.assign(static_cast<std::size_t>(profiles.rows()), 0.0);
  for (int i = 0; i < profiles.rows(); ++i) {
    const int j = assignment[static_cast<std::size_t>(i)];
    if (j < reference.rows()) {
      out.match[static_cast<std::size_t>(i)] = j;
      out.correlation[static_cast<std::size_t>(i)] = correlations(i, j);
    }
  }
  return out;
}

FactorClusterResult cluster_factor_profiles(
    const DenseMatrix& profiles,
    const std::vector<int>& group,
    double min_correlation) {
  const int n = profiles.rows();
  if (static_cast<int>(group.size()) != n) {
    throw std::runtime_error("cluster_factor_profiles: one group per profile row is required");
  }
  FactorClusterResult out;
  if (n == 0) {
    return out;
  }
  // Ownership shares normalized within each factorization.
  DenseMatrix ownership(n, profiles.cols(), 0.0);
  std::vector<int> groups(group);
  std::sort(groups.begin(), groups.end());
  groups.erase(std::unique(groups.begin(), groups.end()), groups.end());
  for (const int g : groups) {
    std::vector<int> rows;
    for (int i = 0; i < n; ++i) {
      if (group[static_cast<std::size_t>(i)] == g) {
        rows.push_back(i);
      }
    }
    DenseMatrix sub(static_cast<int>(rows.size()), profiles.cols(), 0.0);
    for (std::size_t r = 0; r < rows.size(); ++r) {
      for (int j = 0; j < profiles.cols(); ++j) {
        sub(static_cast<int>(r), j) = profiles(rows[r], j);
      }
    }
    const DenseMatrix own = ownership_profiles(sub);
    for (std::size_t r = 0; r < rows.size(); ++r) {
      for (int j = 0; j < profiles.cols(); ++j) {
        ownership(rows[r], j) = own(static_cast<int>(r), j);
      }
    }
  }
  const DenseMatrix corr = row_correlation_matrix(ownership, ownership);

  // Average-linkage agglomeration with Lance-Williams distance updates and
  // cached nearest neighbors.
  const double max_distance = 1.0 - min_correlation;
  const std::size_t nn_ = static_cast<std::size_t>(n);
  std::vector<double> dist(nn_ * nn_, 0.0);
  auto d = [&](int i, int j) -> double& { return dist[static_cast<std::size_t>(i) * nn_ + static_cast<std::size_t>(j)]; };
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      d(i, j) = 1.0 - corr(i, j);
    }
  }
  std::vector<int> size(nn_, 1);
  std::vector<bool> alive(nn_, true);
  std::vector<int> root(nn_);
  for (int i = 0; i < n; ++i) {
    root[static_cast<std::size_t>(i)] = i;
  }
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<int> nn(nn_, -1);
  std::vector<double> nn_dist(nn_, inf);
  auto refresh = [&](int i) {
    nn[static_cast<std::size_t>(i)] = -1;
    nn_dist[static_cast<std::size_t>(i)] = inf;
    for (int j = 0; j < n; ++j) {
      if (j == i || !alive[static_cast<std::size_t>(j)]) {
        continue;
      }
      if (d(i, j) < nn_dist[static_cast<std::size_t>(i)]) {
        nn_dist[static_cast<std::size_t>(i)] = d(i, j);
        nn[static_cast<std::size_t>(i)] = j;
      }
    }
  };
  for (int i = 0; i < n; ++i) {
    refresh(i);
  }
  int clusters = n;
  while (clusters > 1) {
    int a = -1;
    double best = inf;
    for (int i = 0; i < n; ++i) {
      if (alive[static_cast<std::size_t>(i)] && nn_dist[static_cast<std::size_t>(i)] < best) {
        best = nn_dist[static_cast<std::size_t>(i)];
        a = i;
      }
    }
    if (a < 0 || best > max_distance) {
      break;
    }
    int b = nn[static_cast<std::size_t>(a)];
    if (a > b) {
      std::swap(a, b);
    }
    const double wa = static_cast<double>(size[static_cast<std::size_t>(a)]);
    const double wb = static_cast<double>(size[static_cast<std::size_t>(b)]);
    for (int k = 0; k < n; ++k) {
      if (k == a || k == b || !alive[static_cast<std::size_t>(k)]) {
        continue;
      }
      const double merged = (wa * d(a, k) + wb * d(b, k)) / (wa + wb);
      d(a, k) = merged;
      d(k, a) = merged;
    }
    size[static_cast<std::size_t>(a)] += size[static_cast<std::size_t>(b)];
    alive[static_cast<std::size_t>(b)] = false;
    --clusters;
    for (int i = 0; i < n; ++i) {
      if (root[static_cast<std::size_t>(i)] == b) {
        root[static_cast<std::size_t>(i)] = a;
      }
    }
    refresh(a);
    for (int k = 0; k < n; ++k) {
      if (k == a || !alive[static_cast<std::size_t>(k)]) {
        continue;
      }
      const int current = nn[static_cast<std::size_t>(k)];
      if (current == a || current == b) {
        refresh(k);
      } else if (d(k, a) < nn_dist[static_cast<std::size_t>(k)]) {
        nn_dist[static_cast<std::size_t>(k)] = d(k, a);
        nn[static_cast<std::size_t>(k)] = a;
      }
    }
  }

  // Number clusters by first appearance and summarize within-cluster agreement.
  std::vector<int> relabel(nn_, -1);
  out.cluster.assign(nn_, -1);
  int next = 0;
  for (int i = 0; i < n; ++i) {
    const int r = root[static_cast<std::size_t>(i)];
    if (relabel[static_cast<std::size_t>(r)] < 0) {
      relabel[static_cast<std::size_t>(r)] = next++;
    }
    out.cluster[static_cast<std::size_t>(i)] = relabel[static_cast<std::size_t>(r)];
  }
  out.n_clusters = next;
  out.within_correlation.assign(nn_, std::numeric_limits<double>::quiet_NaN());
  for (int i = 0; i < n; ++i) {
    double total = 0.0;
    int count = 0;
    for (int j = 0; j < n; ++j) {
      if (j != i && out.cluster[static_cast<std::size_t>(j)] == out.cluster[static_cast<std::size_t>(i)]) {
        total += corr(i, j);
        ++count;
      }
    }
    if (count > 0) {
      out.within_correlation[static_cast<std::size_t>(i)] = total / static_cast<double>(count);
    }
  }
  return out;
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
