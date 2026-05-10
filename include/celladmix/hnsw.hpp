// Shared hnswlib-backed nearest-neighbor helpers for cell and report UMAPs.

#pragma once

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <queue>
#include <utility>
#include <vector>

#include "knncolle/knncolle.hpp"
#include "subpar/subpar.hpp"
#include "third_party/hnswlib/hnswlib.h"

namespace celladmix {

struct HnswNeighborOptions {
  int m = 16;
  int ef_construction = 200;
  int ef_search = 100;
  unsigned int seed = 100U;
  int exact_threshold = 128;
};

// Convert Eigen column observations into hnswlib's row-major float layout.
inline std::vector<float> hnsw_row_major_from_columns(const Eigen::MatrixXd& data) {
  const int dim = static_cast<int>(data.rows());
  const int nobs = static_cast<int>(data.cols());
  std::vector<float> out(static_cast<std::size_t>(dim * nobs));
  for (int obs = 0; obs < nobs; ++obs) {
    for (int d = 0; d < dim; ++d) {
      out[static_cast<std::size_t>(obs * dim + d)] = static_cast<float>(data(d, obs));
    }
  }
  return out;
}

// Exact fallback for tiny matrices where HNSW setup costs more than it saves.
inline knncolle::NeighborList<int, double> exact_l2_neighbors_from_columns(
    const Eigen::MatrixXd& data,
    int k) {
  const int dim = static_cast<int>(data.rows());
  const int nobs = static_cast<int>(data.cols());
  k = std::max(0, std::min(k, nobs - 1));
  knncolle::NeighborList<int, double> out(static_cast<std::size_t>(nobs));
  if (k <= 0) {
    return out;
  }

  std::vector<std::pair<double, int>> distances;
  distances.reserve(static_cast<std::size_t>(std::max(0, nobs - 1)));
  for (int obs = 0; obs < nobs; ++obs) {
    distances.clear();
    for (int other = 0; other < nobs; ++other) {
      if (other == obs) {
        continue;
      }
      double squared = 0.0;
      for (int d = 0; d < dim; ++d) {
        const double delta = data(d, obs) - data(d, other);
        squared += delta * delta;
      }
      distances.push_back({squared, other});
    }
    if (static_cast<int>(distances.size()) > k) {
      std::partial_sort(
          distances.begin(),
          distances.begin() + k,
          distances.end(),
          [](const auto& lhs, const auto& rhs) {
            if (lhs.first == rhs.first) {
              return lhs.second < rhs.second;
            }
            return lhs.first < rhs.first;
          });
      distances.resize(static_cast<std::size_t>(k));
    } else {
      std::sort(distances.begin(), distances.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.first == rhs.first) {
          return lhs.second < rhs.second;
        }
        return lhs.first < rhs.first;
      });
    }

    auto& current = out[static_cast<std::size_t>(obs)];
    current.reserve(distances.size());
    for (const auto& hit : distances) {
      current.push_back({hit.second, hit.first});
    }
  }
  return out;
}

// Build and query an HNSW L2 index over Eigen column observations.
inline knncolle::NeighborList<int, double> hnsw_l2_neighbors_from_columns(
    const Eigen::MatrixXd& data,
    int k,
    int num_threads,
    HnswNeighborOptions options = {}) {
  const int dim = static_cast<int>(data.rows());
  const int nobs = static_cast<int>(data.cols());
  k = std::max(0, std::min(k, nobs - 1));
  if (nobs == 0 || dim == 0 || k <= 0) {
    return knncolle::NeighborList<int, double>(static_cast<std::size_t>(nobs));
  }
  if (nobs <= std::max(1, options.exact_threshold)) {
    return exact_l2_neighbors_from_columns(data, k);
  }

  const int workers = std::max(1, std::min(num_threads, nobs));
  const int query_k = std::min(k + 1, nobs);
  auto row_major = hnsw_row_major_from_columns(data);

  hnswlib::L2Space space(dim);
  hnswlib::HierarchicalNSW<float> index(
      &space,
      static_cast<std::size_t>(nobs),
      static_cast<std::size_t>(std::max(2, options.m)),
      static_cast<std::size_t>(std::max(options.ef_construction, options.m)),
      static_cast<std::size_t>(options.seed));

  subpar::parallelize_range(workers, nobs, [&](int, int start, int length) {
    for (int obs = start, end = start + length; obs < end; ++obs) {
      index.addPoint(
          row_major.data() + static_cast<std::size_t>(obs * dim),
          static_cast<hnswlib::labeltype>(obs));
    }
  });
  index.setEf(static_cast<std::size_t>(std::max(options.ef_search, k + 1)));

  knncolle::NeighborList<int, double> out(static_cast<std::size_t>(nobs));
  subpar::parallelize_range(workers, nobs, [&](int, int start, int length) {
    for (int obs = start, end = start + length; obs < end; ++obs) {
      auto queue = index.searchKnn(
          row_major.data() + static_cast<std::size_t>(obs * dim),
          static_cast<std::size_t>(query_k));

      std::vector<std::pair<double, int>> reversed;
      reversed.reserve(queue.size());
      while (!queue.empty()) {
        reversed.push_back({
            static_cast<double>(queue.top().first),
            static_cast<int>(queue.top().second)});
        queue.pop();
      }

      auto& current = out[static_cast<std::size_t>(obs)];
      current.reserve(static_cast<std::size_t>(k));
      for (auto it = reversed.rbegin(); it != reversed.rend() && static_cast<int>(current.size()) < k; ++it) {
        if (it->second != obs) {
          current.push_back({it->second, it->first});
        }
      }
    }
  });
  return out;
}

}  // namespace celladmix
