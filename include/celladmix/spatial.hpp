// Exact transcript-level spatial kNN wrapper used by NCV construction and CRF smoothing.

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include "celladmix/types.hpp"

#ifndef NANOFLANN_NO_THREADS
#define NANOFLANN_NO_THREADS
#endif

#include "third_party/nanoflann.hpp"

namespace celladmix {

struct NeighborHit {
  int index = -1;
  double squared_distance = 0.0;
};

namespace detail {

struct TranscriptPointAdaptor {
  const TranscriptTable* table = nullptr;
  const std::vector<int>* indices = nullptr;

  TranscriptPointAdaptor(const TranscriptTable& table_in, const std::vector<int>* indices_in)
      : table(&table_in), indices(indices_in) {}

  inline std::size_t kdtree_get_point_count() const {
    return indices == nullptr ? table->size() : indices->size();
  }

  inline double kdtree_get_pt(const std::size_t idx, const std::size_t dim) const {
    const std::size_t global_idx =
        indices == nullptr ? idx : static_cast<std::size_t>(indices->at(idx));
    if (dim == 0) {
      return table->x[global_idx];
    }
    if (dim == 1) {
      return table->y[global_idx];
    }
    return table->z[global_idx];
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const {
    return false;
  }
};

using TranscriptKdTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, TranscriptPointAdaptor>,
    TranscriptPointAdaptor,
    3,
    int>;

}  // namespace detail

class SpatialKnnIndex {
 public:
  // Build an exact KD-tree over all transcript rows in the table.
  explicit SpatialKnnIndex(const TranscriptTable& table)
      : table_(&table), adaptor_(table, nullptr), tree_(3, adaptor_, nanoflann::KDTreeSingleIndexAdaptorParams(10)) {
    tree_.buildIndex();
  }

  // Build an exact KD-tree over a provided subset of transcript rows.
  SpatialKnnIndex(const TranscriptTable& table, std::vector<int> indices)
      : table_(&table),
        indices_(std::move(indices)),
        adaptor_(table, &indices_),
        tree_(3, adaptor_, nanoflann::KDTreeSingleIndexAdaptorParams(10)) {
    tree_.buildIndex();
  }

  // Return the number of indexed transcript rows.
  std::size_t size() const {
    return indices_.empty() ? table_->size() : indices_.size();
  }

  // Report whether the index contains any transcript rows.
  bool empty() const { return size() == 0; }

  // Query exact nearest transcript neighbors for one transcript row.
  std::vector<NeighborHit> query(int query_index, int k, bool include_self) const {
    if (k <= 0 || empty()) {
      return {};
    }

    const int tree_size = static_cast<int>(size());
    const int requested = std::min(tree_size, k + (include_self ? 0 : 1));
    std::vector<int> internal_indices(static_cast<std::size_t>(requested), -1);
    std::vector<double> squared_distances(static_cast<std::size_t>(requested), 0.0);

    const std::array<double, 3> point = {
        table_->x[static_cast<std::size_t>(query_index)],
        table_->y[static_cast<std::size_t>(query_index)],
        table_->z[static_cast<std::size_t>(query_index)]};

    nanoflann::KNNResultSet<double, int> result_set(requested);
    result_set.init(internal_indices.data(), squared_distances.data());
    tree_.findNeighbors(result_set, point.data(), nanoflann::SearchParameters(0.0f, true));

    std::vector<NeighborHit> out;
    out.reserve(static_cast<std::size_t>(std::min(k, tree_size)));
    for (int i = 0; i < requested; ++i) {
      const int internal = internal_indices[static_cast<std::size_t>(i)];
      const int global_index = indices_.empty()
          ? internal
          : indices_[static_cast<std::size_t>(internal)];
      if (!include_self && global_index == query_index) {
        continue;
      }
      out.push_back({global_index, squared_distances[static_cast<std::size_t>(i)]});
    }

    std::sort(out.begin(), out.end(), [](const NeighborHit& lhs, const NeighborHit& rhs) {
      if (lhs.squared_distance == rhs.squared_distance) {
        return lhs.index < rhs.index;
      }
      return lhs.squared_distance < rhs.squared_distance;
    });

    if (static_cast<int>(out.size()) > k) {
      out.resize(static_cast<std::size_t>(k));
    }
    return out;
  }

  // Query exact nearest transcript neighbors for an arbitrary spatial point.
  std::vector<NeighborHit> query_point(double x, double y, double z, int k) const {
    if (k <= 0 || empty()) {
      return {};
    }

    const int tree_size = static_cast<int>(size());
    const int requested = std::min(tree_size, k);
    std::vector<int> internal_indices(static_cast<std::size_t>(requested), -1);
    std::vector<double> squared_distances(static_cast<std::size_t>(requested), 0.0);

    const std::array<double, 3> point = {x, y, z};
    nanoflann::KNNResultSet<double, int> result_set(requested);
    result_set.init(internal_indices.data(), squared_distances.data());
    tree_.findNeighbors(result_set, point.data(), nanoflann::SearchParameters(0.0f, true));

    std::vector<NeighborHit> out;
    out.reserve(static_cast<std::size_t>(requested));
    for (int i = 0; i < requested; ++i) {
      const int internal = internal_indices[static_cast<std::size_t>(i)];
      const int global_index = indices_.empty()
          ? internal
          : indices_[static_cast<std::size_t>(internal)];
      out.push_back({global_index, squared_distances[static_cast<std::size_t>(i)]});
    }

    std::sort(out.begin(), out.end(), [](const NeighborHit& lhs, const NeighborHit& rhs) {
      if (lhs.squared_distance == rhs.squared_distance) {
        return lhs.index < rhs.index;
      }
      return lhs.squared_distance < rhs.squared_distance;
    });
    return out;
  }

 private:
  const TranscriptTable* table_ = nullptr;
  std::vector<int> indices_;
  detail::TranscriptPointAdaptor adaptor_;
  detail::TranscriptKdTree tree_;
};

}  // namespace celladmix
