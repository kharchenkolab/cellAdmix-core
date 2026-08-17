#include "celladmix/graph.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#include "celladmix/spatial.hpp"
#include "subpar/range.hpp"

namespace celladmix {

// Within-cell graph construction and Potts-style smoothing for transcript labels.

namespace {

// Clamp the requested worker count to the available amount of work.
int effective_threads(int requested, int num_tasks) {
  if (num_tasks <= 0) {
    return 1;
  }
  return std::max(1, subpar::sanitize_num_workers(requested, num_tasks));
}

// Smooth one cell's transcript labels with iterative conditional modes.
std::vector<int> smooth_labels_icm_subset(
    const DenseMatrix& node_scores,
    const std::vector<int>& members,
    const KnnGraph& graph,
    double same_label_ratio,
    int max_iterations) {
  if (members.size() != static_cast<std::size_t>(graph.n_nodes)) {
    throw std::runtime_error("members must match graph nodes");
  }

  const int n_factors = node_scores.cols();
  const double smoothness = std::log(std::max(same_label_ratio, 1.0));
  std::vector<double> log_scores(static_cast<std::size_t>(graph.n_nodes * n_factors), 0.0);
  std::vector<int> labels(static_cast<std::size_t>(graph.n_nodes), 0);

  const auto& score_data = node_scores.data();
  for (int i = 0; i < graph.n_nodes; ++i) {
    const int global = members[static_cast<std::size_t>(i)];
    const std::size_t row_offset = static_cast<std::size_t>(global * n_factors);
    int best = 0;
    double best_value = -std::numeric_limits<double>::infinity();
    for (int factor = 0; factor < n_factors; ++factor) {
      const double value = std::log(std::max(score_data[row_offset + static_cast<std::size_t>(factor)], 1e-12));
      log_scores[static_cast<std::size_t>(i * n_factors + factor)] = value;
      if (value > best_value) {
        best_value = value;
        best = factor;
      }
    }
    labels[static_cast<std::size_t>(i)] = best;
  }

  if (smoothness <= 0.0 || graph.indices.empty()) {
    return labels;
  }

  std::vector<int> neighbor_label_counts(static_cast<std::size_t>(n_factors), 0);
  for (int iteration = 0; iteration < max_iterations; ++iteration) {
    int changed = 0;
    for (int i = 0; i < graph.n_nodes; ++i) {
      std::fill(neighbor_label_counts.begin(), neighbor_label_counts.end(), 0);
      for (int edge = graph.indptr[static_cast<std::size_t>(i)];
           edge < graph.indptr[static_cast<std::size_t>(i + 1)];
           ++edge) {
        const int neighbor = graph.indices[static_cast<std::size_t>(edge)];
        neighbor_label_counts[static_cast<std::size_t>(labels[static_cast<std::size_t>(neighbor)])] += 1;
      }

      int best_label = labels[static_cast<std::size_t>(i)];
      double best_score = -std::numeric_limits<double>::infinity();
      const std::size_t row_offset = static_cast<std::size_t>(i * n_factors);
      for (int factor = 0; factor < n_factors; ++factor) {
        const double score =
            log_scores[row_offset + static_cast<std::size_t>(factor)] +
            smoothness * static_cast<double>(neighbor_label_counts[static_cast<std::size_t>(factor)]);
        if (score > best_score) {
          best_score = score;
          best_label = factor;
        }
      }
      if (best_label != labels[static_cast<std::size_t>(i)]) {
        labels[static_cast<std::size_t>(i)] = best_label;
        changed += 1;
      }
    }
    if (changed == 0) {
      break;
    }
  }

  return labels;
}

}  // namespace

// Build a symmetric transcript KNN graph for a provided set of cell members.
KnnGraph build_cell_knn_graph(const TranscriptTable& table, const std::vector<int>& members, int k) {
  const int n = static_cast<int>(members.size());

  KnnGraph graph;
  graph.n_nodes = n;
  graph.indptr.reserve(static_cast<std::size_t>(n + 1));
  graph.indptr.push_back(0);
  graph.indices.reserve(static_cast<std::size_t>(std::max(0, n * std::min(k, std::max(n - 1, 0)))));

  if (n <= 1 || k <= 0) {
    if (n == 1) {
      graph.indptr.push_back(0);
    }
    return graph;
  }

  std::unordered_map<int, int> local_index;
  local_index.reserve(members.size());
  for (int local_i = 0; local_i < n; ++local_i) {
    local_index.emplace(members[static_cast<std::size_t>(local_i)], local_i);
  }

  std::vector<std::vector<int>> neighbors(static_cast<std::size_t>(n));
  const SpatialKnnIndex knn_index(table, members);
  for (int local_i = 0; local_i < n; ++local_i) {
    const auto hits = knn_index.query(members[static_cast<std::size_t>(local_i)], k, false);
    for (const auto& hit : hits) {
      const int local_j = local_index.at(hit.index);
      neighbors[static_cast<std::size_t>(local_i)].push_back(local_j);
      neighbors[static_cast<std::size_t>(local_j)].push_back(local_i);
    }
  }

  for (int local_i = 0; local_i < n; ++local_i) {
    auto& node_neighbors = neighbors[static_cast<std::size_t>(local_i)];
    std::sort(node_neighbors.begin(), node_neighbors.end());
    node_neighbors.erase(std::unique(node_neighbors.begin(), node_neighbors.end()), node_neighbors.end());
    for (const int local_j : node_neighbors) {
      graph.indices.push_back(local_j);
    }
    graph.indptr.push_back(static_cast<int>(graph.indices.size()));
  }

  return graph;
}

// Convenience overload that builds the within-cell graph from a cell index.
KnnGraph build_cell_knn_graph(const TranscriptTable& table, int cell_index, int k) {
  const auto by_cell = table.transcripts_by_cell();
  return build_cell_knn_graph(table, by_cell.at(static_cast<std::size_t>(cell_index)), k);
}

// Smooth labels on an already-built graph without cell grouping logic.
std::vector<int> smooth_labels_icm(
    const DenseMatrix& node_scores,
    const KnnGraph& graph,
    double same_label_ratio,
    int max_iterations) {
  if (node_scores.rows() != graph.n_nodes) {
    throw std::runtime_error("node_scores rows must match graph nodes");
  }

  const int n_factors = node_scores.cols();
  const double smoothness = std::log(std::max(same_label_ratio, 1.0));
  std::vector<double> log_scores(static_cast<std::size_t>(graph.n_nodes * n_factors), 0.0);
  std::vector<int> labels(static_cast<std::size_t>(graph.n_nodes), 0);
  for (int i = 0; i < graph.n_nodes; ++i) {
    int best = 0;
    double best_value = -std::numeric_limits<double>::infinity();
    const std::size_t row_offset = static_cast<std::size_t>(i * n_factors);
    for (int factor = 0; factor < n_factors; ++factor) {
      const double value = std::log(std::max(node_scores(i, factor), 1e-12));
      log_scores[row_offset + static_cast<std::size_t>(factor)] = value;
      if (value > best_value) {
        best_value = value;
        best = factor;
      }
    }
    labels[static_cast<std::size_t>(i)] = best;
  }

  if (smoothness <= 0.0 || graph.indices.empty()) {
    return labels;
  }

  std::vector<int> neighbor_label_counts(static_cast<std::size_t>(n_factors), 0);
  for (int iteration = 0; iteration < max_iterations; ++iteration) {
    int changed = 0;
    for (int i = 0; i < graph.n_nodes; ++i) {
      std::fill(neighbor_label_counts.begin(), neighbor_label_counts.end(), 0);
      for (int edge = graph.indptr[static_cast<std::size_t>(i)];
           edge < graph.indptr[static_cast<std::size_t>(i + 1)];
           ++edge) {
        const int neighbor = graph.indices[static_cast<std::size_t>(edge)];
        neighbor_label_counts[static_cast<std::size_t>(labels[static_cast<std::size_t>(neighbor)])] += 1;
      }

      int best_label = labels[static_cast<std::size_t>(i)];
      double best_score = -std::numeric_limits<double>::infinity();
      const std::size_t row_offset = static_cast<std::size_t>(i * n_factors);
      for (int factor = 0; factor < n_factors; ++factor) {
        const double score =
            log_scores[row_offset + static_cast<std::size_t>(factor)] +
            smoothness * static_cast<double>(neighbor_label_counts[static_cast<std::size_t>(factor)]);
        if (score > best_score) {
          best_score = score;
          best_label = factor;
        }
      }
      if (best_label != labels[static_cast<std::size_t>(i)]) {
        labels[static_cast<std::size_t>(i)] = best_label;
        changed += 1;
      }
    }
    if (changed == 0) {
      break;
    }
  }

  return labels;
}

// Build every cell's smoothing graph once for reuse across score sets.
CellGraphs build_cell_graphs(
    const TranscriptTable& table,
    int k_neighbors,
    int num_threads) {
  CellGraphs out;
  out.by_cell = table.transcripts_by_cell();
  out.graphs.resize(out.by_cell.size());
  if (k_neighbors <= 0) {
    return out;
  }
  const int n_cells = static_cast<int>(out.by_cell.size());
  auto build_range = [&](int start, int length) {
    for (int cell = start; cell < start + length; ++cell) {
      const auto& members = out.by_cell[static_cast<std::size_t>(cell)];
      if (members.size() > 1) {
        out.graphs[static_cast<std::size_t>(cell)] =
            build_cell_knn_graph(table, members, k_neighbors);
      }
    }
  };
  if (num_threads <= 1 || n_cells <= 1) {
    build_range(0, n_cells);
  } else {
    const int workers = effective_threads(num_threads, n_cells);
    subpar::parallelize_range<true>(workers, n_cells, [&](int, int start, int length) {
      build_range(start, length);
    });
  }
  return out;
}

// Assign transcript factors reusing prebuilt per-cell graphs.
std::vector<int> assign_factors_per_cell(
    const TranscriptTable& table,
    const DenseMatrix& node_scores,
    const CellGraphs& cell_graphs,
    double same_label_ratio,
    int max_iterations,
    int num_threads) {
  if (node_scores.rows() != static_cast<int>(table.size())) {
    throw std::runtime_error("node_scores rows must match number of transcripts");
  }
  if (cell_graphs.by_cell.size() != cell_graphs.graphs.size()) {
    throw std::runtime_error("cell_graphs members and graphs must align");
  }

  std::vector<int> labels(table.size(), 0);
  const bool smoothing_active = max_iterations > 0 && same_label_ratio > 1.0;
  const int n_cells = static_cast<int>(cell_graphs.by_cell.size());
  auto assign_range = [&](int start, int length) {
    for (int cell = start; cell < start + length; ++cell) {
      const auto& members = cell_graphs.by_cell[static_cast<std::size_t>(cell)];
      if (members.empty()) {
        continue;
      }
      const auto& graph = cell_graphs.graphs[static_cast<std::size_t>(cell)];
      if (members.size() == 1 || !smoothing_active || graph.n_nodes == 0) {
        for (const int member : members) {
          labels[static_cast<std::size_t>(member)] = node_scores.row_argmax(member);
        }
        continue;
      }
      const auto cell_labels = smooth_labels_icm_subset(
          node_scores, members, graph, same_label_ratio, max_iterations);
      for (std::size_t local = 0; local < members.size(); ++local) {
        labels[static_cast<std::size_t>(members[local])] = cell_labels[local];
      }
    }
  };
  if (num_threads <= 1 || n_cells <= 1) {
    assign_range(0, n_cells);
  } else {
    const int workers = effective_threads(num_threads, n_cells);
    subpar::parallelize_range<true>(workers, n_cells, [&](int, int start, int length) {
      assign_range(start, length);
    });
  }
  return labels;
}

// Assign transcript factors cell-by-cell and optionally smooth them within each cell.
std::vector<int> assign_factors_per_cell(
    const TranscriptTable& table,
    const DenseMatrix& node_scores,
    int k_neighbors,
    double same_label_ratio,
    int max_iterations,
    int num_threads) {
  if (node_scores.rows() != static_cast<int>(table.size())) {
    throw std::runtime_error("node_scores rows must match number of transcripts");
  }

  const auto by_cell = table.transcripts_by_cell();
  std::vector<int> labels(table.size(), 0);

  if (k_neighbors <= 0 || max_iterations <= 0 || same_label_ratio <= 1.0) {
    if (num_threads <= 1 || table.size() <= 1) {
      for (std::size_t i = 0; i < table.size(); ++i) {
        labels[i] = node_scores.row_argmax(static_cast<int>(i));
      }
    } else {
      const int workers = effective_threads(num_threads, static_cast<int>(table.size()));
      subpar::parallelize_range<true>(workers, static_cast<int>(table.size()), [&](int, int start, int length) {
        for (int i = start; i < start + length; ++i) {
          labels[static_cast<std::size_t>(i)] = node_scores.row_argmax(i);
        }
      });
    }
    return labels;
  }

  if (num_threads <= 1 || by_cell.size() <= 1) {
    for (std::size_t cell = 0; cell < by_cell.size(); ++cell) {
      const auto& members = by_cell[cell];
      if (members.empty()) {
        continue;
      }
      if (members.size() == 1) {
        labels[static_cast<std::size_t>(members.front())] = node_scores.row_argmax(members.front());
        continue;
      }

      const auto graph = build_cell_knn_graph(table, members, k_neighbors);
      const auto cell_labels =
          smooth_labels_icm_subset(node_scores, members, graph, same_label_ratio, max_iterations);
      for (std::size_t local = 0; local < members.size(); ++local) {
        labels[static_cast<std::size_t>(members[local])] = cell_labels[local];
      }
    }
    return labels;
  }

  const int workers = effective_threads(num_threads, static_cast<int>(by_cell.size()));
  subpar::parallelize_range<true>(workers, static_cast<int>(by_cell.size()), [&](int, int start, int length) {
    for (int cell = start; cell < start + length; ++cell) {
      const auto& members = by_cell[static_cast<std::size_t>(cell)];
      if (members.empty()) {
        continue;
      }
      if (members.size() == 1) {
        labels[static_cast<std::size_t>(members.front())] = node_scores.row_argmax(members.front());
        continue;
      }

      const auto graph = build_cell_knn_graph(table, members, k_neighbors);
      const auto cell_labels =
          smooth_labels_icm_subset(node_scores, members, graph, same_label_ratio, max_iterations);
      for (std::size_t local = 0; local < members.size(); ++local) {
        labels[static_cast<std::size_t>(members[local])] = cell_labels[local];
      }
    }
  });
  return labels;
}

}  // namespace celladmix
