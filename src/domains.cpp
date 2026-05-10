#include "celladmix/domains.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "celladmix/spatial.hpp"
#include "subpar/range.hpp"

namespace celladmix {

namespace {

struct DomainNeighbor {
  int index = -1;
  double distance = 0.0;
};

struct KmeansResult {
  std::vector<int> labels;
  DenseMatrix centroids;
  int iterations = 0;
};

int effective_threads(int requested, int tasks) {
  if (tasks <= 0) {
    return 1;
  }
  return std::max(1, subpar::sanitize_num_workers(requested, tasks));
}

void validate_cells(const CellTable& cells) {
  const std::size_t n = cells.cell_ids.size();
  if (n == 0) {
    throw std::runtime_error("domain detection requires at least one cell");
  }
  if (cells.centroid_x.size() != n || cells.centroid_y.size() != n) {
    throw std::runtime_error("cell coordinates must match cell ids for domain detection");
  }
  if (!cells.centroid_z.empty() && cells.centroid_z.size() != n) {
    throw std::runtime_error("cell z coordinates must be empty or match cell ids");
  }
}

int max_neighbor_count(const DomainOptions& options) {
  int out = std::max(0, options.smooth_k);
  for (const int scale : options.scales) {
    out = std::max(out, scale);
  }
  return out;
}

TranscriptTable cells_as_points(const CellTable& cells) {
  TranscriptTable points;
  points.x = cells.centroid_x;
  points.y = cells.centroid_y;
  points.z = cells.centroid_z.empty()
      ? std::vector<double>(cells.cell_ids.size(), 0.0)
      : cells.centroid_z;
  return points;
}

std::vector<std::vector<int>> group_cells_by_sample(const CellTable& cells) {
  const int n = static_cast<int>(cells.cell_ids.size());
  std::unordered_map<std::string, int> sample_index;
  std::vector<std::vector<int>> groups;
  groups.reserve(1);
  for (int i = 0; i < n; ++i) {
    std::string key = "sample";
    if (cells.sample_ids.size() == cells.cell_ids.size() && !cells.sample_ids[static_cast<std::size_t>(i)].empty()) {
      key = cells.sample_ids[static_cast<std::size_t>(i)];
    } else if (cells.fov_ids.size() == cells.cell_ids.size() && !cells.fov_ids[static_cast<std::size_t>(i)].empty()) {
      key = cells.fov_ids[static_cast<std::size_t>(i)];
    }
    auto it = sample_index.find(key);
    if (it == sample_index.end()) {
      const int idx = static_cast<int>(groups.size());
      sample_index.emplace(key, idx);
      groups.emplace_back();
      it = sample_index.find(key);
    }
    groups[static_cast<std::size_t>(it->second)].push_back(i);
  }
  return groups;
}

std::vector<std::vector<DomainNeighbor>> build_cell_knn(
    const CellTable& cells,
    int max_k,
    double max_edge_distance,
    int num_threads) {
  const int n = static_cast<int>(cells.cell_ids.size());
  std::vector<std::vector<DomainNeighbor>> out(static_cast<std::size_t>(n));
  if (n <= 1 || max_k <= 0) {
    return out;
  }

  const TranscriptTable points = cells_as_points(cells);
  const auto groups = group_cells_by_sample(cells);
  const bool has_distance_limit = std::isfinite(max_edge_distance) && max_edge_distance > 0.0;
  const int workers = effective_threads(num_threads, n);

  for (const auto& members : groups) {
    if (members.size() <= 1) {
      continue;
    }
    const int group_k = std::min<int>(max_k, static_cast<int>(members.size()) - 1);
    const SpatialKnnIndex index(points, members);
    subpar::parallelize_range(workers, static_cast<int>(members.size()), [&](int, int start, int length) {
      for (int p = start, end = start + length; p < end; ++p) {
        const int cell = members[static_cast<std::size_t>(p)];
        const auto hits = index.query(cell, group_k, false);
        auto& current = out[static_cast<std::size_t>(cell)];
        current.reserve(hits.size());
        for (const auto& hit : hits) {
          const double distance = std::sqrt(std::max(0.0, hit.squared_distance));
          if (has_distance_limit && distance > max_edge_distance) {
            continue;
          }
          current.push_back({hit.index, distance});
        }
      }
    });
  }
  return out;
}

std::vector<std::vector<int>> undirected_prefix_graph(
    const std::vector<std::vector<DomainNeighbor>>& neighbors,
    int k) {
  const int n = static_cast<int>(neighbors.size());
  std::vector<std::vector<int>> graph(static_cast<std::size_t>(n));
  if (k <= 0) {
    return graph;
  }
  for (int i = 0; i < n; ++i) {
    const auto& current = neighbors[static_cast<std::size_t>(i)];
    const int limit = std::min<int>(k, current.size());
    for (int p = 0; p < limit; ++p) {
      const int j = current[static_cast<std::size_t>(p)].index;
      if (j < 0 || j >= n || j == i) {
        continue;
      }
      graph[static_cast<std::size_t>(i)].push_back(j);
      graph[static_cast<std::size_t>(j)].push_back(i);
    }
  }
  for (auto& current : graph) {
    std::sort(current.begin(), current.end());
    current.erase(std::unique(current.begin(), current.end()), current.end());
  }
  return graph;
}

DenseMatrix build_multiscale_features(
    const std::vector<std::vector<DomainNeighbor>>& neighbors,
    const std::vector<int>& annotation_ids,
    int n_types,
    const DomainOptions& options) {
  const int n = static_cast<int>(annotation_ids.size());
  const int n_scales = static_cast<int>(options.scales.size());
  const int composition_cols = n_scales * n_types;
  const int total_cols = composition_cols + (options.include_center_type ? n_types : 0);
  DenseMatrix features(n, total_cols, 0.0);

  for (int i = 0; i < n; ++i) {
    for (int s = 0; s < n_scales; ++s) {
      const int scale = std::max(0, options.scales[static_cast<std::size_t>(s)]);
      const int offset = s * n_types;
      const auto& current = neighbors[static_cast<std::size_t>(i)];
      const int limit = std::min<int>(scale, current.size());
      double total = 0.0;

      if (options.include_self) {
        const int type = annotation_ids[static_cast<std::size_t>(i)];
        if (type >= 0 && type < n_types) {
          features(i, offset + type) += options.self_weight;
          total += options.self_weight;
        }
      }

      for (int p = 0; p < limit; ++p) {
        const int neighbor = current[static_cast<std::size_t>(p)].index;
        const int type = annotation_ids[static_cast<std::size_t>(neighbor)];
        if (type >= 0 && type < n_types) {
          features(i, offset + type) += 1.0;
          total += 1.0;
        }
      }

      if (total > 0.0) {
        for (int type = 0; type < n_types; ++type) {
          features(i, offset + type) /= total;
        }
      }
    }

    if (options.include_center_type) {
      const int type = annotation_ids[static_cast<std::size_t>(i)];
      if (type >= 0 && type < n_types) {
        features(i, composition_cols + type) = options.center_type_weight;
      }
    }
  }
  return features;
}

void transform_composition_blocks(DenseMatrix& features, int n_types, int n_blocks, const std::string& transform) {
  if (transform == "none") {
    return;
  }
  if (transform != "clr" && transform != "hellinger") {
    throw std::runtime_error("domain feature transform must be one of clr, hellinger, or none");
  }
  constexpr double pseudocount = 1e-3;
  for (int row = 0; row < features.rows(); ++row) {
    for (int block = 0; block < n_blocks; ++block) {
      const int offset = block * n_types;
      if (transform == "hellinger") {
        for (int type = 0; type < n_types; ++type) {
          features(row, offset + type) = std::sqrt(std::max(0.0, features(row, offset + type)));
        }
      } else {
        double mean_log = 0.0;
        for (int type = 0; type < n_types; ++type) {
          mean_log += std::log(features(row, offset + type) + pseudocount);
        }
        mean_log /= static_cast<double>(n_types);
        for (int type = 0; type < n_types; ++type) {
          features(row, offset + type) = std::log(features(row, offset + type) + pseudocount) - mean_log;
        }
      }
    }
  }
}

double squared_distance_row_centroid(const DenseMatrix& features, int row, const DenseMatrix& centroids, int centroid) {
  double out = 0.0;
  for (int col = 0; col < features.cols(); ++col) {
    const double delta = features(row, col) - centroids(centroid, col);
    out += delta * delta;
  }
  return out;
}

DenseMatrix compute_centroids(const DenseMatrix& features, const std::vector<int>& labels, int k) {
  DenseMatrix centroids(k, features.cols(), 0.0);
  std::vector<int> sizes(static_cast<std::size_t>(k), 0);
  for (int row = 0; row < features.rows(); ++row) {
    const int label = labels[static_cast<std::size_t>(row)];
    if (label < 0 || label >= k) {
      continue;
    }
    sizes[static_cast<std::size_t>(label)] += 1;
    for (int col = 0; col < features.cols(); ++col) {
      centroids(label, col) += features(row, col);
    }
  }
  for (int label = 0; label < k; ++label) {
    const double denom = static_cast<double>(std::max(1, sizes[static_cast<std::size_t>(label)]));
    for (int col = 0; col < features.cols(); ++col) {
      centroids(label, col) /= denom;
    }
  }
  return centroids;
}

KmeansResult run_kmeans(const DenseMatrix& features, int k, int max_iterations, unsigned int seed, int num_threads) {
  const int n = features.rows();
  if (n == 0 || k <= 0) {
    throw std::runtime_error("kmeans requires non-empty features and positive k");
  }
  k = std::min(k, n);
  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> uniform_cell(0, n - 1);

  DenseMatrix centroids(k, features.cols(), 0.0);
  std::vector<double> min_distance(static_cast<std::size_t>(n), std::numeric_limits<double>::infinity());
  int first = uniform_cell(rng);
  for (int col = 0; col < features.cols(); ++col) {
    centroids(0, col) = features(first, col);
  }

  for (int center = 1; center < k; ++center) {
    double total = 0.0;
    for (int row = 0; row < n; ++row) {
      min_distance[static_cast<std::size_t>(row)] =
          std::min(min_distance[static_cast<std::size_t>(row)],
                   squared_distance_row_centroid(features, row, centroids, center - 1));
      total += min_distance[static_cast<std::size_t>(row)];
    }
    int chosen = uniform_cell(rng);
    if (total > 0.0 && std::isfinite(total)) {
      std::uniform_real_distribution<double> picker(0.0, total);
      double target = picker(rng);
      for (int row = 0; row < n; ++row) {
        target -= min_distance[static_cast<std::size_t>(row)];
        if (target <= 0.0) {
          chosen = row;
          break;
        }
      }
    }
    for (int col = 0; col < features.cols(); ++col) {
      centroids(center, col) = features(chosen, col);
    }
  }

  std::vector<int> labels(static_cast<std::size_t>(n), 0);
  const int workers = effective_threads(num_threads, n);
  int iterations = 0;
  for (int iter = 0; iter < max_iterations; ++iter) {
    iterations = iter + 1;
    std::vector<int> changed_by_worker(static_cast<std::size_t>(workers), 0);
    subpar::parallelize_range(workers, n, [&](int worker, int start, int length) {
      int changed = 0;
      for (int row = start, end = start + length; row < end; ++row) {
        int best = labels[static_cast<std::size_t>(row)];
        double best_distance = std::numeric_limits<double>::infinity();
        for (int center = 0; center < k; ++center) {
          const double distance = squared_distance_row_centroid(features, row, centroids, center);
          if (distance < best_distance) {
            best_distance = distance;
            best = center;
          }
        }
        if (best != labels[static_cast<std::size_t>(row)]) {
          labels[static_cast<std::size_t>(row)] = best;
          changed += 1;
        }
      }
      changed_by_worker[static_cast<std::size_t>(worker)] = changed;
    });
    const int changed = std::accumulate(changed_by_worker.begin(), changed_by_worker.end(), 0);
    centroids = compute_centroids(features, labels, k);
    if (changed == 0) {
      break;
    }
  }
  return {std::move(labels), std::move(centroids), iterations};
}

std::vector<int> potts_smooth(
    const DenseMatrix& features,
    std::vector<int> labels,
    int k,
    const std::vector<std::vector<int>>& graph,
    double lambda,
    int max_iterations,
    int* iterations_out) {
  int iterations = 0;
  if (lambda <= 0.0 || graph.empty()) {
    if (iterations_out != nullptr) {
      *iterations_out = iterations;
    }
    return labels;
  }
  DenseMatrix centroids = compute_centroids(features, labels, k);
  for (int iter = 0; iter < max_iterations; ++iter) {
    iterations = iter + 1;
    int changed = 0;
    for (int row = 0; row < features.rows(); ++row) {
      int best = labels[static_cast<std::size_t>(row)];
      double best_cost = std::numeric_limits<double>::infinity();
      for (int label = 0; label < k; ++label) {
        double disagreement = 0.0;
        for (const int neighbor : graph[static_cast<std::size_t>(row)]) {
          disagreement += labels[static_cast<std::size_t>(neighbor)] == label ? 0.0 : 1.0;
        }
        const double cost =
            squared_distance_row_centroid(features, row, centroids, label) +
            lambda * disagreement;
        if (cost < best_cost) {
          best_cost = cost;
          best = label;
        }
      }
      if (best != labels[static_cast<std::size_t>(row)]) {
        labels[static_cast<std::size_t>(row)] = best;
        changed += 1;
      }
    }
    centroids = compute_centroids(features, labels, k);
    if (changed == 0) {
      break;
    }
  }
  if (iterations_out != nullptr) {
    *iterations_out = iterations;
  }
  return labels;
}

std::vector<int> connected_components_by_label(
    const std::vector<std::vector<int>>& graph,
    const std::vector<int>& labels,
    std::vector<int>* sizes_out = nullptr) {
  const int n = static_cast<int>(labels.size());
  std::vector<int> component(static_cast<std::size_t>(n), -1);
  std::vector<int> sizes;
  std::vector<int> stack;
  for (int start = 0; start < n; ++start) {
    if (component[static_cast<std::size_t>(start)] >= 0) {
      continue;
    }
    const int comp = static_cast<int>(sizes.size());
    int size = 0;
    stack.clear();
    stack.push_back(start);
    component[static_cast<std::size_t>(start)] = comp;
    while (!stack.empty()) {
      const int current = stack.back();
      stack.pop_back();
      ++size;
      for (const int neighbor : graph[static_cast<std::size_t>(current)]) {
        if (labels[static_cast<std::size_t>(neighbor)] != labels[static_cast<std::size_t>(current)]) {
          continue;
        }
        if (component[static_cast<std::size_t>(neighbor)] >= 0) {
          continue;
        }
        component[static_cast<std::size_t>(neighbor)] = comp;
        stack.push_back(neighbor);
      }
    }
    sizes.push_back(size);
  }
  if (sizes_out != nullptr) {
    *sizes_out = sizes;
  }
  return component;
}

std::vector<int> remove_small_components(
    const std::vector<std::vector<int>>& graph,
    std::vector<int> labels,
    int k,
    int min_size) {
  if (min_size <= 1) {
    return labels;
  }
  std::vector<int> component_sizes;
  const auto components = connected_components_by_label(graph, labels, &component_sizes);
  std::vector<std::vector<int>> members(component_sizes.size());
  for (int i = 0; i < static_cast<int>(labels.size()); ++i) {
    members[static_cast<std::size_t>(components[static_cast<std::size_t>(i)])].push_back(i);
  }
  for (std::size_t comp = 0; comp < members.size(); ++comp) {
    if (component_sizes[comp] >= min_size) {
      continue;
    }
    std::vector<int> label_counts(static_cast<std::size_t>(k), 0);
    for (const int cell : members[comp]) {
      for (const int neighbor : graph[static_cast<std::size_t>(cell)]) {
        if (components[static_cast<std::size_t>(neighbor)] == static_cast<int>(comp)) {
          continue;
        }
        const int label = labels[static_cast<std::size_t>(neighbor)];
        if (label >= 0 && label < k) {
          label_counts[static_cast<std::size_t>(label)] += 1;
        }
      }
    }
    const int best = static_cast<int>(
        std::max_element(label_counts.begin(), label_counts.end()) - label_counts.begin());
    if (label_counts[static_cast<std::size_t>(best)] == 0) {
      continue;
    }
    for (const int cell : members[comp]) {
      labels[static_cast<std::size_t>(cell)] = best;
    }
  }
  return labels;
}

DenseMatrix compute_domain_composition(
    const std::vector<int>& labels,
    const std::vector<int>& annotation_ids,
    int k,
    int n_types,
    std::vector<int>* domain_sizes_out = nullptr) {
  DenseMatrix composition(k, n_types, 0.0);
  std::vector<int> sizes(static_cast<std::size_t>(k), 0);
  for (std::size_t i = 0; i < labels.size(); ++i) {
    const int label = labels[i];
    const int type = annotation_ids[i];
    if (label >= 0 && label < k && type >= 0 && type < n_types) {
      composition(label, type) += 1.0;
      sizes[static_cast<std::size_t>(label)] += 1;
    }
  }
  for (int label = 0; label < k; ++label) {
    const double denom = static_cast<double>(std::max(1, sizes[static_cast<std::size_t>(label)]));
    for (int type = 0; type < n_types; ++type) {
      composition(label, type) /= denom;
    }
  }
  if (domain_sizes_out != nullptr) {
    *domain_sizes_out = sizes;
  }
  return composition;
}

double spatial_coherence(const std::vector<std::vector<int>>& graph, const std::vector<int>& labels) {
  double same = 0.0;
  double total = 0.0;
  for (int i = 0; i < static_cast<int>(graph.size()); ++i) {
    for (const int neighbor : graph[static_cast<std::size_t>(i)]) {
      if (neighbor <= i) {
        continue;
      }
      total += 1.0;
      if (labels[static_cast<std::size_t>(i)] == labels[static_cast<std::size_t>(neighbor)]) {
        same += 1.0;
      }
    }
  }
  return total > 0.0 ? same / total : 1.0;
}

}  // namespace

DomainResult identify_domains(
    const CellTable& cells,
    const std::vector<int>& annotation_ids,
    int n_annotation_types,
    const DomainOptions& options) {
  validate_cells(cells);
  const int n = static_cast<int>(cells.cell_ids.size());
  if (static_cast<int>(annotation_ids.size()) != n) {
    throw std::runtime_error("annotation id vector must match cell count");
  }
  if (n_annotation_types <= 0) {
    throw std::runtime_error("domain detection requires at least one annotation type");
  }
  if (options.scales.empty()) {
    throw std::runtime_error("domain detection requires at least one spatial scale");
  }
  const int k = std::min(std::max(1, options.n_domains), n);
  const int max_k = std::min(std::max(1, max_neighbor_count(options)), std::max(1, n - 1));

  auto neighbors = build_cell_knn(cells, max_k, options.max_edge_distance, options.num_threads);
  DenseMatrix features = build_multiscale_features(neighbors, annotation_ids, n_annotation_types, options);
  transform_composition_blocks(
      features,
      n_annotation_types,
      static_cast<int>(options.scales.size()),
      options.transform);

  auto kmeans = run_kmeans(features, k, options.kmeans_iterations, options.seed, options.num_threads);
  const auto smooth_graph = undirected_prefix_graph(neighbors, std::min(options.smooth_k, max_k));
  int smooth_iterations = 0;
  auto labels = potts_smooth(
      features,
      std::move(kmeans.labels),
      k,
      smooth_graph,
      options.smooth_lambda,
      options.smooth_iterations,
      &smooth_iterations);
  labels = remove_small_components(smooth_graph, std::move(labels), k, options.min_region_size);

  DomainResult out;
  out.labels = std::move(labels);
  out.domain_composition =
      compute_domain_composition(out.labels, annotation_ids, k, n_annotation_types, &out.domain_sizes);
  out.region_ids = connected_components_by_label(smooth_graph, out.labels, &out.region_sizes);
  out.spatial_coherence = spatial_coherence(smooth_graph, out.labels);
  out.boundary_fraction = 1.0 - out.spatial_coherence;
  out.kmeans_iterations = kmeans.iterations;
  out.smooth_iterations = smooth_iterations;
  return out;
}

}  // namespace celladmix
