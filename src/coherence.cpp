#include "celladmix/coherence.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <iterator>
#include <optional>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "celladmix/spatial.hpp"
#include "subpar/range.hpp"

namespace celladmix {

namespace {

struct WeightedKnnGraph {
  int n_nodes = 0;
  std::vector<int> indptr;
  std::vector<int> indices;
  std::vector<double> weights;
};

struct CoherencePreparedData {
  std::vector<std::string> type_names;
  std::vector<int> cell_type_index;
  std::vector<std::vector<int>> members_by_cell;
  std::vector<int> factor_counts_by_type;
  std::vector<int> molecules_by_type;
  std::vector<int> molecules_by_factor;
  std::vector<double> source_log_enrichment;
  std::vector<double> source_probability;
  int n_factors = 0;
};

struct LocalCoherenceMoments {
  int count = 0;
  int active_count = 0;
  double mean_raw_score = 0.0;
  double mean_coherence_score = 0.0;
  double mean_score = 0.0;
  double q75_score = 0.0;
  double mean_edge_weight = 0.0;
  int largest_patch_count = 0;
  double largest_patch_fraction = 0.0;
  double patch_score = 0.0;
};

int effective_threads(int requested, int num_tasks) {
  if (num_tasks <= 0) {
    return 1;
  }
  return std::max(1, subpar::sanitize_num_workers(requested, num_tasks));
}

std::vector<std::string> resolve_cell_types(
    const TranscriptTable& table,
    const CellTable& cells) {
  std::vector<std::string> out(table.num_cells(), "");
  if (cells.cell_ids.size() == table.num_cells() &&
      cells.cell_types.size() == table.num_cells()) {
    out = cells.cell_types;
  } else if (!cells.cell_ids.empty() && cells.cell_types.size() == cells.cell_ids.size()) {
    std::unordered_map<std::string, std::string> by_id;
    by_id.reserve(cells.cell_ids.size());
    for (std::size_t i = 0; i < cells.cell_ids.size(); ++i) {
      by_id.emplace(cells.cell_ids[i], cells.cell_types[i]);
    }
    for (std::size_t cell = 0; cell < table.cells.size(); ++cell) {
      const auto it = by_id.find(table.cells[cell]);
      if (it != by_id.end()) {
        out[cell] = it->second;
      }
    }
  }

  if (std::any_of(out.begin(), out.end(), [](const std::string& value) { return value.empty(); }) &&
      table.cell_types.size() == table.size()) {
    for (std::size_t i = 0; i < table.size(); ++i) {
      auto& value = out[static_cast<std::size_t>(table.cell_index[i])];
      if (value.empty()) {
        value = table.cell_types[i];
      }
    }
  }

  return out;
}

double quantile_sorted(std::vector<double> values, double q) {
  if (values.empty()) {
    return 0.0;
  }
  std::sort(values.begin(), values.end());
  const double pos = std::clamp(q, 0.0, 1.0) * static_cast<double>(values.size() - 1);
  const auto lo = static_cast<std::size_t>(std::floor(pos));
  const auto hi = static_cast<std::size_t>(std::ceil(pos));
  if (lo == hi) {
    return values[lo];
  }
  const double frac = pos - static_cast<double>(lo);
  return values[lo] * (1.0 - frac) + values[hi] * frac;
}

double one_sided_normal_p(double z) {
  if (!std::isfinite(z)) {
    return 1.0;
  }
  if (z <= 0.0) {
    return 1.0;
  }
  return 0.5 * std::erfc(z / std::sqrt(2.0));
}

double paired_wilcoxon_greater(
    const std::vector<double>& real,
    const std::vector<double>& null_scores) {
  if (real.size() != null_scores.size()) {
    throw std::runtime_error("paired Wilcoxon inputs must have equal length");
  }
  std::vector<std::pair<double, double>> abs_and_sign;
  abs_and_sign.reserve(real.size());
  for (std::size_t i = 0; i < real.size(); ++i) {
    const double diff = real[i] - null_scores[i];
    if (std::abs(diff) <= 1e-12) {
      continue;
    }
    abs_and_sign.emplace_back(std::abs(diff), diff > 0.0 ? 1.0 : -1.0);
  }
  if (abs_and_sign.empty()) {
    return 1.0;
  }
  std::sort(abs_and_sign.begin(), abs_and_sign.end(), [](const auto& left, const auto& right) {
    return left.first < right.first;
  });

  double w_plus = 0.0;
  std::size_t i = 0;
  while (i < abs_and_sign.size()) {
    std::size_t j = i + 1U;
    while (j < abs_and_sign.size() && abs_and_sign[j].first == abs_and_sign[i].first) {
      ++j;
    }
    const double rank = (static_cast<double>(i + 1U) + static_cast<double>(j)) / 2.0;
    for (std::size_t k = i; k < j; ++k) {
      if (abs_and_sign[k].second > 0.0) {
        w_plus += rank;
      }
    }
    i = j;
  }

  const double n = static_cast<double>(abs_and_sign.size());
  const double mean = n * (n + 1.0) / 4.0;
  const double variance = n * (n + 1.0) * (2.0 * n + 1.0) / 24.0;
  if (variance <= 0.0) {
    return 1.0;
  }
  const double z = (w_plus - mean - 0.5) / std::sqrt(variance);
  const double p = 0.5 * std::erfc(z / std::sqrt(2.0));
  return std::min(1.0, std::max(0.0, p));
}

void emit_coherence_progress(
    const CoherenceTestOptions& options,
    const std::string& message,
    const std::chrono::steady_clock::time_point& start) {
  if (!options.progress) {
    return;
  }
  options.progress(
      message,
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
}

CoherencePreparedData prepare_coherence_data(
    const TranscriptTable& table,
    const CellTable& cells,
    const std::vector<int>& labels,
    const CoherenceTestOptions& options) {
  if (table.cell_index.empty()) {
    throw std::runtime_error("TranscriptTable must provide cell indices before coherence scoring");
  }
  if (labels.size() != table.size()) {
    throw std::runtime_error("labels must match TranscriptTable size");
  }

  CoherencePreparedData data;
  const auto cell_types = resolve_cell_types(table, cells);
  data.cell_type_index.assign(table.num_cells(), -1);
  std::unordered_map<std::string, int> type_lookup;
  for (std::size_t cell = 0; cell < cell_types.size(); ++cell) {
    const auto& type = cell_types[cell];
    if (type.empty()) {
      continue;
    }
    auto it = type_lookup.find(type);
    if (it == type_lookup.end()) {
      const int idx = static_cast<int>(data.type_names.size());
      it = type_lookup.emplace(type, idx).first;
      data.type_names.push_back(type);
    }
    data.cell_type_index[cell] = it->second;
  }

  data.n_factors = labels.empty() ? 0 : (*std::max_element(labels.begin(), labels.end()) + 1);
  if (data.n_factors <= 0) {
    throw std::runtime_error("Coherence scoring requires non-negative factor labels");
  }

  const int n_types = static_cast<int>(data.type_names.size());
  if (n_types == 0) {
    throw std::runtime_error("Coherence scoring requires cell_type metadata for at least one scored cell");
  }
  data.members_by_cell = table.transcripts_by_cell();
  data.factor_counts_by_type.assign(static_cast<std::size_t>(n_types * data.n_factors), 0);
  data.molecules_by_type.assign(static_cast<std::size_t>(n_types), 0);
  data.molecules_by_factor.assign(static_cast<std::size_t>(data.n_factors), 0);

  for (std::size_t i = 0; i < table.size(); ++i) {
    const int factor = labels[i];
    if (factor < 0 || factor >= data.n_factors) {
      continue;
    }
    const int cell = table.cell_index[i];
    const int type = data.cell_type_index[static_cast<std::size_t>(cell)];
    if (type < 0) {
      continue;
    }
    data.factor_counts_by_type[static_cast<std::size_t>(type * data.n_factors + factor)] += 1;
    data.molecules_by_type[static_cast<std::size_t>(type)] += 1;
    data.molecules_by_factor[static_cast<std::size_t>(factor)] += 1;
  }

  data.source_log_enrichment.assign(static_cast<std::size_t>(n_types * data.n_factors), 0.0);
  data.source_probability.assign(static_cast<std::size_t>(n_types * data.n_factors), 0.0);
  const double total_molecules = static_cast<double>(
      std::accumulate(data.molecules_by_type.begin(), data.molecules_by_type.end(), 0));
  const double eps = std::max(options.source_pseudocount, 1e-9);
  for (int type = 0; type < n_types; ++type) {
    const double type_total = static_cast<double>(data.molecules_by_type[static_cast<std::size_t>(type)]);
    const double other_total = std::max(0.0, total_molecules - type_total);
    for (int factor = 0; factor < data.n_factors; ++factor) {
      const int count = data.factor_counts_by_type[static_cast<std::size_t>(type * data.n_factors + factor)];
      const int factor_total = data.molecules_by_factor[static_cast<std::size_t>(factor)];
      const double type_fraction = (static_cast<double>(count) + eps) /
          (type_total + eps * static_cast<double>(data.n_factors));
      const double other_fraction =
          (static_cast<double>(factor_total - count) + eps) /
          (other_total + eps * static_cast<double>(data.n_factors));
      data.source_log_enrichment[static_cast<std::size_t>(type * data.n_factors + factor)] =
          std::log(type_fraction / std::max(other_fraction, 1e-300));
      data.source_probability[static_cast<std::size_t>(type * data.n_factors + factor)] =
          (static_cast<double>(count) + eps) /
          (static_cast<double>(factor_total) + eps * static_cast<double>(n_types));
    }
  }

  return data;
}

double segment_barrier(
    const Image2D& image,
    double x0,
    double y0,
    double x1,
    double y1,
    int line_samples) {
  const int n = std::max(2, line_samples);
  double best = 0.0;
  for (int i = 0; i < n; ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(n - 1);
    best = std::max(best, image.sample_bilinear(x0 + t * (x1 - x0), y0 + t * (y1 - y0)));
  }
  return best;
}

WeightedKnnGraph build_weighted_cell_knn_graph(
    const TranscriptTable& table,
    const std::vector<int>& members,
    const Image2D* membrane_image,
    const CoherenceTestOptions& options) {
  const int n = static_cast<int>(members.size());
  WeightedKnnGraph graph;
  graph.n_nodes = n;
  graph.indptr.reserve(static_cast<std::size_t>(n + 1));
  graph.indptr.push_back(0);
  if (n <= 1 || options.k_neighbors <= 0) {
    if (n == 1) {
      graph.indptr.push_back(0);
    }
    return graph;
  }

  std::unordered_map<int, int> local_index;
  local_index.reserve(members.size());
  for (int local = 0; local < n; ++local) {
    local_index.emplace(members[static_cast<std::size_t>(local)], local);
  }

  const double max_distance = options.max_neighbor_distance;
  const double max_distance2 = max_distance > 0.0 ? max_distance * max_distance : -1.0;
  const double sigma = options.distance_sigma > 0.0 ? options.distance_sigma : -1.0;
  const double denom = sigma > 0.0 ? 2.0 * sigma * sigma : 1.0;

  std::vector<std::vector<std::pair<int, double>>> neighbors(static_cast<std::size_t>(n));
  const SpatialKnnIndex knn_index(table, members);
  for (int local_i = 0; local_i < n; ++local_i) {
    const int global_i = members[static_cast<std::size_t>(local_i)];
    const auto hits = knn_index.query(global_i, options.k_neighbors, false);
    for (const auto& hit : hits) {
      if (max_distance2 > 0.0 && hit.squared_distance > max_distance2) {
        continue;
      }
      const int local_j = local_index.at(hit.index);
      double weight = sigma > 0.0 ? std::exp(-hit.squared_distance / denom) : 1.0;
      if (membrane_image != nullptr && options.use_membrane_barrier && options.membrane_alpha > 0.0) {
        const double barrier = segment_barrier(
            *membrane_image,
            table.x[static_cast<std::size_t>(global_i)],
            table.y[static_cast<std::size_t>(global_i)],
            table.x[static_cast<std::size_t>(hit.index)],
            table.y[static_cast<std::size_t>(hit.index)],
            options.line_samples);
        weight *= std::exp(-options.membrane_alpha * barrier);
      }
      neighbors[static_cast<std::size_t>(local_i)].push_back({local_j, weight});
      neighbors[static_cast<std::size_t>(local_j)].push_back({local_i, weight});
    }
  }

  for (int local_i = 0; local_i < n; ++local_i) {
    auto& node_neighbors = neighbors[static_cast<std::size_t>(local_i)];
    std::sort(node_neighbors.begin(), node_neighbors.end(),
              [](const auto& left, const auto& right) {
                if (left.first == right.first) return left.second > right.second;
                return left.first < right.first;
              });
    int last = -1;
    double best_weight = 0.0;
    for (const auto& edge : node_neighbors) {
      if (edge.first != last && last >= 0) {
        graph.indices.push_back(last);
        graph.weights.push_back(best_weight);
      }
      if (edge.first != last) {
        last = edge.first;
        best_weight = edge.second;
      } else {
        best_weight = std::max(best_weight, edge.second);
      }
    }
    if (last >= 0) {
      graph.indices.push_back(last);
      graph.weights.push_back(best_weight);
    }
    graph.indptr.push_back(static_cast<int>(graph.indices.size()));
  }

  return graph;
}

LocalCoherenceMoments score_local_subset(
    const TranscriptTable& table,
    const std::vector<int>& members,
    const WeightedKnnGraph& graph,
    const std::vector<int>& subset,
    const std::vector<unsigned char>& selected,
    double source_log_enrichment,
    const std::vector<double>* local_margin_override,
    const std::vector<double>& factor_margin,
    const CoherenceTestOptions& options) {
  LocalCoherenceMoments moments;
  moments.count = static_cast<int>(subset.size());
  if (subset.empty()) {
    return moments;
  }

  std::vector<double> final_scores;
  final_scores.reserve(subset.size());
  std::vector<double> final_score_by_local(static_cast<std::size_t>(graph.n_nodes),
                                           std::numeric_limits<double>::quiet_NaN());
  double raw_total = 0.0;
  double coherence_total = 0.0;
  double score_total = 0.0;
  double edge_weight_total = 0.0;
  int active = 0;
  auto margin_for_local = [&](int local, int global) {
    if (local_margin_override != nullptr &&
        static_cast<int>(local_margin_override->size()) == graph.n_nodes) {
      const double value = (*local_margin_override)[static_cast<std::size_t>(local)];
      if (std::isfinite(value)) {
        return value;
      }
    }
    return factor_margin.empty() ? 0.0 : factor_margin[static_cast<std::size_t>(global)];
  };
  for (const int local_i : subset) {
    const int global_i = members[static_cast<std::size_t>(local_i)];
    const double margin_i = margin_for_local(local_i, global_i);
    const double raw_i = source_log_enrichment + options.beta_margin * margin_i;

    double denom = 0.0;
    double support = 0.0;
    for (int edge = graph.indptr[static_cast<std::size_t>(local_i)];
         edge < graph.indptr[static_cast<std::size_t>(local_i + 1)];
         ++edge) {
      const int local_j = graph.indices[static_cast<std::size_t>(edge)];
      const double weight = graph.weights[static_cast<std::size_t>(edge)];
      denom += weight;
      if (selected[static_cast<std::size_t>(local_j)] != 0) {
        const int global_j = members[static_cast<std::size_t>(local_j)];
        const double margin_j = margin_for_local(local_j, global_j);
        const double raw_j = source_log_enrichment + options.beta_margin * margin_j;
        support += weight * std::max(raw_j, 0.0);
      }
    }
    const double coherence_i = denom > 0.0 ? support / denom : 0.0;
    const double final_i = raw_i + options.lambda_coherence * coherence_i;
    raw_total += raw_i;
    coherence_total += coherence_i;
    score_total += final_i;
    edge_weight_total += denom;
    final_scores.push_back(final_i);
    final_score_by_local[static_cast<std::size_t>(local_i)] = final_i;
    if (final_i > options.score_threshold) {
      ++active;
    }
  }

  const double patch_edge_weight_min = std::max(0.0, options.patch_edge_weight_min);
  std::vector<unsigned char> visited(static_cast<std::size_t>(graph.n_nodes), 0);
  std::vector<int> stack;
  stack.reserve(subset.size());

  for (const int seed : subset) {
    if (seed < 0 || seed >= graph.n_nodes ||
        selected[static_cast<std::size_t>(seed)] == 0) {
      continue;
    }
    int local_patch_count = 1;
    double local_patch_score = final_score_by_local[static_cast<std::size_t>(seed)];
    for (int edge = graph.indptr[static_cast<std::size_t>(seed)];
         edge < graph.indptr[static_cast<std::size_t>(seed + 1)];
         ++edge) {
      if (graph.weights[static_cast<std::size_t>(edge)] < patch_edge_weight_min) {
        continue;
      }
      const int local_j = graph.indices[static_cast<std::size_t>(edge)];
      if (selected[static_cast<std::size_t>(local_j)] == 0) {
        continue;
      }
      const double score_j = final_score_by_local[static_cast<std::size_t>(local_j)];
      if (std::isfinite(score_j)) {
        local_patch_score += score_j;
        ++local_patch_count;
      }
    }
    const double local_mean =
        local_patch_score / static_cast<double>(std::max(1, local_patch_count));
    const double size_weight = std::sqrt(
        static_cast<double>(local_patch_count) /
        static_cast<double>(std::max<int>(1, moments.count)));
    moments.patch_score = std::max(moments.patch_score, local_mean * size_weight);
  }

  for (const int seed : subset) {
    if (seed < 0 || seed >= graph.n_nodes ||
        visited[static_cast<std::size_t>(seed)] != 0 ||
        selected[static_cast<std::size_t>(seed)] == 0) {
      continue;
    }
    visited[static_cast<std::size_t>(seed)] = 1;
    stack.clear();
    stack.push_back(seed);
    int component_count = 0;
    double component_score = 0.0;
    while (!stack.empty()) {
      const int local_i = stack.back();
      stack.pop_back();
      ++component_count;
      const double score_i = final_score_by_local[static_cast<std::size_t>(local_i)];
      if (std::isfinite(score_i)) {
        component_score += score_i;
      }
      for (int edge = graph.indptr[static_cast<std::size_t>(local_i)];
           edge < graph.indptr[static_cast<std::size_t>(local_i + 1)];
           ++edge) {
        if (graph.weights[static_cast<std::size_t>(edge)] < patch_edge_weight_min) {
          continue;
        }
        const int local_j = graph.indices[static_cast<std::size_t>(edge)];
        if (selected[static_cast<std::size_t>(local_j)] == 0 ||
            visited[static_cast<std::size_t>(local_j)] != 0) {
          continue;
        }
        visited[static_cast<std::size_t>(local_j)] = 1;
        stack.push_back(local_j);
      }
    }
    (void)component_score;
    if (component_count > moments.largest_patch_count) {
      moments.largest_patch_count = component_count;
    }
  }

  const double denom = static_cast<double>(std::max<int>(1, moments.count));
  moments.active_count = active;
  moments.mean_raw_score = raw_total / denom;
  moments.mean_coherence_score = coherence_total / denom;
  moments.mean_score = score_total / denom;
  moments.q75_score = quantile_sorted(final_scores, 0.75);
  moments.mean_edge_weight = edge_weight_total / denom;
  moments.largest_patch_fraction =
      static_cast<double>(moments.largest_patch_count) / denom;
  return moments;
}

std::vector<int> quantile_bins(const std::vector<double>& values, int requested_bins) {
  const int bins = std::max(1, requested_bins);
  std::vector<int> out(values.size(), 0);
  if (bins <= 1 || values.empty()) {
    return out;
  }

  std::vector<double> finite;
  finite.reserve(values.size());
  for (const double value : values) {
    if (std::isfinite(value)) {
      finite.push_back(value);
    }
  }
  if (finite.size() < 2) {
    return out;
  }
  std::sort(finite.begin(), finite.end());

  std::vector<double> cuts;
  cuts.reserve(static_cast<std::size_t>(bins - 1));
  for (int bin = 1; bin < bins; ++bin) {
    const double pos =
        static_cast<double>(bin) * static_cast<double>(finite.size() - 1U) /
        static_cast<double>(bins);
    const auto lo = static_cast<std::size_t>(std::floor(pos));
    const auto hi = static_cast<std::size_t>(std::ceil(pos));
    const double frac = pos - static_cast<double>(lo);
    cuts.push_back(finite[lo] * (1.0 - frac) + finite[hi] * frac);
  }

  for (std::size_t i = 0; i < values.size(); ++i) {
    if (!std::isfinite(values[i])) {
      out[i] = 0;
      continue;
    }
    out[i] = static_cast<int>(std::upper_bound(cuts.begin(), cuts.end(), values[i]) - cuts.begin());
  }
  return out;
}

std::vector<int> build_null_strata(
    const TranscriptTable& table,
    const std::vector<int>& members,
    const WeightedKnnGraph& graph,
    const CoherenceTestOptions& options) {
  const int n_local = static_cast<int>(members.size());
  const bool use_nucleus_overlap =
      options.null_match_nucleus && table.overlaps_nucleus.size() == table.size();
  const bool use_nucleus_distance =
      options.null_match_nucleus && table.nucleus_distance.size() == table.size() &&
      options.null_nucleus_distance_bins > 1;
  const bool use_density =
      options.null_match_density && options.null_density_bins > 1;
  if (!use_nucleus_overlap && !use_nucleus_distance && !use_density) {
    return {};
  }

  std::vector<int> strata(static_cast<std::size_t>(n_local), 0);
  int multiplier = 1;
  if (use_nucleus_overlap) {
    for (int local = 0; local < n_local; ++local) {
      const int global = members[static_cast<std::size_t>(local)];
      const int value = table.overlaps_nucleus[static_cast<std::size_t>(global)];
      const int category = value < 0 ? 2 : (value > 0 ? 1 : 0);
      strata[static_cast<std::size_t>(local)] += category * multiplier;
    }
    multiplier *= 3;
  }

  if (use_nucleus_distance) {
    std::vector<double> distances(static_cast<std::size_t>(n_local), 0.0);
    for (int local = 0; local < n_local; ++local) {
      const int global = members[static_cast<std::size_t>(local)];
      distances[static_cast<std::size_t>(local)] =
          table.nucleus_distance[static_cast<std::size_t>(global)];
    }
    const auto bins = quantile_bins(distances, options.null_nucleus_distance_bins);
    for (int local = 0; local < n_local; ++local) {
      strata[static_cast<std::size_t>(local)] += bins[static_cast<std::size_t>(local)] * multiplier;
    }
    multiplier *= std::max(1, options.null_nucleus_distance_bins);
  }

  if (use_density) {
    std::vector<double> density(static_cast<std::size_t>(n_local), 0.0);
    for (int local = 0; local < n_local; ++local) {
      double total_weight = 0.0;
      for (int edge = graph.indptr[static_cast<std::size_t>(local)];
           edge < graph.indptr[static_cast<std::size_t>(local + 1)];
           ++edge) {
        total_weight += graph.weights[static_cast<std::size_t>(edge)];
      }
      density[static_cast<std::size_t>(local)] = total_weight;
    }
    const auto bins = quantile_bins(density, options.null_density_bins);
    for (int local = 0; local < n_local; ++local) {
      strata[static_cast<std::size_t>(local)] += bins[static_cast<std::size_t>(local)] * multiplier;
    }
  }

  return strata;
}

std::vector<int> sample_null_subset(
    const std::vector<int>& local_by_factor,
    int n_local,
    const std::vector<int>* null_strata,
    const std::vector<unsigned char>& observed_selected,
    bool exclude_observed,
    std::mt19937& rng) {
  const int target_count = static_cast<int>(local_by_factor.size());
  if (null_strata != nullptr &&
      static_cast<int>(null_strata->size()) == n_local &&
      target_count > 0) {
    std::vector<int> stratum_order;
    std::unordered_map<int, int> needed_by_stratum;
    needed_by_stratum.reserve(local_by_factor.size());
    for (const int local : local_by_factor) {
      const int stratum = (*null_strata)[static_cast<std::size_t>(local)];
      auto inserted = needed_by_stratum.emplace(stratum, 0);
      if (inserted.second) {
        stratum_order.push_back(stratum);
      }
      inserted.first->second += 1;
    }

    std::vector<unsigned char> already_selected(static_cast<std::size_t>(n_local), 0);
    std::vector<int> out;
    out.reserve(local_by_factor.size());
    auto draw_from_pool = [&](const std::vector<int>& pool, int need) {
      if (need <= 0 || pool.empty()) {
        return 0;
      }
      std::vector<int> candidates;
      candidates.reserve(pool.size());
      for (const int local : pool) {
        if (local >= 0 &&
            local < n_local &&
            already_selected[static_cast<std::size_t>(local)] == 0) {
          candidates.push_back(local);
        }
      }
      const int keep = std::min(need, static_cast<int>(candidates.size()));
      for (int i = 0; i < keep; ++i) {
        std::uniform_int_distribution<int> pick(i, static_cast<int>(candidates.size()) - 1);
        const int j = pick(rng);
        std::swap(candidates[static_cast<std::size_t>(i)], candidates[static_cast<std::size_t>(j)]);
        const int local = candidates[static_cast<std::size_t>(i)];
        already_selected[static_cast<std::size_t>(local)] = 1;
        out.push_back(local);
      }
      return keep;
    };

    int deficit = 0;
    for (const int stratum : stratum_order) {
      int remaining = needed_by_stratum[stratum];
      std::vector<int> pool;
      pool.reserve(static_cast<std::size_t>(n_local));
      for (int local = 0; local < n_local; ++local) {
        if ((*null_strata)[static_cast<std::size_t>(local)] != stratum) {
          continue;
        }
        if (!exclude_observed ||
            observed_selected[static_cast<std::size_t>(local)] == 0) {
          pool.push_back(local);
        }
      }
      remaining -= draw_from_pool(pool, remaining);
      if (remaining > 0 && exclude_observed) {
        pool.clear();
        for (int local = 0; local < n_local; ++local) {
          if ((*null_strata)[static_cast<std::size_t>(local)] == stratum) {
            pool.push_back(local);
          }
        }
        remaining -= draw_from_pool(pool, remaining);
      }
      if (remaining > 0) {
        deficit += remaining;
      }
    }

    if (deficit > 0) {
      std::vector<int> pool;
      pool.reserve(static_cast<std::size_t>(n_local));
      for (int local = 0; local < n_local; ++local) {
        if (!exclude_observed ||
            observed_selected[static_cast<std::size_t>(local)] == 0) {
          pool.push_back(local);
        }
      }
      deficit -= draw_from_pool(pool, deficit);
    }
    if (deficit > 0) {
      std::vector<int> pool(static_cast<std::size_t>(n_local));
      std::iota(pool.begin(), pool.end(), 0);
      (void)draw_from_pool(pool, deficit);
    }
    return out;
  }

  std::vector<int> pool;
  pool.reserve(static_cast<std::size_t>(n_local));
  if (exclude_observed) {
    for (int local = 0; local < n_local; ++local) {
      if (observed_selected[static_cast<std::size_t>(local)] == 0) {
        pool.push_back(local);
      }
    }
  }
  if (static_cast<int>(pool.size()) < target_count) {
    pool.clear();
    pool.reserve(static_cast<std::size_t>(n_local));
    for (int local = 0; local < n_local; ++local) {
      pool.push_back(local);
    }
  }

  const int keep = std::min(target_count, static_cast<int>(pool.size()));
  for (int i = 0; i < keep; ++i) {
    std::uniform_int_distribution<int> pick(i, static_cast<int>(pool.size()) - 1);
    const int j = pick(rng);
    std::swap(pool[static_cast<std::size_t>(i)], pool[static_cast<std::size_t>(j)]);
  }
  pool.resize(static_cast<std::size_t>(keep));
  return pool;
}

std::vector<double> permuted_factor_margins_by_local(
    const std::vector<int>& members,
    const std::vector<int>& observed_subset,
    const std::vector<int>& null_subset,
    const std::vector<int>* null_strata,
    const std::vector<double>& factor_margin,
    std::mt19937& rng) {
  const int n_local = static_cast<int>(members.size());
  if (factor_margin.empty() || observed_subset.empty() || null_subset.empty()) {
    return {};
  }

  const double missing = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> out(static_cast<std::size_t>(n_local), missing);
  auto stratum_for_local = [&](int local) {
    if (null_strata != nullptr &&
        static_cast<int>(null_strata->size()) == n_local) {
      return (*null_strata)[static_cast<std::size_t>(local)];
    }
    return 0;
  };

  std::unordered_map<int, std::vector<double>> margins_by_stratum;
  margins_by_stratum.reserve(observed_subset.size());
  std::vector<double> all_margins;
  all_margins.reserve(observed_subset.size());
  for (const int local : observed_subset) {
    if (local < 0 || local >= n_local) {
      continue;
    }
    const int global = members[static_cast<std::size_t>(local)];
    if (global < 0 || static_cast<std::size_t>(global) >= factor_margin.size()) {
      continue;
    }
    const double margin = factor_margin[static_cast<std::size_t>(global)];
    margins_by_stratum[stratum_for_local(local)].push_back(margin);
    all_margins.push_back(margin);
  }
  if (all_margins.empty()) {
    return {};
  }

  std::unordered_map<int, std::vector<int>> positions_by_stratum;
  positions_by_stratum.reserve(null_subset.size());
  for (const int local : null_subset) {
    if (local >= 0 && local < n_local) {
      positions_by_stratum[stratum_for_local(local)].push_back(local);
    }
  }

  for (auto& entry : positions_by_stratum) {
    auto margin_it = margins_by_stratum.find(entry.first);
    std::vector<double> margins =
        margin_it == margins_by_stratum.end() || margin_it->second.empty()
            ? all_margins
            : margin_it->second;
    std::shuffle(margins.begin(), margins.end(), rng);
    auto& positions = entry.second;
    std::shuffle(positions.begin(), positions.end(), rng);
    for (std::size_t i = 0; i < positions.size(); ++i) {
      const double margin = i < margins.size()
          ? margins[i]
          : margins[static_cast<std::size_t>(
                std::uniform_int_distribution<int>(0, static_cast<int>(margins.size()) - 1)(rng))];
      out[static_cast<std::size_t>(positions[i])] = margin;
    }
  }

  return out;
}

}  // namespace

CoherenceTestResult run_coherence_test(
    const TranscriptTable& table,
    const CellTable& cells,
    const std::vector<int>& labels,
    const std::vector<double>& factor_margin,
    const Image2D* membrane_image,
    const CoherenceTestOptions& options) {
  if (!factor_margin.empty() && factor_margin.size() != table.size()) {
    throw std::runtime_error("factor_margin must be empty or match TranscriptTable size");
  }
  if (options.k_neighbors < 1) {
    throw std::runtime_error("k_neighbors must be positive");
  }
  if (options.min_factor_molecules < 1) {
    throw std::runtime_error("min_factor_molecules must be positive");
  }
  if (options.line_samples < 2) {
    throw std::runtime_error("line_samples must be at least 2");
  }
  if (options.null_iterations < 0) {
    throw std::runtime_error("null_iterations must be non-negative");
  }
  if (options.null_nucleus_distance_bins < 1 || options.null_density_bins < 1) {
    throw std::runtime_error("coherence null bin counts must be positive");
  }
  if (options.null_method != "label_permutation" &&
      options.null_method != "matched_subset") {
    throw std::runtime_error("coherence null_method must be 'label_permutation' or 'matched_subset'");
  }
  if (options.normalize_membrane &&
      !(options.membrane_low_quantile >= 0.0 && options.membrane_low_quantile <= 1.0 &&
        options.membrane_high_quantile >= 0.0 && options.membrane_high_quantile <= 1.0 &&
        options.membrane_high_quantile > options.membrane_low_quantile)) {
    throw std::runtime_error("membrane quantiles must satisfy 0 <= low < high <= 1");
  }

  const auto start = std::chrono::steady_clock::now();
  auto stage_start = start;
  emit_coherence_progress(options, "Starting coherence test: molecules=" + std::to_string(table.size()), start);

  std::optional<Image2D> normalized_membrane;
  const Image2D* scoring_membrane = membrane_image;
  CoherenceTestResult result;
  if (membrane_image != nullptr && options.use_membrane_barrier && options.normalize_membrane) {
    normalized_membrane = *membrane_image;
    const auto scale = normalize_image_quantile(
        *normalized_membrane,
        options.membrane_low_quantile,
        options.membrane_high_quantile);
    scoring_membrane = &(*normalized_membrane);
    result.membrane_normalized = true;
    result.membrane_scale_low = scale.first;
    result.membrane_scale_high = scale.second;
    emit_coherence_progress(
        options,
        "Normalized membrane crop: low=" + std::to_string(scale.first) +
            ", high=" + std::to_string(scale.second),
        stage_start);
  }

  const auto data = prepare_coherence_data(table, cells, labels, options);
  emit_coherence_progress(
      options,
      "Prepared coherence data: cell_types=" + std::to_string(data.type_names.size()) +
          ", factors=" + std::to_string(data.n_factors),
      stage_start);
  result.null_method = options.compute_null && options.null_iterations > 0
      ? options.null_method
      : "none";
  result.null_match_nucleus_overlap_used =
      options.compute_null && options.null_iterations > 0 &&
      options.null_match_nucleus && table.overlaps_nucleus.size() == table.size();
  result.null_match_nucleus_distance_used =
      options.compute_null && options.null_iterations > 0 &&
      options.null_match_nucleus && table.nucleus_distance.size() == table.size() &&
      options.null_nucleus_distance_bins > 1;
  result.null_match_density_used =
      options.compute_null && options.null_iterations > 0 &&
      options.null_match_density && options.null_density_bins > 1;
  result.null_nucleus_distance_bins = result.null_match_nucleus_distance_used
      ? options.null_nucleus_distance_bins
      : 0;
  result.null_density_bins = result.null_match_density_used
      ? options.null_density_bins
      : 0;
  if (options.compute_null && options.null_iterations > 0) {
    emit_coherence_progress(
        options,
        "Coherence null: method=" + result.null_method +
            ", nucleus_overlap=" +
            std::string(result.null_match_nucleus_overlap_used ? "yes" : "no") +
            ", nucleus_distance=" +
            std::string(result.null_match_nucleus_distance_used ? "yes" : "no") +
            ", density=" +
            std::string(result.null_match_density_used ? "yes" : "no"),
        stage_start);
  }

  result.cell_types = data.type_names;
  result.source_log_enrichment = data.source_log_enrichment;
  result.source_probability = data.source_probability;
  result.n_factors = data.n_factors;

  const int n_cells = static_cast<int>(table.num_cells());
  const int n_types = static_cast<int>(data.type_names.size());
  const int workers = effective_threads(options.num_threads, n_cells);
  std::vector<std::vector<CoherenceCellScore>> local_scores(static_cast<std::size_t>(workers));

  stage_start = std::chrono::steady_clock::now();
  subpar::parallelize_range<true>(workers, n_cells, [&](int worker, int start_cell, int length) {
    auto& out = local_scores[static_cast<std::size_t>(worker)];
    for (int cell = start_cell; cell < start_cell + length; ++cell) {
      const auto& members = data.members_by_cell[static_cast<std::size_t>(cell)];
      if (members.empty()) {
        continue;
      }
      const int target_type = data.cell_type_index[static_cast<std::size_t>(cell)];
      if (target_type < 0) {
        continue;
      }

      const auto graph = build_weighted_cell_knn_graph(table, members, scoring_membrane, options);
      const int n_local = static_cast<int>(members.size());
      const auto null_strata = (options.compute_null && options.null_iterations > 0)
          ? build_null_strata(table, members, graph, options)
          : std::vector<int>{};
      const std::vector<int>* null_strata_ptr = null_strata.empty() ? nullptr : &null_strata;
      std::vector<std::vector<int>> local_by_factor(static_cast<std::size_t>(data.n_factors));
      for (int local = 0; local < n_local; ++local) {
        const int global = members[static_cast<std::size_t>(local)];
        const int factor = labels[static_cast<std::size_t>(global)];
        if (factor >= 0 && factor < data.n_factors) {
          local_by_factor[static_cast<std::size_t>(factor)].push_back(local);
        }
      }

      std::vector<unsigned char> observed_selected(static_cast<std::size_t>(n_local), 0);
      std::vector<unsigned char> null_selected(static_cast<std::size_t>(n_local), 0);
      for (int factor = 0; factor < data.n_factors; ++factor) {
        const auto& factor_members = local_by_factor[static_cast<std::size_t>(factor)];
        if (static_cast<int>(factor_members.size()) < options.min_factor_molecules) {
          continue;
        }
        std::fill(observed_selected.begin(), observed_selected.end(), static_cast<unsigned char>(0));
        for (const int local : factor_members) {
          observed_selected[static_cast<std::size_t>(local)] = 1;
        }

        for (int source_type = 0; source_type < n_types; ++source_type) {
          if (!options.include_self_source && source_type == target_type) {
            continue;
          }
          const std::size_t sf = static_cast<std::size_t>(source_type * data.n_factors + factor);
          const double source_log_enrichment = data.source_log_enrichment[sf];
          if (source_log_enrichment < options.min_source_log_enrichment) {
            continue;
          }
          const double source_probability = data.source_probability[sf];

          const auto observed = score_local_subset(
              table,
              members,
              graph,
              factor_members,
              observed_selected,
              source_log_enrichment,
              nullptr,
              factor_margin,
              options);

          double mean_null_score = 0.0;
          double mean_null_patch_score = 0.0;
          if (options.compute_null && options.null_iterations > 0) {
            const unsigned int row_seed = options.seed ^
                static_cast<unsigned int>((cell + 1) * 73856093U) ^
                static_cast<unsigned int>((source_type + 1) * 19349663U) ^
                static_cast<unsigned int>((factor + 1) * 83492791U);
            std::mt19937 rng(row_seed);
            const int iterations = std::max(options.null_iterations, 1);
            const bool use_label_permutation = options.null_method == "label_permutation";
            for (int iter = 0; iter < iterations; ++iter) {
              std::fill(null_selected.begin(), null_selected.end(), static_cast<unsigned char>(0));
              const auto null_subset = sample_null_subset(
                  factor_members,
                  n_local,
                  null_strata_ptr,
                  observed_selected,
                  use_label_permutation ? false : options.null_exclude_factor,
                  rng);
              for (const int local : null_subset) {
                null_selected[static_cast<std::size_t>(local)] = 1;
              }
              const auto null_margin_by_local = use_label_permutation
                  ? permuted_factor_margins_by_local(
                        members,
                        factor_members,
                        null_subset,
                        null_strata_ptr,
                        factor_margin,
                        rng)
                  : std::vector<double>{};
              const std::vector<double>* null_margin_ptr =
                  null_margin_by_local.empty() ? nullptr : &null_margin_by_local;
              const auto null_moments = score_local_subset(
                  table,
                  members,
                  graph,
                  null_subset,
                  null_selected,
                  source_log_enrichment,
                  null_margin_ptr,
                  factor_margin,
                  options);
              mean_null_score += null_moments.mean_score / static_cast<double>(iterations);
              mean_null_patch_score += null_moments.patch_score / static_cast<double>(iterations);
            }
          }

          CoherenceCellScore row;
          row.target_cell = cell;
          row.target_type = target_type;
          row.source_type = source_type;
          row.factor = factor;
          row.factor_count = observed.count;
          row.active_count = observed.active_count;
          const double denom = static_cast<double>(std::max<int>(1, row.factor_count));
          row.mean_raw_score = observed.mean_raw_score;
          row.mean_coherence_score = observed.mean_coherence_score;
          row.mean_score = observed.mean_score;
          row.mean_null_score = mean_null_score;
          row.score_delta = row.mean_score - row.mean_null_score;
          row.q75_score = observed.q75_score;
          row.active_fraction = static_cast<double>(observed.active_count) / denom;
          row.mean_edge_weight = observed.mean_edge_weight;
          row.largest_patch_count = observed.largest_patch_count;
          row.largest_patch_fraction = observed.largest_patch_fraction;
          row.patch_score = observed.patch_score;
          row.mean_null_patch_score = mean_null_patch_score;
          row.patch_score_delta = row.patch_score - row.mean_null_patch_score;
          row.source_log_enrichment = source_log_enrichment;
          row.source_probability = source_probability;
          out.push_back(row);
        }
      }
    }
  });
  emit_coherence_progress(options, "Scored coherence cells", stage_start);

  std::size_t total_scores = 0;
  for (const auto& part : local_scores) {
    total_scores += part.size();
  }
  result.cell_scores.reserve(total_scores);
  for (auto& part : local_scores) {
    result.cell_scores.insert(result.cell_scores.end(),
                              std::make_move_iterator(part.begin()),
                              std::make_move_iterator(part.end()));
  }

  stage_start = std::chrono::steady_clock::now();
  std::sort(result.cell_scores.begin(), result.cell_scores.end(),
            [](const CoherenceCellScore& left, const CoherenceCellScore& right) {
              if (left.target_type != right.target_type) return left.target_type < right.target_type;
              if (left.source_type != right.source_type) return left.source_type < right.source_type;
              if (left.factor != right.factor) return left.factor < right.factor;
              return left.target_cell < right.target_cell;
            });

  struct Accumulator {
    std::vector<int> rows;
  };
  std::unordered_map<std::size_t, Accumulator> groups;
  groups.reserve(result.cell_scores.size());
  for (int row = 0; row < static_cast<int>(result.cell_scores.size()); ++row) {
    const auto& score = result.cell_scores[static_cast<std::size_t>(row)];
    const std::size_t key =
        (static_cast<std::size_t>(score.target_type) * static_cast<std::size_t>(n_types) +
         static_cast<std::size_t>(score.source_type)) *
            static_cast<std::size_t>(data.n_factors) +
        static_cast<std::size_t>(score.factor);
    groups[key].rows.push_back(row);
  }

  result.summaries.reserve(groups.size());
  for (auto& entry : groups) {
    const auto& rows = entry.second.rows;
    if (static_cast<int>(rows.size()) < options.min_cells) {
      continue;
    }
    const auto& first = result.cell_scores[static_cast<std::size_t>(rows.front())];
    int n_molecules = 0;
    int n_active = 0;
    double total_score = 0.0;
    double total_null_score = 0.0;
    double total_patch_score = 0.0;
    double total_null_patch_score = 0.0;
    double total_patch_fraction = 0.0;
    std::vector<double> cell_means;
    std::vector<double> null_means;
    std::vector<double> delta_means;
    std::vector<double> patch_means;
    std::vector<double> null_patch_means;
    cell_means.reserve(rows.size());
    null_means.reserve(rows.size());
    delta_means.reserve(rows.size());
    patch_means.reserve(rows.size());
    null_patch_means.reserve(rows.size());
    for (const int row_i : rows) {
      auto& score = result.cell_scores[static_cast<std::size_t>(row_i)];
      score.used_in_summary = true;
      n_molecules += score.factor_count;
      n_active += score.active_count;
      total_score += score.mean_score;
      total_null_score += score.mean_null_score;
      total_patch_score += score.patch_score;
      total_null_patch_score += score.mean_null_patch_score;
      total_patch_fraction += score.largest_patch_fraction;
      cell_means.push_back(score.mean_score);
      null_means.push_back(score.mean_null_score);
      delta_means.push_back(score.score_delta);
      patch_means.push_back(score.patch_score);
      null_patch_means.push_back(score.mean_null_patch_score);
    }
    const double n = static_cast<double>(rows.size());
    const double mean = total_score / std::max(1.0, n);
    const double null_mean = total_null_score / std::max(1.0, n);
    const double delta_mean = mean - null_mean;
    const double patch_mean = total_patch_score / std::max(1.0, n);
    const double null_patch_mean = total_null_patch_score / std::max(1.0, n);
    const double delta_patch_mean = patch_mean - null_patch_mean;
    const double patch_fraction_mean = total_patch_fraction / std::max(1.0, n);
    const double p = options.compute_null && options.null_iterations > 0
        ? paired_wilcoxon_greater(cell_means, null_means)
        : std::numeric_limits<double>::quiet_NaN();
    const double patch_p = options.compute_null && options.null_iterations > 0
        ? paired_wilcoxon_greater(patch_means, null_patch_means)
        : std::numeric_limits<double>::quiet_NaN();

    CoherenceSummary summary;
    summary.target_type = first.target_type;
    summary.source_type = first.source_type;
    summary.factor = first.factor;
    summary.n_cells = static_cast<int>(rows.size());
    summary.n_molecules = n_molecules;
    summary.n_active_molecules = n_active;
    summary.active_fraction =
        static_cast<double>(n_active) / static_cast<double>(std::max(1, n_molecules));
    summary.mean_score = mean;
    summary.mean_null_score = null_mean;
    summary.mean_delta_score = delta_mean;
    summary.q75_score = quantile_sorted(cell_means, 0.75);
    summary.mean_patch_score = patch_mean;
    summary.mean_null_patch_score = null_patch_mean;
    summary.mean_delta_patch_score = delta_patch_mean;
    summary.mean_largest_patch_fraction = patch_fraction_mean;
    summary.patch_p_value = patch_p;
    if (std::isnan(patch_p)) {
      summary.patch_neg_log10_p = std::numeric_limits<double>::quiet_NaN();
    } else if (patch_p > 0.0) {
      summary.patch_neg_log10_p = -std::log10(patch_p);
    } else {
      summary.patch_neg_log10_p = std::numeric_limits<double>::infinity();
    }
    summary.p_value = p;
    if (std::isnan(p)) {
      summary.neg_log10_p = std::numeric_limits<double>::quiet_NaN();
    } else if (p > 0.0) {
      summary.neg_log10_p = -std::log10(p);
    } else {
      summary.neg_log10_p = std::numeric_limits<double>::infinity();
    }
    summary.source_log_enrichment = first.source_log_enrichment;
    summary.source_probability = first.source_probability;
    result.summaries.push_back(summary);
  }

  std::sort(result.summaries.begin(), result.summaries.end(),
            [](const CoherenceSummary& left, const CoherenceSummary& right) {
              if (left.target_type != right.target_type) return left.target_type < right.target_type;
              if (left.source_type != right.source_type) return left.source_type < right.source_type;
              return left.factor < right.factor;
            });
  emit_coherence_progress(options, "Summarized coherence scores", stage_start);
  emit_coherence_progress(
      options,
      "Finished coherence test: cell_scores=" + std::to_string(result.cell_scores.size()) +
          ", summaries=" + std::to_string(result.summaries.size()) +
          ", compute_null=" + std::string(options.compute_null ? "TRUE" : "FALSE") +
          ", null_iterations=" + std::to_string(options.null_iterations),
      start);
  return result;
}

}  // namespace celladmix
