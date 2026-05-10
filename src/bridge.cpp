#include "celladmix/bridge.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include "celladmix/spatial.hpp"
#include "subpar/range.hpp"

namespace celladmix {

// Bridge-scoring helpers for finding candidate factor transfer between touching cells.

namespace {

int effective_threads(int requested, int num_tasks) {
  if (num_tasks <= 0) {
    return 1;
  }
  return std::max(1, subpar::sanitize_num_workers(requested, num_tasks));
}

double elapsed_seconds(const std::chrono::steady_clock::time_point& start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

std::string format_double(double value, int digits = 3) {
  std::ostringstream out;
  out << std::fixed << std::setprecision(digits) << value;
  return out.str();
}

void emit_bridge_progress(
    const BridgeTestOptions& options,
    const std::string& message,
    const std::chrono::steady_clock::time_point& start) {
  if (options.progress) {
    options.progress(message, elapsed_seconds(start));
  }
}

std::uint64_t pair_key(int first, int second) {
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(first)) << 32) |
      static_cast<std::uint32_t>(second);
}

int key_first(std::uint64_t key) {
  return static_cast<int>(key >> 32);
}

int key_second(std::uint64_t key) {
  return static_cast<int>(key & 0xffffffffU);
}

std::vector<std::string> resolve_cell_types(
    const TranscriptTable& table,
    const CellTable& cells) {
  std::vector<std::string> out(table.num_cells(), "");

  if (cells.cell_types.size() == cells.cell_ids.size() && !cells.cell_types.empty()) {
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

struct BridgePreparedData {
  std::vector<std::string> type_names;
  std::vector<int> cell_type_index;
  std::vector<std::vector<int>> members_by_cell;
  std::vector<int> factor_counts;
  std::vector<double> factor_fractions;
  std::vector<double> centroid_x;
  std::vector<double> centroid_y;
  std::vector<double> centroid_z;
  int n_factors = 0;
  int n_unannotated_cells = 0;
};

double cell_radius_xy(
    const TranscriptTable& table,
    const std::vector<int>& members,
    double cx,
    double cy);

BridgePreparedData prepare_bridge_data(
    const TranscriptTable& table,
    const CellTable& cells,
    const std::vector<int>& labels) {
  if (table.cell_index.empty()) {
    throw std::runtime_error("TranscriptTable must provide cell indices before bridge scoring");
  }
  if (labels.size() != table.size()) {
    throw std::runtime_error("labels must match TranscriptTable size");
  }

  BridgePreparedData data;
  const auto cell_types = resolve_cell_types(table, cells);
  data.cell_type_index.assign(table.num_cells(), -1);
  std::unordered_map<std::string, int> type_lookup;
  for (std::size_t cell = 0; cell < cell_types.size(); ++cell) {
    const auto& type = cell_types[cell];
    if (type.empty()) {
      ++data.n_unannotated_cells;
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
  if (data.type_names.empty()) {
    throw std::runtime_error("Bridge scoring requires at least one annotated cell type");
  }

  data.n_factors = labels.empty() ? 0 : (*std::max_element(labels.begin(), labels.end()) + 1);
  if (data.n_factors <= 0) {
    throw std::runtime_error("Bridge scoring requires non-negative factor labels");
  }

  data.members_by_cell = table.transcripts_by_cell();
  const int n_cells = static_cast<int>(table.num_cells());
  data.factor_counts.assign(static_cast<std::size_t>(n_cells * data.n_factors), 0);
  data.factor_fractions.assign(static_cast<std::size_t>(n_cells * data.n_factors), 0.0);
  for (std::size_t i = 0; i < table.size(); ++i) {
    const int factor = labels[i];
    if (factor < 0 || factor >= data.n_factors) {
      continue;
    }
    const int cell = table.cell_index[i];
    data.factor_counts[static_cast<std::size_t>(cell * data.n_factors + factor)] += 1;
  }
  for (int cell = 0; cell < n_cells; ++cell) {
    const double denom = static_cast<double>(std::max<std::size_t>(
        data.members_by_cell[static_cast<std::size_t>(cell)].size(), static_cast<std::size_t>(1)));
    for (int factor = 0; factor < data.n_factors; ++factor) {
      data.factor_fractions[static_cast<std::size_t>(cell * data.n_factors + factor)] =
          static_cast<double>(data.factor_counts[static_cast<std::size_t>(cell * data.n_factors + factor)]) /
          denom;
    }
  }
  data.centroid_x = infer_cell_centroids_x(table);
  data.centroid_y = infer_cell_centroids_y(table);
  data.centroid_z = infer_cell_centroids_z(table);
  return data;
}

std::vector<BridgeCandidatePair> discover_molecule_global_bridge_candidates(
    const TranscriptTable& table,
    const BridgePreparedData& data,
    const BridgeTestOptions& options) {
  if (options.candidate_k <= 0) {
    return {};
  }
  // TODO: for full Xenium-scale runs, replace this global candidate KD-tree
  // with exact tile-block plus halo discovery. Final two-cell scoring can
  // remain local and exact.
  const SpatialKnnIndex global_index(table);
  const int n_types = static_cast<int>(data.type_names.size());
  const int workers = effective_threads(options.num_threads, static_cast<int>(table.size()));
  std::vector<std::unordered_map<std::uint64_t, int>> local_counts(static_cast<std::size_t>(workers));

  subpar::parallelize_range<true>(workers, static_cast<int>(table.size()), [&](int worker, int start, int length) {
    auto& counts = local_counts[static_cast<std::size_t>(worker)];
    std::vector<int> nearest_cell_by_type(static_cast<std::size_t>(n_types), -1);
    std::vector<int> touched_types;
    touched_types.reserve(static_cast<std::size_t>(n_types));

    for (int row = start; row < start + length; ++row) {
      const int target_cell = table.cell_index[static_cast<std::size_t>(row)];
      if (target_cell < 0 ||
          data.cell_type_index[static_cast<std::size_t>(target_cell)] < 0) {
        continue;
      }
      const auto hits = global_index.query(row, options.candidate_k, false);
      for (const auto& hit : hits) {
        const int source_cell = table.cell_index[static_cast<std::size_t>(hit.index)];
        if (source_cell == target_cell || source_cell < 0) {
          continue;
        }
        const int source_type = data.cell_type_index[static_cast<std::size_t>(source_cell)];
        if (source_type < 0 || nearest_cell_by_type[static_cast<std::size_t>(source_type)] >= 0) {
          continue;
        }
        nearest_cell_by_type[static_cast<std::size_t>(source_type)] = source_cell;
        touched_types.push_back(source_type);
      }

      for (const int source_type : touched_types) {
        const int source_cell = nearest_cell_by_type[static_cast<std::size_t>(source_type)];
        counts[pair_key(target_cell, source_cell)] += 1;
        nearest_cell_by_type[static_cast<std::size_t>(source_type)] = -1;
      }
      touched_types.clear();
    }
  });

  std::unordered_map<std::uint64_t, int> pair_counts;
  for (const auto& local : local_counts) {
    for (const auto& entry : local) {
      pair_counts[entry.first] += entry.second;
    }
  }

  std::unordered_map<std::uint64_t, std::vector<BridgeCandidatePair>> by_type_pair;
  std::unordered_map<std::uint64_t, int> type_pair_contacts;
  by_type_pair.reserve(pair_counts.size());
  for (const auto& entry : pair_counts) {
    const int target_cell = key_first(entry.first);
    const int source_cell = key_second(entry.first);
    const int target_type = data.cell_type_index[static_cast<std::size_t>(target_cell)];
    const int source_type = data.cell_type_index[static_cast<std::size_t>(source_cell)];
    if (target_type < 0 || source_type < 0) {
      continue;
    }
    const auto type_key = pair_key(target_type, source_type);
    by_type_pair[type_key].push_back(
        {target_cell, source_cell, target_type, source_type, entry.second});
    type_pair_contacts[type_key] += entry.second;
  }

  std::vector<BridgeCandidatePair> out;
  for (auto& entry : by_type_pair) {
    if (type_pair_contacts[entry.first] < options.min_type_pair_contacts) {
      continue;
    }
    auto& pairs = entry.second;
    std::sort(pairs.begin(), pairs.end(), [](const auto& lhs, const auto& rhs) {
      if (lhs.contact_count != rhs.contact_count) return lhs.contact_count > rhs.contact_count;
      if (lhs.target_cell != rhs.target_cell) return lhs.target_cell < rhs.target_cell;
      return lhs.source_cell < rhs.source_cell;
    });

    std::unordered_set<int> used_targets;
    std::unordered_set<int> used_sources;
    std::unordered_set<int> used_any_same_type;
    used_targets.reserve(pairs.size());
    used_sources.reserve(pairs.size());
    used_any_same_type.reserve(pairs.size() * 2U);

    for (const auto& candidate : pairs) {
      if (candidate.target_type == candidate.source_type) {
        if (used_any_same_type.count(candidate.target_cell) ||
            used_any_same_type.count(candidate.source_cell)) {
          continue;
        }
        used_any_same_type.insert(candidate.target_cell);
        used_any_same_type.insert(candidate.source_cell);
      } else {
        if (used_targets.count(candidate.target_cell) ||
            used_sources.count(candidate.source_cell)) {
          continue;
        }
        used_targets.insert(candidate.target_cell);
        used_sources.insert(candidate.source_cell);
      }
      out.push_back(candidate);
    }
  }

  std::sort(out.begin(), out.end(), [](const auto& lhs, const auto& rhs) {
    if (lhs.target_type != rhs.target_type) return lhs.target_type < rhs.target_type;
    if (lhs.source_type != rhs.source_type) return lhs.source_type < rhs.source_type;
    if (lhs.contact_count != rhs.contact_count) return lhs.contact_count > rhs.contact_count;
    if (lhs.target_cell != rhs.target_cell) return lhs.target_cell < rhs.target_cell;
    return lhs.source_cell < rhs.source_cell;
  });
  return out;
}

std::vector<BridgeCandidatePair> discover_cell_center_bridge_candidates(
    const TranscriptTable& table,
    const BridgePreparedData& data,
    const BridgeTestOptions& options) {
  if (options.cell_candidate_k <= 0) {
    return {};
  }

  const int n_cells = static_cast<int>(data.members_by_cell.size());
  TranscriptTable centroids;
  centroids.x = data.centroid_x;
  centroids.y = data.centroid_y;
  centroids.z = data.centroid_z;
  const SpatialKnnIndex cell_index(centroids);

  std::vector<double> radius(static_cast<std::size_t>(n_cells), 0.0);
  for (int cell = 0; cell < n_cells; ++cell) {
    radius[static_cast<std::size_t>(cell)] = cell_radius_xy(
        table,
        data.members_by_cell[static_cast<std::size_t>(cell)],
        data.centroid_x[static_cast<std::size_t>(cell)],
        data.centroid_y[static_cast<std::size_t>(cell)]);
  }

  const int workers = effective_threads(options.num_threads, n_cells);
  std::vector<std::unordered_map<std::uint64_t, std::vector<BridgeCandidatePair>>> local_by_type(
      static_cast<std::size_t>(workers));

  subpar::parallelize_range<true>(workers, n_cells, [&](int worker, int start, int length) {
    auto& by_type = local_by_type[static_cast<std::size_t>(worker)];
    for (int target_cell = start; target_cell < start + length; ++target_cell) {
      const int target_type = data.cell_type_index[static_cast<std::size_t>(target_cell)];
      if (target_type < 0 ||
          data.members_by_cell[static_cast<std::size_t>(target_cell)].empty()) {
        continue;
      }
      const auto hits = cell_index.query(target_cell, options.cell_candidate_k, false);
      for (const auto& hit : hits) {
        const int source_cell = hit.index;
        if (source_cell == target_cell || source_cell < 0 ||
            source_cell >= n_cells ||
            data.members_by_cell[static_cast<std::size_t>(source_cell)].empty()) {
          continue;
        }

        const double dx = data.centroid_x[static_cast<std::size_t>(target_cell)] -
            data.centroid_x[static_cast<std::size_t>(source_cell)];
        const double dy = data.centroid_y[static_cast<std::size_t>(target_cell)] -
            data.centroid_y[static_cast<std::size_t>(source_cell)];
        const double distance = std::sqrt(dx * dx + dy * dy);
        const double radius_sum =
            radius[static_cast<std::size_t>(target_cell)] +
            radius[static_cast<std::size_t>(source_cell)];
        if (options.cell_candidate_halo >= 0.0 &&
            distance > radius_sum + options.cell_candidate_halo) {
          continue;
        }

        const double normalized_distance = distance / std::max(radius_sum, 1e-6);
        const int contact_rank = std::max(
            1,
            static_cast<int>(std::round(1000000.0 / (1.0 + normalized_distance))));
        const int source_type = data.cell_type_index[static_cast<std::size_t>(source_cell)];
        if (source_type < 0) {
          continue;
        }
        by_type[pair_key(target_type, source_type)].push_back(
            {target_cell, source_cell, target_type, source_type, contact_rank});
      }
    }
  });

  std::unordered_map<std::uint64_t, std::vector<BridgeCandidatePair>> by_type_pair;
  for (auto& local : local_by_type) {
    for (auto& entry : local) {
      auto& dest = by_type_pair[entry.first];
      dest.insert(dest.end(), entry.second.begin(), entry.second.end());
    }
  }

  const int max_per_type_pair = options.candidate_pairs_per_type_pair > 0
      ? options.candidate_pairs_per_type_pair
      : std::max(options.max_cells_per_type_pair * 5, options.max_cells_per_type_pair);
  std::vector<BridgeCandidatePair> out;
  for (auto& entry : by_type_pair) {
    auto& pairs = entry.second;
    if (static_cast<int>(pairs.size()) < options.min_type_pair_contacts) {
      continue;
    }
    std::sort(pairs.begin(), pairs.end(), [](const auto& lhs, const auto& rhs) {
      if (lhs.contact_count != rhs.contact_count) return lhs.contact_count > rhs.contact_count;
      if (lhs.target_cell != rhs.target_cell) return lhs.target_cell < rhs.target_cell;
      return lhs.source_cell < rhs.source_cell;
    });

    std::unordered_set<int> used_targets;
    std::unordered_set<int> used_sources;
    std::unordered_set<int> used_any_same_type;
    used_targets.reserve(pairs.size());
    used_sources.reserve(pairs.size());
    used_any_same_type.reserve(pairs.size() * 2U);

    int kept = 0;
    for (const auto& candidate : pairs) {
      if (candidate.target_type == candidate.source_type) {
        if (used_any_same_type.count(candidate.target_cell) ||
            used_any_same_type.count(candidate.source_cell)) {
          continue;
        }
        used_any_same_type.insert(candidate.target_cell);
        used_any_same_type.insert(candidate.source_cell);
      } else {
        if (used_targets.count(candidate.target_cell) ||
            used_sources.count(candidate.source_cell)) {
          continue;
        }
        used_targets.insert(candidate.target_cell);
        used_sources.insert(candidate.source_cell);
      }
      out.push_back(candidate);
      ++kept;
      if (max_per_type_pair > 0 && kept >= max_per_type_pair) {
        break;
      }
    }
  }

  std::sort(out.begin(), out.end(), [](const auto& lhs, const auto& rhs) {
    if (lhs.target_type != rhs.target_type) return lhs.target_type < rhs.target_type;
    if (lhs.source_type != rhs.source_type) return lhs.source_type < rhs.source_type;
    if (lhs.contact_count != rhs.contact_count) return lhs.contact_count > rhs.contact_count;
    if (lhs.target_cell != rhs.target_cell) return lhs.target_cell < rhs.target_cell;
    return lhs.source_cell < rhs.source_cell;
  });
  return out;
}

std::vector<BridgeCandidatePair> discover_bridge_candidates(
    const TranscriptTable& table,
    const BridgePreparedData& data,
    const BridgeTestOptions& options) {
  if (options.candidate_mode == "molecule_global") {
    return discover_molecule_global_bridge_candidates(table, data, options);
  }
  if (options.candidate_mode == "cell_center") {
    return discover_cell_center_bridge_candidates(table, data, options);
  }
  throw std::runtime_error("unknown bridge candidate_mode: " + options.candidate_mode);
}

std::vector<BridgePairScore> score_candidate_pair(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    const BridgePreparedData& data,
    const BridgeCandidatePair& candidate,
    const BridgeTestOptions& options) {
  const auto& target_members = data.members_by_cell[static_cast<std::size_t>(candidate.target_cell)];
  const auto& source_members = data.members_by_cell[static_cast<std::size_t>(candidate.source_cell)];
  if (target_members.empty() || source_members.empty()) {
    return {};
  }

  std::vector<int> combined;
  combined.reserve(target_members.size() + source_members.size());
  combined.insert(combined.end(), target_members.begin(), target_members.end());
  combined.insert(combined.end(), source_members.begin(), source_members.end());
  if (combined.size() <= 1U) {
    return {};
  }

  const SpatialKnnIndex pair_index(table, combined);
  std::vector<int> crossing_by_factor(static_cast<std::size_t>(data.n_factors), 0);
  int total_crossing = 0;

  for (const int query : target_members) {
    const auto hits = pair_index.query(query, options.crossing_k, true);
    bool crosses = false;
    for (const auto& hit : hits) {
      if (table.cell_index[static_cast<std::size_t>(hit.index)] == candidate.source_cell) {
        crosses = true;
        break;
      }
    }
    if (!crosses) {
      continue;
    }
    ++total_crossing;
    const int factor = labels[static_cast<std::size_t>(query)];
    if (factor >= 0 && factor < data.n_factors) {
      crossing_by_factor[static_cast<std::size_t>(factor)] += 1;
    }
  }

  std::vector<BridgePairScore> out;
  out.reserve(static_cast<std::size_t>(data.n_factors));
  for (int factor = 0; factor < data.n_factors; ++factor) {
    const int factor_count =
        data.factor_counts[static_cast<std::size_t>(candidate.target_cell * data.n_factors + factor)];
    if (factor_count < options.min_factor_molecules) {
      continue;
    }
    const double target_fraction =
        data.factor_fractions[static_cast<std::size_t>(candidate.target_cell * data.n_factors + factor)];
    const double source_fraction =
        data.factor_fractions[static_cast<std::size_t>(candidate.source_cell * data.n_factors + factor)];
    const double crossing_fraction =
        static_cast<double>(crossing_by_factor[static_cast<std::size_t>(factor)]) /
        static_cast<double>(std::max(factor_count, 1));

    BridgePairScore row;
    row.target_cell = candidate.target_cell;
    row.source_cell = candidate.source_cell;
    row.target_type = candidate.target_type;
    row.source_type = candidate.source_type;
    row.factor = factor;
    row.factor_count = factor_count;
    row.total_crossing_count = total_crossing;
    row.target_fraction = target_fraction;
    row.source_fraction = source_fraction;
    row.crossing_fraction = crossing_fraction;
    row.score = std::abs(source_fraction - target_fraction) * crossing_fraction;
    out.push_back(row);
  }
  return out;
}

std::vector<BridgePairScore> score_bridge_candidates(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    const BridgePreparedData& data,
    const std::vector<BridgeCandidatePair>& candidates,
    const BridgeTestOptions& options) {
  const int workers = effective_threads(options.num_threads, static_cast<int>(candidates.size()));
  std::vector<std::vector<BridgePairScore>> local_scores(static_cast<std::size_t>(workers));
  subpar::parallelize_range<true>(workers, static_cast<int>(candidates.size()), [&](int worker, int start, int length) {
    auto& out = local_scores[static_cast<std::size_t>(worker)];
    for (int i = start; i < start + length; ++i) {
      auto scores = score_candidate_pair(table, labels, data, candidates[static_cast<std::size_t>(i)], options);
      if (options.candidate_mode == "cell_center" &&
          !scores.empty() &&
          scores.front().total_crossing_count <= 0) {
        continue;
      }
      out.insert(out.end(), scores.begin(), scores.end());
    }
  });

  std::size_t total = 0;
  for (const auto& local : local_scores) total += local.size();
  std::vector<BridgePairScore> out;
  out.reserve(total);
  for (auto& local : local_scores) {
    out.insert(out.end(), local.begin(), local.end());
  }
  std::sort(out.begin(), out.end(), [](const auto& lhs, const auto& rhs) {
    if (lhs.target_type != rhs.target_type) return lhs.target_type < rhs.target_type;
    if (lhs.source_type != rhs.source_type) return lhs.source_type < rhs.source_type;
    if (lhs.factor != rhs.factor) return lhs.factor < rhs.factor;
    if (lhs.score != rhs.score) return lhs.score > rhs.score;
    if (lhs.target_cell != rhs.target_cell) return lhs.target_cell < rhs.target_cell;
    return lhs.source_cell < rhs.source_cell;
  });
  return out;
}

double quantile_75(std::vector<double> values) {
  if (values.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  std::sort(values.begin(), values.end());
  const double pos = 0.75 * static_cast<double>(values.size() - 1U);
  const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
  const std::size_t hi = static_cast<std::size_t>(std::ceil(pos));
  if (lo == hi) {
    return values[lo];
  }
  const double frac = pos - static_cast<double>(lo);
  return values[lo] * (1.0 - frac) + values[hi] * frac;
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
  std::sort(abs_and_sign.begin(), abs_and_sign.end(), [](const auto& lhs, const auto& rhs) {
    return lhs.first < rhs.first;
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

double mean_value(const std::vector<double>& values) {
  if (values.empty()) {
    return 0.0;
  }
  return std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
}

double cell_radius_xy(
    const TranscriptTable& table,
    const std::vector<int>& members,
    double cx,
    double cy) {
  double radius = 0.0;
  for (const int row : members) {
    const double dx = table.x[static_cast<std::size_t>(row)] - cx;
    const double dy = table.y[static_cast<std::size_t>(row)] - cy;
    radius = std::max(radius, std::sqrt(dx * dx + dy * dy));
  }
  return radius;
}

struct ShiftedSource {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
};

struct CrossingProfileStats {
  long long calls = 0;
  long long target_queries = 0;
  double table_sec = 0.0;
  double index_sec = 0.0;
  double query_sec = 0.0;

  void add(const CrossingProfileStats& other) {
    calls += other.calls;
    target_queries += other.target_queries;
    table_sec += other.table_sec;
    index_sec += other.index_sec;
    query_sec += other.query_sec;
  }
};

struct NullProfileStats {
  long long pair_iterations = 0;
  long long failed_replacements = 0;
  long long empty_pairs = 0;
  long long replacement_candidates = 0;
  long long replacement_pool_size = 0;
  long long movement_checks = 0;
  long long movement_reached = 0;
  long long movement_unreached = 0;
  long long final_crossing_sum = 0;
  long long target_cache_queries = 0;
  double target_cache_sec = 0.0;
  double replacement_sec = 0.0;
  double initialize_sec = 0.0;
  double movement_sec = 0.0;
  double final_sec = 0.0;
  CrossingProfileStats movement_crossing;
  CrossingProfileStats final_crossing;

  void add(const NullProfileStats& other) {
    pair_iterations += other.pair_iterations;
    failed_replacements += other.failed_replacements;
    empty_pairs += other.empty_pairs;
    replacement_candidates += other.replacement_candidates;
    replacement_pool_size += other.replacement_pool_size;
    movement_checks += other.movement_checks;
    movement_reached += other.movement_reached;
    movement_unreached += other.movement_unreached;
    final_crossing_sum += other.final_crossing_sum;
    target_cache_queries += other.target_cache_queries;
    target_cache_sec += other.target_cache_sec;
    replacement_sec += other.replacement_sec;
    initialize_sec += other.initialize_sec;
    movement_sec += other.movement_sec;
    final_sec += other.final_sec;
    movement_crossing.add(other.movement_crossing);
    final_crossing.add(other.final_crossing);
  }
};

struct ReplacementSample {
  int cell = -1;
  int candidate_count = 0;
  int pool_size = 0;
};

struct TargetCrossingCache {
  std::vector<int> members;
  std::vector<int> factor_labels;
  std::vector<double> target_kth_squared_distance;
};

struct SourcePlacement {
  std::unique_ptr<SpatialKnnIndex> source_index;
  double base_cx = 0.0;
  double base_cy = 0.0;
  double shift_x = 0.0;
  double shift_y = 0.0;
};

int count_shifted_crossing(
    const TranscriptTable& table,
    const std::vector<int>& target_members,
    const std::vector<int>& labels,
    const ShiftedSource& source,
    int factor,
    int crossing_k,
    CrossingProfileStats* profile = nullptr) {
  if (profile != nullptr) {
    ++profile->calls;
  }
  TranscriptTable tmp;
  tmp.x.reserve(target_members.size() + source.x.size());
  tmp.y.reserve(target_members.size() + source.x.size());
  tmp.z.reserve(target_members.size() + source.x.size());
  auto stage_start = std::chrono::steady_clock::now();
  for (const int row : target_members) {
    tmp.x.push_back(table.x[static_cast<std::size_t>(row)]);
    tmp.y.push_back(table.y[static_cast<std::size_t>(row)]);
    tmp.z.push_back(table.z[static_cast<std::size_t>(row)]);
  }
  tmp.x.insert(tmp.x.end(), source.x.begin(), source.x.end());
  tmp.y.insert(tmp.y.end(), source.y.begin(), source.y.end());
  tmp.z.insert(tmp.z.end(), source.z.begin(), source.z.end());
  if (profile != nullptr) {
    profile->table_sec += elapsed_seconds(stage_start);
  }
  if (tmp.x.size() <= 1U || target_members.empty() || source.x.empty()) {
    return 0;
  }

  stage_start = std::chrono::steady_clock::now();
  const SpatialKnnIndex index(tmp);
  if (profile != nullptr) {
    profile->index_sec += elapsed_seconds(stage_start);
  }
  const int source_start = static_cast<int>(target_members.size());
  int crossing = 0;
  stage_start = std::chrono::steady_clock::now();
  for (std::size_t local = 0; local < target_members.size(); ++local) {
    if (factor >= 0 && labels[static_cast<std::size_t>(target_members[local])] != factor) {
      continue;
    }
    const auto hits = index.query(static_cast<int>(local), crossing_k, true);
    bool crosses = false;
    for (const auto& hit : hits) {
      if (hit.index >= source_start) {
        crosses = true;
        break;
      }
    }
    crossing += crosses ? 1 : 0;
  }
  if (profile != nullptr) {
    profile->target_queries += static_cast<long long>(target_members.size());
    profile->query_sec += elapsed_seconds(stage_start);
  }
  return crossing;
}

std::vector<int> count_shifted_crossing_by_factor(
    const TranscriptTable& table,
    const std::vector<int>& target_members,
    const std::vector<int>& labels,
    const ShiftedSource& source,
    int n_factors,
    int crossing_k,
    CrossingProfileStats* profile = nullptr) {
  if (profile != nullptr) {
    ++profile->calls;
  }
  std::vector<int> crossing_by_factor(static_cast<std::size_t>(n_factors), 0);
  TranscriptTable tmp;
  tmp.x.reserve(target_members.size() + source.x.size());
  tmp.y.reserve(target_members.size() + source.x.size());
  tmp.z.reserve(target_members.size() + source.x.size());
  auto stage_start = std::chrono::steady_clock::now();
  for (const int row : target_members) {
    tmp.x.push_back(table.x[static_cast<std::size_t>(row)]);
    tmp.y.push_back(table.y[static_cast<std::size_t>(row)]);
    tmp.z.push_back(table.z[static_cast<std::size_t>(row)]);
  }
  tmp.x.insert(tmp.x.end(), source.x.begin(), source.x.end());
  tmp.y.insert(tmp.y.end(), source.y.begin(), source.y.end());
  tmp.z.insert(tmp.z.end(), source.z.begin(), source.z.end());
  if (profile != nullptr) {
    profile->table_sec += elapsed_seconds(stage_start);
  }
  if (tmp.x.size() <= 1U || target_members.empty() || source.x.empty()) {
    return crossing_by_factor;
  }

  stage_start = std::chrono::steady_clock::now();
  const SpatialKnnIndex index(tmp);
  if (profile != nullptr) {
    profile->index_sec += elapsed_seconds(stage_start);
  }
  const int source_start = static_cast<int>(target_members.size());
  stage_start = std::chrono::steady_clock::now();
  for (std::size_t local = 0; local < target_members.size(); ++local) {
    const auto hits = index.query(static_cast<int>(local), crossing_k, true);
    bool crosses = false;
    for (const auto& hit : hits) {
      if (hit.index >= source_start) {
        crosses = true;
        break;
      }
    }
    if (!crosses) {
      continue;
    }
    const int factor = labels[static_cast<std::size_t>(target_members[local])];
    if (factor >= 0 && factor < n_factors) {
      crossing_by_factor[static_cast<std::size_t>(factor)] += 1;
    }
  }
  if (profile != nullptr) {
    profile->target_queries += static_cast<long long>(target_members.size());
    profile->query_sec += elapsed_seconds(stage_start);
  }
  return crossing_by_factor;
}

TargetCrossingCache build_target_crossing_cache(
    const TranscriptTable& table,
    const std::vector<int>& target_members,
    const std::vector<int>& labels,
    int crossing_k,
    NullProfileStats* profile = nullptr) {
  TargetCrossingCache cache;
  cache.members = target_members;
  cache.factor_labels.reserve(target_members.size());
  cache.target_kth_squared_distance.assign(target_members.size(), std::numeric_limits<double>::infinity());
  for (const int row : target_members) {
    cache.factor_labels.push_back(labels[static_cast<std::size_t>(row)]);
  }

  if (target_members.empty() || crossing_k <= 0 ||
      static_cast<int>(target_members.size()) < crossing_k) {
    return cache;
  }

  const SpatialKnnIndex target_index(table, target_members);
  for (std::size_t local = 0; local < target_members.size(); ++local) {
    const auto hits = target_index.query(target_members[local], crossing_k, true);
    if (profile != nullptr) {
      ++profile->target_cache_queries;
    }
    if (static_cast<int>(hits.size()) >= crossing_k) {
      cache.target_kth_squared_distance[local] = hits.back().squared_distance;
    }
  }
  return cache;
}

SourcePlacement initialize_source_placement(
    const TranscriptTable& table,
    const BridgePreparedData& data,
    const std::vector<int>& source_members,
    int target_cell,
    int source_cell,
    std::mt19937& rng) {
  const double target_cx = data.centroid_x[static_cast<std::size_t>(target_cell)];
  const double target_cy = data.centroid_y[static_cast<std::size_t>(target_cell)];
  const double source_cx = data.centroid_x[static_cast<std::size_t>(source_cell)];
  const double source_cy = data.centroid_y[static_cast<std::size_t>(source_cell)];
  const auto& target_members = data.members_by_cell[static_cast<std::size_t>(target_cell)];
  const double radius =
      cell_radius_xy(table, target_members, target_cx, target_cy) +
      cell_radius_xy(table, source_members, source_cx, source_cy);
  std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * std::acos(-1.0));
  const double angle = angle_dist(rng);
  const double desired_cx = target_cx + radius * std::cos(angle);
  const double desired_cy = target_cy + radius * std::sin(angle);

  SourcePlacement out;
  out.source_index = std::make_unique<SpatialKnnIndex>(table, source_members);
  out.base_cx = source_cx;
  out.base_cy = source_cy;
  out.shift_x = desired_cx - source_cx;
  out.shift_y = desired_cy - source_cy;
  return out;
}

void move_source_toward_target(
    SourcePlacement& source,
    double target_cx,
    double target_cy,
    double step) {
  const double source_cx = source.base_cx + source.shift_x;
  const double source_cy = source.base_cy + source.shift_y;
  source.shift_x += step * (target_cx - source_cx);
  source.shift_y += step * (target_cy - source_cy);
}

bool shifted_source_crosses_target(
    const TranscriptTable& table,
    const TargetCrossingCache& target,
    const SourcePlacement& source,
    std::size_t local,
    CrossingProfileStats* profile) {
  const double threshold = target.target_kth_squared_distance[local];
  if (std::isinf(threshold)) {
    return true;
  }
  const int row = target.members[local];
  const auto hits = source.source_index->query_point(
      table.x[static_cast<std::size_t>(row)] - source.shift_x,
      table.y[static_cast<std::size_t>(row)] - source.shift_y,
      table.z[static_cast<std::size_t>(row)],
      1);
  if (profile != nullptr) {
    ++profile->target_queries;
  }
  return !hits.empty() && hits.front().squared_distance < threshold;
}

int count_shifted_crossing_fast(
    const TranscriptTable& table,
    const TargetCrossingCache& target,
    const SourcePlacement& source,
    int factor,
    CrossingProfileStats* profile = nullptr) {
  if (profile != nullptr) {
    ++profile->calls;
  }
  int crossing = 0;
  const auto query_start = std::chrono::steady_clock::now();
  for (std::size_t local = 0; local < target.members.size(); ++local) {
    if (factor >= 0 && target.factor_labels[local] != factor) {
      continue;
    }
    crossing += shifted_source_crosses_target(table, target, source, local, profile) ? 1 : 0;
  }
  if (profile != nullptr) {
    profile->query_sec += elapsed_seconds(query_start);
  }
  return crossing;
}

std::vector<int> count_shifted_crossing_by_factor_fast(
    const TranscriptTable& table,
    const TargetCrossingCache& target,
    const SourcePlacement& source,
    int n_factors,
    CrossingProfileStats* profile = nullptr) {
  if (profile != nullptr) {
    ++profile->calls;
  }
  std::vector<int> crossing_by_factor(static_cast<std::size_t>(n_factors), 0);
  const auto query_start = std::chrono::steady_clock::now();
  for (std::size_t local = 0; local < target.members.size(); ++local) {
    if (!shifted_source_crosses_target(table, target, source, local, profile)) {
      continue;
    }
    const int factor = target.factor_labels[local];
    if (factor >= 0 && factor < n_factors) {
      crossing_by_factor[static_cast<std::size_t>(factor)] += 1;
    }
  }
  if (profile != nullptr) {
    profile->query_sec += elapsed_seconds(query_start);
  }
  return crossing_by_factor;
}

ShiftedSource initialize_shifted_source(
    const TranscriptTable& table,
    const BridgePreparedData& data,
    const std::vector<int>& source_members,
    int target_cell,
    int source_cell,
    std::mt19937& rng) {
  const double target_cx = data.centroid_x[static_cast<std::size_t>(target_cell)];
  const double target_cy = data.centroid_y[static_cast<std::size_t>(target_cell)];
  const double source_cx = data.centroid_x[static_cast<std::size_t>(source_cell)];
  const double source_cy = data.centroid_y[static_cast<std::size_t>(source_cell)];
  const auto& target_members = data.members_by_cell[static_cast<std::size_t>(target_cell)];
  const double radius =
      cell_radius_xy(table, target_members, target_cx, target_cy) +
      cell_radius_xy(table, source_members, source_cx, source_cy);
  std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * std::acos(-1.0));
  const double angle = angle_dist(rng);
  const double desired_cx = target_cx + radius * std::cos(angle);
  const double desired_cy = target_cy + radius * std::sin(angle);
  const double shift_x = desired_cx - source_cx;
  const double shift_y = desired_cy - source_cy;

  ShiftedSource out;
  out.x.reserve(source_members.size());
  out.y.reserve(source_members.size());
  out.z.reserve(source_members.size());
  for (const int row : source_members) {
    out.x.push_back(table.x[static_cast<std::size_t>(row)] + shift_x);
    out.y.push_back(table.y[static_cast<std::size_t>(row)] + shift_y);
    out.z.push_back(table.z[static_cast<std::size_t>(row)]);
  }
  return out;
}

void move_source_toward_target(
    ShiftedSource& source,
    double target_cx,
    double target_cy,
    double step) {
  if (source.x.empty()) {
    return;
  }
  const double source_cx = mean_value(source.x);
  const double source_cy = mean_value(source.y);
  const double shift_x = step * (target_cx - source_cx);
  const double shift_y = step * (target_cy - source_cy);
  for (std::size_t i = 0; i < source.x.size(); ++i) {
    source.x[i] += shift_x;
    source.y[i] += shift_y;
  }
}

double factor_distance(
    const BridgePreparedData& data,
    int left_cell,
    int right_cell) {
  double total = 0.0;
  for (int factor = 0; factor < data.n_factors; ++factor) {
    const double diff =
        data.factor_fractions[static_cast<std::size_t>(left_cell * data.n_factors + factor)] -
        data.factor_fractions[static_cast<std::size_t>(right_cell * data.n_factors + factor)];
    total += diff * diff;
  }
  return std::sqrt(total);
}

ReplacementSample sample_replacement_source_cell_detailed(
    const BridgePreparedData& data,
    const BridgePairScore& row,
    const BridgeTestOptions& options,
    unsigned int seed) {
  ReplacementSample out;
  std::vector<std::pair<double, int>> candidates;
  const int n_cells = static_cast<int>(data.members_by_cell.size());
  candidates.reserve(static_cast<std::size_t>(n_cells));
  for (int cell = 0; cell < n_cells; ++cell) {
    if (cell == row.source_cell || cell == row.target_cell) {
      continue;
    }
    if (data.cell_type_index[static_cast<std::size_t>(cell)] != row.source_type) {
      continue;
    }
    if (data.members_by_cell[static_cast<std::size_t>(cell)].empty()) {
      continue;
    }
    candidates.emplace_back(factor_distance(data, row.source_cell, cell), cell);
  }
  out.candidate_count = static_cast<int>(candidates.size());
  if (candidates.empty()) {
    return out;
  }
  std::sort(candidates.begin(), candidates.end(), [](const auto& lhs, const auto& rhs) {
    if (lhs.first != rhs.first) return lhs.first < rhs.first;
    return lhs.second < rhs.second;
  });
  out.pool_size = std::max<int>(
      1,
      std::min<int>(
          static_cast<int>(candidates.size()),
          static_cast<int>(std::ceil(candidates.size() * std::max(options.null_pool_fraction, 0.0))) + 1));
  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> pick(0, out.pool_size - 1);
  out.cell = candidates[static_cast<std::size_t>(pick(rng))].second;
  return out;
}

int sample_replacement_source_cell(
    const BridgePreparedData& data,
    const BridgePairScore& row,
    const BridgeTestOptions& options,
    unsigned int seed) {
  return sample_replacement_source_cell_detailed(data, row, options, seed).cell;
}

double compute_null_score(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    const BridgePreparedData& data,
    const BridgePairScore& row,
    const BridgeTestOptions& options,
    unsigned int seed) {
  const int replacement = sample_replacement_source_cell(data, row, options, seed);
  if (replacement < 0) {
    return 0.0;
  }
  const auto& target_members = data.members_by_cell[static_cast<std::size_t>(row.target_cell)];
  const auto& source_members = data.members_by_cell[static_cast<std::size_t>(replacement)];
  if (target_members.empty() || source_members.empty()) {
    return 0.0;
  }

  std::mt19937 rng(seed ^ 0x9e3779b9U);
  ShiftedSource shifted = initialize_shifted_source(
      table,
      data,
      source_members,
      row.target_cell,
      replacement,
      rng);
  const double target_cx = data.centroid_x[static_cast<std::size_t>(row.target_cell)];
  const double target_cy = data.centroid_y[static_cast<std::size_t>(row.target_cell)];
  for (int iter = 0; iter < std::max(options.null_max_iterations, 0); ++iter) {
    const int crossing = count_shifted_crossing(
        table,
        target_members,
        labels,
        shifted,
        -1,
        options.crossing_k);
    if (crossing >= row.total_crossing_count) {
      break;
    }
    move_source_toward_target(shifted, target_cx, target_cy, options.null_step);
  }

  const int factor_count =
      data.factor_counts[static_cast<std::size_t>(row.target_cell * data.n_factors + row.factor)];
  if (factor_count <= 0) {
    return 0.0;
  }
  const int crossing = count_shifted_crossing(
      table,
      target_members,
      labels,
      shifted,
      row.factor,
      options.crossing_k);
  const double crossing_fraction =
      static_cast<double>(crossing) / static_cast<double>(factor_count);
  const double target_fraction =
      data.factor_fractions[static_cast<std::size_t>(row.target_cell * data.n_factors + row.factor)];
  const double source_fraction =
      data.factor_fractions[static_cast<std::size_t>(replacement * data.n_factors + row.factor)];
  return std::abs(source_fraction - target_fraction) * crossing_fraction;
}

std::vector<double> compute_null_scores_for_pair(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    const BridgePreparedData& data,
    const BridgePairScore& row,
    const BridgeTestOptions& options,
    unsigned int seed,
    const TargetCrossingCache* target_cache = nullptr,
    NullProfileStats* profile = nullptr) {
  if (profile != nullptr) {
    ++profile->pair_iterations;
  }
  std::vector<double> out(static_cast<std::size_t>(data.n_factors), 0.0);
  auto stage_start = std::chrono::steady_clock::now();
  const auto replacement_sample = sample_replacement_source_cell_detailed(data, row, options, seed);
  if (profile != nullptr) {
    profile->replacement_sec += elapsed_seconds(stage_start);
    profile->replacement_candidates += replacement_sample.candidate_count;
    profile->replacement_pool_size += replacement_sample.pool_size;
  }
  const int replacement = replacement_sample.cell;
  if (replacement < 0) {
    if (profile != nullptr) {
      ++profile->failed_replacements;
    }
    return out;
  }
  const auto& target_members = data.members_by_cell[static_cast<std::size_t>(row.target_cell)];
  const auto& source_members = data.members_by_cell[static_cast<std::size_t>(replacement)];
  if (target_members.empty() || source_members.empty()) {
    if (profile != nullptr) {
      ++profile->empty_pairs;
    }
    return out;
  }

  std::mt19937 rng(seed ^ 0x9e3779b9U);
  const double target_cx = data.centroid_x[static_cast<std::size_t>(row.target_cell)];
  const double target_cy = data.centroid_y[static_cast<std::size_t>(row.target_cell)];
  int movement_checks = 0;
  int final_crossing = 0;
  bool reached = false;

  std::vector<int> crossing_by_factor;
  if (options.fast_null_crossing && target_cache != nullptr) {
    stage_start = std::chrono::steady_clock::now();
    SourcePlacement source = initialize_source_placement(
        table,
        data,
        source_members,
        row.target_cell,
        replacement,
        rng);
    if (profile != nullptr) {
      profile->initialize_sec += elapsed_seconds(stage_start);
    }

    stage_start = std::chrono::steady_clock::now();
    for (int iter = 0; iter < std::max(options.null_max_iterations, 0); ++iter) {
      final_crossing = count_shifted_crossing_fast(
          table,
          *target_cache,
          source,
          -1,
          profile == nullptr ? nullptr : &profile->movement_crossing);
      ++movement_checks;
      if (final_crossing >= row.total_crossing_count) {
        reached = true;
        break;
      }
      move_source_toward_target(source, target_cx, target_cy, options.null_step);
    }
    if (profile != nullptr) {
      profile->movement_sec += elapsed_seconds(stage_start);
      profile->movement_checks += movement_checks;
      profile->final_crossing_sum += final_crossing;
      if (reached) {
        ++profile->movement_reached;
      } else {
        ++profile->movement_unreached;
      }
    }

    stage_start = std::chrono::steady_clock::now();
    crossing_by_factor = count_shifted_crossing_by_factor_fast(
        table,
        *target_cache,
        source,
        data.n_factors,
        profile == nullptr ? nullptr : &profile->final_crossing);
    if (profile != nullptr) {
      profile->final_sec += elapsed_seconds(stage_start);
    }
  } else {
    stage_start = std::chrono::steady_clock::now();
    ShiftedSource shifted = initialize_shifted_source(
        table,
        data,
        source_members,
        row.target_cell,
        replacement,
        rng);
    if (profile != nullptr) {
      profile->initialize_sec += elapsed_seconds(stage_start);
    }

    stage_start = std::chrono::steady_clock::now();
    for (int iter = 0; iter < std::max(options.null_max_iterations, 0); ++iter) {
      final_crossing = count_shifted_crossing(
          table,
          target_members,
          labels,
          shifted,
          -1,
          options.crossing_k,
          profile == nullptr ? nullptr : &profile->movement_crossing);
      ++movement_checks;
      if (final_crossing >= row.total_crossing_count) {
        reached = true;
        break;
      }
      move_source_toward_target(shifted, target_cx, target_cy, options.null_step);
    }
    if (profile != nullptr) {
      profile->movement_sec += elapsed_seconds(stage_start);
      profile->movement_checks += movement_checks;
      profile->final_crossing_sum += final_crossing;
      if (reached) {
        ++profile->movement_reached;
      } else {
        ++profile->movement_unreached;
      }
    }

    stage_start = std::chrono::steady_clock::now();
    crossing_by_factor = count_shifted_crossing_by_factor(
        table,
        target_members,
        labels,
        shifted,
        data.n_factors,
        options.crossing_k,
        profile == nullptr ? nullptr : &profile->final_crossing);
    if (profile != nullptr) {
      profile->final_sec += elapsed_seconds(stage_start);
    }
  }
  for (int factor = 0; factor < data.n_factors; ++factor) {
    const int factor_count =
        data.factor_counts[static_cast<std::size_t>(row.target_cell * data.n_factors + factor)];
    if (factor_count <= 0) {
      continue;
    }
    const double crossing_fraction =
        static_cast<double>(crossing_by_factor[static_cast<std::size_t>(factor)]) /
        static_cast<double>(factor_count);
    const double target_fraction =
        data.factor_fractions[static_cast<std::size_t>(row.target_cell * data.n_factors + factor)];
    const double source_fraction =
        data.factor_fractions[static_cast<std::size_t>(replacement * data.n_factors + factor)];
    out[static_cast<std::size_t>(factor)] =
        std::abs(source_fraction - target_fraction) * crossing_fraction;
  }
  return out;
}

std::vector<std::vector<int>> select_summary_groups(
    std::vector<BridgePairScore>& pair_scores,
    const BridgeTestOptions& options) {
  std::map<std::tuple<int, int, int>, std::vector<int>> by_group;
  for (int i = 0; i < static_cast<int>(pair_scores.size()); ++i) {
    const auto& row = pair_scores[static_cast<std::size_t>(i)];
    by_group[std::make_tuple(row.target_type, row.source_type, row.factor)].push_back(i);
  }

  std::vector<std::vector<int>> groups;
  groups.reserve(by_group.size());
  for (const auto& entry : by_group) {
    auto selected = entry.second;
    if (static_cast<int>(selected.size()) < options.min_pairs) {
      continue;
    }
    std::sort(selected.begin(), selected.end(), [&](int lhs, int rhs) {
      const auto& lrow = pair_scores[static_cast<std::size_t>(lhs)];
      const auto& rrow = pair_scores[static_cast<std::size_t>(rhs)];
      if (lrow.score != rrow.score) return lrow.score > rrow.score;
      if (lrow.target_cell != rrow.target_cell) return lrow.target_cell < rrow.target_cell;
      return lrow.source_cell < rrow.source_cell;
    });
    if (options.max_cells_per_type_pair > 0 &&
        static_cast<int>(selected.size()) > options.max_cells_per_type_pair) {
      const auto [target_type, source_type, factor] = entry.first;
      const unsigned int group_seed = options.seed ^
          static_cast<unsigned int>((target_type + 1) * 73856093U) ^
          static_cast<unsigned int>((source_type + 1) * 19349663U) ^
          static_cast<unsigned int>((factor + 1) * 83492791U);
      std::mt19937 rng(group_seed);
      std::shuffle(selected.begin(), selected.end(), rng);
      selected.resize(static_cast<std::size_t>(options.max_cells_per_type_pair));
      std::sort(selected.begin(), selected.end());
    }
    for (const int idx : selected) {
      pair_scores[static_cast<std::size_t>(idx)].used_in_summary = true;
    }
    groups.push_back(std::move(selected));
  }
  return groups;
}

std::vector<double> precompute_null_scores(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    const BridgePreparedData& data,
    const std::vector<BridgePairScore>& pair_scores,
    const std::vector<std::vector<int>>& groups,
    const BridgeTestOptions& options) {
  std::vector<double> null_scores(pair_scores.size(), 0.0);
  if (!options.compute_null || options.null_iterations <= 0) {
    return null_scores;
  }

  const auto setup_start = std::chrono::steady_clock::now();
  std::unordered_map<std::uint64_t, int> pair_slot;
  std::vector<std::vector<int>> indices_by_pair;
  long long selected_rows = 0;
  for (const auto& group : groups) {
    for (const int idx : group) {
      ++selected_rows;
      const auto& row = pair_scores[static_cast<std::size_t>(idx)];
      const auto key = pair_key(row.target_cell, row.source_cell);
      auto it = pair_slot.find(key);
      if (it == pair_slot.end()) {
        const int slot = static_cast<int>(indices_by_pair.size());
        it = pair_slot.emplace(key, slot).first;
        indices_by_pair.emplace_back();
      }
      indices_by_pair[static_cast<std::size_t>(it->second)].push_back(idx);
    }
  }

  const int iterations = std::max(options.null_iterations, 1);
  const int workers = effective_threads(options.num_threads, static_cast<int>(indices_by_pair.size()));
  emit_bridge_progress(
      options,
      "Bridge null setup: groups=" + std::to_string(groups.size()) +
          ", selected_rows=" + std::to_string(selected_rows) +
          ", unique_pairs=" + std::to_string(indices_by_pair.size()) +
          ", iterations=" + std::to_string(iterations) +
          ", workers=" + std::to_string(workers) +
          ", fast_crossing=" + std::string(options.fast_null_crossing ? "TRUE" : "FALSE"),
      setup_start);

  const auto null_start = std::chrono::steady_clock::now();
  std::vector<NullProfileStats> worker_stats(static_cast<std::size_t>(workers));
  subpar::parallelize_range<true>(workers, static_cast<int>(indices_by_pair.size()), [&](int worker, int start, int length) {
    auto& stats = worker_stats[static_cast<std::size_t>(worker)];
    for (int slot = start; slot < start + length; ++slot) {
      const auto& indices = indices_by_pair[static_cast<std::size_t>(slot)];
      if (indices.empty()) {
        continue;
      }
      const auto& first = pair_scores[static_cast<std::size_t>(indices.front())];
      TargetCrossingCache target_cache;
      const TargetCrossingCache* target_cache_ptr = nullptr;
      if (options.fast_null_crossing) {
        const auto cache_start = std::chrono::steady_clock::now();
        const auto& target_members = data.members_by_cell[static_cast<std::size_t>(first.target_cell)];
        target_cache = build_target_crossing_cache(
            table,
            target_members,
            labels,
            options.crossing_k,
            &stats);
        stats.target_cache_sec += elapsed_seconds(cache_start);
        target_cache_ptr = &target_cache;
      }
      for (int iter = 0; iter < iterations; ++iter) {
        const unsigned int seed = options.seed ^
            static_cast<unsigned int>((first.target_cell + 1) * 73856093U) ^
            static_cast<unsigned int>((first.source_cell + 1) * 19349663U) ^
            static_cast<unsigned int>((iter + 1) * 83492791U);
        const auto factor_scores = compute_null_scores_for_pair(
            table,
            labels,
            data,
            first,
            options,
            seed,
            target_cache_ptr,
            &stats);
        for (const int idx : indices) {
          const auto& row = pair_scores[static_cast<std::size_t>(idx)];
          null_scores[static_cast<std::size_t>(idx)] +=
              factor_scores[static_cast<std::size_t>(row.factor)] / static_cast<double>(iterations);
        }
      }
    }
  });

  NullProfileStats stats;
  for (const auto& worker_stats_entry : worker_stats) {
    stats.add(worker_stats_entry);
  }
  const double wall_sec = elapsed_seconds(null_start);
  const double pair_iterations = static_cast<double>(std::max<long long>(stats.pair_iterations, 1));
  const double movement_checks = static_cast<double>(std::max<long long>(stats.movement_checks, 1));
  emit_bridge_progress(
      options,
      "Bridge null counters: pair_iterations=" + std::to_string(stats.pair_iterations) +
          ", failed_replacements=" + std::to_string(stats.failed_replacements) +
          ", empty_pairs=" + std::to_string(stats.empty_pairs) +
          ", movement_checks=" + std::to_string(stats.movement_checks) +
          ", avg_checks_per_pair_iter=" + format_double(stats.movement_checks / pair_iterations, 2) +
          ", movement_reached=" + std::to_string(stats.movement_reached) +
          ", movement_unreached=" + std::to_string(stats.movement_unreached) +
          ", avg_replacement_candidates=" + format_double(stats.replacement_candidates / pair_iterations, 1) +
          ", avg_replacement_pool=" + format_double(stats.replacement_pool_size / pair_iterations, 1) +
          ", avg_final_crossing=" + format_double(stats.final_crossing_sum / pair_iterations, 1),
      null_start);
  emit_bridge_progress(
      options,
      "Bridge null timing: wall=" + format_double(wall_sec) +
          "s, target_cache_sum=" + format_double(stats.target_cache_sec) +
          "s, replacement_sum=" + format_double(stats.replacement_sec) +
          "s, initialize_sum=" + format_double(stats.initialize_sec) +
          "s, movement_sum=" + format_double(stats.movement_sec) +
          "s, final_score_sum=" + format_double(stats.final_sec) + "s",
      null_start);
  emit_bridge_progress(
      options,
      "Bridge null target cache: queries=" + std::to_string(stats.target_cache_queries),
      null_start);
  emit_bridge_progress(
      options,
      "Bridge null movement kNN: calls=" + std::to_string(stats.movement_crossing.calls) +
          ", target_queries=" + std::to_string(stats.movement_crossing.target_queries) +
          ", table_sum=" + format_double(stats.movement_crossing.table_sec) +
          "s, index_sum=" + format_double(stats.movement_crossing.index_sec) +
          "s, query_sum=" + format_double(stats.movement_crossing.query_sec) + "s",
      null_start);
  emit_bridge_progress(
      options,
      "Bridge null final kNN: calls=" + std::to_string(stats.final_crossing.calls) +
          ", target_queries=" + std::to_string(stats.final_crossing.target_queries) +
          ", table_sum=" + format_double(stats.final_crossing.table_sec) +
          "s, index_sum=" + format_double(stats.final_crossing.index_sec) +
          "s, query_sum=" + format_double(stats.final_crossing.query_sec) + "s",
      null_start);

  return null_scores;
}

std::vector<BridgeSummary> summarize_bridge_scores(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    const BridgePreparedData& data,
    std::vector<BridgePairScore>& pair_scores,
    const BridgeTestOptions& options) {
  auto stage_start = std::chrono::steady_clock::now();
  auto groups = select_summary_groups(pair_scores, options);
  long long selected_rows = 0;
  for (const auto& group : groups) {
    selected_rows += static_cast<long long>(group.size());
  }
  emit_bridge_progress(
      options,
      "Selected bridge summary groups: groups=" + std::to_string(groups.size()) +
          ", selected_rows=" + std::to_string(selected_rows),
      stage_start);
  stage_start = std::chrono::steady_clock::now();
  const auto precomputed_null_scores = precompute_null_scores(
      table,
      labels,
      data,
      pair_scores,
      groups,
      options);
  emit_bridge_progress(options, "Computed bridge null scores", stage_start);
  stage_start = std::chrono::steady_clock::now();
  const int workers = effective_threads(options.num_threads, static_cast<int>(groups.size()));
  std::vector<std::vector<BridgeSummary>> local(static_cast<std::size_t>(workers));

  subpar::parallelize_range<true>(workers, static_cast<int>(groups.size()), [&](int worker, int start, int length) {
    auto& out = local[static_cast<std::size_t>(worker)];
    for (int group_idx = start; group_idx < start + length; ++group_idx) {
      const auto& indices = groups[static_cast<std::size_t>(group_idx)];
      if (indices.empty()) {
        continue;
      }
      const auto& first = pair_scores[static_cast<std::size_t>(indices.front())];
      std::vector<double> real;
      std::vector<double> null_scores;
      real.reserve(indices.size());
      null_scores.reserve(indices.size());

      for (const int idx : indices) {
        const auto& row = pair_scores[static_cast<std::size_t>(idx)];
        real.push_back(row.score);
        if (options.compute_null && options.null_iterations > 0) {
          null_scores.push_back(precomputed_null_scores[static_cast<std::size_t>(idx)]);
        } else {
          null_scores.push_back(0.0);
        }
      }

      BridgeSummary summary;
      summary.target_type = first.target_type;
      summary.source_type = first.source_type;
      summary.factor = first.factor;
      summary.n_pairs = static_cast<int>(real.size());
      summary.mean_score = mean_value(real);
      summary.mean_null_score = mean_value(null_scores);
      summary.q75_score = quantile_75(real);
      summary.p_value = options.compute_null && options.null_iterations > 0
          ? paired_wilcoxon_greater(real, null_scores)
          : std::numeric_limits<double>::quiet_NaN();
      if (std::isnan(summary.p_value)) {
        summary.neg_log10_p = std::numeric_limits<double>::quiet_NaN();
      } else if (summary.p_value > 0.0) {
        summary.neg_log10_p = -std::log10(summary.p_value);
      } else {
        summary.neg_log10_p = std::numeric_limits<double>::infinity();
      }
      out.push_back(summary);
    }
  });

  std::vector<BridgeSummary> out;
  for (auto& chunk : local) {
    out.insert(out.end(), chunk.begin(), chunk.end());
  }
  std::sort(out.begin(), out.end(), [](const auto& lhs, const auto& rhs) {
    if (lhs.target_type != rhs.target_type) return lhs.target_type < rhs.target_type;
    if (lhs.source_type != rhs.source_type) return lhs.source_type < rhs.source_type;
    return lhs.factor < rhs.factor;
  });
  emit_bridge_progress(
      options,
      "Summarized bridge scores: summaries=" + std::to_string(out.size()),
      stage_start);
  return out;
}

}  // namespace

// Enumerate cell pairs that have transcript-level spatial contacts across boundaries.
std::vector<CellPair> discover_adjacent_cell_pairs(const TranscriptTable& table, int k) {
  std::set<std::pair<int, int>> pairs;
  const SpatialKnnIndex global_index(table);
  for (std::size_t i = 0; i < table.size(); ++i) {
    const auto hits = global_index.query(static_cast<int>(i), k, false);
    for (const auto& hit : hits) {
      const int other = hit.index;
      const int cell_a = table.cell_index[i];
      const int cell_b = table.cell_index[static_cast<std::size_t>(other)];
      if (cell_a == cell_b) {
        continue;
      }
      pairs.insert(std::minmax(cell_a, cell_b));
    }
  }

  std::vector<CellPair> out;
  out.reserve(pairs.size());
  for (const auto& pair : pairs) {
    CellPair cell_pair;
    cell_pair.cell_a = pair.first;
    cell_pair.cell_b = pair.second;
    out.push_back(cell_pair);
  }
  return out;
}

// Compute the fraction of a cell's transcripts assigned to one factor.
double factor_fraction(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    int cell_index,
    int factor) {
  int total = 0;
  int match = 0;
  for (std::size_t i = 0; i < table.size(); ++i) {
    if (table.cell_index[i] != cell_index) {
      continue;
    }
    total += 1;
    if (labels[i] == factor) {
      match += 1;
    }
  }
  return total == 0 ? 0.0 : static_cast<double>(match) / static_cast<double>(total);
}

// Measure how often factor-labeled transcripts in cell A neighbor transcripts from cell B.
double crossing_fraction(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    int cell_a,
    int cell_b,
    int factor,
    int k) {
  std::vector<int> combined;
  std::vector<int> factor_members_a;
  combined.reserve(table.size());
  for (std::size_t i = 0; i < table.size(); ++i) {
    if (table.cell_index[i] == cell_a || table.cell_index[i] == cell_b) {
      combined.push_back(static_cast<int>(i));
      if (table.cell_index[i] == cell_a && labels[i] == factor) {
        factor_members_a.push_back(static_cast<int>(i));
      }
    }
  }
  if (factor_members_a.empty()) {
    return 0.0;
  }

  const SpatialKnnIndex pair_index(table, combined);
  int crossing = 0;
  for (const int query : factor_members_a) {
    const auto hits = pair_index.query(query, k, false);
    bool crosses = false;
    for (const auto& hit : hits) {
      const int candidate = hit.index;
      if (table.cell_index[static_cast<std::size_t>(candidate)] == cell_b) {
        crosses = true;
        break;
      }
    }
    crossing += crosses ? 1 : 0;
  }

  return static_cast<double>(crossing) / static_cast<double>(factor_members_a.size());
}

// Combine abundance shift and boundary crossing into one bridge evidence score.
BridgeEvidence compute_bridge_evidence(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    int cell_a,
    int cell_b,
    int factor,
    int k) {
  const double fraction_a = factor_fraction(table, labels, cell_a, factor);
  const double fraction_b = factor_fraction(table, labels, cell_b, factor);
  const double cross = crossing_fraction(table, labels, cell_a, cell_b, factor, k);
  BridgeEvidence evidence;
  evidence.pair.cell_a = cell_a;
  evidence.pair.cell_b = cell_b;
  evidence.factor = factor;
  evidence.fraction_a = fraction_a;
  evidence.fraction_b = fraction_b;
  evidence.crossing_fraction = cross;
  evidence.score = std::abs(fraction_b - fraction_a) * cross;
  return evidence;
}

BridgeTestResult run_bridge_test(
    const TranscriptTable& table,
    const CellTable& cells,
    const std::vector<int>& labels,
    const BridgeTestOptions& options) {
  if (options.candidate_mode != "molecule_global" && options.candidate_mode != "cell_center") {
    throw std::runtime_error("candidate_mode must be 'molecule_global' or 'cell_center'");
  }
  if (options.crossing_k <= 0) {
    throw std::runtime_error("crossing_k must be positive");
  }
  if (options.min_factor_molecules < 1) {
    throw std::runtime_error("min_factor_molecules must be positive");
  }
  if (options.min_pairs < 1) {
    throw std::runtime_error("min_pairs must be positive");
  }
  const auto bridge_start = std::chrono::steady_clock::now();
  emit_bridge_progress(
      options,
      "Starting bridge test: molecules=" + std::to_string(table.size()) +
          ", cells=" + std::to_string(table.num_cells()) +
          ", candidate_mode=" + options.candidate_mode +
          ", compute_null=" + std::string(options.compute_null ? "TRUE" : "FALSE"),
      bridge_start);
  BridgeTestResult result;
  auto stage_start = std::chrono::steady_clock::now();
  auto prepared = prepare_bridge_data(table, cells, labels);
  result.cell_types = prepared.type_names;
  emit_bridge_progress(
      options,
      "Prepared bridge data: cell_types=" + std::to_string(prepared.type_names.size()) +
          ", factors=" + std::to_string(prepared.n_factors) +
          ", skipped_unannotated_cells=" + std::to_string(prepared.n_unannotated_cells),
      stage_start);
  stage_start = std::chrono::steady_clock::now();
  result.candidates = discover_bridge_candidates(table, prepared, options);
  emit_bridge_progress(
      options,
      "Discovered bridge candidates: pairs=" + std::to_string(result.candidates.size()),
      stage_start);
  stage_start = std::chrono::steady_clock::now();
  result.pair_scores = score_bridge_candidates(table, labels, prepared, result.candidates, options);
  emit_bridge_progress(
      options,
      "Scored bridge candidate pairs: score_rows=" + std::to_string(result.pair_scores.size()),
      stage_start);
  stage_start = std::chrono::steady_clock::now();
  result.summaries = summarize_bridge_scores(table, labels, prepared, result.pair_scores, options);
  emit_bridge_progress(options, "Completed bridge summary stage", stage_start);
  emit_bridge_progress(
      options,
      "Finished bridge test: summaries=" + std::to_string(result.summaries.size()),
      bridge_start);
  return result;
}

}  // namespace celladmix
