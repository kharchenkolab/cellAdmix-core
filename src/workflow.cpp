#include "celladmix/workflow.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_set>

namespace celladmix {

// Dataset subsetting and summary helpers shared across loading and reporting paths.

namespace {

// Decide whether one transcript row survives the requested crop and mask filters.
bool keep_transcript_row(
    std::size_t i,
    const std::optional<std::string>& crop_id,
    const std::vector<std::string>* crop_ids,
    const std::vector<bool>* keep_mask,
    std::size_t total_size) {
  if (crop_id.has_value()) {
    if (crop_ids == nullptr || crop_ids->size() != total_size) {
      throw std::runtime_error("crop_id selection requires crop ids matching TranscriptTable size");
    }
    if ((*crop_ids)[i] != *crop_id) {
      return false;
    }
  }
  if (keep_mask != nullptr) {
    if (keep_mask->size() != total_size) {
      throw std::runtime_error("keep_mask must match TranscriptTable size");
    }
    if (!(*keep_mask)[i]) {
      return false;
    }
  }
  return true;
}

// Pre-size subset output vectors to avoid repeated reallocation during filtering.
void reserve_subset_transcript_capacity(TranscriptTable& out, const TranscriptTable& in, std::size_t n) {
  out.gene_key.reserve(n);
  out.orig_cell_id.reserve(n);
  out.x.reserve(n);
  out.y.reserve(n);
  out.z.reserve(n);
  if (!in.transcript_ids.empty()) out.transcript_ids.reserve(n);
  if (!in.sample_ids.empty()) out.sample_ids.reserve(n);
  if (!in.fov_ids.empty()) out.fov_ids.reserve(n);
  if (!in.cell_types.empty()) out.cell_types.reserve(n);
  if (!in.orig_nucleus_ids.empty()) out.orig_nucleus_ids.reserve(n);
  if (!in.qv.empty()) out.qv.reserve(n);
  if (!in.overlaps_nucleus.empty()) out.overlaps_nucleus.reserve(n);
  if (!in.nucleus_distance.empty()) out.nucleus_distance.reserve(n);
}

}  // namespace

// Subset a transcript table by crop and/or explicit keep mask.
TranscriptTable subset_transcripts(
    const TranscriptTable& table,
    const std::optional<std::string>& crop_id,
    const std::vector<std::string>* crop_ids,
    const std::vector<bool>* keep_mask) {
  TranscriptTable out;
  reserve_subset_transcript_capacity(out, table, table.size());

  for (std::size_t i = 0; i < table.size(); ++i) {
    if (!keep_transcript_row(i, crop_id, crop_ids, keep_mask, table.size())) {
      continue;
    }
    out.gene_key.push_back(table.gene_key[i]);
    out.orig_cell_id.push_back(table.orig_cell_id[i]);
    out.x.push_back(table.x[i]);
    out.y.push_back(table.y[i]);
    out.z.push_back(table.z[i]);
    if (!table.transcript_ids.empty()) out.transcript_ids.push_back(table.transcript_ids[i]);
    if (!table.sample_ids.empty()) out.sample_ids.push_back(table.sample_ids[i]);
    if (!table.fov_ids.empty()) out.fov_ids.push_back(table.fov_ids[i]);
    if (!table.cell_types.empty()) out.cell_types.push_back(table.cell_types[i]);
    if (!table.orig_nucleus_ids.empty()) out.orig_nucleus_ids.push_back(table.orig_nucleus_ids[i]);
    if (!table.qv.empty()) out.qv.push_back(table.qv[i]);
    if (!table.overlaps_nucleus.empty()) out.overlaps_nucleus.push_back(table.overlaps_nucleus[i]);
    if (!table.nucleus_distance.empty()) out.nucleus_distance.push_back(table.nucleus_distance[i]);
  }

  out.finalize();
  return out;
}

// Subset a cell table by crop and/or explicit cell-id whitelist.
CellTable subset_cells(
    const CellTable& cells,
    const std::optional<std::string>& crop_id,
    const std::vector<std::string>* crop_ids,
    const std::vector<std::string>* keep_cell_ids) {
  std::unordered_set<std::string> keep_cells;
  if (keep_cell_ids != nullptr) {
    keep_cells.reserve(keep_cell_ids->size());
    keep_cells.insert(keep_cell_ids->begin(), keep_cell_ids->end());
  }

  CellTable out;
  out.cell_ids.reserve(cells.size());
  out.sample_ids.reserve(cells.sample_ids.size());
  out.fov_ids.reserve(cells.fov_ids.size());
  out.cell_types.reserve(cells.cell_types.size());
  out.centroid_x.reserve(cells.centroid_x.size());
  out.centroid_y.reserve(cells.centroid_y.size());
  out.centroid_z.reserve(cells.centroid_z.size());

  for (std::size_t i = 0; i < cells.size(); ++i) {
    if (crop_id.has_value()) {
      if (crop_ids == nullptr || crop_ids->size() != cells.size()) {
        throw std::runtime_error("crop_id selection requires crop ids matching CellTable size");
      }
      if ((*crop_ids)[i] != *crop_id) {
        continue;
      }
    }
    if (keep_cell_ids != nullptr && keep_cells.find(cells.cell_ids[i]) == keep_cells.end()) {
      continue;
    }
    out.cell_ids.push_back(cells.cell_ids[i]);
    if (!cells.sample_ids.empty()) out.sample_ids.push_back(cells.sample_ids[i]);
    if (!cells.fov_ids.empty()) out.fov_ids.push_back(cells.fov_ids[i]);
    if (!cells.cell_types.empty()) out.cell_types.push_back(cells.cell_types[i]);
    out.centroid_x.push_back(cells.centroid_x[i]);
    out.centroid_y.push_back(cells.centroid_y[i]);
    out.centroid_z.push_back(cells.centroid_z[i]);
  }

  return out;
}

// Reconstruct a minimal cell table from transcript-level annotations.
CellTable infer_cell_table(const TranscriptTable& table) {
  CellTable out;
  out.cell_ids = table.cells;
  out.centroid_x = infer_cell_centroids_x(table);
  out.centroid_y = infer_cell_centroids_y(table);
  out.centroid_z = infer_cell_centroids_z(table);
  out.sample_ids.assign(table.num_cells(), "");
  out.fov_ids.assign(table.num_cells(), "");
  out.cell_types.assign(table.num_cells(), "");

  for (std::size_t i = 0; i < table.size(); ++i) {
    const std::size_t cell = static_cast<std::size_t>(table.cell_index[i]);
    if (!table.sample_ids.empty() && out.sample_ids[cell].empty()) {
      out.sample_ids[cell] = table.sample_ids[i];
    }
    if (!table.fov_ids.empty() && out.fov_ids[cell].empty()) {
      out.fov_ids[cell] = table.fov_ids[i];
    }
    if (!table.cell_types.empty() && out.cell_types[cell].empty()) {
      out.cell_types[cell] = table.cell_types[i];
    }
  }

  if (table.sample_ids.empty()) out.sample_ids.clear();
  if (table.fov_ids.empty()) out.fov_ids.clear();
  if (table.cell_types.empty()) out.cell_types.clear();
  return out;
}

// Summarize transcript labels into per-cell dominant factors and fractions.
std::vector<CellFactorSummary> summarize_cells(
    const TranscriptTable& table,
    const std::vector<int>* labels) {
  if (labels != nullptr && labels->size() != table.size()) {
    throw std::runtime_error("labels must match TranscriptTable size");
  }

  int n_factors = 0;
  if (labels != nullptr && !labels->empty()) {
    n_factors = *std::max_element(labels->begin(), labels->end()) + 1;
  }

  std::vector<CellFactorSummary> summaries(table.num_cells());
  const auto centroid_x = infer_cell_centroids_x(table);
  const auto centroid_y = infer_cell_centroids_y(table);
  const auto centroid_z = infer_cell_centroids_z(table);

  for (std::size_t cell = 0; cell < table.num_cells(); ++cell) {
    auto& summary = summaries[cell];
    summary.cell_id = table.cells[cell];
    summary.centroid_x = centroid_x[cell];
    summary.centroid_y = centroid_y[cell];
    summary.centroid_z = centroid_z[cell];
    summary.factor_fractions.assign(static_cast<std::size_t>(n_factors), 0.0);
  }

  for (std::size_t i = 0; i < table.size(); ++i) {
    auto& summary = summaries[static_cast<std::size_t>(table.cell_index[i])];
    summary.transcript_count += 1;
    if (!table.cell_types.empty() && summary.cell_type.empty()) {
      summary.cell_type = table.cell_types[i];
    }
    if (labels != nullptr) {
      const int factor = (*labels)[i];
      if (factor >= 0) {
        summary.factor_fractions[static_cast<std::size_t>(factor)] += 1.0;
      }
    }
  }

  if (labels != nullptr) {
    for (auto& summary : summaries) {
      if (summary.transcript_count <= 0 || summary.factor_fractions.empty()) {
        continue;
      }
      int best_factor = 0;
      double best_count = -1.0;
      for (std::size_t factor = 0; factor < summary.factor_fractions.size(); ++factor) {
        if (summary.factor_fractions[factor] > best_count) {
          best_count = summary.factor_fractions[factor];
          best_factor = static_cast<int>(factor);
        }
        summary.factor_fractions[factor] /=
            static_cast<double>(summary.transcript_count);
      }
      summary.dominant_factor = best_factor;
      summary.dominant_fraction = summary.factor_fractions[static_cast<std::size_t>(best_factor)];
    }
  }

  return summaries;
}

// Score all neighboring cell pairs for one candidate bridge factor.
std::vector<BridgeEvidence> score_bridge_evidence_all(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    int factor,
    int adjacency_k,
    int crossing_k) {
  if (labels.size() != table.size()) {
    throw std::runtime_error("labels must match TranscriptTable size");
  }
  const auto pairs = discover_adjacent_cell_pairs(table, adjacency_k);
  std::vector<BridgeEvidence> out;
  out.reserve(pairs.size());
  for (const auto& pair : pairs) {
    out.push_back(
        compute_bridge_evidence(table, labels, pair.cell_a, pair.cell_b, factor, crossing_k));
  }
  std::sort(out.begin(), out.end(), [](const auto& lhs, const auto& rhs) {
    if (lhs.score == rhs.score) {
      if (lhs.pair.cell_a == rhs.pair.cell_a) {
        return lhs.pair.cell_b < rhs.pair.cell_b;
      }
      return lhs.pair.cell_a < rhs.pair.cell_a;
    }
    return lhs.score > rhs.score;
  });
  return out;
}

}  // namespace celladmix
