// Core transcript and cell table types plus light derived summaries.

#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace celladmix {

// Single source of truth for the core package version.
inline constexpr const char* kCelladmixVersion = "0.0.1";

struct CellTable {
  std::vector<std::string> cell_ids;
  std::vector<std::string> sample_ids;
  std::vector<std::string> fov_ids;
  std::vector<std::string> cell_types;
  std::vector<double> centroid_x;
  std::vector<double> centroid_y;
  std::vector<double> centroid_z;

  // Return the number of cells represented in the table.
  std::size_t size() const { return cell_ids.size(); }
};

// Transcript-centric table used by the core algorithms after preprocessing.
struct TranscriptTable {
  std::vector<std::string> gene_key;
  std::vector<std::string> orig_cell_id;
  std::vector<std::string> transcript_ids;
  std::vector<std::string> sample_ids;
  std::vector<std::string> fov_ids;
  std::vector<std::string> cell_types;
  std::vector<std::string> orig_nucleus_ids;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> qv;
  std::vector<int> overlaps_nucleus;
  std::vector<double> nucleus_distance;

  std::vector<std::string> genes;
  std::vector<std::string> cells;
  std::vector<int> gene_index;
  std::vector<int> cell_index;

  // Return the number of transcript rows.
  std::size_t size() const { return x.size(); }
  // Return the number of unique genes after finalize().
  std::size_t num_genes() const { return genes.size(); }
  // Return the number of unique cells after finalize().
  std::size_t num_cells() const { return cells.size(); }

  // Report whether transcript-level cell-type annotations are available.
  bool has_cell_types() const { return cell_types.size() == size(); }
  // Report whether Xenium-style nucleus compartment hints are available.
  bool has_nucleus_compartment() const {
    return overlaps_nucleus.size() == size() || nucleus_distance.size() == size();
  }

  // Finalize derived gene/cell dictionaries and integer indices.
  void finalize() {
    const std::size_t n = size();
    if (gene_key.size() != n || orig_cell_id.size() != n || y.size() != n) {
      throw std::runtime_error("TranscriptTable fields must have consistent lengths");
    }
    if (z.empty()) {
      z.assign(n, 0.0);
    }
    if (z.size() != n) {
      throw std::runtime_error("TranscriptTable z must be empty or match x length");
    }
    if (!transcript_ids.empty() && transcript_ids.size() != n) {
      throw std::runtime_error("TranscriptTable transcript_ids must be empty or match x length");
    }
    if (!sample_ids.empty() && sample_ids.size() != n) {
      throw std::runtime_error("TranscriptTable sample_ids must be empty or match x length");
    }
    if (!fov_ids.empty() && fov_ids.size() != n) {
      throw std::runtime_error("TranscriptTable fov_ids must be empty or match x length");
    }
    if (!cell_types.empty() && cell_types.size() != n) {
      throw std::runtime_error("TranscriptTable cell_types must be empty or match x length");
    }
    if (!orig_nucleus_ids.empty() && orig_nucleus_ids.size() != n) {
      throw std::runtime_error("TranscriptTable orig_nucleus_ids must be empty or match x length");
    }
    if (!qv.empty() && qv.size() != n) {
      throw std::runtime_error("TranscriptTable qv must be empty or match x length");
    }
    if (!overlaps_nucleus.empty() && overlaps_nucleus.size() != n) {
      throw std::runtime_error("TranscriptTable overlaps_nucleus must be empty or match x length");
    }
    if (!nucleus_distance.empty() && nucleus_distance.size() != n) {
      throw std::runtime_error("TranscriptTable nucleus_distance must be empty or match x length");
    }

    genes.clear();
    cells.clear();
    gene_index.resize(n);
    cell_index.resize(n);

    std::unordered_map<std::string, int> gene_map;
    std::unordered_map<std::string, int> cell_map;
    gene_map.reserve(n);
    cell_map.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
      const auto& gene = gene_key[i];
      const auto& cell = orig_cell_id[i];

      auto gene_it = gene_map.find(gene);
      if (gene_it == gene_map.end()) {
        const int idx = static_cast<int>(genes.size());
        gene_map.emplace(gene, idx);
        genes.push_back(gene);
        gene_index[i] = idx;
      } else {
        gene_index[i] = gene_it->second;
      }

      auto cell_it = cell_map.find(cell);
      if (cell_it == cell_map.end()) {
        const int idx = static_cast<int>(cells.size());
        cell_map.emplace(cell, idx);
        cells.push_back(cell);
        cell_index[i] = idx;
      } else {
        cell_index[i] = cell_it->second;
      }
    }
  }

  // Group transcript row indices by inferred cell index.
  std::vector<std::vector<int>> transcripts_by_cell() const {
    std::vector<std::vector<int>> out(num_cells());
    for (std::size_t i = 0; i < size(); ++i) {
      out.at(static_cast<std::size_t>(cell_index.at(i))).push_back(static_cast<int>(i));
    }
    return out;
  }
};

// Infer x centroids by averaging transcript coordinates per cell.
inline std::vector<double> infer_cell_centroids_x(const TranscriptTable& table) {
  std::vector<double> sums(table.num_cells(), 0.0);
  std::vector<int> counts(table.num_cells(), 0);
  for (std::size_t i = 0; i < table.size(); ++i) {
    sums[static_cast<std::size_t>(table.cell_index[i])] += table.x[i];
    counts[static_cast<std::size_t>(table.cell_index[i])] += 1;
  }
  for (std::size_t i = 0; i < sums.size(); ++i) {
    sums[i] /= static_cast<double>(counts[i] == 0 ? 1 : counts[i]);
  }
  return sums;
}

// Infer y centroids by averaging transcript coordinates per cell.
inline std::vector<double> infer_cell_centroids_y(const TranscriptTable& table) {
  std::vector<double> sums(table.num_cells(), 0.0);
  std::vector<int> counts(table.num_cells(), 0);
  for (std::size_t i = 0; i < table.size(); ++i) {
    sums[static_cast<std::size_t>(table.cell_index[i])] += table.y[i];
    counts[static_cast<std::size_t>(table.cell_index[i])] += 1;
  }
  for (std::size_t i = 0; i < sums.size(); ++i) {
    sums[i] /= static_cast<double>(counts[i] == 0 ? 1 : counts[i]);
  }
  return sums;
}

// Infer z centroids by averaging transcript coordinates per cell.
inline std::vector<double> infer_cell_centroids_z(const TranscriptTable& table) {
  std::vector<double> sums(table.num_cells(), 0.0);
  std::vector<int> counts(table.num_cells(), 0);
  for (std::size_t i = 0; i < table.size(); ++i) {
    sums[static_cast<std::size_t>(table.cell_index[i])] += table.z[i];
    counts[static_cast<std::size_t>(table.cell_index[i])] += 1;
  }
  for (std::size_t i = 0; i < sums.size(); ++i) {
    sums[i] /= static_cast<double>(counts[i] == 0 ? 1 : counts[i]);
  }
  return sums;
}

}  // namespace celladmix
