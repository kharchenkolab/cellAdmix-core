// Higher-level table manipulation and summary helpers shared across workflows.

#pragma once

#include <optional>
#include <string>
#include <vector>

#include "celladmix/bridge.hpp"
#include "celladmix/types.hpp"
#include "celladmix/xenium.hpp"

namespace celladmix {

// Subset a transcript table by crop id and/or keep mask.
TranscriptTable subset_transcripts(
    const TranscriptTable& table,
    const std::optional<std::string>& crop_id = std::nullopt,
    const std::vector<std::string>* crop_ids = nullptr,
    const std::vector<bool>* keep_mask = nullptr);

// Subset a cell table by crop id and/or explicit cell ids.
CellTable subset_cells(
    const CellTable& cells,
    const std::optional<std::string>& crop_id = std::nullopt,
    const std::vector<std::string>* crop_ids = nullptr,
    const std::vector<std::string>* keep_cell_ids = nullptr);

// Infer a minimal cell table directly from a transcript table.
CellTable infer_cell_table(const TranscriptTable& table);

// Per-cell aggregate summary derived from transcript labels.
struct CellFactorSummary {
  std::string cell_id;
  std::string cell_type;
  double centroid_x = 0.0;
  double centroid_y = 0.0;
  double centroid_z = 0.0;
  int transcript_count = 0;
  int dominant_factor = -1;
  double dominant_fraction = 0.0;
  std::vector<double> factor_fractions;
};

// Summarize transcript labels as per-cell dominant factors and fractions.
std::vector<CellFactorSummary> summarize_cells(
    const TranscriptTable& table,
    const std::vector<int>* labels = nullptr);

// Score bridge evidence for all adjacent cell pairs for one factor.
std::vector<BridgeEvidence> score_bridge_evidence_all(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    int factor,
    int adjacency_k = 5,
    int crossing_k = 10);

}  // namespace celladmix
