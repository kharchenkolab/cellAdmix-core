// Bridge-scoring helpers for quantifying factor-specific evidence across
// ordered adjacent cell pairs and cell-type pairs.

#pragma once

#include <functional>
#include <string>
#include <vector>

#include "celladmix/types.hpp"

namespace celladmix {

struct CellPair {
  int cell_a = -1;
  int cell_b = -1;
};

// Summary of bridge evidence for one factor across one ordered cell pair.
struct BridgeEvidence {
  CellPair pair;
  int factor = -1;
  double fraction_a = 0.0;
  double fraction_b = 0.0;
  double crossing_fraction = 0.0;
  double score = 0.0;
};

// Controls the original-style bridge test over ordered target/source cell-type
// pairs. The candidate step uses molecule-level physical kNN, while pair
// scoring uses exact kNN inside each two-cell local neighborhood.
struct BridgeTestOptions {
  std::string candidate_mode = "molecule_global";
  int candidate_k = 10;
  int cell_candidate_k = 50;
  int candidate_pairs_per_type_pair = 1000;
  int crossing_k = 20;
  int min_type_pair_contacts = 5;
  int min_factor_molecules = 5;
  int min_pairs = 5;
  int max_cells_per_type_pair = 200;
  int null_iterations = 3;
  int null_max_iterations = 50;
  double cell_candidate_halo = -1.0;
  double null_step = 0.02;
  double null_pool_fraction = 0.1;
  int num_threads = 1;
  unsigned int seed = 1U;
  bool compute_null = true;
  bool fast_null_crossing = true;
  std::function<void(const std::string&, double)> progress;
};

// One ordered target/source cell pair selected from molecule-level adjacency.
struct BridgeCandidatePair {
  int target_cell = -1;
  int source_cell = -1;
  int target_type = -1;
  int source_type = -1;
  int contact_count = 0;
};

// Per-factor score for one ordered target/source cell pair.
struct BridgePairScore {
  int target_cell = -1;
  int source_cell = -1;
  int target_type = -1;
  int source_type = -1;
  int factor = -1;
  int factor_count = 0;
  int total_crossing_count = 0;
  double target_fraction = 0.0;
  double source_fraction = 0.0;
  double crossing_fraction = 0.0;
  double score = 0.0;
  bool used_in_summary = false;
};

// Statistical summary for one ordered target/source cell-type pair and factor.
struct BridgeSummary {
  int target_type = -1;
  int source_type = -1;
  int factor = -1;
  int n_pairs = 0;
  double mean_score = 0.0;
  double mean_null_score = 0.0;
  double q75_score = 0.0;
  double p_value = 1.0;
  double neg_log10_p = 0.0;
};

struct BridgeTestResult {
  std::vector<std::string> cell_types;
  std::vector<BridgeCandidatePair> candidates;
  std::vector<BridgePairScore> pair_scores;
  std::vector<BridgeSummary> summaries;
};

// Find unique neighboring cell pairs from transcript-level spatial proximity.
std::vector<CellPair> discover_adjacent_cell_pairs(const TranscriptTable& table, int k = 5);

// Measure the fraction of one cell assigned to a given factor label.
double factor_fraction(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    int cell_index,
    int factor);

// Measure how often factor-positive transcripts in cell_a cross into cell_b.
double crossing_fraction(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    int cell_a,
    int cell_b,
    int factor,
    int k = 10);

// Combine within-cell fractions and crossing evidence into one bridge score.
BridgeEvidence compute_bridge_evidence(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    int cell_a,
    int cell_b,
    int factor,
    int k = 10);

// Run the original-style bridge test over all ordered cell-type pairs and
// factors using fitted molecule factor labels.
BridgeTestResult run_bridge_test(
    const TranscriptTable& table,
    const CellTable& cells,
    const std::vector<int>& labels,
    const BridgeTestOptions& options = {});

}  // namespace celladmix
