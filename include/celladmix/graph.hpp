// Within-cell kNN graph construction and Potts-style label smoothing.

#pragma once

#include <vector>

#include "celladmix/matrix.hpp"
#include "celladmix/types.hpp"

namespace celladmix {

struct KnnGraph {
  int n_nodes = 0;
  std::vector<int> indptr;
  std::vector<int> indices;
};

// Build a transcript-level kNN graph for a specific list of transcript members.
KnnGraph build_cell_knn_graph(const TranscriptTable& table, const std::vector<int>& members, int k);

// Build a transcript-level kNN graph for one inferred cell in the table.
KnnGraph build_cell_knn_graph(const TranscriptTable& table, int cell_index, int k);

// Smooth node labels on a fixed graph using iterative conditional modes.
std::vector<int> smooth_labels_icm(
    const DenseMatrix& node_scores,
    const KnnGraph& graph,
    double same_label_ratio = 5.0,
    int max_iterations = 20);

// Apply per-cell smoothing to dense transcript factor scores across a table.
std::vector<int> assign_factors_per_cell(
    const TranscriptTable& table,
    const DenseMatrix& node_scores,
    int k_neighbors = 10,
    double same_label_ratio = 5.0,
    int max_iterations = 20,
    int num_threads = 1);

}  // namespace celladmix
