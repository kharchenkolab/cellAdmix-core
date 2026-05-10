// Lightweight tissue-domain detection from cell positions and annotations.
//
// The first implementation builds multiscale local annotation-composition
// features, initializes domains with k-means, smooths labels with a Potts-style
// spatial prior, and extracts connected regions.

#pragma once

#include <string>
#include <vector>

#include "celladmix/matrix.hpp"
#include "celladmix/types.hpp"

namespace celladmix {

struct DomainOptions {
  std::vector<int> scales = {15, 50, 150};
  std::string transform = "clr";
  bool include_self = true;
  double self_weight = 1.0;
  bool include_center_type = false;
  double center_type_weight = 0.25;
  int n_domains = 8;
  int smooth_k = 20;
  double smooth_lambda = 0.5;
  double max_edge_distance = -1.0;
  int min_region_size = 25;
  int kmeans_iterations = 50;
  int smooth_iterations = 10;
  unsigned int seed = 1;
  int num_threads = 1;
};

struct DomainResult {
  std::vector<int> labels;
  std::vector<int> region_ids;
  std::vector<int> domain_sizes;
  std::vector<int> region_sizes;
  DenseMatrix domain_composition;
  double spatial_coherence = 0.0;
  double boundary_fraction = 1.0;
  int kmeans_iterations = 0;
  int smooth_iterations = 0;
};

// Identify spatial domains from cell coordinates and integer annotation ids.
DomainResult identify_domains(
    const CellTable& cells,
    const std::vector<int>& annotation_ids,
    int n_annotation_types,
    const DomainOptions& options = {});

}  // namespace celladmix
