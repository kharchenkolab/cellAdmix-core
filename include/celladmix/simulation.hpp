// Synthetic-data generators used by tests and quick local debugging.

#pragma once

#include <string>
#include <vector>

#include "celladmix/membrane.hpp"
#include "celladmix/types.hpp"

namespace celladmix {

struct SimulationParams {
  int transcripts_per_cell = 120;
  int admixture_per_target_cell = 20;
  unsigned int seed = 1;
};

// Small synthetic Xenium-like dataset with known admixture truth.
struct SimulatedDataset {
  TranscriptTable transcripts;
  CellTable cells;
  Image2D membrane_image;
  std::vector<bool> true_admixture;
  std::vector<std::string> malignant_markers;
  std::vector<std::string> fibro_markers;
  int malignant_cell = -1;
  int fibro_cell = -1;
};

// Simulate a four-cell NSCLC/fibroblast admixture toy dataset.
SimulatedDataset simulate_nsclc_admixture(const SimulationParams& params = {});

}  // namespace celladmix
