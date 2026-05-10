#include "celladmix/simulation.hpp"

#include <random>
#include <stdexcept>

namespace celladmix {

// Synthetic NSCLC-like transcript simulation used for testing and examples.

namespace {

struct Rect {
  double x0;
  double x1;
  double y0;
  double y1;
};

// Draw a continuous coordinate from a uniform interval.
double uniform(std::mt19937& rng, double lo, double hi) {
  std::uniform_real_distribution<double> dist(lo, hi);
  return dist(rng);
}

// Sample one gene according to a multinomial expression profile.
std::string sample_gene(
    std::mt19937& rng,
    const std::vector<std::string>& genes,
    const std::vector<double>& probabilities) {
  std::discrete_distribution<int> dist(probabilities.begin(), probabilities.end());
  return genes[static_cast<std::size_t>(dist(rng))];
}

// Add one simulated cell and all of its transcripts to the synthetic dataset.
void add_cell_transcripts(
    TranscriptTable& table,
    CellTable& cells,
    std::vector<bool>& true_admixture,
    std::mt19937& rng,
    const std::string& cell_id,
    const std::string& cell_type,
    const Rect& rect,
    int n_transcripts,
    const std::vector<std::string>& genes,
    const std::vector<double>& probabilities,
    bool mark_admixture = false) {
  const auto cx = 0.5 * (rect.x0 + rect.x1);
  const auto cy = 0.5 * (rect.y0 + rect.y1);

  cells.cell_ids.push_back(cell_id);
  cells.cell_types.push_back(cell_type);
  cells.centroid_x.push_back(cx);
  cells.centroid_y.push_back(cy);
  cells.centroid_z.push_back(0.0);

  for (int i = 0; i < n_transcripts; ++i) {
    table.gene_key.push_back(sample_gene(rng, genes, probabilities));
    table.orig_cell_id.push_back(cell_id);
    table.cell_types.push_back(cell_type);
    table.x.push_back(uniform(rng, rect.x0, rect.x1));
    table.y.push_back(uniform(rng, rect.y0, rect.y1));
    table.z.push_back(0.0);
    table.transcript_ids.push_back(cell_id + "_tx_" + std::to_string(i));
    true_admixture.push_back(mark_admixture);
  }
}

}  // namespace

// Build a small admixture toy dataset with malignant and fibroblast compartments.
SimulatedDataset simulate_nsclc_admixture(const SimulationParams& params) {
  SimulatedDataset out;
  out.malignant_markers = {"KRT19", "KRT8", "KRT17"};
  out.fibro_markers = {"COL1A1", "COL3A1", "DCN"};
  const std::vector<std::string> genes = {
      "KRT19", "KRT8", "KRT17", "COL1A1", "COL3A1", "DCN", "ACTB", "GAPDH"};

  const std::vector<double> malignant_probs = {20, 18, 16, 1, 1, 1, 6, 6};
  const std::vector<double> fibro_probs = {1, 1, 1, 20, 18, 16, 6, 6};

  std::mt19937 rng(params.seed);

  const Rect malignant_rect{10.0, 40.0, 10.0, 50.0};
  const Rect fibro_rect{45.0, 75.0, 10.0, 50.0};
  const Rect malignant_rect_2{10.0, 40.0, 55.0, 95.0};
  const Rect fibro_rect_2{45.0, 75.0, 55.0, 95.0};

  add_cell_transcripts(
      out.transcripts,
      out.cells,
      out.true_admixture,
      rng,
      "malignant_1",
      "malignant",
      malignant_rect,
      params.transcripts_per_cell,
      genes,
      malignant_probs);
  add_cell_transcripts(
      out.transcripts,
      out.cells,
      out.true_admixture,
      rng,
      "fibroblast_1",
      "fibroblast",
      fibro_rect,
      params.transcripts_per_cell,
      genes,
      fibro_probs);
  add_cell_transcripts(
      out.transcripts,
      out.cells,
      out.true_admixture,
      rng,
      "malignant_2",
      "malignant",
      malignant_rect_2,
      params.transcripts_per_cell,
      genes,
      malignant_probs);
  add_cell_transcripts(
      out.transcripts,
      out.cells,
      out.true_admixture,
      rng,
      "fibroblast_2",
      "fibroblast",
      fibro_rect_2,
      params.transcripts_per_cell,
      genes,
      fibro_probs);

  for (int i = 0; i < params.admixture_per_target_cell; ++i) {
    out.transcripts.gene_key.push_back(sample_gene(rng, genes, malignant_probs));
    out.transcripts.orig_cell_id.push_back("fibroblast_1");
    out.transcripts.cell_types.push_back("fibroblast");
    out.transcripts.x.push_back(uniform(rng, 44.5, 46.5));
    out.transcripts.y.push_back(uniform(rng, 12.0, 48.0));
    out.transcripts.z.push_back(0.0);
    out.transcripts.transcript_ids.push_back("fibroblast_1_admix_" + std::to_string(i));
    out.true_admixture.push_back(true);
  }

  for (int i = 0; i < params.admixture_per_target_cell; ++i) {
    out.transcripts.gene_key.push_back(sample_gene(rng, genes, malignant_probs));
    out.transcripts.orig_cell_id.push_back("fibroblast_2");
    out.transcripts.cell_types.push_back("fibroblast");
    out.transcripts.x.push_back(uniform(rng, 44.5, 46.5));
    out.transcripts.y.push_back(uniform(rng, 57.0, 93.0));
    out.transcripts.z.push_back(0.0);
    out.transcripts.transcript_ids.push_back("fibroblast_2_admix_" + std::to_string(i));
    out.true_admixture.push_back(true);
  }

  out.transcripts.finalize();

  out.malignant_cell = 0;
  out.fibro_cell = 1;

  out.membrane_image.width = 100;
  out.membrane_image.height = 110;
  out.membrane_image.pixel_size = 1.0;
  out.membrane_image.values.assign(
      static_cast<std::size_t>(out.membrane_image.width * out.membrane_image.height),
      0.05);
  for (int y = 0; y < out.membrane_image.height; ++y) {
    out.membrane_image.values[static_cast<std::size_t>(y * out.membrane_image.width + 48)] = 1.0;
    out.membrane_image.values[static_cast<std::size_t>(y * out.membrane_image.width + 49)] = 1.0;
  }

  return out;
}

}  // namespace celladmix
