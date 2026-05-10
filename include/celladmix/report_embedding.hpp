// Developer/report-only helpers for building sampled training-molecule UMAP
// sidecars from persisted runs.

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "celladmix/pipeline.hpp"
#include "celladmix/types.hpp"

namespace celladmix {

struct TrainingMoleculeUmapData {
  std::vector<int> training_rank;
  std::vector<std::int64_t> obs_id;
  std::vector<std::string> cell_id;
  std::vector<std::string> gene;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<int> factor_label;
  std::vector<double> factor_margin;
  std::vector<int> dominant_component;
  std::vector<double> dominant_component_margin;
  std::vector<int> projected_component;
  std::vector<double> projected_component_margin;
  std::vector<double> umap_x;
  std::vector<double> umap_y;
  DenseMatrix component_weights;
  DenseMatrix projected_component_weights;
};

// Resolve the canonical parquet path for the training-molecule report sidecar.
std::string training_molecule_umap_parquet_path(const std::string& run_path_or_dir);

// Build and write the sampled training-molecule UMAP parquet sidecar.
void write_training_molecule_umap_report(
    const std::string& run_root_dir,
    const TranscriptTable& table,
    const BasicPipelineResult& fit,
    int ncv_k,
    int umap_neighbors = 15,
    int umap_epochs = 200,
    unsigned int seed = 1U,
    double normalization_scale = 5000.0,
    int pca_dims = 30,
    int num_threads = 1,
    const std::vector<std::int64_t>* obs_id_override = nullptr);

// Load the sampled training-molecule UMAP sidecar into memory.
TrainingMoleculeUmapData load_training_molecule_umap_report(const std::string& run_path_or_dir);

}  // namespace celladmix
