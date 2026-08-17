// Persisted run layout, manifest metadata, and parquet read/write helpers.

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "celladmix/clustering.hpp"
#include "celladmix/pipeline.hpp"
#include "celladmix/types.hpp"
#include "celladmix/xenium.hpp"

namespace celladmix {

struct RunSourceFile {
  std::string role;
  std::string path;
  bool exists = false;
};

// Source-level metadata describing the raw input used for a run.
struct RunSourceInfo {
  std::string type;
  std::string path;
  bool used_parquet = false;
  std::vector<RunSourceFile> files;
};

// Canonical file layout for a persisted run directory.
struct RunPaths {
  std::string root_dir;
  std::string run_json;
  std::string factors_parquet;
  std::string cells_parquet;
  std::string molecules_parquet;
  std::string training_rows_parquet;
  std::string scores_dir;
  std::string corrected_dir;
};

// Storage controls for persisted parquet outputs.
struct RunStorageOptions {
  double tile_size = 100.0;
  int parquet_row_group_size = 65536;
};

// Compact NMF restart diagnostics stored with persisted runs. Stability
// values are matched gene-ownership correlations (see nmf_stability.hpp);
// stability_metric records which definition produced them.
struct NmfRunDiagnostics {
  double final_objective = 0.0;
  unsigned int selected_seed = 1;
  int selected_run = 0;
  std::vector<double> candidate_final_objectives;
  std::vector<double> candidate_best_match_correlations;
  std::vector<double> selected_factor_stability;
  double candidate_final_objective_mean = 0.0;
  double candidate_final_objective_sd = 0.0;
  double candidate_best_match_correlation_mean = 0.0;
  std::string stability_metric = "ownership_matched";
  int stability_comparison_runs = 0;
  int stable_factor_count = 0;
  double stability_threshold = 0.0;
};

struct SparseNmfResult;

// Copy multirun diagnostics from a fit result into manifest form.
NmfRunDiagnostics nmf_diagnostics_from_fit(const SparseNmfResult& fit);

// Top-level manifest describing one persisted run.
struct RunManifest {
  std::string format_version = "0.2";
  std::string package_version;
  std::string annotation_hash;
  std::string run_type = "fit";
  RunSourceInfo source;
  RunPaths paths;
  BasicPipelineOptions pipeline_options;
  RunStorageOptions storage_options;
  std::optional<std::string> analysis_crop;
  std::optional<std::string> parent_run;
  std::vector<std::string> genes;
  std::vector<double> nmf_gene_weights;
  NmfRunDiagnostics nmf_diagnostics;
  std::vector<std::string> crop_ids;
  double nmf_transform_target_row_sum = 0.0;
  std::size_t n_transcripts = 0;
  std::size_t n_cells = 0;
  std::size_t n_factors = 0;
  std::size_t n_training_rows = 0;
  bool has_z = false;
  bool has_qv = false;
  bool has_transcript_id = false;
  bool has_cell_type = false;
  bool has_sample_id = false;
  bool has_fov_id = false;
  bool has_nucleus_id = false;
  bool has_overlaps_nucleus = false;
  bool has_nucleus_distance = false;
};

// Filtering options for loading subsets back from persisted run artifacts.
struct RunLoadOptions {
  std::optional<std::string> crop_id;
  std::optional<CropBox> region;
  int sample_n = -1;
  unsigned int seed = 1U;
  // When >= 0, molecule labels are replaced by this ensemble member's
  // labels (see ensure_ensemble_labels), matched per molecule by obs id.
  int ensemble_member = -1;
};

// In-memory representation of a run reloaded from parquet outputs.
struct RunData {
  RunManifest manifest;
  CellTable cells;
  TranscriptTable transcripts;
  std::vector<std::int64_t> obs_ids;
  std::vector<int> labels;
  std::vector<double> factor_margin;
  std::vector<std::string> transcript_crop_ids;
};

// Minimal fitted-run data needed by bridge scoring. This avoids loading gene
// strings, factor margins, and other molecule columns that bridge tests do not
// use.
struct BridgeRunData {
  RunManifest manifest;
  CellTable cells;
  TranscriptTable transcripts;
  std::vector<int> labels;
};

// Minimal run subset needed to rebuild sampled-training NCV report sidecars.
struct TrainingRunData {
  RunManifest manifest;
  CellTable cells;
  TranscriptTable transcripts;
  std::vector<std::int64_t> obs_ids;
  std::vector<int> labels;
  std::vector<double> factor_margin;
  std::vector<int> training_query_indices;
  std::vector<std::int64_t> training_obs_ids;
};

// Derive the standard file layout for a run rooted at one directory.
RunPaths make_run_paths(const std::string& root_dir);

// Return the per-cell correction summary path for a corrected run.
std::string correction_summary_parquet_path(const std::string& run_path_or_dir);

// Persist one basic fit run to disk and return its manifest.
RunManifest write_basic_run(
    const std::string& root_dir,
    const RunSourceInfo& source,
    const TranscriptTable& table,
    const BasicPipelineResult& fit,
    const BasicPipelineOptions& pipeline_options,
    const RunStorageOptions& storage_options = {},
    const std::optional<std::string>& analysis_crop = std::nullopt,
    const std::vector<std::string>* transcript_crop_ids = nullptr,
    const CellTable* source_cells = nullptr);

// Serialize a run manifest to its JSON metadata file.
void write_manifest_json(const RunManifest& manifest);

// Read run metadata from a run directory or its run.json path.
RunManifest read_run_manifest(const std::string& path_or_dir);

// Load the persisted cell summary table for a run.
CellTable load_run_cells(const std::string& path_or_dir);

// Load the persisted list of observation ids used for NMF training rows.
std::vector<std::int64_t> load_run_training_obs_ids(const std::string& path_or_dir);

// Reload persisted transcript, cell, and label data from a run.
RunData load_run_data(
    const std::string& path_or_dir,
    const RunLoadOptions& options = {});

// Aggregate persisted run molecules into sparse cell-by-gene counts.
CellCountMatrix collect_run_counts(
    const std::string& path_or_dir,
    const RunLoadOptions& options = {});

// Load only numeric molecule columns required for bridge scoring.
BridgeRunData load_run_bridge_data(
    const std::string& path_or_dir,
    const RunLoadOptions& options = {});

// Load only molecules from cells containing sampled training observations.
TrainingRunData load_run_training_cell_data(const std::string& path_or_dir);

// Persist a corrected child run derived from a parent run and keep mask.
RunManifest write_corrected_run(
    const std::string& root_dir,
    const std::string& parent_run_path_or_dir,
    const RunData& run_data,
    const std::vector<bool>& keep_mask);

// --- Ensemble member pool (per-restart factorizations and labelings) ---

// Path of the persisted restart-H pool inside a run directory.
std::string ensemble_h_parquet_path(const std::string& run_path_or_dir);
// Path of one member's persisted molecule labels.
std::string ensemble_labels_parquet_path(
    const std::string& run_path_or_dir,
    int member);

// Persist the restart member pool as full-gene-space H matrices.
void write_ensemble_h_parquet(
    const std::string& path,
    const std::vector<DenseMatrix>& candidate_h,
    int row_group_size);

// Load the restart member pool; empty when the run has none persisted.
std::vector<DenseMatrix> load_ensemble_h(const std::string& run_path_or_dir);

// Number of members whose molecule labels are already persisted.
int ensemble_member_count(const std::string& run_path_or_dir);

// Compute (or reuse) molecule labels for every restart in the member pool
// and persist them aligned row-for-row with molecules.parquet. Returns the
// member count (0 when the run carries no member pool).
int ensure_ensemble_labels(
    const std::string& run_path_or_dir,
    int num_threads);

// One member's molecule labels keyed by obs id (ascending, for lookup).
struct EnsembleMemberLabels {
  std::vector<std::int64_t> obs_ids;
  std::vector<int> labels;
  int label_for(std::int64_t obs_id) const;
};
EnsembleMemberLabels load_ensemble_member_labels(
    const std::string& run_path_or_dir,
    int member);

// Cell-by-factor molecule fractions under one member's labeling; rows
// follow cells.parquet order, columns are factors.
DenseMatrix ensemble_member_cell_fractions(
    const std::string& run_path_or_dir,
    int member);

}  // namespace celladmix
