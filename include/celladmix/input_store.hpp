// Normalized source-store access for clustering and fit entrypoints.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "celladmix/clustering.hpp"
#include "celladmix/tabular.hpp"
#include "celladmix/types.hpp"
#include "celladmix/xenium.hpp"

namespace celladmix {

struct InputStoreManifest {
  std::string format_version = "0.2";
  std::string root_dir;
  std::string source_type;
  std::string store_mode;
  std::string source_path;
  std::string source_fingerprint;
  std::string filter_signature;
  bool has_molecules = false;
  bool has_molecule_rows = false;
  bool has_cell_gene_counts = false;
  bool has_cell_offsets = false;
  bool has_spatial_tiles = false;
  bool has_labels = false;
  bool has_factor_scores = false;
  std::size_t n_molecules = 0;
  std::size_t n_cells = 0;
  std::size_t n_genes = 0;
  bool has_z = true;
  bool has_qv = false;
  bool has_transcript_id = false;
  bool has_cell_type = false;
  bool has_sample_id = false;
  bool has_fov_id = false;
  bool has_overlaps_nucleus = false;
  bool has_nucleus_distance = false;
};

struct InputStorePaths {
  std::string root_dir;
  std::string manifest_json;
  std::string genes_parquet;
  std::string cells_parquet;
  std::string counts_parquet;
  std::string molecules_parquet;
  std::string cell_offsets_parquet;
};

struct InputStoreBuildOptions {
  std::string store_dir;
  bool materialize_molecules = true;
  bool force = false;
  int num_threads = 1;
  int parquet_row_group_size = 65536;
};

struct InputStoreData {
  InputStoreManifest manifest;
  TranscriptTable transcripts;
  CellTable cells;
  std::vector<std::string> transcript_crop_ids;
};

struct InputStoreCellOffset {
  int cell_idx = 0;
  std::int64_t start = -1;
  std::int64_t end = -1;
};

struct InputStoreMoleculeBlock {
  std::vector<std::int64_t> row_index;
  std::vector<std::int64_t> obs_id;
  std::vector<int> cell_idx;
  std::vector<int> gene_idx;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> qv;
  std::vector<int> overlaps_nucleus;
  std::vector<double> nucleus_distance;
  std::vector<std::string> crop_id;

  std::size_t size() const { return x.size(); }
};

InputStorePaths make_input_store_paths(const std::string& root_dir);

InputStoreManifest read_input_store_manifest(const std::string& store_dir);

InputStoreManifest build_xenium_input_store(
    const std::string& manifest_path_or_dir,
    const XeniumLoadOptions& load_options,
    const InputStoreBuildOptions& store_options);

InputStoreManifest build_tabular_input_store(
    const TabularSourceSpec& source,
    const TabularLoadOptions& load_options,
    const InputStoreBuildOptions& store_options);

CellCountMatrix load_input_store_counts(const std::string& store_dir);

std::vector<InputStoreCellOffset> load_input_store_cell_offsets(const std::string& store_dir);

InputStoreMoleculeBlock load_input_store_molecule_range(
    const std::string& store_dir,
    std::int64_t start,
    std::int64_t end);

// Load molecule coordinates and encoded ids without expanding per-row strings.
InputStoreData load_input_store_numeric_data(const std::string& store_dir);

InputStoreData load_input_store_data(const std::string& store_dir);

// Load only the cell table from a normalized input store.
CellTable load_input_store_cells(const std::string& store_dir);

}  // namespace celladmix
