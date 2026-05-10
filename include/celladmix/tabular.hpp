// Generic tabular molecule loading from CSV/Parquet plus optional labeled TIFF
// segmentation assignment.

#pragma once

#include <string>
#include <vector>

#include "celladmix/types.hpp"
#include "celladmix/xenium.hpp"

namespace celladmix {

// Column mapping and sidecar files used to interpret a generic molecules table.
struct TabularSourceSpec {
  std::string molecules_path;
  std::string x_col = "x";
  std::string y_col = "y";
  std::string z_col;
  std::string gene_col = "gene";
  std::string qv_col;
  std::string cell_id_col;
  std::string cell_type_col;
  std::string cell_metadata_path;
  std::string cell_metadata_cell_id_col = "cell";
  std::string cell_metadata_cell_type_col;
  std::string sample_id_col;
  std::string fov_id_col;
  std::string sample_id;
  std::string fov_id;
  std::string segmentation_mask_path;
};

// Crop and QV filters applied while loading generic tabular molecules.
struct TabularLoadOptions {
  std::vector<CropBox> crops;
  double min_qv = -1.0;
  bool keep_unassigned = false;
};

// In-memory tabular dataset after cell assignment and basic filtering.
struct TabularBundleData {
  TranscriptTable transcripts;
  std::vector<std::string> transcript_crop_ids;
  bool used_parquet = false;
};

// Load one CSV/Parquet molecules table into the core transcript representation.
TabularBundleData load_tabular_bundle(
    const TabularSourceSpec& source,
    const TabularLoadOptions& options = {});

}  // namespace celladmix
