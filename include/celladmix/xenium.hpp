// Xenium manifest parsing and transcript/cell loading utilities.

#pragma once

#include <algorithm>
#include <cctype>
#include <functional>
#include <string>
#include <vector>

#include "celladmix/types.hpp"

namespace celladmix {

struct CropBox {
  std::string crop_id;
  double xmin = 0.0;
  double xmax = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;
  double zmin = 0.0;
  double zmax = 0.0;
  bool has_z = false;
};

// Subset of Xenium manifest fields needed by the core loader.
struct XeniumManifest {
  std::string run_name;
  double pixel_size = 1.0;
  double z_step_size = 0.0;
  std::string morphology_filepath;
  std::string morphology_focus_filepath;
  std::string transcripts_csv_path = "transcripts.csv.gz";
  std::string transcripts_parquet_path = "transcripts.parquet";
};

// Controls for crop-aware Xenium loading and filtering.
struct XeniumLoadOptions {
  std::vector<CropBox> crops;
  std::vector<std::string> cell_filter;
  std::vector<std::string> gene_filter;
  double min_qv = -1.0;
  bool keep_unassigned = false;
  bool keep_non_gene = false;
  bool read_cells = true;
  bool prefer_parquet = true;
  std::function<void(const std::string&, double)> progress;
};

// Raw Xenium data bundle loaded into transcript and cell tables.
struct XeniumBundleData {
  XeniumManifest manifest;
  TranscriptTable transcripts;
  CellTable cells;
  std::vector<std::string> transcript_crop_ids;
  std::vector<std::string> cell_crop_ids;
  bool used_parquet = false;
};

// Parse a Xenium experiment manifest from disk.
XeniumManifest read_xenium_manifest(const std::string& path);
// Load a Xenium bundle directory or manifest into in-memory tables.
XeniumBundleData load_xenium_bundle(
    const std::string& manifest_path_or_dir,
    const XeniumLoadOptions& options = {});

// Classify Xenium decoded features for biological-expression workflows.
// Modern Xenium parquet files may expose `is_gene` and `codeword_category`;
// older bundles can still be filtered by the stable 10x control feature names.
inline std::string xenium_lower_ascii(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
    return static_cast<char>(std::tolower(ch));
  });
  return value;
}

inline bool xenium_ascii_starts_with(const std::string& value, const char* prefix) {
  std::size_t i = 0;
  for (; prefix[i] != '\0'; ++i) {
    if (i >= value.size()) {
      return false;
    }
    const auto lhs = static_cast<unsigned char>(value[i]);
    const auto rhs = static_cast<unsigned char>(prefix[i]);
    if (std::tolower(lhs) != std::tolower(rhs)) {
      return false;
    }
  }
  return true;
}

inline bool xenium_parse_bool(const std::string& value, bool default_value) {
  const std::string lower = xenium_lower_ascii(value);
  if (lower == "true" || lower == "t" || lower == "1" || lower == "yes" || lower == "y") {
    return true;
  }
  if (lower == "false" || lower == "f" || lower == "0" || lower == "no" || lower == "n") {
    return false;
  }
  return default_value;
}

inline bool xenium_feature_name_is_control(const std::string& feature_name) {
  return xenium_ascii_starts_with(feature_name, "DeprecatedCodeword") ||
         xenium_ascii_starts_with(feature_name, "UnassignedCodeword") ||
         xenium_ascii_starts_with(feature_name, "NegControl") ||
         xenium_ascii_starts_with(feature_name, "NegativeControl") ||
         xenium_ascii_starts_with(feature_name, "Negative_Control") ||
         xenium_ascii_starts_with(feature_name, "GenomicControl") ||
         xenium_ascii_starts_with(feature_name, "Genomic_Control") ||
         xenium_ascii_starts_with(feature_name, "Intergenic_Region") ||
         xenium_ascii_starts_with(feature_name, "Blank");
}

inline bool xenium_category_is_biological_gene(const std::string& codeword_category) {
  if (codeword_category.empty()) {
    return true;
  }
  const std::string lower = xenium_lower_ascii(codeword_category);
  return lower == "gene" ||
         lower == "predesigned_gene" ||
         lower == "custom_gene" ||
         lower == "gene_expression" ||
         lower == "gene expression";
}

inline bool xenium_is_biological_feature(
    const std::string& feature_name,
    const std::string& codeword_category = {},
    const std::string& is_gene = {}) {
  if (feature_name.empty()) {
    return false;
  }
  if (!codeword_category.empty() && !xenium_category_is_biological_gene(codeword_category)) {
    return false;
  }
  if (!is_gene.empty() && !xenium_parse_bool(is_gene, true)) {
    return false;
  }
  return !xenium_feature_name_is_control(feature_name);
}

}  // namespace celladmix
