#include "celladmix/tabular.hpp"
#include "celladmix/tabular_mask.hpp"

#include <arrow/api.h>
#include <arrow/compute/api.h>
#include <arrow/csv/api.h>
#include <arrow/io/api.h>
#include <arrow/result.h>
#include <arrow/util/bit_util.h>
#include <parquet/arrow/reader.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace celladmix {
namespace {

// Arrow- and TIFF-backed helpers for generic CSV/Parquet molecule loading.

#define CELLADMIX_ARROW_CHECK_OK(expr)                                                \
  do {                                                                                \
    auto _status = (expr);                                                            \
    if (!_status.ok()) {                                                              \
      throw std::runtime_error(std::string("Arrow error: ") + _status.ToString());    \
    }                                                                                 \
  } while (0)

template <typename T>
T arrow_unwrap(arrow::Result<T>&& result) {
  if (!result.ok()) {
    throw std::runtime_error(std::string("Arrow error: ") + result.status().ToString());
  }
  return std::move(*result);
}

std::string file_extension(const std::string& path) {
  const auto dot = path.rfind('.');
  if (dot == std::string::npos) {
    return {};
  }
  std::string ext = path.substr(dot + 1);
  std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char ch) {
    return static_cast<char>(std::tolower(ch));
  });
  return ext;
}

std::shared_ptr<arrow::Table> read_csv_arrow(const std::string& path) {
  auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path));
  auto read_options = arrow::csv::ReadOptions::Defaults();
  auto parse_options = arrow::csv::ParseOptions::Defaults();
  auto convert_options = arrow::csv::ConvertOptions::Defaults();
  if (file_extension(path) == "tsv") {
    parse_options.delimiter = '\t';
  }
  auto reader = arrow_unwrap(
      arrow::csv::TableReader::Make(
          arrow::io::default_io_context(),
          input,
          read_options,
          parse_options,
          convert_options));
  return arrow_unwrap(reader->Read());
}

std::shared_ptr<arrow::Table> read_parquet_arrow(const std::string& path) {
  auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path));
  parquet::arrow::FileReaderBuilder builder;
  CELLADMIX_ARROW_CHECK_OK(builder.Open(input));
  auto reader = arrow_unwrap(builder.Build());
  std::shared_ptr<arrow::Table> table;
  CELLADMIX_ARROW_CHECK_OK(reader->ReadTable(&table));
  return table;
}

std::shared_ptr<arrow::Table> read_table(const std::string& path) {
  const std::string ext = file_extension(path);
  if (ext == "csv" || ext == "tsv") {
    return read_csv_arrow(path);
  }
  if (ext == "parquet" || ext == "pq") {
    return read_parquet_arrow(path);
  }
  throw std::runtime_error("Unsupported tabular file format for '" + path + "'");
}

int find_column_index(const std::shared_ptr<arrow::Table>& table, const std::string& name) {
  const int idx = table->schema()->GetFieldIndex(name);
  if (idx < 0) {
    throw std::runtime_error("Column '" + name + "' not found in tabular input");
  }
  return idx;
}

bool has_column(const std::shared_ptr<arrow::Table>& table, const std::string& name) {
  return table->schema()->GetFieldIndex(name) >= 0;
}

std::vector<double> extract_double_column(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& name) {
  const int idx = find_column_index(table, name);
  const auto chunked = table->column(idx);
  std::vector<double> out(static_cast<std::size_t>(chunked->length()));
  int64_t offset = 0;
  for (int c = 0; c < chunked->num_chunks(); ++c) {
    const auto chunk = chunked->chunk(c);
    const int64_t len = chunk->length();
    if (chunk->type_id() == arrow::Type::DOUBLE) {
      const auto arr = std::static_pointer_cast<arrow::DoubleArray>(chunk);
      for (int64_t i = 0; i < len; ++i) out[static_cast<std::size_t>(offset + i)] = arr->Value(i);
    } else if (chunk->type_id() == arrow::Type::FLOAT) {
      const auto arr = std::static_pointer_cast<arrow::FloatArray>(chunk);
      for (int64_t i = 0; i < len; ++i) out[static_cast<std::size_t>(offset + i)] = static_cast<double>(arr->Value(i));
    } else if (chunk->type_id() == arrow::Type::INT64) {
      const auto arr = std::static_pointer_cast<arrow::Int64Array>(chunk);
      for (int64_t i = 0; i < len; ++i) out[static_cast<std::size_t>(offset + i)] = static_cast<double>(arr->Value(i));
    } else if (chunk->type_id() == arrow::Type::INT32) {
      const auto arr = std::static_pointer_cast<arrow::Int32Array>(chunk);
      for (int64_t i = 0; i < len; ++i) out[static_cast<std::size_t>(offset + i)] = static_cast<double>(arr->Value(i));
    } else {
      auto casted = arrow_unwrap(arrow::compute::Cast(*chunk, arrow::float64()));
      const auto arr = std::static_pointer_cast<arrow::DoubleArray>(casted);
      for (int64_t i = 0; i < len; ++i) out[static_cast<std::size_t>(offset + i)] = arr->Value(i);
    }
    offset += len;
  }
  return out;
}

std::vector<std::string> extract_string_column(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& name) {
  const int idx = find_column_index(table, name);
  const auto chunked = table->column(idx);
  std::vector<std::string> out(static_cast<std::size_t>(chunked->length()));
  int64_t offset = 0;
  for (int c = 0; c < chunked->num_chunks(); ++c) {
    const auto chunk = chunked->chunk(c);
    const int64_t len = chunk->length();
    if (chunk->type_id() == arrow::Type::STRING) {
      const auto arr = std::static_pointer_cast<arrow::StringArray>(chunk);
      for (int64_t i = 0; i < len; ++i) out[static_cast<std::size_t>(offset + i)] = arr->GetString(i);
    } else if (chunk->type_id() == arrow::Type::LARGE_STRING) {
      const auto arr = std::static_pointer_cast<arrow::LargeStringArray>(chunk);
      for (int64_t i = 0; i < len; ++i) out[static_cast<std::size_t>(offset + i)] = arr->GetString(i);
    } else if (chunk->type_id() == arrow::Type::DICTIONARY) {
      const auto dict_arr = std::static_pointer_cast<arrow::DictionaryArray>(chunk);
      const auto dict = dict_arr->dictionary();
      if (dict->type_id() == arrow::Type::STRING) {
        const auto dict_values = std::static_pointer_cast<arrow::StringArray>(dict);
        auto indices_cast = arrow_unwrap(
            arrow::compute::Cast(*dict_arr->indices(), arrow::int64()));
        const auto idx_values = std::static_pointer_cast<arrow::Int64Array>(indices_cast);
        for (int64_t i = 0; i < len; ++i) {
          if (dict_arr->IsNull(i)) {
            out[static_cast<std::size_t>(offset + i)] = "";
          } else {
            out[static_cast<std::size_t>(offset + i)] =
                dict_values->GetString(static_cast<int64_t>(idx_values->Value(i)));
          }
        }
      } else {
        for (int64_t i = 0; i < len; ++i) {
          out[static_cast<std::size_t>(offset + i)] = arrow_unwrap(chunk->GetScalar(i))->ToString();
        }
      }
    } else {
      for (int64_t i = 0; i < len; ++i) {
        out[static_cast<std::size_t>(offset + i)] = arrow_unwrap(chunk->GetScalar(i))->ToString();
      }
    }
    offset += len;
  }
  return out;
}

struct UnionBounds {
  bool active = false;
  double xmin = 0.0;
  double xmax = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;
  double zmin = 0.0;
  double zmax = 0.0;
  bool has_z = false;
};

UnionBounds union_bounds(const std::vector<CropBox>& crops) {
  UnionBounds out;
  if (crops.empty()) {
    return out;
  }
  out.active = true;
  out.xmin = crops.front().xmin;
  out.xmax = crops.front().xmax;
  out.ymin = crops.front().ymin;
  out.ymax = crops.front().ymax;
  out.zmin = crops.front().zmin;
  out.zmax = crops.front().zmax;
  out.has_z = crops.front().has_z;
  for (std::size_t i = 1; i < crops.size(); ++i) {
    out.xmin = std::min(out.xmin, crops[i].xmin);
    out.xmax = std::max(out.xmax, crops[i].xmax);
    out.ymin = std::min(out.ymin, crops[i].ymin);
    out.ymax = std::max(out.ymax, crops[i].ymax);
    out.zmin = std::min(out.zmin, crops[i].zmin);
    out.zmax = std::max(out.zmax, crops[i].zmax);
    out.has_z = out.has_z || crops[i].has_z;
  }
  return out;
}

std::string assign_crop_id(double x, double y, double z, const std::vector<CropBox>& crops) {
  if (crops.empty()) {
    return {};
  }
  for (const auto& crop : crops) {
    if (x < crop.xmin || x > crop.xmax || y < crop.ymin || y > crop.ymax) {
      continue;
    }
    if (crop.has_z && (z < crop.zmin || z > crop.zmax)) {
      continue;
    }
    return crop.crop_id;
  }
  return {};
}

bool passes_union_bounds(double x, double y, double z, const UnionBounds& bounds) {
  if (!std::isfinite(x) || !std::isfinite(y)) {
    return false;
  }
  if (!bounds.active) {
    return true;
  }
  if (x < bounds.xmin || x > bounds.xmax || y < bounds.ymin || y > bounds.ymax) {
    return false;
  }
  if (bounds.has_z && (!std::isfinite(z) || z < bounds.zmin || z > bounds.zmax)) {
    return false;
  }
  return true;
}

bool is_unassigned_cell_id(const std::string& value) {
  return value.empty() || value == "0" || value == "NA" || value == "NaN" || value == "null";
}

}  // namespace

TabularBundleData load_tabular_bundle(
    const TabularSourceSpec& source,
    const TabularLoadOptions& options) {
  if (source.molecules_path.empty()) {
    throw std::runtime_error("molecules_path must be provided for tabular input");
  }
  if (source.gene_col.empty() || source.x_col.empty() || source.y_col.empty()) {
    throw std::runtime_error("x_col, y_col, and gene_col must be provided for tabular input");
  }
  if (!source.cell_id_col.empty() && !source.segmentation_mask_path.empty()) {
    throw std::runtime_error("Specify either cell_id_col or segmentation_mask_path, not both");
  }
  if (source.cell_id_col.empty() && source.segmentation_mask_path.empty()) {
    throw std::runtime_error("Tabular input requires either cell_id_col or segmentation_mask_path");
  }

  const auto table_in = read_table(source.molecules_path);
  const bool used_parquet = file_extension(source.molecules_path) == "parquet" ||
                            file_extension(source.molecules_path) == "pq";

  std::vector<double> x = extract_double_column(table_in, source.x_col);
  std::vector<double> y = extract_double_column(table_in, source.y_col);
  std::vector<double> z;
  if (!source.z_col.empty() && has_column(table_in, source.z_col)) {
    z = extract_double_column(table_in, source.z_col);
  }
  if (z.empty()) {
    z.assign(x.size(), 0.0);
  }
  std::vector<std::string> gene = extract_string_column(table_in, source.gene_col);
  std::vector<double> qv;
  if (!source.qv_col.empty() && has_column(table_in, source.qv_col)) {
    qv = extract_double_column(table_in, source.qv_col);
  }
  std::vector<std::string> cell_id;
  if (!source.cell_id_col.empty()) {
    cell_id = extract_string_column(table_in, source.cell_id_col);
  } else {
    const auto mask = read_labeled_tiff_mask(source.segmentation_mask_path);
    cell_id = assign_cell_ids_from_mask(x, y, mask);
  }

  std::vector<std::string> cell_type;
  if (!source.cell_type_col.empty() && has_column(table_in, source.cell_type_col)) {
    cell_type = extract_string_column(table_in, source.cell_type_col);
  }
  std::vector<std::string> sample_id;
  if (!source.sample_id_col.empty() && has_column(table_in, source.sample_id_col)) {
    sample_id = extract_string_column(table_in, source.sample_id_col);
  } else if (!source.sample_id.empty()) {
    sample_id.assign(x.size(), source.sample_id);
  }
  std::vector<std::string> fov_id;
  if (!source.fov_id_col.empty() && has_column(table_in, source.fov_id_col)) {
    fov_id = extract_string_column(table_in, source.fov_id_col);
  } else if (!source.fov_id.empty()) {
    fov_id.assign(x.size(), source.fov_id);
  }
  std::vector<std::string> transcript_id;
  if (has_column(table_in, "transcript_id")) {
    transcript_id = extract_string_column(table_in, "transcript_id");
  }

  const UnionBounds bounds = union_bounds(options.crops);
  std::size_t n_unassigned_rows = 0;
  std::unordered_set<std::string> assigned_cell_ids;
  assigned_cell_ids.reserve(cell_id.size() / 16 + 1);
  TabularBundleData out;
  out.used_parquet = used_parquet;
  out.transcripts.gene_key.reserve(x.size());
  out.transcripts.orig_cell_id.reserve(x.size());
  out.transcripts.x.reserve(x.size());
  out.transcripts.y.reserve(x.size());
  out.transcripts.z.reserve(x.size());
  out.transcript_crop_ids.reserve(x.size());
  if (!qv.empty()) out.transcripts.qv.reserve(x.size());
  if (!cell_type.empty()) out.transcripts.cell_types.reserve(x.size());
  if (!sample_id.empty()) out.transcripts.sample_ids.reserve(x.size());
  if (!fov_id.empty()) out.transcripts.fov_ids.reserve(x.size());
  if (!transcript_id.empty()) out.transcripts.transcript_ids.reserve(x.size());

  for (std::size_t i = 0; i < x.size(); ++i) {
    if (!passes_union_bounds(x[i], y[i], z[i], bounds)) {
      continue;
    }
    if (!qv.empty() && options.min_qv >= 0.0 &&
        (!std::isfinite(qv[i]) || qv[i] < options.min_qv)) {
      continue;
    }
    const std::string cell = cell_id[i];
    if (is_unassigned_cell_id(cell)) {
      ++n_unassigned_rows;
      if (!options.keep_unassigned) {
        continue;
      }
    }
    const std::string crop_id = assign_crop_id(x[i], y[i], z[i], options.crops);
    if (!options.crops.empty() && crop_id.empty()) {
      continue;
    }

    out.transcripts.gene_key.push_back(gene[i]);
    out.transcripts.orig_cell_id.push_back(is_unassigned_cell_id(cell) ? "__unassigned__" : cell);
    if (!is_unassigned_cell_id(cell)) {
      assigned_cell_ids.insert(cell);
    }
    out.transcripts.x.push_back(x[i]);
    out.transcripts.y.push_back(y[i]);
    out.transcripts.z.push_back(z[i]);
    out.transcript_crop_ids.push_back(crop_id);
    if (!qv.empty()) out.transcripts.qv.push_back(qv[i]);
    if (!cell_type.empty()) out.transcripts.cell_types.push_back(cell_type[i]);
    if (!sample_id.empty()) out.transcripts.sample_ids.push_back(sample_id[i]);
    if (!fov_id.empty()) out.transcripts.fov_ids.push_back(fov_id[i]);
    if (!transcript_id.empty()) out.transcripts.transcript_ids.push_back(transcript_id[i]);
  }

  out.transcripts.finalize();
  const std::size_t total_rows = x.size();
  const std::size_t assigned_rows = total_rows >= n_unassigned_rows ? (total_rows - n_unassigned_rows) : 0;
  const double assigned_fraction = total_rows == 0
      ? 0.0
      : static_cast<double>(assigned_rows) / static_cast<double>(total_rows);
  std::cerr << "[INFO] Loaded tabular source: " << total_rows
            << " rows total, " << assigned_rows
            << " assigned to " << assigned_cell_ids.size()
            << " cells (" << std::round(assigned_fraction * 1000.0) / 10.0
            << "%), " << out.transcripts.size()
            << " kept after filters\n";
  return out;
}

}  // namespace celladmix
