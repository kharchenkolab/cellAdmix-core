// Normalized disk-backed input store for source-independent downstream access.

#include "celladmix/input_store.hpp"
#include "celladmix/tabular_mask.hpp"

#include <arrow/api.h>
#include <arrow/csv/api.h>
#include <arrow/io/compressed.h>
#include <arrow/io/api.h>
#include <arrow/util/compression.h>
#include <arrow/util/bit_util.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/file_reader.h>
#include <parquet/statistics.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "celladmix/workflow.hpp"

namespace celladmix {
namespace {

template <class T>
T arrow_unwrap(arrow::Result<T>&& result, const char* context) {
  if (!result.ok()) {
    throw std::runtime_error(std::string(context) + ": " + result.status().ToString());
  }
  return std::move(*result);
}

void arrow_check(const arrow::Status& status, const char* context) {
  if (!status.ok()) {
    throw std::runtime_error(std::string(context) + ": " + status.ToString());
  }
}

void ensure_directory(const std::filesystem::path& path) {
  std::error_code ec;
  std::filesystem::create_directories(path, ec);
  if (ec) {
    throw std::runtime_error("Could not create directory '" + path.string() + "': " + ec.message());
  }
}

std::string normalize_path_for_store(const std::string& path) {
  if (path.empty()) {
    return {};
  }
  std::error_code ec;
  const auto absolute = std::filesystem::absolute(std::filesystem::path(path), ec);
  if (ec) {
    return path;
  }
  const auto canonical = std::filesystem::weakly_canonical(absolute, ec);
  return ec ? absolute.string() : canonical.string();
}

bool file_exists(const std::string& path) {
  return !path.empty() && std::filesystem::exists(std::filesystem::path(path));
}

std::string dirname_of(const std::string& path) {
  return std::filesystem::path(path).parent_path().string();
}

std::string join_path(const std::string& base, const std::string& child) {
  if (child.empty()) {
    return base;
  }
  const auto child_path = std::filesystem::path(child);
  if (child_path.is_absolute()) {
    return child_path.string();
  }
  return (std::filesystem::path(base) / child_path).string();
}

std::string resolve_xenium_manifest_path(const std::string& path_or_dir) {
  const auto path = std::filesystem::path(path_or_dir);
  if (std::filesystem::is_directory(path)) {
    return (path / "experiment.xenium").string();
  }
  return path.string();
}

std::string file_fingerprint(const std::string& path) {
  const auto fs_path = std::filesystem::path(path);
  if (!std::filesystem::exists(fs_path)) {
    return normalize_path_for_store(path) + "|missing";
  }
  std::error_code ec;
  const auto size = std::filesystem::file_size(fs_path, ec);
  const auto size_text = ec ? std::string("size?") : std::to_string(size);
  ec.clear();
  const auto mtime = std::filesystem::last_write_time(fs_path, ec);
  const auto mtime_ticks = std::chrono::duration_cast<std::chrono::nanoseconds>(
      mtime.time_since_epoch()).count();
  const auto mtime_text = ec ? std::string("mtime?") : std::to_string(mtime_ticks);
  return normalize_path_for_store(path) + "|" + size_text + "|" + mtime_text;
}

std::string combined_fingerprint(const std::vector<std::string>& paths) {
  std::ostringstream out;
  for (const auto& path : paths) {
    if (path.empty()) {
      continue;
    }
    out << file_fingerprint(path) << "\n";
  }
  return out.str();
}

std::string file_extension_lower(const std::string& path) {
  std::string ext = std::filesystem::path(path).extension().string();
  if (!ext.empty() && ext.front() == '.') {
    ext.erase(ext.begin());
  }
  std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char ch) {
    return static_cast<char>(std::tolower(ch));
  });
  return ext;
}

bool has_suffix(const std::string& value, const std::string& suffix) {
  return value.size() >= suffix.size() &&
      value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool is_gzip_path(const std::string& path) {
  std::string lower = path;
  std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char ch) {
    return static_cast<char>(std::tolower(ch));
  });
  return has_suffix(lower, ".gz");
}

std::string logical_file_extension_lower(const std::string& path) {
  std::string lower = path;
  std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char ch) {
    return static_cast<char>(std::tolower(ch));
  });
  if (has_suffix(lower, ".gz")) {
    lower.resize(lower.size() - 3);
  }
  return file_extension_lower(lower);
}

bool is_parquet_extension(const std::string& ext) {
  return ext == "parquet" || ext == "pq";
}

bool is_delimited_text_extension(const std::string& ext) {
  return ext == "csv" || ext == "tsv";
}

std::shared_ptr<arrow::io::InputStream> open_delimited_input_stream(
    const std::string& path,
    const char* context) {
  auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path), context);
  if (!is_gzip_path(path)) {
    return input;
  }
  auto codec = arrow_unwrap(
      arrow::util::Codec::Create(arrow::Compression::GZIP),
      "Create gzip codec");
  return arrow_unwrap(
      arrow::io::CompressedInputStream::Make(codec.get(), input),
      "Open gzip-compressed tabular input");
}

std::string crop_signature(const std::vector<CropBox>& crops) {
  nlohmann::json out = nlohmann::json::array();
  for (const auto& crop : crops) {
    out.push_back({
        {"crop_id", crop.crop_id},
        {"xmin", crop.xmin},
        {"xmax", crop.xmax},
        {"ymin", crop.ymin},
        {"ymax", crop.ymax},
        {"zmin", crop.zmin},
        {"zmax", crop.zmax},
        {"has_z", crop.has_z}});
  }
  return out.dump();
}

std::string xenium_filter_signature(const XeniumLoadOptions& options) {
  return nlohmann::json{
      {"crops", nlohmann::json::parse(crop_signature(options.crops))},
      {"cell_filter", options.cell_filter},
      {"gene_filter", options.gene_filter},
      {"min_qv", options.min_qv},
      {"keep_unassigned", options.keep_unassigned},
      {"keep_non_gene", options.keep_non_gene},
      {"read_cells", options.read_cells},
      {"prefer_parquet", options.prefer_parquet}}.dump();
}

std::string tabular_filter_signature(
    const TabularSourceSpec& source,
    const TabularLoadOptions& options) {
  return nlohmann::json{
      {"crops", nlohmann::json::parse(crop_signature(options.crops))},
      {"min_qv", options.min_qv},
      {"keep_unassigned", options.keep_unassigned},
      {"x_col", source.x_col},
      {"y_col", source.y_col},
      {"z_col", source.z_col},
      {"gene_col", source.gene_col},
      {"qv_col", source.qv_col},
      {"cell_id_col", source.cell_id_col},
      {"segmentation_mask_path", normalize_path_for_store(source.segmentation_mask_path)},
      {"cell_type_col", source.cell_type_col},
      {"cell_metadata_path", normalize_path_for_store(source.cell_metadata_path)},
      {"cell_metadata_cell_id_col", source.cell_metadata_cell_id_col},
      {"cell_metadata_cell_type_col", source.cell_metadata_cell_type_col},
      {"sample_id_col", source.sample_id_col},
      {"fov_id_col", source.fov_id_col},
      {"sample_id", source.sample_id},
      {"fov_id", source.fov_id}}.dump();
}

std::unordered_set<std::string> make_string_filter_set(const std::vector<std::string>& values) {
  std::unordered_set<std::string> out;
  out.reserve(values.size());
  for (const auto& value : values) {
    if (!value.empty()) {
      out.insert(value);
    }
  }
  return out;
}

bool passes_string_filter(
    const std::string& value,
    const std::unordered_set<std::string>& filter) {
  return filter.empty() || filter.find(value) != filter.end();
}

bool compatible_existing_store(
    const std::string& store_dir,
    const std::string& source_type,
    const std::string& source_path,
    const std::string& source_fingerprint,
    const std::string& filter_signature,
    bool require_molecules,
    InputStoreManifest* manifest_out = nullptr) {
  try {
    const auto manifest = read_input_store_manifest(store_dir);
    if (manifest_out != nullptr) {
      *manifest_out = manifest;
    }
    return manifest.format_version == InputStoreManifest{}.format_version &&
        manifest.source_type == source_type &&
        manifest.source_path == source_path &&
        manifest.source_fingerprint == source_fingerprint &&
        manifest.filter_signature == filter_signature &&
        manifest.has_cell_gene_counts &&
        (!require_molecules || manifest.has_molecules);
  } catch (...) {
    return false;
  }
}

std::shared_ptr<arrow::Array> build_string_array(const std::vector<std::string>& values) {
  arrow::StringBuilder builder;
  arrow_check(builder.Reserve(static_cast<int64_t>(values.size())), "Reserve string builder");
  for (const auto& value : values) {
    arrow_check(builder.Append(value), "Append string");
  }
  return arrow_unwrap(builder.Finish(), "Finish string array");
}

std::shared_ptr<arrow::Array> build_double_array(const std::vector<double>& values) {
  arrow::DoubleBuilder builder;
  arrow_check(builder.Reserve(static_cast<int64_t>(values.size())), "Reserve double builder");
  for (double value : values) {
    arrow_check(builder.Append(value), "Append double");
  }
  return arrow_unwrap(builder.Finish(), "Finish double array");
}

std::shared_ptr<arrow::Array> build_int32_array(const std::vector<int>& values) {
  arrow::Int32Builder builder;
  arrow_check(builder.Reserve(static_cast<int64_t>(values.size())), "Reserve int32 builder");
  for (int value : values) {
    arrow_check(builder.Append(value), "Append int32");
  }
  return arrow_unwrap(builder.Finish(), "Finish int32 array");
}

std::shared_ptr<arrow::Array> build_int64_array(const std::vector<std::int64_t>& values) {
  arrow::Int64Builder builder;
  arrow_check(builder.Reserve(static_cast<int64_t>(values.size())), "Reserve int64 builder");
  for (std::int64_t value : values) {
    arrow_check(builder.Append(value), "Append int64");
  }
  return arrow_unwrap(builder.Finish(), "Finish int64 array");
}

void write_parquet_table(
    const std::shared_ptr<arrow::Table>& table,
    const std::filesystem::path& path,
    int row_group_size) {
  ensure_directory(path.parent_path());
  auto sink = arrow_unwrap(arrow::io::FileOutputStream::Open(path.string()), "Open parquet output");
  auto writer = arrow_unwrap(
      parquet::arrow::FileWriter::Open(*table->schema(), arrow::default_memory_pool(), sink),
      "Open parquet writer");
  arrow_check(
      writer->WriteTable(
          *table,
          std::max<int64_t>(
              1,
              std::min<int64_t>(
                  table->num_rows(),
                  static_cast<int64_t>(std::max(row_group_size, 1))))),
      "Write parquet table");
  arrow_check(writer->Close(), "Close parquet writer");
  arrow_check(sink->Close(), "Close parquet sink");
}

std::shared_ptr<arrow::Table> read_parquet_table(const std::filesystem::path& path) {
  auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path.string()), "Open parquet input");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open parquet reader");
  auto reader = arrow_unwrap(builder.Build(), "Build parquet reader");
  reader->set_use_threads(true);
  std::shared_ptr<arrow::Table> table;
  arrow_check(reader->ReadTable(&table), "Read parquet table");
  return table;
}

std::shared_ptr<arrow::Table> read_delimited_table_as_strings(
    const std::filesystem::path& path,
    const std::string& id_col,
    const std::string& label_col) {
  auto input = open_delimited_input_stream(path.string(), "Open tabular metadata input");
  auto read_options = arrow::csv::ReadOptions::Defaults();
  auto parse_options = arrow::csv::ParseOptions::Defaults();
  auto convert_options = arrow::csv::ConvertOptions::Defaults();
  read_options.block_size = 1 << 20;
  const auto ext = logical_file_extension_lower(path.string());
  if (ext == "tsv") {
    parse_options.delimiter = '\t';
  }
  convert_options.column_types[id_col] = arrow::utf8();
  convert_options.column_types[label_col] = arrow::utf8();
  auto reader = arrow_unwrap(
      arrow::csv::TableReader::Make(
          arrow::io::default_io_context(),
          input,
          read_options,
          parse_options,
          convert_options),
      "Create tabular metadata reader");
  return arrow_unwrap(reader->Read(), "Read tabular metadata table");
}

int column_index(const std::shared_ptr<arrow::Table>& table, const std::string& name) {
  const int idx = table->schema()->GetFieldIndex(name);
  if (idx < 0) {
    throw std::runtime_error("Column '" + name + "' not found in store parquet");
  }
  return idx;
}

bool has_column(const std::shared_ptr<arrow::Table>& table, const std::string& name) {
  return table->schema()->GetFieldIndex(name) >= 0;
}

std::vector<int> read_int32_column(const std::shared_ptr<arrow::Table>& table, const std::string& name) {
  const auto column = table->column(column_index(table, name));
  std::vector<int> out(static_cast<std::size_t>(column->length()));
  int64_t offset = 0;
  for (const auto& chunk : column->chunks()) {
    if (chunk->type_id() == arrow::Type::INT32) {
      const auto arr = std::static_pointer_cast<arrow::Int32Array>(chunk);
      for (int64_t i = 0; i < arr->length(); ++i) out[static_cast<std::size_t>(offset + i)] = arr->Value(i);
    } else if (chunk->type_id() == arrow::Type::INT64) {
      const auto arr = std::static_pointer_cast<arrow::Int64Array>(chunk);
      for (int64_t i = 0; i < arr->length(); ++i) out[static_cast<std::size_t>(offset + i)] = static_cast<int>(arr->Value(i));
    } else {
      throw std::runtime_error("Column '" + name + "' must be int32/int64");
    }
    offset += chunk->length();
  }
  return out;
}

std::vector<std::int64_t> read_int64_column(const std::shared_ptr<arrow::Table>& table, const std::string& name) {
  const auto column = table->column(column_index(table, name));
  std::vector<std::int64_t> out(static_cast<std::size_t>(column->length()));
  int64_t offset = 0;
  for (const auto& chunk : column->chunks()) {
    if (chunk->type_id() == arrow::Type::INT64) {
      const auto arr = std::static_pointer_cast<arrow::Int64Array>(chunk);
      for (int64_t i = 0; i < arr->length(); ++i) out[static_cast<std::size_t>(offset + i)] = arr->Value(i);
    } else if (chunk->type_id() == arrow::Type::INT32) {
      const auto arr = std::static_pointer_cast<arrow::Int32Array>(chunk);
      for (int64_t i = 0; i < arr->length(); ++i) out[static_cast<std::size_t>(offset + i)] = arr->Value(i);
    } else {
      throw std::runtime_error("Column '" + name + "' must be int32/int64");
    }
    offset += chunk->length();
  }
  return out;
}

std::vector<double> read_double_column(const std::shared_ptr<arrow::Table>& table, const std::string& name) {
  const auto column = table->column(column_index(table, name));
  std::vector<double> out(static_cast<std::size_t>(column->length()));
  int64_t offset = 0;
  for (const auto& chunk : column->chunks()) {
    if (chunk->type_id() == arrow::Type::DOUBLE) {
      const auto arr = std::static_pointer_cast<arrow::DoubleArray>(chunk);
      for (int64_t i = 0; i < arr->length(); ++i) out[static_cast<std::size_t>(offset + i)] = arr->Value(i);
    } else if (chunk->type_id() == arrow::Type::FLOAT) {
      const auto arr = std::static_pointer_cast<arrow::FloatArray>(chunk);
      for (int64_t i = 0; i < arr->length(); ++i) out[static_cast<std::size_t>(offset + i)] = arr->Value(i);
    } else {
      throw std::runtime_error("Column '" + name + "' must be float/double");
    }
    offset += chunk->length();
  }
  return out;
}

std::vector<std::string> read_string_column(const std::shared_ptr<arrow::Table>& table, const std::string& name) {
  const auto column = table->column(column_index(table, name));
  std::vector<std::string> out;
  out.reserve(static_cast<std::size_t>(column->length()));
  for (const auto& chunk : column->chunks()) {
    if (chunk->type_id() == arrow::Type::STRING) {
      const auto arr = std::static_pointer_cast<arrow::StringArray>(chunk);
      for (int64_t i = 0; i < arr->length(); ++i) out.push_back(arr->IsNull(i) ? "" : arr->GetString(i));
    } else if (chunk->type_id() == arrow::Type::LARGE_STRING) {
      const auto arr = std::static_pointer_cast<arrow::LargeStringArray>(chunk);
      for (int64_t i = 0; i < arr->length(); ++i) out.push_back(arr->IsNull(i) ? "" : arr->GetString(i));
    } else {
      throw std::runtime_error("Column '" + name + "' must be string");
    }
  }
  return out;
}

std::unordered_map<std::string, std::string> read_cell_type_metadata(
    const TabularSourceSpec& source) {
  std::unordered_map<std::string, std::string> out;
  if (source.cell_metadata_path.empty() || source.cell_metadata_cell_type_col.empty()) {
    return out;
  }
  if (source.cell_metadata_cell_id_col.empty()) {
    throw std::runtime_error("cell_metadata_cell_id_col must be provided when cell metadata is supplied");
  }

  const auto path = std::filesystem::path(source.cell_metadata_path);
  const std::string ext = logical_file_extension_lower(path.string());
  std::shared_ptr<arrow::Table> table;
  if (is_parquet_extension(ext)) {
    table = read_parquet_table(path);
  } else if (is_delimited_text_extension(ext)) {
    table = read_delimited_table_as_strings(
        path,
        source.cell_metadata_cell_id_col,
        source.cell_metadata_cell_type_col);
  } else {
    throw std::runtime_error(
        "Unsupported cell metadata file format for '" + source.cell_metadata_path + "'");
  }

  if (!has_column(table, source.cell_metadata_cell_id_col)) {
    throw std::runtime_error(
        "Column '" + source.cell_metadata_cell_id_col + "' not found in cell metadata sidecar");
  }
  if (!has_column(table, source.cell_metadata_cell_type_col)) {
    throw std::runtime_error(
        "Column '" + source.cell_metadata_cell_type_col + "' not found in cell metadata sidecar");
  }

  const auto cell_ids = read_string_column(table, source.cell_metadata_cell_id_col);
  const auto cell_types = read_string_column(table, source.cell_metadata_cell_type_col);
  if (cell_ids.size() != cell_types.size()) {
    throw std::runtime_error("Cell metadata sidecar has inconsistent column lengths");
  }

  out.reserve(cell_ids.size());
  for (std::size_t i = 0; i < cell_ids.size(); ++i) {
    if (cell_ids[i].empty()) {
      continue;
    }
    const auto inserted = out.emplace(cell_ids[i], cell_types[i]);
    if (!inserted.second && inserted.first->second != cell_types[i]) {
      throw std::runtime_error(
          "Cell metadata sidecar contains conflicting labels for cell '" + cell_ids[i] + "'");
    }
  }
  return out;
}

std::string cell_type_for_cell(
    const std::unordered_map<std::string, std::string>& cell_type_by_cell,
    const std::string& cell_id) {
  const auto it = cell_type_by_cell.find(cell_id);
  return it == cell_type_by_cell.end() ? std::string() : it->second;
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

bool is_unassigned_cell_id(const std::string& cell_id) {
  return cell_id.empty() ||
      cell_id == "0" ||
      cell_id == "NA" ||
      cell_id == "NaN" ||
      cell_id == "null" ||
      cell_id == "UNASSIGNED" ||
      cell_id == "unassigned" ||
      cell_id == "cell_0";
}

int parse_bool_or_missing(const std::string& value) {
  if (value.empty()) {
    return -1;
  }
  std::string lower;
  lower.reserve(value.size());
  for (const char ch : value) {
    lower.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
  }
  if (lower == "true" || lower == "t" || lower == "yes" || lower == "y") {
    return 1;
  }
  if (lower == "false" || lower == "f" || lower == "no" || lower == "n") {
    return 0;
  }
  try {
    return std::stod(value) != 0.0 ? 1 : 0;
  } catch (...) {
    return -1;
  }
}

struct NumericArrayView {
  enum class Kind { Invalid, Double, Float, Int64, Int32, Int16, Int8, UInt64, UInt32, UInt16, UInt8, Fallback };

  const uint8_t* null_bitmap = nullptr;
  int64_t offset = 0;
  Kind kind = Kind::Invalid;
  const arrow::DoubleArray* double_arr = nullptr;
  const arrow::FloatArray* float_arr = nullptr;
  const arrow::Int64Array* int64_arr = nullptr;
  const arrow::Int32Array* int32_arr = nullptr;
  const arrow::Int16Array* int16_arr = nullptr;
  const arrow::Int8Array* int8_arr = nullptr;
  const arrow::UInt64Array* uint64_arr = nullptr;
  const arrow::UInt32Array* uint32_arr = nullptr;
  const arrow::UInt16Array* uint16_arr = nullptr;
  const arrow::UInt8Array* uint8_arr = nullptr;
  std::shared_ptr<arrow::Array> fallback;

  NumericArrayView() = default;
  explicit NumericArrayView(const std::shared_ptr<arrow::Array>& arr) { reset(arr); }

  void reset(const std::shared_ptr<arrow::Array>& arr) {
    null_bitmap = nullptr;
    offset = 0;
    kind = Kind::Invalid;
    double_arr = nullptr;
    float_arr = nullptr;
    int64_arr = nullptr;
    int32_arr = nullptr;
    int16_arr = nullptr;
    int8_arr = nullptr;
    uint64_arr = nullptr;
    uint32_arr = nullptr;
    uint16_arr = nullptr;
    uint8_arr = nullptr;
    fallback.reset();
    if (!arr) return;
    null_bitmap = arr->null_bitmap_data();
    offset = arr->offset();
    switch (arr->type_id()) {
      case arrow::Type::DOUBLE: double_arr = static_cast<const arrow::DoubleArray*>(arr.get()); kind = Kind::Double; break;
      case arrow::Type::FLOAT: float_arr = static_cast<const arrow::FloatArray*>(arr.get()); kind = Kind::Float; break;
      case arrow::Type::INT64: int64_arr = static_cast<const arrow::Int64Array*>(arr.get()); kind = Kind::Int64; break;
      case arrow::Type::INT32: int32_arr = static_cast<const arrow::Int32Array*>(arr.get()); kind = Kind::Int32; break;
      case arrow::Type::INT16: int16_arr = static_cast<const arrow::Int16Array*>(arr.get()); kind = Kind::Int16; break;
      case arrow::Type::INT8: int8_arr = static_cast<const arrow::Int8Array*>(arr.get()); kind = Kind::Int8; break;
      case arrow::Type::UINT64: uint64_arr = static_cast<const arrow::UInt64Array*>(arr.get()); kind = Kind::UInt64; break;
      case arrow::Type::UINT32: uint32_arr = static_cast<const arrow::UInt32Array*>(arr.get()); kind = Kind::UInt32; break;
      case arrow::Type::UINT16: uint16_arr = static_cast<const arrow::UInt16Array*>(arr.get()); kind = Kind::UInt16; break;
      case arrow::Type::UINT8: uint8_arr = static_cast<const arrow::UInt8Array*>(arr.get()); kind = Kind::UInt8; break;
      default: fallback = arr; kind = Kind::Fallback; break;
    }
  }

  bool is_valid(int64_t i) const {
    if (kind == Kind::Invalid) return false;
    return null_bitmap == nullptr || arrow::bit_util::GetBit(null_bitmap, offset + i);
  }

  double value(int64_t i) const {
    if (!is_valid(i)) return std::numeric_limits<double>::quiet_NaN();
    switch (kind) {
      case Kind::Double: return double_arr->Value(i);
      case Kind::Float: return static_cast<double>(float_arr->Value(i));
      case Kind::Int64: return static_cast<double>(int64_arr->Value(i));
      case Kind::Int32: return static_cast<double>(int32_arr->Value(i));
      case Kind::Int16: return static_cast<double>(int16_arr->Value(i));
      case Kind::Int8: return static_cast<double>(int8_arr->Value(i));
      case Kind::UInt64: return static_cast<double>(uint64_arr->Value(i));
      case Kind::UInt32: return static_cast<double>(uint32_arr->Value(i));
      case Kind::UInt16: return static_cast<double>(uint16_arr->Value(i));
      case Kind::UInt8: return static_cast<double>(uint8_arr->Value(i));
      case Kind::Fallback: {
        auto scalar = arrow_unwrap(fallback->GetScalar(i), "Read numeric scalar");
        return std::stod(scalar->ToString());
      }
      default: break;
    }
    return std::numeric_limits<double>::quiet_NaN();
  }

  int64_t int64_value(int64_t i) const {
    if (!is_valid(i)) return -1;
    switch (kind) {
      case Kind::Int64: return int64_arr->Value(i);
      case Kind::Int32: return int32_arr->Value(i);
      case Kind::Int16: return int16_arr->Value(i);
      case Kind::Int8: return int8_arr->Value(i);
      case Kind::UInt64: return static_cast<int64_t>(uint64_arr->Value(i));
      case Kind::UInt32: return static_cast<int64_t>(uint32_arr->Value(i));
      case Kind::UInt16: return static_cast<int64_t>(uint16_arr->Value(i));
      case Kind::UInt8: return static_cast<int64_t>(uint8_arr->Value(i));
      case Kind::Double: return static_cast<int64_t>(std::llround(double_arr->Value(i)));
      case Kind::Float: return static_cast<int64_t>(std::llround(float_arr->Value(i)));
      case Kind::Fallback: {
        auto scalar = arrow_unwrap(fallback->GetScalar(i), "Read integer scalar");
        return static_cast<int64_t>(std::stoll(scalar->ToString()));
      }
      default: break;
    }
    return -1;
  }

  std::string string_value(int64_t i) const {
    if (!is_valid(i)) return {};
    switch (kind) {
      case Kind::Int64: return std::to_string(int64_arr->Value(i));
      case Kind::Int32: return std::to_string(int32_arr->Value(i));
      case Kind::Int16: return std::to_string(int16_arr->Value(i));
      case Kind::Int8: return std::to_string(static_cast<int>(int8_arr->Value(i)));
      case Kind::UInt64: return std::to_string(uint64_arr->Value(i));
      case Kind::UInt32: return std::to_string(uint32_arr->Value(i));
      case Kind::UInt16: return std::to_string(uint16_arr->Value(i));
      case Kind::UInt8: return std::to_string(static_cast<unsigned int>(uint8_arr->Value(i)));
      case Kind::Double: {
        std::ostringstream out;
        out << std::setprecision(std::numeric_limits<double>::max_digits10) << double_arr->Value(i);
        return out.str();
      }
      case Kind::Float: {
        std::ostringstream out;
        out << std::setprecision(std::numeric_limits<float>::max_digits10) << float_arr->Value(i);
        return out.str();
      }
      case Kind::Fallback: {
        auto scalar = arrow_unwrap(fallback->GetScalar(i), "Read string scalar");
        return scalar->ToString();
      }
      default: break;
    }
    return {};
  }
};

struct StringArrayView {
  const arrow::StringArray* string_arr = nullptr;
  const arrow::LargeStringArray* large_string_arr = nullptr;
  const arrow::DictionaryArray* dict_arr = nullptr;
  const arrow::StringArray* dict_values = nullptr;
  const arrow::LargeStringArray* dict_large_values = nullptr;
  NumericArrayView dict_indices;
  NumericArrayView numeric_arr;
  NumericArrayView dict_numeric_values;
  bool numeric_mode = false;
  bool dict_numeric_mode = false;
  std::shared_ptr<arrow::Array> fallback;
  std::shared_ptr<arrow::Array> dict_fallback_values;

  StringArrayView() = default;
  explicit StringArrayView(const std::shared_ptr<arrow::Array>& arr) { reset(arr); }

  void reset(const std::shared_ptr<arrow::Array>& arr) {
    string_arr = nullptr;
    large_string_arr = nullptr;
    dict_arr = nullptr;
    dict_values = nullptr;
    dict_large_values = nullptr;
    dict_indices = NumericArrayView();
    numeric_arr = NumericArrayView();
    dict_numeric_values = NumericArrayView();
    numeric_mode = false;
    dict_numeric_mode = false;
    fallback.reset();
    dict_fallback_values.reset();
    if (!arr) return;
    if (arr->type_id() == arrow::Type::STRING) {
      string_arr = static_cast<const arrow::StringArray*>(arr.get());
      return;
    }
    if (arr->type_id() == arrow::Type::LARGE_STRING) {
      large_string_arr = static_cast<const arrow::LargeStringArray*>(arr.get());
      return;
    }
    if (arr->type_id() == arrow::Type::DICTIONARY) {
      dict_arr = static_cast<const arrow::DictionaryArray*>(arr.get());
      const auto dict = dict_arr->dictionary();
      if (dict->type_id() == arrow::Type::STRING) {
        dict_values = static_cast<const arrow::StringArray*>(dict.get());
      } else if (dict->type_id() == arrow::Type::LARGE_STRING) {
        dict_large_values = static_cast<const arrow::LargeStringArray*>(dict.get());
      } else if (
          dict->type_id() == arrow::Type::DOUBLE ||
          dict->type_id() == arrow::Type::FLOAT ||
          dict->type_id() == arrow::Type::INT64 ||
          dict->type_id() == arrow::Type::INT32 ||
          dict->type_id() == arrow::Type::INT16 ||
          dict->type_id() == arrow::Type::INT8 ||
          dict->type_id() == arrow::Type::UINT64 ||
          dict->type_id() == arrow::Type::UINT32 ||
          dict->type_id() == arrow::Type::UINT16 ||
          dict->type_id() == arrow::Type::UINT8) {
        dict_numeric_values.reset(dict);
        dict_numeric_mode = true;
      } else {
        dict_fallback_values = dict;
      }
      dict_indices.reset(dict_arr->indices());
      return;
    }
    if (
        arr->type_id() == arrow::Type::DOUBLE ||
        arr->type_id() == arrow::Type::FLOAT ||
        arr->type_id() == arrow::Type::INT64 ||
        arr->type_id() == arrow::Type::INT32 ||
        arr->type_id() == arrow::Type::INT16 ||
        arr->type_id() == arrow::Type::INT8 ||
        arr->type_id() == arrow::Type::UINT64 ||
        arr->type_id() == arrow::Type::UINT32 ||
        arr->type_id() == arrow::Type::UINT16 ||
        arr->type_id() == arrow::Type::UINT8) {
      numeric_arr.reset(arr);
      numeric_mode = true;
      return;
    }
    fallback = arr;
  }

  std::string value(int64_t i) const {
    if (string_arr != nullptr) return string_arr->IsNull(i) ? std::string() : string_arr->GetString(i);
    if (large_string_arr != nullptr) return large_string_arr->IsNull(i) ? std::string() : large_string_arr->GetString(i);
    if (dict_arr != nullptr && !dict_arr->IsNull(i)) {
      const int64_t idx = dict_indices.int64_value(i);
      if (dict_values != nullptr) return dict_values->GetString(idx);
      if (dict_large_values != nullptr) return dict_large_values->GetString(idx);
      if (dict_numeric_mode) return dict_numeric_values.string_value(idx);
      if (dict_fallback_values != nullptr) {
        auto scalar = arrow_unwrap(dict_fallback_values->GetScalar(idx), "Read dictionary scalar");
        return scalar->ToString();
      }
      return {};
    }
    if (numeric_mode) return numeric_arr.string_value(i);
    if (fallback != nullptr) {
      auto scalar = arrow_unwrap(fallback->GetScalar(i), "Read string scalar");
      return scalar->ToString();
    }
    return {};
  }
};

std::size_t find_first_field(
    const std::shared_ptr<arrow::Schema>& schema,
    const std::vector<std::string>& candidates,
    const std::string& label) {
  for (const auto& candidate : candidates) {
    const int idx = schema->GetFieldIndex(candidate);
    if (idx >= 0) return static_cast<std::size_t>(idx);
  }
  throw std::runtime_error("Missing required Xenium parquet column: " + label);
}

std::optional<std::size_t> find_optional_field(
    const std::shared_ptr<arrow::Schema>& schema,
    const std::vector<std::string>& candidates) {
  for (const auto& candidate : candidates) {
    const int idx = schema->GetFieldIndex(candidate);
    if (idx >= 0) return static_cast<std::size_t>(idx);
  }
  return std::nullopt;
}

std::optional<std::pair<double, double>> parquet_numeric_minmax(
    const parquet::ColumnChunkMetaData& column) {
  if (!column.is_stats_set()) return std::nullopt;
  const auto stats = column.statistics();
  if (!stats || !stats->HasMinMax()) return std::nullopt;
  switch (stats->physical_type()) {
    case parquet::Type::DOUBLE: {
      const auto typed = std::static_pointer_cast<parquet::DoubleStatistics>(stats);
      return std::make_pair(static_cast<double>(typed->min()), static_cast<double>(typed->max()));
    }
    case parquet::Type::FLOAT: {
      const auto typed = std::static_pointer_cast<parquet::FloatStatistics>(stats);
      return std::make_pair(static_cast<double>(typed->min()), static_cast<double>(typed->max()));
    }
    case parquet::Type::INT64: {
      const auto typed = std::static_pointer_cast<parquet::Int64Statistics>(stats);
      return std::make_pair(static_cast<double>(typed->min()), static_cast<double>(typed->max()));
    }
    case parquet::Type::INT32: {
      const auto typed = std::static_pointer_cast<parquet::Int32Statistics>(stats);
      return std::make_pair(static_cast<double>(typed->min()), static_cast<double>(typed->max()));
    }
    default: break;
  }
  return std::nullopt;
}

bool row_group_numeric_range_overlaps(
    const parquet::RowGroupMetaData& row_group,
    int column_index,
    double min_value,
    double max_value) {
  if (column_index < 0 || column_index >= row_group.num_columns()) return true;
  const auto column = row_group.ColumnChunk(column_index);
  if (!column) return true;
  const auto minmax = parquet_numeric_minmax(*column);
  if (!minmax.has_value() || !std::isfinite(minmax->first) || !std::isfinite(minmax->second)) return true;
  return !(minmax->second < min_value || minmax->first > max_value);
}

bool row_group_numeric_max_at_least(
    const parquet::RowGroupMetaData& row_group,
    int column_index,
    double min_value) {
  if (column_index < 0 || column_index >= row_group.num_columns()) return true;
  const auto column = row_group.ColumnChunk(column_index);
  if (!column) return true;
  const auto minmax = parquet_numeric_minmax(*column);
  if (!minmax.has_value() || !std::isfinite(minmax->second)) return true;
  return minmax->second >= min_value;
}

std::vector<int> select_row_groups_for_numeric_filters(
    const std::string& path,
    std::size_t x_col,
    std::size_t y_col,
    const std::optional<std::size_t>& z_col,
    const std::optional<std::size_t>& qv_col,
    const UnionBounds& bounds,
    double min_qv) {
  if (bounds.has_z && !z_col.has_value() && (0.0 < bounds.zmin || 0.0 > bounds.zmax)) {
    return {};
  }
  auto parquet_reader = parquet::ParquetFileReader::OpenFile(path, false);
  const auto metadata = parquet_reader->metadata();
  std::vector<int> row_groups;
  row_groups.reserve(static_cast<std::size_t>(metadata->num_row_groups()));
  for (int rg = 0; rg < metadata->num_row_groups(); ++rg) {
    const auto row_group = metadata->RowGroup(rg);
    bool keep = true;
    if (bounds.active) {
      keep = keep &&
          row_group_numeric_range_overlaps(*row_group, static_cast<int>(x_col), bounds.xmin, bounds.xmax) &&
          row_group_numeric_range_overlaps(*row_group, static_cast<int>(y_col), bounds.ymin, bounds.ymax);
      if (keep && bounds.has_z && z_col.has_value()) {
        keep = row_group_numeric_range_overlaps(*row_group, static_cast<int>(*z_col), bounds.zmin, bounds.zmax);
      }
    }
    if (keep && min_qv >= 0.0 && qv_col.has_value()) {
      keep = row_group_numeric_max_at_least(*row_group, static_cast<int>(*qv_col), min_qv);
    }
    if (keep) row_groups.push_back(rg);
  }
  return row_groups;
}

void annotate_transcripts_from_cell_table(TranscriptTable& table, const CellTable& cells) {
  if (cells.cell_ids.empty()) return;
  std::unordered_map<std::string, std::size_t> cell_lookup;
  cell_lookup.reserve(cells.cell_ids.size());
  for (std::size_t i = 0; i < cells.cell_ids.size(); ++i) {
    cell_lookup.emplace(cells.cell_ids[i], i);
  }
  const bool add_cell_types = table.cell_types.empty() && !cells.cell_types.empty();
  const bool add_sample_ids = table.sample_ids.empty() && !cells.sample_ids.empty();
  const bool add_fov_ids = table.fov_ids.empty() && !cells.fov_ids.empty();
  if (add_cell_types) table.cell_types.assign(table.size(), "");
  if (add_sample_ids) table.sample_ids.assign(table.size(), "");
  if (add_fov_ids) table.fov_ids.assign(table.size(), "");
  for (std::size_t i = 0; i < table.size(); ++i) {
    const auto it = cell_lookup.find(table.orig_cell_id[i]);
    if (it == cell_lookup.end()) continue;
    const std::size_t cell = it->second;
    if (add_cell_types) table.cell_types[i] = cells.cell_types[cell];
    if (add_sample_ids) table.sample_ids[i] = cells.sample_ids[cell];
    if (add_fov_ids) table.fov_ids[i] = cells.fov_ids[cell];
  }
}

std::vector<std::string> crop_ids_by_cell(
    const TranscriptTable& table,
    const std::vector<std::string>* transcript_crop_ids) {
  std::vector<std::string> out(table.num_cells(), "");
  if (transcript_crop_ids == nullptr || transcript_crop_ids->size() != table.size()) {
    return out;
  }
  for (std::size_t i = 0; i < table.size(); ++i) {
    auto& current = out[static_cast<std::size_t>(table.cell_index[i])];
    if (current.empty()) {
      current = (*transcript_crop_ids)[i];
    }
  }
  return out;
}

CellCountMatrix build_counts_from_table(
    const TranscriptTable& table,
    const CellTable& cells,
    const std::vector<std::string>* transcript_crop_ids) {
  CellCountMatrix out;
  out.genes = table.genes;
  out.cells = cells;
  if (out.cells.cell_ids.empty()) {
    out.cells = infer_cell_table(table);
  }
  const auto by_cell = table.transcripts_by_cell();
  const auto cell_crop_ids = crop_ids_by_cell(table, transcript_crop_ids);
  out.crop_ids = cell_crop_ids;
  out.transcript_counts.reserve(by_cell.size());
  out.indptr.reserve(by_cell.size() + 1);
  out.indptr.push_back(0);
  for (const auto& members : by_cell) {
    std::unordered_map<int, int> gene_counts;
    gene_counts.reserve(members.size());
    for (const int row : members) {
      gene_counts[table.gene_index[static_cast<std::size_t>(row)]] += 1;
    }
    std::vector<std::pair<int, int>> entries;
    entries.reserve(gene_counts.size());
    for (const auto& kv : gene_counts) entries.emplace_back(kv.first, kv.second);
    std::sort(entries.begin(), entries.end());
    for (const auto& entry : entries) {
      out.indices.push_back(entry.first);
      out.values.push_back(static_cast<double>(entry.second));
    }
    out.transcript_counts.push_back(static_cast<int>(members.size()));
    out.indptr.push_back(static_cast<int>(out.indices.size()));
  }
  return out;
}

void validate_counts_upgrade_compatibility(
    const CellCountMatrix& existing,
    const CellCountMatrix& candidate,
    const std::string& context) {
  if (existing.genes != candidate.genes) {
    throw std::runtime_error(context + ": full-store gene dictionary differs from existing counts store");
  }
  if (existing.cells.cell_ids != candidate.cells.cell_ids) {
    throw std::runtime_error(context + ": full-store cell dictionary differs from existing counts store");
  }
  if (existing.transcript_counts != candidate.transcript_counts) {
    throw std::runtime_error(context + ": full-store transcript counts differ from existing counts store");
  }
  if (existing.indptr != candidate.indptr ||
      existing.indices != candidate.indices ||
      existing.values != candidate.values) {
    throw std::runtime_error(context + ": full-store sparse cell-gene counts differ from existing counts store");
  }
  if (!existing.crop_ids.empty() && !candidate.crop_ids.empty() && existing.crop_ids != candidate.crop_ids) {
    throw std::runtime_error(context + ": full-store crop ids differ from existing counts store");
  }
}

std::size_t find_first_field(
    const std::shared_ptr<arrow::Schema>& schema,
    const std::vector<std::string>& candidates,
    const std::string& label);

std::optional<std::size_t> find_optional_field(
    const std::shared_ptr<arrow::Schema>& schema,
    const std::vector<std::string>& candidates);

void enrich_counts_from_cells_parquet(CellCountMatrix& counts, const std::string& path) {
  if (!file_exists(path) || counts.cells.cell_ids.empty()) {
    return;
  }
  const auto table = read_parquet_table(path);
  const auto schema = table->schema();
  const auto cell_id_idx = find_first_field(schema, {"cell_id", "cell"}, "cell id");
  const auto x_idx = find_first_field(schema, {"x_centroid", "x"}, "cell x");
  const auto y_idx = find_first_field(schema, {"y_centroid", "y"}, "cell y");
  const auto z_idx = find_optional_field(schema, {"z_centroid", "z"});
  const auto cell_type_idx = find_optional_field(schema, {"cell_type", "celltype"});
  const auto sample_id_idx = find_optional_field(schema, {"sample_id", "sample"});
  const auto fov_id_idx = find_optional_field(schema, {"fov_name", "fov_id", "fovn"});

  const auto field_name = [&](std::size_t idx) {
    return schema->field(static_cast<int>(idx))->name();
  };
  const auto cell_ids = read_string_column(table, field_name(cell_id_idx));
  const auto centroid_x = read_double_column(table, field_name(x_idx));
  const auto centroid_y = read_double_column(table, field_name(y_idx));
  const auto centroid_z = z_idx.has_value()
      ? read_double_column(table, field_name(*z_idx))
      : std::vector<double>(cell_ids.size(), 0.0);
  const auto cell_types = cell_type_idx.has_value()
      ? read_string_column(table, field_name(*cell_type_idx))
      : std::vector<std::string>();
  const auto sample_ids = sample_id_idx.has_value()
      ? read_string_column(table, field_name(*sample_id_idx))
      : std::vector<std::string>();
  const auto fov_ids = fov_id_idx.has_value()
      ? read_string_column(table, field_name(*fov_id_idx))
      : std::vector<std::string>();

  std::unordered_map<std::string, std::size_t> row_by_cell;
  row_by_cell.reserve(cell_ids.size());
  for (std::size_t i = 0; i < cell_ids.size(); ++i) {
    row_by_cell.emplace(cell_ids[i], i);
  }

  if (cell_type_idx.has_value()) counts.cells.cell_types.assign(counts.cells.cell_ids.size(), "");
  if (sample_id_idx.has_value()) counts.cells.sample_ids.assign(counts.cells.cell_ids.size(), "");
  if (fov_id_idx.has_value()) counts.cells.fov_ids.assign(counts.cells.cell_ids.size(), "");
  for (std::size_t cell = 0; cell < counts.cells.cell_ids.size(); ++cell) {
    const auto it = row_by_cell.find(counts.cells.cell_ids[cell]);
    if (it == row_by_cell.end()) {
      continue;
    }
    const std::size_t row = it->second;
    if (row < centroid_x.size()) counts.cells.centroid_x[cell] = centroid_x[row];
    if (row < centroid_y.size()) counts.cells.centroid_y[cell] = centroid_y[row];
    if (row < centroid_z.size()) counts.cells.centroid_z[cell] = centroid_z[row];
    if (cell_type_idx.has_value() && row < cell_types.size()) counts.cells.cell_types[cell] = cell_types[row];
    if (sample_id_idx.has_value() && row < sample_ids.size()) counts.cells.sample_ids[cell] = sample_ids[row];
    if (fov_id_idx.has_value() && row < fov_ids.size()) counts.cells.fov_ids[cell] = fov_ids[row];
  }
}

struct XeniumParquetCountsResult {
  CellCountMatrix counts;
  bool has_z = true;
  bool has_qv = false;
  bool has_transcript_id = false;
  std::int64_t scanned_rows = 0;
  std::int64_t kept_rows = 0;
};

struct EncodedMoleculeStore {
  CellCountMatrix counts;
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
  bool has_z = true;
  bool has_qv = false;
  bool has_transcript_id = false;
  bool has_overlaps_nucleus = false;
  bool has_nucleus_distance = false;
  std::int64_t scanned_rows = 0;
  std::int64_t kept_rows = 0;
};

class CellGeneCountAccumulator {
 public:
  int add_gene(const std::string& gene) {
    const auto it = gene_lookup_.find(gene);
    if (it != gene_lookup_.end()) {
      return it->second;
    }
    const int id = static_cast<int>(counts_.genes.size());
    gene_lookup_.emplace(gene, id);
    counts_.genes.push_back(gene);
    return id;
  }

  void ensure_cell_metadata_vectors(
      const std::string& cell_type,
      const std::string& sample_id,
      const std::string& fov_id) {
    if (!cell_type.empty() && counts_.cells.cell_types.empty()) {
      counts_.cells.cell_types.assign(counts_.cells.cell_ids.size(), "");
    }
    if (!sample_id.empty() && counts_.cells.sample_ids.empty()) {
      counts_.cells.sample_ids.assign(counts_.cells.cell_ids.size(), "");
    }
    if (!fov_id.empty() && counts_.cells.fov_ids.empty()) {
      counts_.cells.fov_ids.assign(counts_.cells.cell_ids.size(), "");
    }
  }

  int add_cell(
      const std::string& cell_id,
      const std::string& crop_id,
      const std::string& cell_type = "",
      const std::string& sample_id = "",
      const std::string& fov_id = "") {
    const auto it = cell_lookup_.find(cell_id);
    if (it != cell_lookup_.end()) {
      const auto cell = static_cast<std::size_t>(it->second);
      ensure_cell_metadata_vectors(cell_type, sample_id, fov_id);
      if (!cell_type.empty() && counts_.cells.cell_types[cell].empty()) {
        counts_.cells.cell_types[cell] = cell_type;
      }
      if (!sample_id.empty() && counts_.cells.sample_ids[cell].empty()) {
        counts_.cells.sample_ids[cell] = sample_id;
      }
      if (!fov_id.empty() && counts_.cells.fov_ids[cell].empty()) {
        counts_.cells.fov_ids[cell] = fov_id;
      }
      return it->second;
    }
    ensure_cell_metadata_vectors(cell_type, sample_id, fov_id);
    const int id = static_cast<int>(counts_.cells.cell_ids.size());
    cell_lookup_.emplace(cell_id, id);
    counts_.cells.cell_ids.push_back(cell_id);
    counts_.cells.centroid_x.push_back(0.0);
    counts_.cells.centroid_y.push_back(0.0);
    counts_.cells.centroid_z.push_back(0.0);
    counts_.transcript_counts.push_back(0);
    counts_.crop_ids.push_back(crop_id);
    sum_x_.push_back(0.0);
    sum_y_.push_back(0.0);
    sum_z_.push_back(0.0);
    cell_gene_counts_.emplace_back();
    if (!counts_.cells.cell_types.empty()) counts_.cells.cell_types.push_back(cell_type);
    if (!counts_.cells.sample_ids.empty()) counts_.cells.sample_ids.push_back(sample_id);
    if (!counts_.cells.fov_ids.empty()) counts_.cells.fov_ids.push_back(fov_id);
    return id;
  }

  void add_transcript(
      const std::string& cell_id,
      const std::string& gene,
      double x,
      double y,
      double z,
      const std::string& crop_id,
      const std::string& cell_type = "",
      const std::string& sample_id = "",
      const std::string& fov_id = "") {
    const int gene_id = add_gene(gene);
    const int cell_id_int = add_cell(cell_id, crop_id, cell_type, sample_id, fov_id);
    const auto cell = static_cast<std::size_t>(cell_id_int);
    cell_gene_counts_[cell][gene_id] += 1;
    counts_.transcript_counts[cell] += 1;
    sum_x_[cell] += x;
    sum_y_[cell] += y;
    sum_z_[cell] += z;
    if (counts_.crop_ids[cell].empty()) {
      counts_.crop_ids[cell] = crop_id;
    }
  }

  CellCountMatrix finalize() {
    counts_.indptr.reserve(counts_.cells.cell_ids.size() + 1);
    counts_.indptr.push_back(0);
    for (std::size_t cell = 0; cell < counts_.cells.cell_ids.size(); ++cell) {
      const int n = std::max(counts_.transcript_counts[cell], 1);
      counts_.cells.centroid_x[cell] = sum_x_[cell] / static_cast<double>(n);
      counts_.cells.centroid_y[cell] = sum_y_[cell] / static_cast<double>(n);
      counts_.cells.centroid_z[cell] = sum_z_[cell] / static_cast<double>(n);
      std::vector<std::pair<int, int>> entries;
      entries.reserve(cell_gene_counts_[cell].size());
      for (const auto& kv : cell_gene_counts_[cell]) {
        entries.emplace_back(kv.first, kv.second);
      }
      std::sort(entries.begin(), entries.end());
      for (const auto& entry : entries) {
        counts_.indices.push_back(entry.first);
        counts_.values.push_back(static_cast<double>(entry.second));
      }
      counts_.indptr.push_back(static_cast<int>(counts_.indices.size()));
    }
    return std::move(counts_);
  }

 private:
  CellCountMatrix counts_;
  std::unordered_map<std::string, int> gene_lookup_;
  std::unordered_map<std::string, int> cell_lookup_;
  std::vector<double> sum_x_;
  std::vector<double> sum_y_;
  std::vector<double> sum_z_;
  std::vector<std::unordered_map<int, int>> cell_gene_counts_;
};

class EncodedMoleculeAccumulator {
 public:
  explicit EncodedMoleculeAccumulator(
      bool store_qv,
      bool store_overlaps_nucleus,
      bool store_nucleus_distance)
      : store_qv_(store_qv),
        store_overlaps_nucleus_(store_overlaps_nucleus),
        store_nucleus_distance_(store_nucleus_distance) {}

  int add_gene(const std::string& gene) {
    const auto it = gene_lookup_.find(gene);
    if (it != gene_lookup_.end()) {
      return it->second;
    }
    const int id = static_cast<int>(out_.counts.genes.size());
    gene_lookup_.emplace(gene, id);
    out_.counts.genes.push_back(gene);
    return id;
  }

  void ensure_cell_metadata_vectors(
      const std::string& cell_type,
      const std::string& sample_id,
      const std::string& fov_id) {
    if (!cell_type.empty() && out_.counts.cells.cell_types.empty()) {
      out_.counts.cells.cell_types.assign(out_.counts.cells.cell_ids.size(), "");
    }
    if (!sample_id.empty() && out_.counts.cells.sample_ids.empty()) {
      out_.counts.cells.sample_ids.assign(out_.counts.cells.cell_ids.size(), "");
    }
    if (!fov_id.empty() && out_.counts.cells.fov_ids.empty()) {
      out_.counts.cells.fov_ids.assign(out_.counts.cells.cell_ids.size(), "");
    }
  }

  int add_cell(
      const std::string& cell_id,
      const std::string& crop_id,
      const std::string& cell_type = "",
      const std::string& sample_id = "",
      const std::string& fov_id = "") {
    const auto it = cell_lookup_.find(cell_id);
    if (it != cell_lookup_.end()) {
      const auto cell = static_cast<std::size_t>(it->second);
      ensure_cell_metadata_vectors(cell_type, sample_id, fov_id);
      if (!cell_type.empty() && out_.counts.cells.cell_types[cell].empty()) {
        out_.counts.cells.cell_types[cell] = cell_type;
      }
      if (!sample_id.empty() && out_.counts.cells.sample_ids[cell].empty()) {
        out_.counts.cells.sample_ids[cell] = sample_id;
      }
      if (!fov_id.empty() && out_.counts.cells.fov_ids[cell].empty()) {
        out_.counts.cells.fov_ids[cell] = fov_id;
      }
      return it->second;
    }
    ensure_cell_metadata_vectors(cell_type, sample_id, fov_id);
    const int id = static_cast<int>(out_.counts.cells.cell_ids.size());
    cell_lookup_.emplace(cell_id, id);
    out_.counts.cells.cell_ids.push_back(cell_id);
    out_.counts.cells.centroid_x.push_back(0.0);
    out_.counts.cells.centroid_y.push_back(0.0);
    out_.counts.cells.centroid_z.push_back(0.0);
    out_.counts.transcript_counts.push_back(0);
    out_.counts.crop_ids.push_back(crop_id);
    sum_x_.push_back(0.0);
    sum_y_.push_back(0.0);
    sum_z_.push_back(0.0);
    cell_gene_counts_.emplace_back();
    if (!out_.counts.cells.cell_types.empty()) out_.counts.cells.cell_types.push_back(cell_type);
    if (!out_.counts.cells.sample_ids.empty()) out_.counts.cells.sample_ids.push_back(sample_id);
    if (!out_.counts.cells.fov_ids.empty()) out_.counts.cells.fov_ids.push_back(fov_id);
    return id;
  }

  void add_transcript(
      const std::string& cell_id,
      const std::string& gene,
      double x,
      double y,
      double z,
      double qv,
      int overlaps_nucleus,
      double nucleus_distance,
      const std::string& crop_id,
      const std::string& cell_type = "",
      const std::string& sample_id = "",
      const std::string& fov_id = "") {
    const int gene_id = add_gene(gene);
    const int cell_id_int = add_cell(cell_id, crop_id, cell_type, sample_id, fov_id);
    const auto cell = static_cast<std::size_t>(cell_id_int);
    cell_gene_counts_[cell][gene_id] += 1;
    out_.counts.transcript_counts[cell] += 1;
    sum_x_[cell] += x;
    sum_y_[cell] += y;
    sum_z_[cell] += z;
    if (out_.counts.crop_ids[cell].empty()) {
      out_.counts.crop_ids[cell] = crop_id;
    }

    out_.obs_id.push_back(static_cast<std::int64_t>(out_.obs_id.size()));
    out_.cell_idx.push_back(cell_id_int);
    out_.gene_idx.push_back(gene_id);
    out_.x.push_back(x);
    out_.y.push_back(y);
    out_.z.push_back(z);
    if (store_qv_) {
      out_.qv.push_back(qv);
    }
    if (store_overlaps_nucleus_) {
      out_.overlaps_nucleus.push_back(overlaps_nucleus);
    }
    if (store_nucleus_distance_) {
      out_.nucleus_distance.push_back(nucleus_distance);
    }
    out_.crop_id.push_back(crop_id);
  }

  EncodedMoleculeStore finalize() {
    out_.counts.indptr.reserve(out_.counts.cells.cell_ids.size() + 1);
    out_.counts.indptr.push_back(0);
    for (std::size_t cell = 0; cell < out_.counts.cells.cell_ids.size(); ++cell) {
      const int n = std::max(out_.counts.transcript_counts[cell], 1);
      out_.counts.cells.centroid_x[cell] = sum_x_[cell] / static_cast<double>(n);
      out_.counts.cells.centroid_y[cell] = sum_y_[cell] / static_cast<double>(n);
      out_.counts.cells.centroid_z[cell] = sum_z_[cell] / static_cast<double>(n);
      std::vector<std::pair<int, int>> entries;
      entries.reserve(cell_gene_counts_[cell].size());
      for (const auto& kv : cell_gene_counts_[cell]) {
        entries.emplace_back(kv.first, kv.second);
      }
      std::sort(entries.begin(), entries.end());
      for (const auto& entry : entries) {
        out_.counts.indices.push_back(entry.first);
        out_.counts.values.push_back(static_cast<double>(entry.second));
      }
      out_.counts.indptr.push_back(static_cast<int>(out_.counts.indices.size()));
    }
    out_.has_qv = store_qv_;
    out_.has_overlaps_nucleus = store_overlaps_nucleus_;
    out_.has_nucleus_distance = store_nucleus_distance_;
    return std::move(out_);
  }

 private:
  bool store_qv_ = false;
  bool store_overlaps_nucleus_ = false;
  bool store_nucleus_distance_ = false;
  EncodedMoleculeStore out_;
  std::unordered_map<std::string, int> gene_lookup_;
  std::unordered_map<std::string, int> cell_lookup_;
  std::vector<double> sum_x_;
  std::vector<double> sum_y_;
  std::vector<double> sum_z_;
  std::vector<std::unordered_map<int, int>> cell_gene_counts_;
};

XeniumParquetCountsResult build_xenium_counts_from_parquet(
    const std::string& path,
    const XeniumLoadOptions& options) {
  auto stage_start = std::chrono::steady_clock::now();
  auto emit_progress = [&](const std::string& message) {
    if (!options.progress) return;
    const auto now = std::chrono::steady_clock::now();
    options.progress(message, std::chrono::duration<double>(now - stage_start).count());
    stage_start = now;
  };

  emit_progress("Opening Xenium transcripts parquet for counts-only store: " + path);
  auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path), "Open Xenium transcript parquet");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open Xenium transcript parquet reader");
  auto reader = arrow_unwrap(builder.Build(), "Build Xenium transcript parquet reader");
  reader->set_use_threads(true);
  reader->set_batch_size(65536);

  std::shared_ptr<arrow::Schema> schema;
  arrow_check(reader->GetSchema(&schema), "Read Xenium transcript parquet schema");
  const auto gene_col = find_first_field(schema, {"feature_name", "gene", "gene_name", "gene_id"}, "gene");
  const auto cell_col = find_first_field(schema, {"cell_id", "cell"}, "cell id");
  const auto x_col = find_first_field(schema, {"x_location", "x"}, "x");
  const auto y_col = find_first_field(schema, {"y_location", "y"}, "y");
  const auto z_col = find_optional_field(schema, {"z_location", "z"});
  const auto qv_col = find_optional_field(schema, {"qv"});
  const auto txid_col = find_optional_field(schema, {"transcript_id"});
  const auto codeword_category_col = find_optional_field(schema, {"codeword_category"});
  const auto is_gene_col = find_optional_field(schema, {"is_gene"});
  const bool use_codeword_category_filter = !options.keep_non_gene && codeword_category_col.has_value();
  const bool use_is_gene_filter = !options.keep_non_gene && !use_codeword_category_filter && is_gene_col.has_value();
  const UnionBounds bounds = union_bounds(options.crops);

  std::vector<int> projected = {static_cast<int>(x_col), static_cast<int>(y_col)};
  if (z_col.has_value()) projected.push_back(static_cast<int>(*z_col));
  if (qv_col.has_value()) projected.push_back(static_cast<int>(*qv_col));
  projected.push_back(static_cast<int>(gene_col));
  projected.push_back(static_cast<int>(cell_col));
  if (use_codeword_category_filter) projected.push_back(static_cast<int>(*codeword_category_col));
  if (use_is_gene_filter) projected.push_back(static_cast<int>(*is_gene_col));

  const auto row_groups = select_row_groups_for_numeric_filters(
      path, x_col, y_col, z_col, qv_col, bounds, options.min_qv);
  emit_progress(
      "Selected " + std::to_string(row_groups.size()) +
      " transcript row groups for counts-only store");

  XeniumParquetCountsResult out;
  out.has_z = z_col.has_value();
  out.has_qv = qv_col.has_value();
  out.has_transcript_id = txid_col.has_value();
  if (row_groups.empty()) {
    return out;
  }

  auto batch_reader = arrow_unwrap(
      reader->GetRecordBatchReader(row_groups, projected),
      "Create Xenium transcript batch reader");
  std::shared_ptr<arrow::RecordBatch> batch;
  CellGeneCountAccumulator accumulator;
  const auto cell_filter = make_string_filter_set(options.cell_filter);
  const auto gene_filter = make_string_filter_set(options.gene_filter);
  int batch_count = 0;
  constexpr std::int64_t progress_row_interval = 1000000;
  std::int64_t next_progress_rows = progress_row_interval;
  while (true) {
    arrow_check(batch_reader->ReadNext(&batch), "Read Xenium transcript batch");
    if (!batch) break;
    ++batch_count;
    out.scanned_rows += batch->num_rows();

    int col = 0;
    NumericArrayView x_view(batch->column(col++));
    NumericArrayView y_view(batch->column(col++));
    std::optional<NumericArrayView> z_view;
    if (z_col.has_value()) z_view.emplace(batch->column(col++));
    std::optional<NumericArrayView> qv_view;
    if (qv_col.has_value()) qv_view.emplace(batch->column(col++));
    StringArrayView gene_view(batch->column(col++));
    StringArrayView cell_view(batch->column(col++));
    std::optional<StringArrayView> codeword_category_view;
    if (use_codeword_category_filter) codeword_category_view.emplace(batch->column(col++));
    std::optional<StringArrayView> is_gene_view;
    if (use_is_gene_filter) is_gene_view.emplace(batch->column(col++));

    for (int64_t row = 0; row < batch->num_rows(); ++row) {
      const double x = x_view.value(row);
      const double y = y_view.value(row);
      const double z = z_view.has_value() ? z_view->value(row) : 0.0;
      if (!passes_union_bounds(x, y, z, bounds)) {
        continue;
      }
      const double qv = qv_view.has_value() ? qv_view->value(row) : 0.0;
      if (options.min_qv >= 0.0 && qv_view.has_value() &&
          (!std::isfinite(qv) || qv < options.min_qv)) {
        continue;
      }
      const std::string crop_id = assign_crop_id(x, y, z, options.crops);
      if (!options.crops.empty() && crop_id.empty()) {
        continue;
      }
      const std::string cell_id = cell_view.value(row);
      if (!options.keep_unassigned && is_unassigned_cell_id(cell_id)) {
        continue;
      }
      if (!passes_string_filter(cell_id, cell_filter)) {
        continue;
      }
      const std::string gene = gene_view.value(row);
      if (gene.empty()) {
        continue;
      }
      if (!options.keep_non_gene &&
          !xenium_is_biological_feature(
              gene,
              codeword_category_view.has_value() ? codeword_category_view->value(row) : std::string(),
              is_gene_view.has_value() ? is_gene_view->value(row) : std::string())) {
        continue;
      }
      if (!passes_string_filter(gene, gene_filter)) {
        continue;
      }
      accumulator.add_transcript(cell_id, gene, x, y, z, crop_id);
      ++out.kept_rows;
    }

    if (options.progress && (out.scanned_rows >= next_progress_rows || batch_count == 1)) {
      emit_progress(
          "Scanned counts-only transcript batches: " + std::to_string(batch_count) +
          ", rows=" + std::to_string(out.scanned_rows) +
          ", kept=" + std::to_string(out.kept_rows));
      while (next_progress_rows <= out.scanned_rows) {
        next_progress_rows += progress_row_interval;
      }
    }
  }
  out.counts = accumulator.finalize();
  emit_progress(
      "Finished counts-only transcript scan: rows=" + std::to_string(out.scanned_rows) +
      ", kept=" + std::to_string(out.kept_rows) +
      ", cells=" + std::to_string(out.counts.cells.cell_ids.size()) +
      ", genes=" + std::to_string(out.counts.genes.size()));
  return out;
}

EncodedMoleculeStore build_xenium_encoded_store_from_parquet(
    const std::string& path,
    const XeniumLoadOptions& options) {
  auto stage_start = std::chrono::steady_clock::now();
  auto emit_progress = [&](const std::string& message) {
    if (!options.progress) return;
    const auto now = std::chrono::steady_clock::now();
    options.progress(message, std::chrono::duration<double>(now - stage_start).count());
    stage_start = now;
  };

  emit_progress("Opening Xenium transcripts parquet for full input store: " + path);
  auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path), "Open Xenium transcript parquet");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open Xenium transcript parquet reader");
  auto reader = arrow_unwrap(builder.Build(), "Build Xenium transcript parquet reader");
  reader->set_use_threads(true);
  reader->set_batch_size(65536);

  std::shared_ptr<arrow::Schema> schema;
  arrow_check(reader->GetSchema(&schema), "Read Xenium transcript parquet schema");
  const auto gene_col = find_first_field(schema, {"feature_name", "gene", "gene_name", "gene_id"}, "gene");
  const auto cell_col = find_first_field(schema, {"cell_id", "cell"}, "cell id");
  const auto x_col = find_first_field(schema, {"x_location", "x"}, "x");
  const auto y_col = find_first_field(schema, {"y_location", "y"}, "y");
  const auto z_col = find_optional_field(schema, {"z_location", "z"});
  const auto qv_col = find_optional_field(schema, {"qv"});
  const auto overlaps_col = find_optional_field(schema, {"overlaps_nucleus"});
  const auto nucleus_distance_col = find_optional_field(schema, {"nucleus_distance"});
  const auto codeword_category_col = find_optional_field(schema, {"codeword_category"});
  const auto is_gene_col = find_optional_field(schema, {"is_gene"});
  const bool use_codeword_category_filter = !options.keep_non_gene && codeword_category_col.has_value();
  const bool use_is_gene_filter = !options.keep_non_gene && !use_codeword_category_filter && is_gene_col.has_value();
  if (overlaps_col.has_value() || nucleus_distance_col.has_value()) {
    std::string message = "Full input store will preserve Xenium nucleus columns:";
    if (overlaps_col.has_value()) message += " overlaps_nucleus";
    if (nucleus_distance_col.has_value()) message += " nucleus_distance";
    emit_progress(message);
  }
  const UnionBounds bounds = union_bounds(options.crops);

  std::vector<int> projected = {static_cast<int>(x_col), static_cast<int>(y_col)};
  if (z_col.has_value()) projected.push_back(static_cast<int>(*z_col));
  if (qv_col.has_value()) projected.push_back(static_cast<int>(*qv_col));
  if (overlaps_col.has_value()) projected.push_back(static_cast<int>(*overlaps_col));
  if (nucleus_distance_col.has_value()) projected.push_back(static_cast<int>(*nucleus_distance_col));
  projected.push_back(static_cast<int>(gene_col));
  projected.push_back(static_cast<int>(cell_col));
  if (use_codeword_category_filter) projected.push_back(static_cast<int>(*codeword_category_col));
  if (use_is_gene_filter) projected.push_back(static_cast<int>(*is_gene_col));

  const auto row_groups = select_row_groups_for_numeric_filters(
      path, x_col, y_col, z_col, qv_col, bounds, options.min_qv);
  emit_progress(
      "Selected " + std::to_string(row_groups.size()) +
      " transcript row groups for full input store");

  EncodedMoleculeStore out;
  out.has_z = z_col.has_value();
  out.has_qv = qv_col.has_value();
  out.has_overlaps_nucleus = overlaps_col.has_value();
  out.has_nucleus_distance = nucleus_distance_col.has_value();
  // The high-performance store keeps a stable integer obs_id and avoids
  // materializing per-molecule transcript_id strings.
  out.has_transcript_id = false;
  if (row_groups.empty()) {
    return out;
  }

  auto batch_reader = arrow_unwrap(
      reader->GetRecordBatchReader(row_groups, projected),
      "Create Xenium transcript batch reader");
  std::shared_ptr<arrow::RecordBatch> batch;
  EncodedMoleculeAccumulator accumulator(
      qv_col.has_value(),
      overlaps_col.has_value(),
      nucleus_distance_col.has_value());
  const auto cell_filter = make_string_filter_set(options.cell_filter);
  const auto gene_filter = make_string_filter_set(options.gene_filter);
  int batch_count = 0;
  constexpr std::int64_t progress_row_interval = 1000000;
  std::int64_t next_progress_rows = progress_row_interval;
  while (true) {
    arrow_check(batch_reader->ReadNext(&batch), "Read Xenium transcript batch");
    if (!batch) break;
    ++batch_count;
    out.scanned_rows += batch->num_rows();

    int col = 0;
    NumericArrayView x_view(batch->column(col++));
    NumericArrayView y_view(batch->column(col++));
    std::optional<NumericArrayView> z_view;
    if (z_col.has_value()) z_view.emplace(batch->column(col++));
    std::optional<NumericArrayView> qv_view;
    if (qv_col.has_value()) qv_view.emplace(batch->column(col++));
    std::optional<StringArrayView> overlaps_view;
    if (overlaps_col.has_value()) overlaps_view.emplace(batch->column(col++));
    std::optional<NumericArrayView> nucleus_distance_view;
    if (nucleus_distance_col.has_value()) nucleus_distance_view.emplace(batch->column(col++));
    StringArrayView gene_view(batch->column(col++));
    StringArrayView cell_view(batch->column(col++));
    std::optional<StringArrayView> codeword_category_view;
    if (use_codeword_category_filter) codeword_category_view.emplace(batch->column(col++));
    std::optional<StringArrayView> is_gene_view;
    if (use_is_gene_filter) is_gene_view.emplace(batch->column(col++));

    for (int64_t row = 0; row < batch->num_rows(); ++row) {
      const double x = x_view.value(row);
      const double y = y_view.value(row);
      const double z = z_view.has_value() ? z_view->value(row) : 0.0;
      if (!passes_union_bounds(x, y, z, bounds)) {
        continue;
      }
      const double qv = qv_view.has_value() ? qv_view->value(row) : 0.0;
      if (options.min_qv >= 0.0 && qv_view.has_value() &&
          (!std::isfinite(qv) || qv < options.min_qv)) {
        continue;
      }
      const std::string crop_id = assign_crop_id(x, y, z, options.crops);
      if (!options.crops.empty() && crop_id.empty()) {
        continue;
      }
      const std::string cell_id = cell_view.value(row);
      if (!options.keep_unassigned && is_unassigned_cell_id(cell_id)) {
        continue;
      }
      if (!passes_string_filter(cell_id, cell_filter)) {
        continue;
      }
      const std::string gene = gene_view.value(row);
      if (gene.empty()) {
        continue;
      }
      if (!options.keep_non_gene &&
          !xenium_is_biological_feature(
              gene,
              codeword_category_view.has_value() ? codeword_category_view->value(row) : std::string(),
              is_gene_view.has_value() ? is_gene_view->value(row) : std::string())) {
        continue;
      }
      if (!passes_string_filter(gene, gene_filter)) {
        continue;
      }
      accumulator.add_transcript(
          cell_id,
          gene,
          x,
          y,
          z,
          qv,
          overlaps_view.has_value() ? parse_bool_or_missing(overlaps_view->value(row)) : -1,
          nucleus_distance_view.has_value()
              ? nucleus_distance_view->value(row)
              : std::numeric_limits<double>::quiet_NaN(),
          crop_id);
      ++out.kept_rows;
    }

    if (options.progress && (out.scanned_rows >= next_progress_rows || batch_count == 1)) {
      emit_progress(
          "Scanned full-store transcript batches: " + std::to_string(batch_count) +
          ", rows=" + std::to_string(out.scanned_rows) +
          ", kept=" + std::to_string(out.kept_rows));
      while (next_progress_rows <= out.scanned_rows) {
        next_progress_rows += progress_row_interval;
      }
    }
  }

  auto encoded = accumulator.finalize();
  encoded.has_z = out.has_z;
  encoded.has_qv = out.has_qv;
  encoded.has_transcript_id = out.has_transcript_id;
  encoded.has_overlaps_nucleus = out.has_overlaps_nucleus;
  encoded.has_nucleus_distance = out.has_nucleus_distance;
  encoded.scanned_rows = out.scanned_rows;
  encoded.kept_rows = out.kept_rows;
  emit_progress(
      "Finished full-store transcript scan: rows=" + std::to_string(encoded.scanned_rows) +
      ", kept=" + std::to_string(encoded.kept_rows) +
      ", cells=" + std::to_string(encoded.counts.cells.cell_ids.size()) +
      ", genes=" + std::to_string(encoded.counts.genes.size()));
  return encoded;
}

struct TabularStreamSource {
  std::shared_ptr<arrow::RecordBatchReader> reader;
  bool used_parquet = false;
  bool has_z = false;
  bool has_qv = false;
  bool has_cell_type = false;
  bool has_sample_id = false;
  bool has_fov_id = false;
};

int required_named_field(
    const std::shared_ptr<arrow::Schema>& schema,
    const std::string& name,
    const std::string& label) {
  const int idx = schema->GetFieldIndex(name);
  if (idx < 0) {
    throw std::runtime_error("Column '" + name + "' not found in tabular input for " + label);
  }
  return idx;
}

std::optional<int> optional_named_field(
    const std::shared_ptr<arrow::Schema>& schema,
    const std::string& name) {
  if (name.empty()) {
    return std::nullopt;
  }
  const int idx = schema->GetFieldIndex(name);
  if (idx < 0) {
    return std::nullopt;
  }
  return idx;
}

void add_projected_field(std::vector<int>& projected, int idx) {
  if (idx < 0) {
    return;
  }
  if (std::find(projected.begin(), projected.end(), idx) == projected.end()) {
    projected.push_back(idx);
  }
}

bool can_stream_tabular_input(const TabularSourceSpec& source) {
  const bool has_cell_id_col = !source.cell_id_col.empty();
  const bool has_mask = !source.segmentation_mask_path.empty();
  if (has_cell_id_col == has_mask) {
    return false;
  }
  const std::string ext = logical_file_extension_lower(source.molecules_path);
  return is_parquet_extension(ext) || is_delimited_text_extension(ext);
}

TabularStreamSource make_tabular_stream_source(
    const TabularSourceSpec& source,
    const TabularLoadOptions& options) {
  const bool has_cell_id_col = !source.cell_id_col.empty();
  const bool has_mask = !source.segmentation_mask_path.empty();
  if (has_cell_id_col == has_mask) {
    throw std::runtime_error(has_cell_id_col
        ? "Specify either cell_id_col or segmentation_mask_path, not both"
        : "Tabular input requires either cell_id_col or segmentation_mask_path");
  }
  if (source.gene_col.empty() || source.x_col.empty() || source.y_col.empty()) {
    throw std::runtime_error("x_col, y_col, and gene_col must be provided for streaming tabular stores");
  }

  const std::string ext = logical_file_extension_lower(source.molecules_path);
  const UnionBounds bounds = union_bounds(options.crops);
  TabularStreamSource out;

  if (is_parquet_extension(ext)) {
    out.used_parquet = true;
    auto input = arrow_unwrap(arrow::io::ReadableFile::Open(source.molecules_path), "Open tabular parquet input");
    parquet::arrow::FileReaderBuilder builder;
    arrow_check(builder.Open(input), "Open tabular parquet reader");
    auto reader = arrow_unwrap(builder.Build(), "Build tabular parquet reader");
    reader->set_use_threads(true);
    reader->set_batch_size(65536);

    std::shared_ptr<arrow::Schema> schema;
    arrow_check(reader->GetSchema(&schema), "Read tabular parquet schema");
    const int x_col = required_named_field(schema, source.x_col, "x");
    const int y_col = required_named_field(schema, source.y_col, "y");
    const int gene_col = required_named_field(schema, source.gene_col, "gene");
    const auto cell_col = has_cell_id_col
        ? std::optional<int>(required_named_field(schema, source.cell_id_col, "cell id"))
        : std::nullopt;
    const auto z_col = optional_named_field(schema, source.z_col);
    const auto qv_col = optional_named_field(schema, source.qv_col);
    const auto cell_type_col = optional_named_field(schema, source.cell_type_col);
    const auto sample_id_col = optional_named_field(schema, source.sample_id_col);
    const auto fov_id_col = optional_named_field(schema, source.fov_id_col);

    std::vector<int> projected;
    add_projected_field(projected, x_col);
    add_projected_field(projected, y_col);
    if (z_col.has_value()) add_projected_field(projected, *z_col);
    if (qv_col.has_value()) add_projected_field(projected, *qv_col);
    add_projected_field(projected, gene_col);
    if (cell_col.has_value()) add_projected_field(projected, *cell_col);
    if (cell_type_col.has_value()) add_projected_field(projected, *cell_type_col);
    if (sample_id_col.has_value()) add_projected_field(projected, *sample_id_col);
    if (fov_id_col.has_value()) add_projected_field(projected, *fov_id_col);

    const auto row_groups = select_row_groups_for_numeric_filters(
        source.molecules_path,
        static_cast<std::size_t>(x_col),
        static_cast<std::size_t>(y_col),
        z_col.has_value() ? std::optional<std::size_t>(static_cast<std::size_t>(*z_col)) : std::nullopt,
        qv_col.has_value() ? std::optional<std::size_t>(static_cast<std::size_t>(*qv_col)) : std::nullopt,
        bounds,
        options.min_qv);
    out.has_z = z_col.has_value();
    out.has_qv = qv_col.has_value();
    out.has_cell_type = cell_type_col.has_value();
    out.has_sample_id = sample_id_col.has_value() || !source.sample_id.empty();
    out.has_fov_id = fov_id_col.has_value() || !source.fov_id.empty();
    if (row_groups.empty()) {
      return out;
    }
    out.reader = arrow_unwrap(
        reader->GetRecordBatchReader(row_groups, projected),
        "Create tabular parquet batch reader");
    return out;
  }

  if (is_delimited_text_extension(ext)) {
    auto input = open_delimited_input_stream(source.molecules_path, "Open tabular delimited input");
    auto read_options = arrow::csv::ReadOptions::Defaults();
    auto parse_options = arrow::csv::ParseOptions::Defaults();
    auto convert_options = arrow::csv::ConvertOptions::Defaults();
    read_options.block_size = 1 << 20;
    if (ext == "tsv") {
      parse_options.delimiter = '\t';
    }
    convert_options.column_types[source.x_col] = arrow::float64();
    convert_options.column_types[source.y_col] = arrow::float64();
    if (!source.z_col.empty()) convert_options.column_types[source.z_col] = arrow::float64();
    if (!source.qv_col.empty()) convert_options.column_types[source.qv_col] = arrow::float64();
    convert_options.column_types[source.gene_col] = arrow::utf8();
    if (has_cell_id_col) convert_options.column_types[source.cell_id_col] = arrow::utf8();
    if (!source.cell_type_col.empty()) convert_options.column_types[source.cell_type_col] = arrow::utf8();
    if (!source.sample_id_col.empty()) convert_options.column_types[source.sample_id_col] = arrow::utf8();
    if (!source.fov_id_col.empty()) convert_options.column_types[source.fov_id_col] = arrow::utf8();
    auto csv_reader = arrow_unwrap(
        arrow::csv::StreamingReader::Make(
            arrow::io::default_io_context(),
            input,
            read_options,
            parse_options,
            convert_options),
        "Create tabular delimited batch reader");
    const auto schema = csv_reader->schema();
    required_named_field(schema, source.x_col, "x");
    required_named_field(schema, source.y_col, "y");
    required_named_field(schema, source.gene_col, "gene");
    if (has_cell_id_col) {
      required_named_field(schema, source.cell_id_col, "cell id");
    }
    out.has_z = optional_named_field(schema, source.z_col).has_value();
    out.has_qv = optional_named_field(schema, source.qv_col).has_value();
    out.has_cell_type = optional_named_field(schema, source.cell_type_col).has_value();
    out.has_sample_id = optional_named_field(schema, source.sample_id_col).has_value() || !source.sample_id.empty();
    out.has_fov_id = optional_named_field(schema, source.fov_id_col).has_value() || !source.fov_id.empty();
    out.reader = csv_reader;
    return out;
  }

  throw std::runtime_error("Unsupported tabular file format for '" + source.molecules_path + "'");
}

struct TabularStreamCountsResult {
  CellCountMatrix counts;
  bool has_z = false;
  bool has_qv = false;
  bool has_transcript_id = false;
  std::int64_t scanned_rows = 0;
  std::int64_t kept_rows = 0;
};

TabularStreamCountsResult build_tabular_counts_from_stream(
    const TabularSourceSpec& source,
    const TabularLoadOptions& options) {
  const auto stream = make_tabular_stream_source(source, options);
  TabularStreamCountsResult out;
  out.has_z = stream.has_z;
  out.has_qv = stream.has_qv;
  out.has_transcript_id = false;
  if (!stream.reader) {
    return out;
  }

  const auto schema = stream.reader->schema();
  const int x_col = required_named_field(schema, source.x_col, "x");
  const int y_col = required_named_field(schema, source.y_col, "y");
  const int gene_col = required_named_field(schema, source.gene_col, "gene");
  const auto cell_col = !source.cell_id_col.empty()
      ? std::optional<int>(required_named_field(schema, source.cell_id_col, "cell id"))
      : std::nullopt;
  const auto z_col = optional_named_field(schema, source.z_col);
  const auto qv_col = optional_named_field(schema, source.qv_col);
  const auto cell_type_col = optional_named_field(schema, source.cell_type_col);
  const auto sample_id_col = optional_named_field(schema, source.sample_id_col);
  const auto fov_id_col = optional_named_field(schema, source.fov_id_col);
  if (cell_type_col.has_value() && !source.cell_metadata_path.empty() &&
      !source.cell_metadata_cell_type_col.empty()) {
    throw std::runtime_error("Specify either molecule-level cell_type_col or cell metadata sidecar, not both");
  }
  const auto cell_type_by_cell = read_cell_type_metadata(source);
  const UnionBounds bounds = union_bounds(options.crops);
  std::optional<TabularTiffMask> mask;
  if (!source.segmentation_mask_path.empty()) {
    mask = read_labeled_tiff_mask(source.segmentation_mask_path);
  }

  CellGeneCountAccumulator accumulator;
  std::shared_ptr<arrow::RecordBatch> batch;
  while (true) {
    arrow_check(stream.reader->ReadNext(&batch), "Read tabular count batch");
    if (!batch) break;
    out.scanned_rows += batch->num_rows();

    NumericArrayView x_view(batch->column(x_col));
    NumericArrayView y_view(batch->column(y_col));
    StringArrayView gene_view(batch->column(gene_col));
    std::optional<StringArrayView> cell_view;
    if (cell_col.has_value()) cell_view.emplace(batch->column(*cell_col));
    std::optional<NumericArrayView> z_view;
    if (z_col.has_value()) z_view.emplace(batch->column(*z_col));
    std::optional<NumericArrayView> qv_view;
    if (qv_col.has_value()) qv_view.emplace(batch->column(*qv_col));
    std::optional<StringArrayView> cell_type_view;
    if (cell_type_col.has_value()) cell_type_view.emplace(batch->column(*cell_type_col));
    std::optional<StringArrayView> sample_id_view;
    if (sample_id_col.has_value()) sample_id_view.emplace(batch->column(*sample_id_col));
    std::optional<StringArrayView> fov_id_view;
    if (fov_id_col.has_value()) fov_id_view.emplace(batch->column(*fov_id_col));

    for (int64_t row = 0; row < batch->num_rows(); ++row) {
      const double x = x_view.value(row);
      const double y = y_view.value(row);
      const double z = z_view.has_value() ? z_view->value(row) : 0.0;
      if (!passes_union_bounds(x, y, z, bounds)) {
        continue;
      }
      const double qv = qv_view.has_value() ? qv_view->value(row) : 0.0;
      if (options.min_qv >= 0.0 && qv_view.has_value() &&
          (!std::isfinite(qv) || qv < options.min_qv)) {
        continue;
      }
      const std::string crop_id = assign_crop_id(x, y, z, options.crops);
      if (!options.crops.empty() && crop_id.empty()) {
        continue;
      }
      std::string cell_id = cell_view.has_value()
          ? cell_view->value(row)
          : assign_cell_id_from_mask(x, y, *mask);
      if (is_unassigned_cell_id(cell_id)) {
        if (!options.keep_unassigned) {
          continue;
        }
        cell_id = "__unassigned__";
      }
      const std::string sample_id =
          sample_id_view.has_value() ? sample_id_view->value(row) : source.sample_id;
      const std::string fov_id =
          fov_id_view.has_value() ? fov_id_view->value(row) : source.fov_id;
      const std::string cell_type = cell_type_view.has_value()
          ? cell_type_view->value(row)
          : cell_type_for_cell(cell_type_by_cell, cell_id);
      accumulator.add_transcript(
          cell_id,
          gene_view.value(row),
          x,
          y,
          z,
          crop_id,
          cell_type,
          sample_id,
          fov_id);
      ++out.kept_rows;
    }
  }
  out.counts = accumulator.finalize();
  return out;
}

EncodedMoleculeStore build_tabular_encoded_store_from_stream(
    const TabularSourceSpec& source,
    const TabularLoadOptions& options) {
  const auto stream = make_tabular_stream_source(source, options);
  EncodedMoleculeStore out;
  out.has_z = stream.has_z;
  out.has_qv = stream.has_qv;
  out.has_transcript_id = false;
  if (!stream.reader) {
    return out;
  }

  const auto schema = stream.reader->schema();
  const int x_col = required_named_field(schema, source.x_col, "x");
  const int y_col = required_named_field(schema, source.y_col, "y");
  const int gene_col = required_named_field(schema, source.gene_col, "gene");
  const auto cell_col = !source.cell_id_col.empty()
      ? std::optional<int>(required_named_field(schema, source.cell_id_col, "cell id"))
      : std::nullopt;
  const auto z_col = optional_named_field(schema, source.z_col);
  const auto qv_col = optional_named_field(schema, source.qv_col);
  const auto cell_type_col = optional_named_field(schema, source.cell_type_col);
  const auto sample_id_col = optional_named_field(schema, source.sample_id_col);
  const auto fov_id_col = optional_named_field(schema, source.fov_id_col);
  if (cell_type_col.has_value() && !source.cell_metadata_path.empty() &&
      !source.cell_metadata_cell_type_col.empty()) {
    throw std::runtime_error("Specify either molecule-level cell_type_col or cell metadata sidecar, not both");
  }
  const auto cell_type_by_cell = read_cell_type_metadata(source);
  const UnionBounds bounds = union_bounds(options.crops);
  std::optional<TabularTiffMask> mask;
  if (!source.segmentation_mask_path.empty()) {
    mask = read_labeled_tiff_mask(source.segmentation_mask_path);
  }

  EncodedMoleculeAccumulator accumulator(
      qv_col.has_value(),
      /*store_overlaps_nucleus=*/false,
      /*store_nucleus_distance=*/false);
  std::shared_ptr<arrow::RecordBatch> batch;
  while (true) {
    arrow_check(stream.reader->ReadNext(&batch), "Read tabular molecule batch");
    if (!batch) break;
    out.scanned_rows += batch->num_rows();

    NumericArrayView x_view(batch->column(x_col));
    NumericArrayView y_view(batch->column(y_col));
    StringArrayView gene_view(batch->column(gene_col));
    std::optional<StringArrayView> cell_view;
    if (cell_col.has_value()) cell_view.emplace(batch->column(*cell_col));
    std::optional<NumericArrayView> z_view;
    if (z_col.has_value()) z_view.emplace(batch->column(*z_col));
    std::optional<NumericArrayView> qv_view;
    if (qv_col.has_value()) qv_view.emplace(batch->column(*qv_col));
    std::optional<StringArrayView> cell_type_view;
    if (cell_type_col.has_value()) cell_type_view.emplace(batch->column(*cell_type_col));
    std::optional<StringArrayView> sample_id_view;
    if (sample_id_col.has_value()) sample_id_view.emplace(batch->column(*sample_id_col));
    std::optional<StringArrayView> fov_id_view;
    if (fov_id_col.has_value()) fov_id_view.emplace(batch->column(*fov_id_col));

    for (int64_t row = 0; row < batch->num_rows(); ++row) {
      const double x = x_view.value(row);
      const double y = y_view.value(row);
      const double z = z_view.has_value() ? z_view->value(row) : 0.0;
      if (!passes_union_bounds(x, y, z, bounds)) {
        continue;
      }
      const double qv = qv_view.has_value() ? qv_view->value(row) : 0.0;
      if (options.min_qv >= 0.0 && qv_view.has_value() &&
          (!std::isfinite(qv) || qv < options.min_qv)) {
        continue;
      }
      const std::string crop_id = assign_crop_id(x, y, z, options.crops);
      if (!options.crops.empty() && crop_id.empty()) {
        continue;
      }
      std::string cell_id = cell_view.has_value()
          ? cell_view->value(row)
          : assign_cell_id_from_mask(x, y, *mask);
      if (is_unassigned_cell_id(cell_id)) {
        if (!options.keep_unassigned) {
          continue;
        }
        cell_id = "__unassigned__";
      }
      const std::string sample_id =
          sample_id_view.has_value() ? sample_id_view->value(row) : source.sample_id;
      const std::string fov_id =
          fov_id_view.has_value() ? fov_id_view->value(row) : source.fov_id;
      const std::string cell_type = cell_type_view.has_value()
          ? cell_type_view->value(row)
          : cell_type_for_cell(cell_type_by_cell, cell_id);
      accumulator.add_transcript(
          cell_id,
          gene_view.value(row),
          x,
          y,
          z,
          qv,
          -1,
          std::numeric_limits<double>::quiet_NaN(),
          crop_id,
          cell_type,
          sample_id,
          fov_id);
      ++out.kept_rows;
    }
  }

  auto encoded = accumulator.finalize();
  encoded.has_z = out.has_z;
  encoded.has_qv = out.has_qv;
  encoded.has_transcript_id = out.has_transcript_id;
  encoded.scanned_rows = out.scanned_rows;
  encoded.kept_rows = out.kept_rows;
  return encoded;
}

void write_genes(const std::filesystem::path& path, const std::vector<std::string>& genes, int row_group_size) {
  std::vector<int> gene_idx(genes.size());
  std::iota(gene_idx.begin(), gene_idx.end(), 0);
  const auto table = arrow::Table::Make(
      arrow::schema({arrow::field("gene_idx", arrow::int32()), arrow::field("gene", arrow::utf8())}),
      {build_int32_array(gene_idx), build_string_array(genes)});
  write_parquet_table(table, path, row_group_size);
}

void write_cells(
    const std::filesystem::path& path,
    const CellCountMatrix& counts,
    int row_group_size) {
  const int n = static_cast<int>(counts.cells.cell_ids.size());
  std::vector<int> cell_idx(static_cast<std::size_t>(n));
  std::iota(cell_idx.begin(), cell_idx.end(), 0);
  std::vector<std::shared_ptr<arrow::Field>> fields = {
      arrow::field("cell_idx", arrow::int32()),
      arrow::field("cell_id", arrow::utf8()),
      arrow::field("transcript_count", arrow::int32()),
      arrow::field("centroid_x", arrow::float64()),
      arrow::field("centroid_y", arrow::float64()),
      arrow::field("centroid_z", arrow::float64()),
      arrow::field("crop_id", arrow::utf8())};
  std::vector<std::shared_ptr<arrow::Array>> arrays = {
      build_int32_array(cell_idx),
      build_string_array(counts.cells.cell_ids),
      build_int32_array(counts.transcript_counts),
      build_double_array(counts.cells.centroid_x),
      build_double_array(counts.cells.centroid_y),
      build_double_array(counts.cells.centroid_z),
      build_string_array(counts.crop_ids)};
  if (!counts.cells.cell_types.empty()) {
    fields.push_back(arrow::field("cell_type", arrow::utf8()));
    arrays.push_back(build_string_array(counts.cells.cell_types));
  }
  if (!counts.cells.sample_ids.empty()) {
    fields.push_back(arrow::field("sample_id", arrow::utf8()));
    arrays.push_back(build_string_array(counts.cells.sample_ids));
  }
  if (!counts.cells.fov_ids.empty()) {
    fields.push_back(arrow::field("fov_id", arrow::utf8()));
    arrays.push_back(build_string_array(counts.cells.fov_ids));
  }
  write_parquet_table(arrow::Table::Make(std::make_shared<arrow::Schema>(fields), arrays), path, row_group_size);
}

void write_counts(const std::filesystem::path& path, const CellCountMatrix& counts, int row_group_size) {
  std::vector<int> cell_idx;
  std::vector<int> gene_idx;
  std::vector<int> count;
  cell_idx.reserve(counts.values.size());
  gene_idx.reserve(counts.values.size());
  count.reserve(counts.values.size());
  for (int cell = 0; cell < static_cast<int>(counts.transcript_counts.size()); ++cell) {
    for (int p = counts.indptr[static_cast<std::size_t>(cell)];
         p < counts.indptr[static_cast<std::size_t>(cell + 1)];
         ++p) {
      cell_idx.push_back(cell);
      gene_idx.push_back(counts.indices[static_cast<std::size_t>(p)]);
      count.push_back(static_cast<int>(counts.values[static_cast<std::size_t>(p)]));
    }
  }
  const auto table = arrow::Table::Make(
      arrow::schema({
          arrow::field("cell_idx", arrow::int32()),
          arrow::field("gene_idx", arrow::int32()),
          arrow::field("count", arrow::int32())}),
      {build_int32_array(cell_idx), build_int32_array(gene_idx), build_int32_array(count)});
  write_parquet_table(table, path, row_group_size);
}

void write_molecules(
    const std::filesystem::path& path,
    const std::filesystem::path& offsets_path,
    const TranscriptTable& table,
    const std::vector<std::string>* transcript_crop_ids,
    int row_group_size) {
  std::vector<std::size_t> order(table.size());
  std::iota(order.begin(), order.end(), 0U);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
    if (table.cell_index[lhs] != table.cell_index[rhs]) return table.cell_index[lhs] < table.cell_index[rhs];
    return lhs < rhs;
  });

  std::vector<std::int64_t> obs_id;
  std::vector<int> cell_idx;
  std::vector<int> gene_idx;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> qv;
  std::vector<int> overlaps_nucleus;
  std::vector<double> nucleus_distance;
  std::vector<std::string> transcript_id;
  std::vector<std::string> crop_id;
  obs_id.reserve(order.size());
  cell_idx.reserve(order.size());
  gene_idx.reserve(order.size());
  x.reserve(order.size());
  y.reserve(order.size());
  z.reserve(order.size());
  if (!table.qv.empty()) qv.reserve(order.size());
  if (!table.overlaps_nucleus.empty()) overlaps_nucleus.reserve(order.size());
  if (!table.nucleus_distance.empty()) nucleus_distance.reserve(order.size());
  if (!table.transcript_ids.empty()) transcript_id.reserve(order.size());
  if (transcript_crop_ids != nullptr && transcript_crop_ids->size() == table.size()) crop_id.reserve(order.size());

  std::vector<int> offset_cell_idx(table.num_cells());
  std::vector<std::int64_t> offset_start(table.num_cells(), -1);
  std::vector<std::int64_t> offset_end(table.num_cells(), -1);
  std::iota(offset_cell_idx.begin(), offset_cell_idx.end(), 0);

  for (std::size_t sorted = 0; sorted < order.size(); ++sorted) {
    const std::size_t row = order[sorted];
    const int cell = table.cell_index[row];
    if (offset_start[static_cast<std::size_t>(cell)] < 0) {
      offset_start[static_cast<std::size_t>(cell)] = static_cast<std::int64_t>(sorted);
    }
    offset_end[static_cast<std::size_t>(cell)] = static_cast<std::int64_t>(sorted + 1);
    obs_id.push_back(static_cast<std::int64_t>(row));
    cell_idx.push_back(cell);
    gene_idx.push_back(table.gene_index[row]);
    x.push_back(table.x[row]);
    y.push_back(table.y[row]);
    z.push_back(table.z[row]);
    if (!table.qv.empty()) qv.push_back(table.qv[row]);
    if (!table.overlaps_nucleus.empty()) overlaps_nucleus.push_back(table.overlaps_nucleus[row]);
    if (!table.nucleus_distance.empty()) nucleus_distance.push_back(table.nucleus_distance[row]);
    if (!table.transcript_ids.empty()) transcript_id.push_back(table.transcript_ids[row]);
    if (transcript_crop_ids != nullptr && transcript_crop_ids->size() == table.size()) {
      crop_id.push_back((*transcript_crop_ids)[row]);
    }
  }

  std::vector<std::shared_ptr<arrow::Field>> fields = {
      arrow::field("obs_id", arrow::int64()),
      arrow::field("cell_idx", arrow::int32()),
      arrow::field("gene_idx", arrow::int32()),
      arrow::field("x", arrow::float64()),
      arrow::field("y", arrow::float64()),
      arrow::field("z", arrow::float64())};
  std::vector<std::shared_ptr<arrow::Array>> arrays = {
      build_int64_array(obs_id),
      build_int32_array(cell_idx),
      build_int32_array(gene_idx),
      build_double_array(x),
      build_double_array(y),
      build_double_array(z)};
  if (!qv.empty()) {
    fields.push_back(arrow::field("qv", arrow::float64()));
    arrays.push_back(build_double_array(qv));
  }
  if (!overlaps_nucleus.empty()) {
    fields.push_back(arrow::field("overlaps_nucleus", arrow::int32()));
    arrays.push_back(build_int32_array(overlaps_nucleus));
  }
  if (!nucleus_distance.empty()) {
    fields.push_back(arrow::field("nucleus_distance", arrow::float64()));
    arrays.push_back(build_double_array(nucleus_distance));
  }
  if (!transcript_id.empty()) {
    fields.push_back(arrow::field("transcript_id", arrow::utf8()));
    arrays.push_back(build_string_array(transcript_id));
  }
  if (!crop_id.empty()) {
    fields.push_back(arrow::field("crop_id", arrow::utf8()));
    arrays.push_back(build_string_array(crop_id));
  }
  write_parquet_table(arrow::Table::Make(std::make_shared<arrow::Schema>(fields), arrays), path, row_group_size);

  const auto offsets = arrow::Table::Make(
      arrow::schema({
          arrow::field("cell_idx", arrow::int32()),
          arrow::field("start", arrow::int64()),
          arrow::field("end", arrow::int64())}),
      {build_int32_array(offset_cell_idx), build_int64_array(offset_start), build_int64_array(offset_end)});
  write_parquet_table(offsets, offsets_path, row_group_size);
}

void write_encoded_molecules(
    const std::filesystem::path& path,
    const std::filesystem::path& offsets_path,
    const EncodedMoleculeStore& store,
    int row_group_size) {
  const std::size_t n = store.obs_id.size();
  const bool has_crop_id =
      std::any_of(store.crop_id.begin(), store.crop_id.end(), [](const std::string& value) {
        return !value.empty();
      });
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), 0U);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
    if (store.cell_idx[lhs] != store.cell_idx[rhs]) return store.cell_idx[lhs] < store.cell_idx[rhs];
    return store.obs_id[lhs] < store.obs_id[rhs];
  });

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
  obs_id.reserve(n);
  cell_idx.reserve(n);
  gene_idx.reserve(n);
  x.reserve(n);
  y.reserve(n);
  z.reserve(n);
  if (!store.qv.empty()) qv.reserve(n);
  if (!store.overlaps_nucleus.empty()) overlaps_nucleus.reserve(n);
  if (!store.nucleus_distance.empty()) nucleus_distance.reserve(n);
  if (has_crop_id) crop_id.reserve(n);

  std::vector<int> offset_cell_idx(store.counts.cells.cell_ids.size());
  std::vector<std::int64_t> offset_start(store.counts.cells.cell_ids.size(), -1);
  std::vector<std::int64_t> offset_end(store.counts.cells.cell_ids.size(), -1);
  std::iota(offset_cell_idx.begin(), offset_cell_idx.end(), 0);

  for (std::size_t sorted = 0; sorted < order.size(); ++sorted) {
    const std::size_t row = order[sorted];
    const int cell = store.cell_idx[row];
    if (offset_start[static_cast<std::size_t>(cell)] < 0) {
      offset_start[static_cast<std::size_t>(cell)] = static_cast<std::int64_t>(sorted);
    }
    offset_end[static_cast<std::size_t>(cell)] = static_cast<std::int64_t>(sorted + 1);
    obs_id.push_back(store.obs_id[row]);
    cell_idx.push_back(cell);
    gene_idx.push_back(store.gene_idx[row]);
    x.push_back(store.x[row]);
    y.push_back(store.y[row]);
    z.push_back(store.z[row]);
    if (!store.qv.empty()) qv.push_back(store.qv[row]);
    if (!store.overlaps_nucleus.empty()) overlaps_nucleus.push_back(store.overlaps_nucleus[row]);
    if (!store.nucleus_distance.empty()) nucleus_distance.push_back(store.nucleus_distance[row]);
    if (has_crop_id) crop_id.push_back(store.crop_id[row]);
  }

  std::vector<std::shared_ptr<arrow::Field>> fields = {
      arrow::field("obs_id", arrow::int64()),
      arrow::field("cell_idx", arrow::int32()),
      arrow::field("gene_idx", arrow::int32()),
      arrow::field("x", arrow::float64()),
      arrow::field("y", arrow::float64()),
      arrow::field("z", arrow::float64())};
  std::vector<std::shared_ptr<arrow::Array>> arrays = {
      build_int64_array(obs_id),
      build_int32_array(cell_idx),
      build_int32_array(gene_idx),
      build_double_array(x),
      build_double_array(y),
      build_double_array(z)};
  if (!qv.empty()) {
    fields.push_back(arrow::field("qv", arrow::float64()));
    arrays.push_back(build_double_array(qv));
  }
  if (!overlaps_nucleus.empty()) {
    fields.push_back(arrow::field("overlaps_nucleus", arrow::int32()));
    arrays.push_back(build_int32_array(overlaps_nucleus));
  }
  if (!nucleus_distance.empty()) {
    fields.push_back(arrow::field("nucleus_distance", arrow::float64()));
    arrays.push_back(build_double_array(nucleus_distance));
  }
  if (has_crop_id) {
    fields.push_back(arrow::field("crop_id", arrow::utf8()));
    arrays.push_back(build_string_array(crop_id));
  }
  write_parquet_table(arrow::Table::Make(std::make_shared<arrow::Schema>(fields), arrays), path, row_group_size);

  const auto offsets = arrow::Table::Make(
      arrow::schema({
          arrow::field("cell_idx", arrow::int32()),
          arrow::field("start", arrow::int64()),
          arrow::field("end", arrow::int64())}),
      {build_int32_array(offset_cell_idx), build_int64_array(offset_start), build_int64_array(offset_end)});
  write_parquet_table(offsets, offsets_path, row_group_size);
}

nlohmann::json manifest_to_json(const InputStoreManifest& manifest) {
  auto normalized = manifest;
  normalized.has_molecule_rows = normalized.has_molecule_rows || normalized.has_molecules;
  normalized.has_molecules = normalized.has_molecules || normalized.has_molecule_rows;
  if (normalized.store_mode.empty()) {
    normalized.store_mode = normalized.has_molecule_rows ? "full" : "counts";
  }
  return nlohmann::json{
      {"format_version", normalized.format_version},
      {"source_type", normalized.source_type},
      {"store_mode", normalized.store_mode},
      {"source_path", normalized.source_path},
      {"source_fingerprint", normalized.source_fingerprint},
      {"filter_signature", normalized.filter_signature},
      {"has_molecules", normalized.has_molecules},
      {"has_molecule_rows", normalized.has_molecule_rows},
      {"has_cell_gene_counts", normalized.has_cell_gene_counts},
      {"has_cell_offsets", normalized.has_cell_offsets},
      {"has_spatial_tiles", normalized.has_spatial_tiles},
      {"has_labels", normalized.has_labels},
      {"has_factor_scores", normalized.has_factor_scores},
      {"n_molecules", normalized.n_molecules},
      {"n_cells", normalized.n_cells},
      {"n_genes", normalized.n_genes},
      {"has_z", normalized.has_z},
      {"has_qv", normalized.has_qv},
      {"has_transcript_id", normalized.has_transcript_id},
      {"has_cell_type", normalized.has_cell_type},
      {"has_sample_id", normalized.has_sample_id},
      {"has_fov_id", normalized.has_fov_id},
      {"has_overlaps_nucleus", normalized.has_overlaps_nucleus},
      {"has_nucleus_distance", normalized.has_nucleus_distance}};
}

InputStoreManifest manifest_from_json(const nlohmann::json& json, const std::filesystem::path& root) {
  InputStoreManifest manifest;
  manifest.root_dir = root.string();
  manifest.format_version = json.value("format_version", "0.1");
  manifest.source_type = json.value("source_type", "");
  manifest.store_mode = json.value("store_mode", "");
  manifest.source_path = json.value("source_path", "");
  manifest.source_fingerprint = json.value("source_fingerprint", "");
  manifest.filter_signature = json.value("filter_signature", "");
  manifest.has_molecules = json.value("has_molecules", false);
  manifest.has_molecule_rows = json.value("has_molecule_rows", manifest.has_molecules);
  manifest.has_cell_gene_counts = json.value("has_cell_gene_counts", false);
  manifest.has_cell_offsets = json.value("has_cell_offsets", false);
  manifest.has_spatial_tiles = json.value("has_spatial_tiles", false);
  manifest.has_labels = json.value("has_labels", false);
  manifest.has_factor_scores = json.value("has_factor_scores", false);
  manifest.n_molecules = json.value("n_molecules", static_cast<std::size_t>(0));
  manifest.n_cells = json.value("n_cells", static_cast<std::size_t>(0));
  manifest.n_genes = json.value("n_genes", static_cast<std::size_t>(0));
  manifest.has_z = json.value("has_z", true);
  manifest.has_qv = json.value("has_qv", false);
  manifest.has_transcript_id = json.value("has_transcript_id", false);
  manifest.has_cell_type = json.value("has_cell_type", false);
  manifest.has_sample_id = json.value("has_sample_id", false);
  manifest.has_fov_id = json.value("has_fov_id", false);
  manifest.has_overlaps_nucleus = json.value("has_overlaps_nucleus", false);
  manifest.has_nucleus_distance = json.value("has_nucleus_distance", false);
  manifest.has_molecule_rows = manifest.has_molecule_rows || manifest.has_molecules;
  manifest.has_molecules = manifest.has_molecules || manifest.has_molecule_rows;
  if (manifest.store_mode.empty()) {
    manifest.store_mode = manifest.has_molecule_rows ? "full" : "counts";
  }
  return manifest;
}

void write_manifest(const InputStoreManifest& manifest) {
  const auto paths = make_input_store_paths(manifest.root_dir);
  std::ofstream out(paths.manifest_json);
  if (!out) {
    throw std::runtime_error("Could not open input store manifest for writing: " + paths.manifest_json);
  }
  out << manifest_to_json(manifest).dump(2) << "\n";
}

InputStoreManifest write_store_from_table(
    const std::string& source_type,
    const std::string& source_path,
    const std::string& source_fingerprint,
    const std::string& filter_signature,
    const TranscriptTable& table,
    const CellTable& cells,
    const std::vector<std::string>* transcript_crop_ids,
    const InputStoreBuildOptions& options) {
  const auto paths = make_input_store_paths(options.store_dir);
  ensure_directory(paths.root_dir);
  const auto counts = build_counts_from_table(table, cells, transcript_crop_ids);
  write_genes(paths.genes_parquet, counts.genes, options.parquet_row_group_size);
  write_cells(paths.cells_parquet, counts, options.parquet_row_group_size);
  write_counts(paths.counts_parquet, counts, options.parquet_row_group_size);
  if (options.materialize_molecules) {
    write_molecules(
        paths.molecules_parquet,
        paths.cell_offsets_parquet,
        table,
        transcript_crop_ids,
        options.parquet_row_group_size);
  }

  InputStoreManifest manifest;
  manifest.root_dir = paths.root_dir;
  manifest.source_type = source_type;
  manifest.store_mode = options.materialize_molecules ? "full" : "counts";
  manifest.source_path = source_path;
  manifest.source_fingerprint = source_fingerprint;
  manifest.filter_signature = filter_signature;
  manifest.has_molecules = options.materialize_molecules;
  manifest.has_molecule_rows = options.materialize_molecules;
  manifest.has_cell_gene_counts = true;
  manifest.has_cell_offsets = options.materialize_molecules;
  manifest.n_molecules = table.size();
  manifest.n_cells = table.num_cells();
  manifest.n_genes = table.num_genes();
  manifest.has_z = true;
  manifest.has_qv = !table.qv.empty();
  manifest.has_transcript_id = !table.transcript_ids.empty();
  manifest.has_cell_type = !counts.cells.cell_types.empty();
  manifest.has_sample_id = !counts.cells.sample_ids.empty();
  manifest.has_fov_id = !counts.cells.fov_ids.empty();
  manifest.has_overlaps_nucleus = options.materialize_molecules && !table.overlaps_nucleus.empty();
  manifest.has_nucleus_distance = options.materialize_molecules && !table.nucleus_distance.empty();
  write_manifest(manifest);
  return manifest;
}

InputStoreManifest write_counts_only_store(
    const std::string& source_type,
    const std::string& source_path,
    const std::string& source_fingerprint,
    const std::string& filter_signature,
    const CellCountMatrix& counts,
    std::size_t n_molecules,
    bool has_z,
    bool has_qv,
    bool has_transcript_id,
    const InputStoreBuildOptions& options) {
  const auto paths = make_input_store_paths(options.store_dir);
  ensure_directory(paths.root_dir);
  write_genes(paths.genes_parquet, counts.genes, options.parquet_row_group_size);
  write_cells(paths.cells_parquet, counts, options.parquet_row_group_size);
  write_counts(paths.counts_parquet, counts, options.parquet_row_group_size);

  InputStoreManifest manifest;
  manifest.root_dir = paths.root_dir;
  manifest.source_type = source_type;
  manifest.store_mode = "counts";
  manifest.source_path = source_path;
  manifest.source_fingerprint = source_fingerprint;
  manifest.filter_signature = filter_signature;
  manifest.has_molecules = false;
  manifest.has_molecule_rows = false;
  manifest.has_cell_gene_counts = true;
  manifest.has_cell_offsets = false;
  manifest.n_molecules = n_molecules;
  manifest.n_cells = counts.cells.cell_ids.size();
  manifest.n_genes = counts.genes.size();
  manifest.has_z = has_z;
  manifest.has_qv = has_qv;
  manifest.has_transcript_id = has_transcript_id;
  manifest.has_cell_type = !counts.cells.cell_types.empty();
  manifest.has_sample_id = !counts.cells.sample_ids.empty();
  manifest.has_fov_id = !counts.cells.fov_ids.empty();
  manifest.has_overlaps_nucleus = false;
  manifest.has_nucleus_distance = false;
  write_manifest(manifest);
  return manifest;
}

InputStoreManifest write_encoded_full_store(
    const std::string& source_type,
    const std::string& source_path,
    const std::string& source_fingerprint,
    const std::string& filter_signature,
    const EncodedMoleculeStore& encoded,
    const InputStoreBuildOptions& options) {
  const auto paths = make_input_store_paths(options.store_dir);
  ensure_directory(paths.root_dir);
  write_genes(paths.genes_parquet, encoded.counts.genes, options.parquet_row_group_size);
  write_cells(paths.cells_parquet, encoded.counts, options.parquet_row_group_size);
  write_counts(paths.counts_parquet, encoded.counts, options.parquet_row_group_size);
  write_encoded_molecules(
      paths.molecules_parquet,
      paths.cell_offsets_parquet,
      encoded,
      options.parquet_row_group_size);

  InputStoreManifest manifest;
  manifest.root_dir = paths.root_dir;
  manifest.source_type = source_type;
  manifest.store_mode = "full";
  manifest.source_path = source_path;
  manifest.source_fingerprint = source_fingerprint;
  manifest.filter_signature = filter_signature;
  manifest.has_molecules = true;
  manifest.has_molecule_rows = true;
  manifest.has_cell_gene_counts = true;
  manifest.has_cell_offsets = true;
  manifest.n_molecules = encoded.obs_id.size();
  manifest.n_cells = encoded.counts.cells.cell_ids.size();
  manifest.n_genes = encoded.counts.genes.size();
  manifest.has_z = encoded.has_z;
  manifest.has_qv = encoded.has_qv;
  manifest.has_transcript_id = encoded.has_transcript_id;
  manifest.has_cell_type = !encoded.counts.cells.cell_types.empty();
  manifest.has_sample_id = !encoded.counts.cells.sample_ids.empty();
  manifest.has_fov_id = !encoded.counts.cells.fov_ids.empty();
  manifest.has_overlaps_nucleus = encoded.has_overlaps_nucleus;
  manifest.has_nucleus_distance = encoded.has_nucleus_distance;
  write_manifest(manifest);
  return manifest;
}

}  // namespace

InputStorePaths make_input_store_paths(const std::string& root_dir) {
  const auto root = std::filesystem::path(root_dir);
  return {
      root.string(),
      (root / "store.json").string(),
      (root / "genes.parquet").string(),
      (root / "cells.parquet").string(),
      (root / "cell_gene_counts.parquet").string(),
      (root / "molecules.parquet").string(),
      (root / "cell_offsets.parquet").string()};
}

InputStoreManifest read_input_store_manifest(const std::string& store_dir) {
  const auto paths = make_input_store_paths(store_dir);
  std::ifstream in(paths.manifest_json);
  if (!in) {
    throw std::runtime_error("Could not open input store manifest: " + paths.manifest_json);
  }
  nlohmann::json json;
  in >> json;
  return manifest_from_json(json, std::filesystem::path(store_dir));
}

InputStoreManifest build_xenium_input_store(
    const std::string& manifest_path_or_dir,
    const XeniumLoadOptions& load_options,
    const InputStoreBuildOptions& store_options) {
  const std::string manifest_path = normalize_path_for_store(resolve_xenium_manifest_path(manifest_path_or_dir));
  const auto manifest = read_xenium_manifest(manifest_path);
  const std::string dataset_dir = dirname_of(manifest_path);
  const std::string transcripts_parquet = join_path(dataset_dir, manifest.transcripts_parquet_path);
  const std::string transcripts_csv = join_path(dataset_dir, manifest.transcripts_csv_path);
  const bool use_parquet = load_options.prefer_parquet && file_exists(transcripts_parquet);
  std::vector<std::string> source_files = {manifest_path, use_parquet ? transcripts_parquet : transcripts_csv};
  const std::string cells_parquet = (std::filesystem::path(dataset_dir) / "cells.parquet").string();
  const std::string cells_csv = (std::filesystem::path(dataset_dir) / "cells.csv.gz").string();
  if (load_options.read_cells) {
    if (file_exists(cells_parquet)) {
      source_files.push_back(cells_parquet);
    } else if (file_exists(cells_csv)) {
      source_files.push_back(cells_csv);
    }
  }
  const std::string source_fingerprint = combined_fingerprint(source_files);
  const std::string filter_signature = xenium_filter_signature(load_options);

  InputStoreManifest existing;
  if (!store_options.force &&
      compatible_existing_store(
          store_options.store_dir,
          "xenium",
          manifest_path,
          source_fingerprint,
          filter_signature,
          store_options.materialize_molecules,
          &existing)) {
    if (load_options.progress) {
      load_options.progress("Reused compatible Xenium input store", 0.0);
    }
    return existing;
  }

  std::optional<CellCountMatrix> upgrade_counts_reference;
  if (store_options.materialize_molecules && !store_options.force &&
      compatible_existing_store(
          store_options.store_dir,
          "xenium",
          manifest_path,
          source_fingerprint,
          filter_signature,
          false,
          &existing) &&
      !existing.has_molecule_rows) {
    upgrade_counts_reference = load_input_store_counts(store_options.store_dir);
  }

  if (!store_options.materialize_molecules && use_parquet) {
    auto counts = build_xenium_counts_from_parquet(transcripts_parquet, load_options);
    if (load_options.read_cells) {
      enrich_counts_from_cells_parquet(counts.counts, cells_parquet);
    }
    return write_counts_only_store(
        "xenium",
        manifest_path,
        source_fingerprint,
        filter_signature,
        counts.counts,
        static_cast<std::size_t>(std::max<std::int64_t>(counts.kept_rows, 0)),
        counts.has_z,
        counts.has_qv,
        counts.has_transcript_id,
        store_options);
  }

  if (store_options.materialize_molecules && use_parquet) {
    auto encoded = build_xenium_encoded_store_from_parquet(transcripts_parquet, load_options);
    if (load_options.read_cells) {
      enrich_counts_from_cells_parquet(encoded.counts, cells_parquet);
    }
    if (upgrade_counts_reference.has_value()) {
      validate_counts_upgrade_compatibility(
          *upgrade_counts_reference,
          encoded.counts,
          "Xenium counts-to-full store upgrade");
    }
    return write_encoded_full_store(
        "xenium",
        manifest_path,
        source_fingerprint,
        filter_signature,
        encoded,
        store_options);
  }

  auto bundle = load_xenium_bundle(manifest_path, load_options);
  if (!bundle.cells.cell_ids.empty()) {
    annotate_transcripts_from_cell_table(bundle.transcripts, bundle.cells);
  }
  const CellTable cells = bundle.cells.cell_ids.empty() ? infer_cell_table(bundle.transcripts) : bundle.cells;
  if (upgrade_counts_reference.has_value()) {
    const auto candidate_counts = build_counts_from_table(
        bundle.transcripts,
        cells,
        bundle.transcript_crop_ids.empty() ? nullptr : &bundle.transcript_crop_ids);
    validate_counts_upgrade_compatibility(
        *upgrade_counts_reference,
        candidate_counts,
        "Xenium counts-to-full store upgrade");
  }
  return write_store_from_table(
      "xenium",
      manifest_path,
      source_fingerprint,
      filter_signature,
      bundle.transcripts,
      cells,
      bundle.transcript_crop_ids.empty() ? nullptr : &bundle.transcript_crop_ids,
      store_options);
}

InputStoreManifest build_tabular_input_store(
    const TabularSourceSpec& source,
    const TabularLoadOptions& load_options,
    const InputStoreBuildOptions& store_options) {
  const std::string source_path = normalize_path_for_store(source.molecules_path);
  const std::string source_fingerprint = combined_fingerprint({
      source.molecules_path,
      source.segmentation_mask_path,
      source.cell_metadata_path});
  const std::string filter_signature = tabular_filter_signature(source, load_options);
  InputStoreManifest existing;
  if (!store_options.force &&
      compatible_existing_store(
          store_options.store_dir,
          "tabular",
          source_path,
          source_fingerprint,
          filter_signature,
          store_options.materialize_molecules,
          &existing)) {
    return existing;
  }

  std::optional<CellCountMatrix> upgrade_counts_reference;
  if (store_options.materialize_molecules && !store_options.force &&
      compatible_existing_store(
          store_options.store_dir,
          "tabular",
          source_path,
          source_fingerprint,
          filter_signature,
          false,
          &existing) &&
      !existing.has_molecule_rows) {
    upgrade_counts_reference = load_input_store_counts(store_options.store_dir);
  }

  if (can_stream_tabular_input(source)) {
    if (!store_options.materialize_molecules) {
      auto counts = build_tabular_counts_from_stream(source, load_options);
      return write_counts_only_store(
          "tabular",
          source_path,
          source_fingerprint,
          filter_signature,
          counts.counts,
          static_cast<std::size_t>(std::max<std::int64_t>(counts.kept_rows, 0)),
          counts.has_z,
          counts.has_qv,
          counts.has_transcript_id,
          store_options);
    }

    auto encoded = build_tabular_encoded_store_from_stream(source, load_options);
    if (upgrade_counts_reference.has_value()) {
      validate_counts_upgrade_compatibility(
          *upgrade_counts_reference,
          encoded.counts,
          "Tabular counts-to-full store upgrade");
    }
    return write_encoded_full_store(
        "tabular",
        source_path,
        source_fingerprint,
        filter_signature,
        encoded,
        store_options);
  }

  auto bundle = load_tabular_bundle(source, load_options);
  const CellTable cells = infer_cell_table(bundle.transcripts);
  if (upgrade_counts_reference.has_value()) {
    const auto candidate_counts = build_counts_from_table(
        bundle.transcripts,
        cells,
        bundle.transcript_crop_ids.empty() ? nullptr : &bundle.transcript_crop_ids);
    validate_counts_upgrade_compatibility(
        *upgrade_counts_reference,
        candidate_counts,
        "Tabular counts-to-full store upgrade");
  }
  return write_store_from_table(
      "tabular",
      source_path,
      source_fingerprint,
      filter_signature,
      bundle.transcripts,
      cells,
      bundle.transcript_crop_ids.empty() ? nullptr : &bundle.transcript_crop_ids,
      store_options);
}

std::vector<std::string> load_input_store_gene_names(const InputStorePaths& paths) {
  const auto genes = read_parquet_table(paths.genes_parquet);
  return read_string_column(genes, "gene");
}

CellTable load_input_store_cell_table(
    const InputStorePaths& paths,
    std::vector<int>* transcript_counts,
    std::vector<std::string>* crop_ids) {
  const auto cells = read_parquet_table(paths.cells_parquet);
  CellTable out;
  out.cell_ids = read_string_column(cells, "cell_id");
  if (transcript_counts != nullptr) {
    *transcript_counts = read_int32_column(cells, "transcript_count");
  }
  out.centroid_x = read_double_column(cells, "centroid_x");
  out.centroid_y = read_double_column(cells, "centroid_y");
  out.centroid_z = read_double_column(cells, "centroid_z");
  if (crop_ids != nullptr) {
    *crop_ids = read_string_column(cells, "crop_id");
  }
  if (has_column(cells, "cell_type")) out.cell_types = read_string_column(cells, "cell_type");
  if (has_column(cells, "sample_id")) out.sample_ids = read_string_column(cells, "sample_id");
  if (has_column(cells, "fov_id")) out.fov_ids = read_string_column(cells, "fov_id");
  return out;
}

CellTable load_input_store_cells(const std::string& store_dir) {
  return load_input_store_cell_table(make_input_store_paths(store_dir), nullptr, nullptr);
}

std::vector<InputStoreCellOffset> load_input_store_cell_offsets(const std::string& store_dir) {
  const auto paths = make_input_store_paths(store_dir);
  const auto offsets = read_parquet_table(paths.cell_offsets_parquet);
  const auto cell_idx = read_int32_column(offsets, "cell_idx");
  const auto start = read_int64_column(offsets, "start");
  const auto end = read_int64_column(offsets, "end");
  if (cell_idx.size() != start.size() || cell_idx.size() != end.size()) {
    throw std::runtime_error("cell_offsets.parquet columns must have equal length");
  }
  std::vector<InputStoreCellOffset> out(cell_idx.size());
  for (std::size_t i = 0; i < cell_idx.size(); ++i) {
    out[i] = InputStoreCellOffset{cell_idx[i], start[i], end[i]};
  }
  return out;
}

std::vector<int> row_groups_overlapping_row_range(
    const std::string& parquet_path,
    std::int64_t start,
    std::int64_t end,
    std::vector<std::int64_t>* row_group_starts = nullptr) {
  auto parquet_reader = parquet::ParquetFileReader::OpenFile(parquet_path, false);
  const auto metadata = parquet_reader->metadata();
  std::vector<int> out;
  std::int64_t cursor = 0;
  if (row_group_starts != nullptr) {
    row_group_starts->assign(static_cast<std::size_t>(metadata->num_row_groups()), 0);
  }
  for (int rg = 0; rg < metadata->num_row_groups(); ++rg) {
    if (row_group_starts != nullptr) {
      (*row_group_starts)[static_cast<std::size_t>(rg)] = cursor;
    }
    const auto rows = metadata->RowGroup(rg)->num_rows();
    const std::int64_t rg_start = cursor;
    const std::int64_t rg_end = cursor + rows;
    if (rg_end > start && rg_start < end) {
      out.push_back(rg);
    }
    cursor = rg_end;
  }
  return out;
}

InputStoreMoleculeBlock load_input_store_molecule_range(
    const std::string& store_dir,
    std::int64_t start,
    std::int64_t end) {
  if (start < 0 || end < start) {
    throw std::runtime_error("invalid input-store molecule row range");
  }
  InputStoreMoleculeBlock out;
  if (end == start) {
    return out;
  }

  const auto paths = make_input_store_paths(store_dir);
  auto input = arrow_unwrap(
      arrow::io::ReadableFile::Open(paths.molecules_parquet),
      "Open input-store molecules parquet");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open input-store molecules parquet reader");
  auto reader = arrow_unwrap(builder.Build(), "Build input-store molecules parquet reader");
  reader->set_use_threads(true);
  reader->set_batch_size(65536);

  std::shared_ptr<arrow::Schema> schema;
  arrow_check(reader->GetSchema(&schema), "Read input-store molecules schema");
  const int obs_col = schema->GetFieldIndex("obs_id");
  const int cell_col = schema->GetFieldIndex("cell_idx");
  const int gene_col = schema->GetFieldIndex("gene_idx");
  const int x_col = schema->GetFieldIndex("x");
  const int y_col = schema->GetFieldIndex("y");
  const int z_col = schema->GetFieldIndex("z");
  if (obs_col < 0 || cell_col < 0 || gene_col < 0 || x_col < 0 || y_col < 0 || z_col < 0) {
    throw std::runtime_error("input-store molecules parquet is missing required columns");
  }
  const int qv_col = schema->GetFieldIndex("qv");
  const int overlaps_col = schema->GetFieldIndex("overlaps_nucleus");
  const int nucleus_distance_col = schema->GetFieldIndex("nucleus_distance");
  const int crop_col = schema->GetFieldIndex("crop_id");
  std::vector<int> projected = {obs_col, cell_col, gene_col, x_col, y_col, z_col};
  if (qv_col >= 0) projected.push_back(qv_col);
  if (overlaps_col >= 0) projected.push_back(overlaps_col);
  if (nucleus_distance_col >= 0) projected.push_back(nucleus_distance_col);
  if (crop_col >= 0) projected.push_back(crop_col);

  std::vector<std::int64_t> row_group_starts;
  const auto row_groups =
      row_groups_overlapping_row_range(paths.molecules_parquet, start, end, &row_group_starts);
  out.row_index.reserve(static_cast<std::size_t>(end - start));
  out.obs_id.reserve(static_cast<std::size_t>(end - start));
  out.cell_idx.reserve(static_cast<std::size_t>(end - start));
  out.gene_idx.reserve(static_cast<std::size_t>(end - start));
  out.x.reserve(static_cast<std::size_t>(end - start));
  out.y.reserve(static_cast<std::size_t>(end - start));
  out.z.reserve(static_cast<std::size_t>(end - start));
  if (qv_col >= 0) out.qv.reserve(static_cast<std::size_t>(end - start));
  if (overlaps_col >= 0) out.overlaps_nucleus.reserve(static_cast<std::size_t>(end - start));
  if (nucleus_distance_col >= 0) out.nucleus_distance.reserve(static_cast<std::size_t>(end - start));
  if (crop_col >= 0) out.crop_id.reserve(static_cast<std::size_t>(end - start));

  for (const int rg : row_groups) {
    auto batch_reader = arrow_unwrap(
        reader->GetRecordBatchReader(std::vector<int>{rg}, projected),
        "Create input-store molecule batch reader");
    std::shared_ptr<arrow::RecordBatch> batch;
    std::int64_t batch_start = row_group_starts[static_cast<std::size_t>(rg)];
    while (true) {
      arrow_check(batch_reader->ReadNext(&batch), "Read input-store molecule batch");
      if (!batch) break;

      int col = 0;
      NumericArrayView obs_view(batch->column(col++));
      NumericArrayView cell_view(batch->column(col++));
      NumericArrayView gene_view(batch->column(col++));
      NumericArrayView x_view(batch->column(col++));
      NumericArrayView y_view(batch->column(col++));
      NumericArrayView z_view(batch->column(col++));
      std::optional<NumericArrayView> qv_view;
      if (qv_col >= 0) qv_view.emplace(batch->column(col++));
      std::optional<NumericArrayView> overlaps_view;
      if (overlaps_col >= 0) overlaps_view.emplace(batch->column(col++));
      std::optional<NumericArrayView> nucleus_distance_view;
      if (nucleus_distance_col >= 0) nucleus_distance_view.emplace(batch->column(col++));
      std::optional<StringArrayView> crop_view;
      if (crop_col >= 0) crop_view.emplace(batch->column(col++));

      for (std::int64_t row = 0; row < batch->num_rows(); ++row) {
        const std::int64_t global = batch_start + row;
        if (global < start || global >= end) {
          continue;
        }
        out.row_index.push_back(global);
        out.obs_id.push_back(static_cast<std::int64_t>(obs_view.value(row)));
        out.cell_idx.push_back(static_cast<int>(cell_view.value(row)));
        out.gene_idx.push_back(static_cast<int>(gene_view.value(row)));
        out.x.push_back(x_view.value(row));
        out.y.push_back(y_view.value(row));
        out.z.push_back(z_view.value(row));
        if (qv_view.has_value()) out.qv.push_back(qv_view->value(row));
        if (overlaps_view.has_value()) out.overlaps_nucleus.push_back(static_cast<int>(overlaps_view->value(row)));
        if (nucleus_distance_view.has_value()) out.nucleus_distance.push_back(nucleus_distance_view->value(row));
        if (crop_view.has_value()) out.crop_id.push_back(crop_view->value(row));
      }
      batch_start += batch->num_rows();
    }
  }
  return out;
}

CellCountMatrix load_input_store_counts(const std::string& store_dir) {
  const auto paths = make_input_store_paths(store_dir);
  CellCountMatrix out;
  out.genes = load_input_store_gene_names(paths);
  out.cells = load_input_store_cell_table(paths, &out.transcript_counts, &out.crop_ids);

  const auto counts = read_parquet_table(paths.counts_parquet);
  const auto count_cell = read_int32_column(counts, "cell_idx");
  const auto count_gene = read_int32_column(counts, "gene_idx");
  const auto count_value = read_int32_column(counts, "count");
  out.indptr.assign(out.cells.cell_ids.size() + 1, 0);
  for (const int cell : count_cell) {
    out.indptr[static_cast<std::size_t>(cell + 1)] += 1;
  }
  for (std::size_t i = 1; i < out.indptr.size(); ++i) {
    out.indptr[i] += out.indptr[i - 1];
  }
  out.indices.assign(count_cell.size(), 0);
  out.values.assign(count_cell.size(), 0.0);
  std::vector<int> write_pos = out.indptr;
  for (std::size_t i = 0; i < count_cell.size(); ++i) {
    const int pos = write_pos[static_cast<std::size_t>(count_cell[i])]++;
    out.indices[static_cast<std::size_t>(pos)] = count_gene[i];
    out.values[static_cast<std::size_t>(pos)] = static_cast<double>(count_value[i]);
  }
  return out;
}

InputStoreData load_input_store_numeric_data(const std::string& store_dir) {
  InputStoreData out;
  out.manifest = read_input_store_manifest(store_dir);
  const auto paths = make_input_store_paths(store_dir);
  if (!out.manifest.has_molecules) {
    throw std::runtime_error("input store does not contain materialized molecules");
  }

  const auto gene_names = load_input_store_gene_names(paths);
  out.cells = load_input_store_cell_table(paths, nullptr, nullptr);
  const auto molecules = read_parquet_table(paths.molecules_parquet);
  const auto cell_idx = read_int32_column(molecules, "cell_idx");
  const auto gene_idx = read_int32_column(molecules, "gene_idx");

  out.transcripts.x = read_double_column(molecules, "x");
  out.transcripts.y = read_double_column(molecules, "y");
  out.transcripts.z = read_double_column(molecules, "z");
  if (has_column(molecules, "qv")) out.transcripts.qv = read_double_column(molecules, "qv");
  if (has_column(molecules, "overlaps_nucleus")) {
    out.transcripts.overlaps_nucleus = read_int32_column(molecules, "overlaps_nucleus");
  }
  if (has_column(molecules, "nucleus_distance")) {
    out.transcripts.nucleus_distance = read_double_column(molecules, "nucleus_distance");
  }
  if (has_column(molecules, "transcript_id")) {
    out.transcripts.transcript_ids = read_string_column(molecules, "transcript_id");
  }
  if (has_column(molecules, "crop_id")) out.transcript_crop_ids = read_string_column(molecules, "crop_id");

  out.transcripts.genes = gene_names;
  out.transcripts.cells = out.cells.cell_ids;
  out.transcripts.gene_index = gene_idx;
  out.transcripts.cell_index = cell_idx;
  return out;
}

InputStoreData load_input_store_data(const std::string& store_dir) {
  InputStoreData out;
  out.manifest = read_input_store_manifest(store_dir);
  const auto paths = make_input_store_paths(store_dir);
  if (!out.manifest.has_molecules) {
    throw std::runtime_error("input store does not contain materialized molecules");
  }
  const auto gene_names = load_input_store_gene_names(paths);
  out.cells = load_input_store_cell_table(paths, nullptr, nullptr);
  const auto molecules = read_parquet_table(paths.molecules_parquet);
  const auto obs_id = read_int64_column(molecules, "obs_id");
  const auto cell_idx = read_int32_column(molecules, "cell_idx");
  const auto gene_idx = read_int32_column(molecules, "gene_idx");
  out.transcripts.x = read_double_column(molecules, "x");
  out.transcripts.y = read_double_column(molecules, "y");
  out.transcripts.z = read_double_column(molecules, "z");
  if (has_column(molecules, "qv")) out.transcripts.qv = read_double_column(molecules, "qv");
  if (has_column(molecules, "overlaps_nucleus")) out.transcripts.overlaps_nucleus = read_int32_column(molecules, "overlaps_nucleus");
  if (has_column(molecules, "nucleus_distance")) out.transcripts.nucleus_distance = read_double_column(molecules, "nucleus_distance");
  if (has_column(molecules, "transcript_id")) out.transcripts.transcript_ids = read_string_column(molecules, "transcript_id");
  if (has_column(molecules, "crop_id")) out.transcript_crop_ids = read_string_column(molecules, "crop_id");

  out.transcripts.genes = gene_names;
  out.transcripts.cells = out.cells.cell_ids;
  out.transcripts.gene_index = gene_idx;
  out.transcripts.cell_index = cell_idx;
  out.transcripts.gene_key.reserve(gene_idx.size());
  out.transcripts.orig_cell_id.reserve(cell_idx.size());
  for (std::size_t i = 0; i < gene_idx.size(); ++i) {
    out.transcripts.gene_key.push_back(gene_names[static_cast<std::size_t>(gene_idx[i])]);
    out.transcripts.orig_cell_id.push_back(out.cells.cell_ids[static_cast<std::size_t>(cell_idx[i])]);
  }
  if (!out.cells.cell_types.empty()) {
    out.transcripts.cell_types.reserve(cell_idx.size());
    for (const int cell : cell_idx) out.transcripts.cell_types.push_back(out.cells.cell_types[static_cast<std::size_t>(cell)]);
  }
  if (!out.cells.sample_ids.empty()) {
    out.transcripts.sample_ids.reserve(cell_idx.size());
    for (const int cell : cell_idx) out.transcripts.sample_ids.push_back(out.cells.sample_ids[static_cast<std::size_t>(cell)]);
  }
  if (!out.cells.fov_ids.empty()) {
    out.transcripts.fov_ids.reserve(cell_idx.size());
    for (const int cell : cell_idx) out.transcripts.fov_ids.push_back(out.cells.fov_ids[static_cast<std::size_t>(cell)]);
  }
  (void) obs_id;
  return out;
}

}  // namespace celladmix
