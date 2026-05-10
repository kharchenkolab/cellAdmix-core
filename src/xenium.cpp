#include "celladmix/xenium.hpp"

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <arrow/result.h>
#include <arrow/util/bit_util.h>
#include <parquet/arrow/reader.h>
#include <parquet/file_reader.h>
#include <parquet/statistics.h>
#include <zlib.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace celladmix {

// Xenium manifest parsing and transcript/cell loading from CSV or parquet bundles.

namespace {

#define CELLADMIX_ARROW_CHECK_OK(expr)                                                \
  do {                                                                                \
    auto _status = (expr);                                                            \
    if (!_status.ok()) {                                                              \
      throw std::runtime_error(std::string("Arrow error: ") + _status.ToString());    \
    }                                                                                 \
  } while (0)

// Unwrap an Arrow result or raise a generic Arrow error.
template <typename T>
T arrow_unwrap(arrow::Result<T>&& result) {
  if (!result.ok()) {
    throw std::runtime_error(std::string("Arrow error: ") + result.status().ToString());
  }
  return std::move(*result);
}

// Emit optional loader progress with elapsed seconds since the previous update.
void emit_load_progress(
    const XeniumLoadOptions& options,
    const std::string& message,
    std::chrono::steady_clock::time_point& stage_start) {
  if (!options.progress) {
    return;
  }
  const auto now = std::chrono::steady_clock::now();
  options.progress(message, std::chrono::duration<double>(now - stage_start).count());
  stage_start = now;
}

// Return the parent directory of a manifest or data file path.
std::string dirname_of(const std::string& path) {
  return std::filesystem::path(path).parent_path().string();
}

// Join a dataset directory with a relative file path from the manifest.
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

// Check whether one expected Xenium file is present.
bool file_exists(const std::string& path) {
  return std::filesystem::exists(std::filesystem::path(path));
}

// Read one small manifest-like text file into memory.
std::string slurp(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("Failed to open Xenium manifest: " + path);
  }
  std::ostringstream buffer;
  buffer << input.rdbuf();
  return buffer.str();
}

// Extract a string field from the manifest JSON with a regex fallback parser.
std::string extract_string(const std::string& text, const std::string& key) {
  const std::regex re("\"" + key + "\"\\s*:\\s*\"([^\"]*)\"");
  std::smatch match;
  if (std::regex_search(text, match, re)) {
    return match[1];
  }
  return {};
}

// Extract a numeric field from the manifest JSON with a regex fallback parser.
double extract_double(const std::string& text, const std::string& key, double default_value) {
  const std::regex re("\"" + key + "\"\\s*:\\s*([0-9]+(?:\\.[0-9]+)?)");
  std::smatch match;
  if (std::regex_search(text, match, re)) {
    return std::stod(match[1]);
  }
  return default_value;
}

// Resolve either a Xenium directory or explicit manifest path to experiment.xenium.
std::string resolve_manifest_path(const std::string& path_or_dir) {
  const auto path = std::filesystem::path(path_or_dir);
  if (std::filesystem::is_directory(path)) {
    return (path / "experiment.xenium").string();
  }
  return path.string();
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

// Compute one union bounding box covering all requested crops.
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

// Assign one transcript or cell centroid to the first matching crop.
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

// Quickly reject coordinates that are outside the union of all requested crops.
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

// Split one CSV line while respecting quoted commas and escaped quotes.
std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> fields;
  std::string current;
  bool in_quotes = false;
  for (std::size_t i = 0; i < line.size(); ++i) {
    const char ch = line[i];
    if (ch == '"') {
      if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
        current.push_back('"');
        i += 1;
      } else {
        in_quotes = !in_quotes;
      }
    } else if (ch == ',' && !in_quotes) {
      fields.push_back(current);
      current.clear();
    } else {
      current.push_back(ch);
    }
  }
  fields.push_back(current);
  return fields;
}

class LineReader {
 public:
  explicit LineReader(const std::string& path) {
    gz_ = path.size() >= 3 && path.substr(path.size() - 3) == ".gz";
    if (gz_) {
      gz_file_ = gzopen(path.c_str(), "rb");
      if (gz_file_ == nullptr) {
        throw std::runtime_error("Failed to open gzip file: " + path);
      }
    } else {
      input_.open(path);
      if (!input_) {
        throw std::runtime_error("Failed to open file: " + path);
      }
    }
  }

  ~LineReader() {
    if (gz_file_ != nullptr) {
      gzclose(gz_file_);
    }
  }

  bool next_line(std::string& line) {
    if (gz_) {
      return next_gz_line(line);
    }
    return static_cast<bool>(std::getline(input_, line));
  }

 private:
  bool next_gz_line(std::string& line) {
    line.clear();
    if (gz_file_ == nullptr) {
      return false;
    }
    char buffer[1 << 16];
    while (true) {
      char* result = gzgets(gz_file_, buffer, static_cast<int>(sizeof(buffer)));
      if (result == nullptr) {
        return !line.empty();
      }
      line.append(buffer);
      if (!line.empty() && line.back() == '\n') {
        line.pop_back();
        if (!line.empty() && line.back() == '\r') {
          line.pop_back();
        }
        return true;
      }
      if (gzeof(gz_file_) != 0) {
        return !line.empty();
      }
    }
  }

  bool gz_ = false;
  gzFile gz_file_ = nullptr;
  std::ifstream input_;
};

// Look up a required CSV/parquet header using a list of acceptable aliases.
std::size_t find_required_header(
    const std::unordered_map<std::string, std::size_t>& header_index,
    const std::vector<std::string>& candidates,
    const std::string& label) {
  for (const auto& candidate : candidates) {
    const auto it = header_index.find(candidate);
    if (it != header_index.end()) {
      return it->second;
    }
  }
  throw std::runtime_error("Missing required Xenium column: " + label);
}

// Look up an optional CSV/parquet header using a list of acceptable aliases.
std::optional<std::size_t> find_optional_header(
    const std::unordered_map<std::string, std::size_t>& header_index,
    const std::vector<std::string>& candidates) {
  for (const auto& candidate : candidates) {
    const auto it = header_index.find(candidate);
    if (it != header_index.end()) {
      return it->second;
    }
  }
  return std::nullopt;
}

// Parse one numeric CSV field with a fallback default.
double parse_double_or_default(const std::string& value, double default_value) {
  if (value.empty()) {
    return default_value;
  }
  return std::stod(value);
}

// Parse Xenium boolean-like fields such as overlaps_nucleus.
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

// Typed numeric view over Arrow arrays used during parquet scans.
struct NumericArrayView {
  enum class Kind {
    Invalid,
    Double,
    Float,
    Int64,
    Int32,
    Int16,
    Int8,
    UInt64,
    UInt32,
    UInt16,
    UInt8,
    Fallback
  };

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
    if (!arr) {
      return;
    }
    null_bitmap = arr->null_bitmap_data();
    offset = arr->offset();
    switch (arr->type_id()) {
      case arrow::Type::DOUBLE:
        double_arr = static_cast<const arrow::DoubleArray*>(arr.get());
        kind = Kind::Double;
        break;
      case arrow::Type::FLOAT:
        float_arr = static_cast<const arrow::FloatArray*>(arr.get());
        kind = Kind::Float;
        break;
      case arrow::Type::INT64:
        int64_arr = static_cast<const arrow::Int64Array*>(arr.get());
        kind = Kind::Int64;
        break;
      case arrow::Type::INT32:
        int32_arr = static_cast<const arrow::Int32Array*>(arr.get());
        kind = Kind::Int32;
        break;
      case arrow::Type::INT16:
        int16_arr = static_cast<const arrow::Int16Array*>(arr.get());
        kind = Kind::Int16;
        break;
      case arrow::Type::INT8:
        int8_arr = static_cast<const arrow::Int8Array*>(arr.get());
        kind = Kind::Int8;
        break;
      case arrow::Type::UINT64:
        uint64_arr = static_cast<const arrow::UInt64Array*>(arr.get());
        kind = Kind::UInt64;
        break;
      case arrow::Type::UINT32:
        uint32_arr = static_cast<const arrow::UInt32Array*>(arr.get());
        kind = Kind::UInt32;
        break;
      case arrow::Type::UINT16:
        uint16_arr = static_cast<const arrow::UInt16Array*>(arr.get());
        kind = Kind::UInt16;
        break;
      case arrow::Type::UINT8:
        uint8_arr = static_cast<const arrow::UInt8Array*>(arr.get());
        kind = Kind::UInt8;
        break;
      default:
        fallback = arr;
        kind = Kind::Fallback;
        break;
    }
  }

  bool is_valid(int64_t i) const {
    if (kind == Kind::Invalid) {
      return false;
    }
    return null_bitmap == nullptr || arrow::bit_util::GetBit(null_bitmap, offset + i);
  }

  double value(int64_t i) const {
    if (!is_valid(i)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
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
        auto scalar = arrow_unwrap(fallback->GetScalar(i));
        return std::stod(scalar->ToString());
      }
      default: break;
    }
    return std::numeric_limits<double>::quiet_NaN();
  }

  int64_t int64_value(int64_t i) const {
    if (!is_valid(i)) {
      return -1;
    }
    switch (kind) {
      case Kind::Int64: return int64_arr->Value(i);
      case Kind::Int32: return static_cast<int64_t>(int32_arr->Value(i));
      case Kind::Int16: return static_cast<int64_t>(int16_arr->Value(i));
      case Kind::Int8: return static_cast<int64_t>(int8_arr->Value(i));
      case Kind::UInt64: return static_cast<int64_t>(uint64_arr->Value(i));
      case Kind::UInt32: return static_cast<int64_t>(uint32_arr->Value(i));
      case Kind::UInt16: return static_cast<int64_t>(uint16_arr->Value(i));
      case Kind::UInt8: return static_cast<int64_t>(uint8_arr->Value(i));
      case Kind::Double: return static_cast<int64_t>(std::llround(double_arr->Value(i)));
      case Kind::Float: return static_cast<int64_t>(std::llround(float_arr->Value(i)));
      case Kind::Fallback: {
        auto scalar = arrow_unwrap(fallback->GetScalar(i));
        return static_cast<int64_t>(std::stoll(scalar->ToString()));
      }
      default: break;
    }
    return -1;
  }

  std::string string_value(int64_t i) const {
    if (!is_valid(i)) {
      return {};
    }
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
        out << std::setprecision(std::numeric_limits<double>::max_digits10)
            << double_arr->Value(i);
        return out.str();
      }
      case Kind::Float: {
        std::ostringstream out;
        out << std::setprecision(std::numeric_limits<float>::max_digits10)
            << float_arr->Value(i);
        return out.str();
      }
      case Kind::Fallback: {
        auto scalar = arrow_unwrap(fallback->GetScalar(i));
        return scalar->ToString();
      }
      default: break;
    }
    return {};
  }
};

// Typed string view over Arrow arrays used during parquet scans.
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
    if (!arr) {
      return;
    }
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
    if (string_arr != nullptr) {
      if (string_arr->IsNull(i)) {
        return {};
      }
      return string_arr->GetString(i);
    }
    if (large_string_arr != nullptr) {
      if (large_string_arr->IsNull(i)) {
        return {};
      }
      return large_string_arr->GetString(i);
    }
    if (dict_arr != nullptr && !dict_arr->IsNull(i)) {
      const int64_t idx = dict_indices.int64_value(i);
      if (dict_values != nullptr) return dict_values->GetString(idx);
      if (dict_large_values != nullptr) return dict_large_values->GetString(idx);
      if (dict_numeric_mode) return dict_numeric_values.string_value(idx);
      if (dict_fallback_values != nullptr) {
        auto scalar = arrow_unwrap(dict_fallback_values->GetScalar(idx));
        return scalar->ToString();
      }
      return {};
    }
    if (numeric_mode) {
      return numeric_arr.string_value(i);
    }
    if (fallback != nullptr) {
      auto scalar = arrow_unwrap(fallback->GetScalar(i));
      return scalar->ToString();
    }
    return {};
  }
};

// Find the first required field name present in a parquet schema.
std::size_t find_first_field(
    const std::shared_ptr<arrow::Schema>& schema,
    const std::vector<std::string>& candidates,
    const std::string& label) {
  for (const auto& candidate : candidates) {
    const int idx = schema->GetFieldIndex(candidate);
    if (idx >= 0) {
      return static_cast<std::size_t>(idx);
    }
  }
  throw std::runtime_error("Missing required Xenium parquet column: " + label);
}

// Find the first optional field name present in a parquet schema.
std::optional<std::size_t> find_optional_field(
    const std::shared_ptr<arrow::Schema>& schema,
    const std::vector<std::string>& candidates) {
  for (const auto& candidate : candidates) {
    const int idx = schema->GetFieldIndex(candidate);
    if (idx >= 0) {
      return static_cast<std::size_t>(idx);
    }
  }
  return std::nullopt;
}

// Detect the common unassigned-cell sentinels used in Xenium output.
bool is_unassigned_cell_id(const std::string& cell_id) {
  return cell_id.empty() || cell_id == "UNASSIGNED" || cell_id == "0" || cell_id == "cell_0";
}

std::optional<std::pair<double, double>> parquet_numeric_minmax(
    const parquet::ColumnChunkMetaData& column) {
  if (!column.is_stats_set()) {
    return std::nullopt;
  }
  const auto stats = column.statistics();
  if (!stats || !stats->HasMinMax()) {
    return std::nullopt;
  }

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
    default:
      break;
  }
  return std::nullopt;
}

bool row_group_numeric_range_overlaps(
    const parquet::RowGroupMetaData& row_group,
    int column_index,
    double min_value,
    double max_value) {
  if (column_index < 0 || column_index >= row_group.num_columns()) {
    return true;
  }
  const auto column = row_group.ColumnChunk(column_index);
  if (!column) {
    return true;
  }
  const auto minmax = parquet_numeric_minmax(*column);
  if (!minmax.has_value()) {
    return true;
  }
  if (!std::isfinite(minmax->first) || !std::isfinite(minmax->second)) {
    return true;
  }
  return !(minmax->second < min_value || minmax->first > max_value);
}

bool row_group_numeric_max_at_least(
    const parquet::RowGroupMetaData& row_group,
    int column_index,
    double min_value) {
  if (column_index < 0 || column_index >= row_group.num_columns()) {
    return true;
  }
  const auto column = row_group.ColumnChunk(column_index);
  if (!column) {
    return true;
  }
  const auto minmax = parquet_numeric_minmax(*column);
  if (!minmax.has_value()) {
    return true;
  }
  if (!std::isfinite(minmax->second)) {
    return true;
  }
  return minmax->second >= min_value;
}

// Select parquet row groups that overlap crop bounds and optional quality filters.
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
        keep = row_group_numeric_range_overlaps(
            *row_group, static_cast<int>(*z_col), bounds.zmin, bounds.zmax);
      }
    }

    if (keep && min_qv >= 0.0 && qv_col.has_value()) {
      keep = row_group_numeric_max_at_least(*row_group, static_cast<int>(*qv_col), min_qv);
    }

    if (keep) {
      row_groups.push_back(rg);
    }
  }

  return row_groups;
}

// Drop cell-table rows for cells that have no loaded transcripts after filtering.
void keep_loaded_cells_only(XeniumBundleData& data) {
  if (data.cells.cell_ids.empty() || data.transcripts.orig_cell_id.empty()) {
    return;
  }
  std::unordered_map<std::string, std::size_t> keep;
  keep.reserve(data.transcripts.orig_cell_id.size());
  for (const auto& cell_id : data.transcripts.orig_cell_id) {
    keep.emplace(cell_id, keep.size());
  }

  CellTable filtered;
  std::vector<std::string> crop_ids;
  filtered.cell_ids.reserve(data.cells.cell_ids.size());
  filtered.sample_ids.reserve(data.cells.sample_ids.size());
  filtered.fov_ids.reserve(data.cells.fov_ids.size());
  filtered.cell_types.reserve(data.cells.cell_types.size());
  filtered.centroid_x.reserve(data.cells.centroid_x.size());
  filtered.centroid_y.reserve(data.cells.centroid_y.size());
  filtered.centroid_z.reserve(data.cells.centroid_z.size());
  crop_ids.reserve(data.cell_crop_ids.size());

  for (std::size_t i = 0; i < data.cells.cell_ids.size(); ++i) {
    if (keep.find(data.cells.cell_ids[i]) == keep.end()) {
      continue;
    }
    filtered.cell_ids.push_back(data.cells.cell_ids[i]);
    if (i < data.cells.sample_ids.size()) filtered.sample_ids.push_back(data.cells.sample_ids[i]);
    if (i < data.cells.fov_ids.size()) filtered.fov_ids.push_back(data.cells.fov_ids[i]);
    if (i < data.cells.cell_types.size()) filtered.cell_types.push_back(data.cells.cell_types[i]);
    if (i < data.cells.centroid_x.size()) filtered.centroid_x.push_back(data.cells.centroid_x[i]);
    if (i < data.cells.centroid_y.size()) filtered.centroid_y.push_back(data.cells.centroid_y[i]);
    if (i < data.cells.centroid_z.size()) filtered.centroid_z.push_back(data.cells.centroid_z[i]);
    if (i < data.cell_crop_ids.size()) crop_ids.push_back(data.cell_crop_ids[i]);
  }

  data.cells = std::move(filtered);
  data.cell_crop_ids = std::move(crop_ids);
}

// Load transcript rows from Xenium CSV while applying crop and quality filters.
void load_xenium_transcripts_csv(
    const std::string& path,
    const XeniumLoadOptions& options,
    XeniumBundleData& data) {
  LineReader reader(path);
  std::string line;
  if (!reader.next_line(line)) {
    throw std::runtime_error("Empty Xenium transcripts CSV: " + path);
  }
  const auto header = split_csv_line(line);
  std::unordered_map<std::string, std::size_t> header_index;
  for (std::size_t i = 0; i < header.size(); ++i) {
    header_index.emplace(header[i], i);
  }

  const auto gene_col =
      find_required_header(header_index, {"feature_name", "gene", "gene_name", "gene_id"}, "gene");
  const auto cell_col = find_required_header(header_index, {"cell_id", "cell"}, "cell id");
  const auto x_col = find_required_header(header_index, {"x_location", "x"}, "x");
  const auto y_col = find_required_header(header_index, {"y_location", "y"}, "y");
  const auto z_col = find_optional_header(header_index, {"z_location", "z"});
  const auto qv_col = find_optional_header(header_index, {"qv"});
  const auto txid_col = find_optional_header(header_index, {"transcript_id"});
  const auto fov_col = find_optional_header(header_index, {"fov_name", "fov_id", "fovn"});
  const auto overlaps_col = find_optional_header(header_index, {"overlaps_nucleus"});
  const auto nucleus_distance_col = find_optional_header(header_index, {"nucleus_distance"});
  const auto codeword_category_col = find_optional_header(header_index, {"codeword_category"});
  const auto is_gene_col = find_optional_header(header_index, {"is_gene"});
  const bool use_codeword_category_filter = !options.keep_non_gene && codeword_category_col.has_value();
  const bool use_is_gene_filter = !options.keep_non_gene && !use_codeword_category_filter && is_gene_col.has_value();
  const UnionBounds bounds = union_bounds(options.crops);

  while (reader.next_line(line)) {
    if (line.empty()) {
      continue;
    }
    const auto fields = split_csv_line(line);
    if (fields.size() <= std::max({gene_col, cell_col, x_col, y_col})) {
      continue;
    }

    const std::string cell_id = fields[cell_col];
    if (!options.keep_unassigned &&
        (cell_id.empty() || cell_id == "UNASSIGNED" || cell_id == "0" || cell_id == "cell_0")) {
      continue;
    }

    const double x = parse_double_or_default(fields[x_col], 0.0);
    const double y = parse_double_or_default(fields[y_col], 0.0);
    const double z = z_col.has_value() ? parse_double_or_default(fields[*z_col], 0.0) : 0.0;
    if (!passes_union_bounds(x, y, z, bounds)) {
      continue;
    }

    const double qv = qv_col.has_value() ? parse_double_or_default(fields[*qv_col], 0.0) : 0.0;
    if (options.min_qv >= 0.0 && qv_col.has_value() &&
        (!std::isfinite(qv) || qv < options.min_qv)) {
      continue;
    }

    const std::string crop_id = assign_crop_id(x, y, z, options.crops);
    if (!options.crops.empty() && crop_id.empty()) {
      continue;
    }

    const std::string gene = fields[gene_col];
    if (gene.empty()) {
      continue;
    }
    if (!options.keep_non_gene &&
        !xenium_is_biological_feature(
            gene,
            use_codeword_category_filter && fields.size() > *codeword_category_col
                ? fields[*codeword_category_col]
                : std::string(),
            use_is_gene_filter && fields.size() > *is_gene_col
                ? fields[*is_gene_col]
                : std::string())) {
      continue;
    }
    data.transcripts.gene_key.push_back(gene);
    data.transcripts.orig_cell_id.push_back(cell_id);
    data.transcripts.x.push_back(x);
    data.transcripts.y.push_back(y);
    data.transcripts.z.push_back(z);
    if (txid_col.has_value()) data.transcripts.transcript_ids.push_back(fields[*txid_col]);
    if (fov_col.has_value()) data.transcripts.fov_ids.push_back(fields[*fov_col]);
    if (qv_col.has_value()) data.transcripts.qv.push_back(qv);
    if (overlaps_col.has_value()) {
      data.transcripts.overlaps_nucleus.push_back(
          fields.size() > *overlaps_col ? parse_bool_or_missing(fields[*overlaps_col]) : -1);
    }
    if (nucleus_distance_col.has_value()) {
      data.transcripts.nucleus_distance.push_back(
          fields.size() > *nucleus_distance_col
              ? parse_double_or_default(fields[*nucleus_distance_col],
                                        std::numeric_limits<double>::quiet_NaN())
              : std::numeric_limits<double>::quiet_NaN());
    }
    data.transcript_crop_ids.push_back(crop_id);
  }
}

// Load cell centroids from Xenium CSV while applying crop filters.
void load_xenium_cells_csv(
    const std::string& path,
    const XeniumLoadOptions& options,
    XeniumBundleData& data) {
  if (!file_exists(path)) {
    return;
  }
  LineReader reader(path);
  std::string line;
  if (!reader.next_line(line)) {
    return;
  }
  const auto header = split_csv_line(line);
  std::unordered_map<std::string, std::size_t> header_index;
  for (std::size_t i = 0; i < header.size(); ++i) {
    header_index.emplace(header[i], i);
  }

  const auto cell_id_col = find_required_header(header_index, {"cell_id", "cell"}, "cell id");
  const auto x_col = find_required_header(header_index, {"x_centroid", "x"}, "cell x");
  const auto y_col = find_required_header(header_index, {"y_centroid", "y"}, "cell y");
  const auto z_col = find_optional_header(header_index, {"z_centroid", "z"});
  const auto fov_col = find_optional_header(header_index, {"fov_name", "fov_id", "fovn"});
  const auto cell_type_col = find_optional_header(header_index, {"cell_type", "celltype"});
  const UnionBounds bounds = union_bounds(options.crops);

  while (reader.next_line(line)) {
    if (line.empty()) {
      continue;
    }
    const auto fields = split_csv_line(line);
    if (fields.size() <= std::max({cell_id_col, x_col, y_col})) {
      continue;
    }
    const double x = parse_double_or_default(fields[x_col], 0.0);
    const double y = parse_double_or_default(fields[y_col], 0.0);
    const double z = z_col.has_value() ? parse_double_or_default(fields[*z_col], 0.0) : 0.0;
    if (!passes_union_bounds(x, y, z, bounds)) {
      continue;
    }
    const std::string crop_id = assign_crop_id(x, y, z, options.crops);
    if (!options.crops.empty() && crop_id.empty()) {
      continue;
    }

    if (fields[cell_id_col].empty()) {
      continue;
    }
    data.cells.cell_ids.push_back(fields[cell_id_col]);
    data.cells.centroid_x.push_back(x);
    data.cells.centroid_y.push_back(y);
    data.cells.centroid_z.push_back(z);
    if (fov_col.has_value()) data.cells.fov_ids.push_back(fields[*fov_col]);
    if (cell_type_col.has_value()) data.cells.cell_types.push_back(fields[*cell_type_col]);
    data.cell_crop_ids.push_back(crop_id);
  }
}

// Load transcript rows from Xenium parquet while pruning row groups by bounds.
void load_xenium_transcripts_parquet(
    const std::string& path,
    const XeniumLoadOptions& options,
    XeniumBundleData& data) {
  auto stage_start = std::chrono::steady_clock::now();
  emit_load_progress(options, "Opening Xenium transcripts parquet: " + path, stage_start);
  auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path));
  parquet::arrow::FileReaderBuilder builder;
  CELLADMIX_ARROW_CHECK_OK(builder.Open(input));
  auto reader = arrow_unwrap(builder.Build());
  reader->set_use_threads(true);
  reader->set_batch_size(65536);

  std::shared_ptr<arrow::Schema> schema;
  CELLADMIX_ARROW_CHECK_OK(reader->GetSchema(&schema));
  emit_load_progress(options, "Read transcripts parquet schema", stage_start);

  const auto gene_col =
      find_first_field(schema, {"feature_name", "gene", "gene_name", "gene_id"}, "gene");
  const auto cell_col = find_first_field(schema, {"cell_id", "cell"}, "cell id");
  const auto x_col = find_first_field(schema, {"x_location", "x"}, "x");
  const auto y_col = find_first_field(schema, {"y_location", "y"}, "y");
  const auto z_col = find_optional_field(schema, {"z_location", "z"});
  const auto qv_col = find_optional_field(schema, {"qv"});
  const auto txid_col = find_optional_field(schema, {"transcript_id"});
  const auto fov_col = find_optional_field(schema, {"fov_name", "fov_id", "fovn"});
  const auto overlaps_col = find_optional_field(schema, {"overlaps_nucleus"});
  const auto nucleus_distance_col = find_optional_field(schema, {"nucleus_distance"});
  const auto codeword_category_col = find_optional_field(schema, {"codeword_category"});
  const auto is_gene_col = find_optional_field(schema, {"is_gene"});
  const bool use_codeword_category_filter = !options.keep_non_gene && codeword_category_col.has_value();
  const bool use_is_gene_filter = !options.keep_non_gene && !use_codeword_category_filter && is_gene_col.has_value();
  if (overlaps_col.has_value() || nucleus_distance_col.has_value()) {
    std::string fields = "Reading Xenium nucleus columns:";
    if (overlaps_col.has_value()) fields += " overlaps_nucleus";
    if (nucleus_distance_col.has_value()) fields += " nucleus_distance";
    emit_load_progress(options, fields, stage_start);
  }

  const UnionBounds bounds = union_bounds(options.crops);
  std::vector<int> projected = {
      static_cast<int>(x_col),
      static_cast<int>(y_col)};
  if (z_col.has_value()) projected.push_back(static_cast<int>(*z_col));
  if (qv_col.has_value()) projected.push_back(static_cast<int>(*qv_col));
  projected.push_back(static_cast<int>(gene_col));
  projected.push_back(static_cast<int>(cell_col));
  if (txid_col.has_value()) projected.push_back(static_cast<int>(*txid_col));
  if (fov_col.has_value()) projected.push_back(static_cast<int>(*fov_col));
  if (overlaps_col.has_value()) projected.push_back(static_cast<int>(*overlaps_col));
  if (nucleus_distance_col.has_value()) projected.push_back(static_cast<int>(*nucleus_distance_col));
  if (use_codeword_category_filter) projected.push_back(static_cast<int>(*codeword_category_col));
  if (use_is_gene_filter) projected.push_back(static_cast<int>(*is_gene_col));

  const auto row_groups = select_row_groups_for_numeric_filters(
      path, x_col, y_col, z_col, qv_col, bounds, options.min_qv);
  if (row_groups.empty()) {
    emit_load_progress(options, "Selected 0 transcript row groups after metadata pruning", stage_start);
    return;
  }
  emit_load_progress(
      options,
      "Selected " + std::to_string(row_groups.size()) +
          " transcript row groups after metadata pruning",
      stage_start);

  auto batch_reader = arrow_unwrap(reader->GetRecordBatchReader(row_groups, projected));
  std::shared_ptr<arrow::RecordBatch> batch;
  std::int64_t scanned_rows = 0;
  std::int64_t kept_transcripts = 0;
  int batch_count = 0;
  constexpr std::int64_t progress_row_interval = 1000000;
  std::int64_t next_progress_rows = progress_row_interval;
  while (true) {
    CELLADMIX_ARROW_CHECK_OK(batch_reader->ReadNext(&batch));
    if (!batch) {
      break;
    }
    ++batch_count;
    scanned_rows += batch->num_rows();

    int col = 0;
    NumericArrayView x_view(batch->column(col++));
    NumericArrayView y_view(batch->column(col++));
    std::optional<NumericArrayView> z_view;
    if (z_col.has_value()) z_view.emplace(batch->column(col++));
    std::optional<NumericArrayView> qv_view;
    if (qv_col.has_value()) qv_view.emplace(batch->column(col++));
    const int gene_batch_col = col++;
    const int cell_batch_col = col++;
    const int txid_batch_col = txid_col.has_value() ? col++ : -1;
    const int fov_batch_col = fov_col.has_value() ? col++ : -1;
    const int overlaps_batch_col = overlaps_col.has_value() ? col++ : -1;
    const int nucleus_distance_batch_col = nucleus_distance_col.has_value() ? col++ : -1;
    const int codeword_category_batch_col = use_codeword_category_filter ? col++ : -1;
    const int is_gene_batch_col = use_is_gene_filter ? col++ : -1;

    std::vector<int64_t> kept_rows;
    std::vector<double> kept_x;
    std::vector<double> kept_y;
    std::vector<double> kept_z;
    std::vector<double> kept_qv;
    std::vector<std::string> kept_crop_ids;
    kept_rows.reserve(static_cast<std::size_t>(batch->num_rows()));
    kept_x.reserve(static_cast<std::size_t>(batch->num_rows()));
    kept_y.reserve(static_cast<std::size_t>(batch->num_rows()));
    kept_z.reserve(static_cast<std::size_t>(batch->num_rows()));
    if (qv_view.has_value()) {
      kept_qv.reserve(static_cast<std::size_t>(batch->num_rows()));
    }
    kept_crop_ids.reserve(static_cast<std::size_t>(batch->num_rows()));

    for (int64_t i = 0; i < batch->num_rows(); ++i) {
      const double x = x_view.value(i);
      const double y = y_view.value(i);
      const double z = z_view.has_value() ? z_view->value(i) : 0.0;
      if (!passes_union_bounds(x, y, z, bounds)) {
        continue;
      }

      const double qv = qv_view.has_value() ? qv_view->value(i) : 0.0;
      if (options.min_qv >= 0.0 && qv_view.has_value() &&
          (!std::isfinite(qv) || qv < options.min_qv)) {
        continue;
      }

      const std::string crop_id = assign_crop_id(x, y, z, options.crops);
      if (!options.crops.empty() && crop_id.empty()) {
        continue;
      }

      kept_rows.push_back(i);
      kept_x.push_back(x);
      kept_y.push_back(y);
      kept_z.push_back(z);
      if (qv_view.has_value()) kept_qv.push_back(qv);
      kept_crop_ids.push_back(std::move(crop_id));
    }

    if (kept_rows.empty()) {
      if (options.progress &&
          (scanned_rows >= next_progress_rows || batch_count == 1)) {
        emit_load_progress(
            options,
            "Scanned transcript parquet batches: " + std::to_string(batch_count) +
                ", rows=" + std::to_string(scanned_rows) +
                ", kept=" + std::to_string(kept_transcripts),
            stage_start);
        while (next_progress_rows <= scanned_rows) {
          next_progress_rows += progress_row_interval;
        }
      }
      continue;
    }

    data.transcripts.gene_key.reserve(data.transcripts.gene_key.size() + kept_rows.size());
    data.transcripts.orig_cell_id.reserve(data.transcripts.orig_cell_id.size() + kept_rows.size());
    data.transcripts.x.reserve(data.transcripts.x.size() + kept_rows.size());
    data.transcripts.y.reserve(data.transcripts.y.size() + kept_rows.size());
    data.transcripts.z.reserve(data.transcripts.z.size() + kept_rows.size());
    if (qv_view.has_value()) data.transcripts.qv.reserve(data.transcripts.qv.size() + kept_rows.size());
    if (txid_col.has_value()) {
      data.transcripts.transcript_ids.reserve(
          data.transcripts.transcript_ids.size() + kept_rows.size());
    }
    if (fov_col.has_value()) {
      data.transcripts.fov_ids.reserve(data.transcripts.fov_ids.size() + kept_rows.size());
    }
    if (overlaps_col.has_value()) {
      data.transcripts.overlaps_nucleus.reserve(
          data.transcripts.overlaps_nucleus.size() + kept_rows.size());
    }
    if (nucleus_distance_col.has_value()) {
      data.transcripts.nucleus_distance.reserve(
          data.transcripts.nucleus_distance.size() + kept_rows.size());
    }
    data.transcript_crop_ids.reserve(data.transcript_crop_ids.size() + kept_rows.size());

    StringArrayView gene_view(batch->column(gene_batch_col));
    StringArrayView cell_view(batch->column(cell_batch_col));
    std::optional<StringArrayView> txid_view;
    if (txid_batch_col >= 0) txid_view.emplace(batch->column(txid_batch_col));
    std::optional<StringArrayView> fov_view;
    if (fov_batch_col >= 0) fov_view.emplace(batch->column(fov_batch_col));
    std::optional<StringArrayView> overlaps_view;
    if (overlaps_batch_col >= 0) overlaps_view.emplace(batch->column(overlaps_batch_col));
    std::optional<NumericArrayView> nucleus_distance_view;
    if (nucleus_distance_batch_col >= 0) {
      nucleus_distance_view.emplace(batch->column(nucleus_distance_batch_col));
    }
    std::optional<StringArrayView> codeword_category_view;
    if (codeword_category_batch_col >= 0) {
      codeword_category_view.emplace(batch->column(codeword_category_batch_col));
    }
    std::optional<StringArrayView> is_gene_view;
    if (is_gene_batch_col >= 0) {
      is_gene_view.emplace(batch->column(is_gene_batch_col));
    }

    const std::size_t before_append = data.transcripts.gene_key.size();
    for (std::size_t kept_idx = 0; kept_idx < kept_rows.size(); ++kept_idx) {
      const int64_t row = kept_rows[kept_idx];
      const std::string cell_id = cell_view.value(row);
      if (!options.keep_unassigned && is_unassigned_cell_id(cell_id)) {
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

      data.transcripts.gene_key.push_back(gene);
      data.transcripts.orig_cell_id.push_back(cell_id);
      data.transcripts.x.push_back(kept_x[kept_idx]);
      data.transcripts.y.push_back(kept_y[kept_idx]);
      data.transcripts.z.push_back(kept_z[kept_idx]);
      if (txid_view.has_value()) data.transcripts.transcript_ids.push_back(txid_view->value(row));
      if (fov_view.has_value()) data.transcripts.fov_ids.push_back(fov_view->value(row));
      if (qv_view.has_value()) data.transcripts.qv.push_back(kept_qv[kept_idx]);
      if (overlaps_view.has_value()) {
        data.transcripts.overlaps_nucleus.push_back(
            parse_bool_or_missing(overlaps_view->value(row)));
      }
      if (nucleus_distance_view.has_value()) {
        data.transcripts.nucleus_distance.push_back(nucleus_distance_view->value(row));
      }
      data.transcript_crop_ids.push_back(kept_crop_ids[kept_idx]);
    }
    kept_transcripts += static_cast<std::int64_t>(data.transcripts.gene_key.size() - before_append);
    if (options.progress &&
        (scanned_rows >= next_progress_rows || batch_count == 1)) {
      emit_load_progress(
          options,
          "Scanned transcript parquet batches: " + std::to_string(batch_count) +
              ", rows=" + std::to_string(scanned_rows) +
              ", kept=" + std::to_string(kept_transcripts),
          stage_start);
      while (next_progress_rows <= scanned_rows) {
        next_progress_rows += progress_row_interval;
      }
    }
  }
  emit_load_progress(
      options,
      "Finished transcript parquet load: rows=" + std::to_string(scanned_rows) +
          ", kept=" + std::to_string(kept_transcripts),
      stage_start);
}

// Load cell centroids from Xenium parquet while pruning row groups by bounds.
void load_xenium_cells_parquet(
    const std::string& path,
    const XeniumLoadOptions& options,
    XeniumBundleData& data) {
  auto stage_start = std::chrono::steady_clock::now();
  emit_load_progress(options, "Opening Xenium cells parquet: " + path, stage_start);
  auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path));
  parquet::arrow::FileReaderBuilder builder;
  CELLADMIX_ARROW_CHECK_OK(builder.Open(input));
  auto reader = arrow_unwrap(builder.Build());
  reader->set_use_threads(true);
  reader->set_batch_size(65536);

  std::shared_ptr<arrow::Schema> schema;
  CELLADMIX_ARROW_CHECK_OK(reader->GetSchema(&schema));
  emit_load_progress(options, "Read cells parquet schema", stage_start);

  const auto cell_id_col = find_first_field(schema, {"cell_id", "cell"}, "cell id");
  const auto x_col = find_first_field(schema, {"x_centroid", "x"}, "cell x");
  const auto y_col = find_first_field(schema, {"y_centroid", "y"}, "cell y");
  const auto z_col = find_optional_field(schema, {"z_centroid", "z"});
  const auto fov_col = find_optional_field(schema, {"fov_name", "fov_id", "fovn"});
  const auto cell_type_col = find_optional_field(schema, {"cell_type", "celltype"});

  const UnionBounds bounds = union_bounds(options.crops);
  std::vector<int> projected = {
      static_cast<int>(x_col),
      static_cast<int>(y_col)};
  if (z_col.has_value()) projected.push_back(static_cast<int>(*z_col));
  projected.push_back(static_cast<int>(cell_id_col));
  if (fov_col.has_value()) projected.push_back(static_cast<int>(*fov_col));
  if (cell_type_col.has_value()) projected.push_back(static_cast<int>(*cell_type_col));

  const auto row_groups =
      select_row_groups_for_numeric_filters(path, x_col, y_col, z_col, std::nullopt, bounds, -1.0);
  if (row_groups.empty()) {
    emit_load_progress(options, "Selected 0 cell row groups after metadata pruning", stage_start);
    return;
  }
  emit_load_progress(
      options,
      "Selected " + std::to_string(row_groups.size()) +
          " cell row groups after metadata pruning",
      stage_start);

  auto batch_reader = arrow_unwrap(reader->GetRecordBatchReader(row_groups, projected));
  std::shared_ptr<arrow::RecordBatch> batch;
  std::int64_t scanned_rows = 0;
  std::int64_t kept_cells = 0;
  int batch_count = 0;
  constexpr std::int64_t progress_row_interval = 500000;
  std::int64_t next_progress_rows = progress_row_interval;
  while (true) {
    CELLADMIX_ARROW_CHECK_OK(batch_reader->ReadNext(&batch));
    if (!batch) {
      break;
    }
    ++batch_count;
    scanned_rows += batch->num_rows();

    int col = 0;
    NumericArrayView x_view(batch->column(col++));
    NumericArrayView y_view(batch->column(col++));
    std::optional<NumericArrayView> z_view;
    if (z_col.has_value()) z_view.emplace(batch->column(col++));
    const int cell_id_batch_col = col++;
    const int fov_batch_col = fov_col.has_value() ? col++ : -1;
    const int cell_type_batch_col = cell_type_col.has_value() ? col++ : -1;

    std::vector<int64_t> kept_rows;
    std::vector<double> kept_x;
    std::vector<double> kept_y;
    std::vector<double> kept_z;
    std::vector<std::string> kept_crop_ids;
    kept_rows.reserve(static_cast<std::size_t>(batch->num_rows()));
    kept_x.reserve(static_cast<std::size_t>(batch->num_rows()));
    kept_y.reserve(static_cast<std::size_t>(batch->num_rows()));
    kept_z.reserve(static_cast<std::size_t>(batch->num_rows()));
    kept_crop_ids.reserve(static_cast<std::size_t>(batch->num_rows()));

    for (int64_t i = 0; i < batch->num_rows(); ++i) {
      const double x = x_view.value(i);
      const double y = y_view.value(i);
      const double z = z_view.has_value() ? z_view->value(i) : 0.0;
      if (!passes_union_bounds(x, y, z, bounds)) {
        continue;
      }
      const std::string crop_id = assign_crop_id(x, y, z, options.crops);
      if (!options.crops.empty() && crop_id.empty()) {
        continue;
      }

      kept_rows.push_back(i);
      kept_x.push_back(x);
      kept_y.push_back(y);
      kept_z.push_back(z);
      kept_crop_ids.push_back(std::move(crop_id));
    }

    if (kept_rows.empty()) {
      if (options.progress &&
          (scanned_rows >= next_progress_rows || batch_count == 1)) {
        emit_load_progress(
            options,
            "Scanned cell parquet batches: " + std::to_string(batch_count) +
                ", rows=" + std::to_string(scanned_rows) +
                ", kept=" + std::to_string(kept_cells),
            stage_start);
        while (next_progress_rows <= scanned_rows) {
          next_progress_rows += progress_row_interval;
        }
      }
      continue;
    }

    data.cells.cell_ids.reserve(data.cells.cell_ids.size() + kept_rows.size());
    data.cells.centroid_x.reserve(data.cells.centroid_x.size() + kept_rows.size());
    data.cells.centroid_y.reserve(data.cells.centroid_y.size() + kept_rows.size());
    data.cells.centroid_z.reserve(data.cells.centroid_z.size() + kept_rows.size());
    if (fov_col.has_value()) data.cells.fov_ids.reserve(data.cells.fov_ids.size() + kept_rows.size());
    if (cell_type_col.has_value()) {
      data.cells.cell_types.reserve(data.cells.cell_types.size() + kept_rows.size());
    }
    data.cell_crop_ids.reserve(data.cell_crop_ids.size() + kept_rows.size());

    StringArrayView cell_id_view(batch->column(cell_id_batch_col));
    std::optional<StringArrayView> fov_view;
    if (fov_batch_col >= 0) fov_view.emplace(batch->column(fov_batch_col));
    std::optional<StringArrayView> cell_type_view;
    if (cell_type_batch_col >= 0) cell_type_view.emplace(batch->column(cell_type_batch_col));

    const std::size_t before_append = data.cells.cell_ids.size();
    for (std::size_t kept_idx = 0; kept_idx < kept_rows.size(); ++kept_idx) {
      const int64_t row = kept_rows[kept_idx];
      const std::string cell_id = cell_id_view.value(row);
      if (cell_id.empty()) {
        continue;
      }
      data.cells.cell_ids.push_back(cell_id);
      data.cells.centroid_x.push_back(kept_x[kept_idx]);
      data.cells.centroid_y.push_back(kept_y[kept_idx]);
      data.cells.centroid_z.push_back(kept_z[kept_idx]);
      if (fov_view.has_value()) data.cells.fov_ids.push_back(fov_view->value(row));
      if (cell_type_view.has_value()) data.cells.cell_types.push_back(cell_type_view->value(row));
      data.cell_crop_ids.push_back(kept_crop_ids[kept_idx]);
    }
    kept_cells += static_cast<std::int64_t>(data.cells.cell_ids.size() - before_append);
    if (options.progress &&
        (scanned_rows >= next_progress_rows || batch_count == 1)) {
      emit_load_progress(
          options,
          "Scanned cell parquet batches: " + std::to_string(batch_count) +
              ", rows=" + std::to_string(scanned_rows) +
              ", kept=" + std::to_string(kept_cells),
          stage_start);
      while (next_progress_rows <= scanned_rows) {
        next_progress_rows += progress_row_interval;
      }
    }
  }
  emit_load_progress(
      options,
      "Finished cell parquet load: rows=" + std::to_string(scanned_rows) +
          ", kept=" + std::to_string(kept_cells),
      stage_start);
}

}  // namespace

// Parse the lightweight Xenium experiment manifest.
XeniumManifest read_xenium_manifest(const std::string& path) {
  const auto text = slurp(path);

  XeniumManifest manifest;
  manifest.run_name = extract_string(text, "run_name");
  manifest.pixel_size = extract_double(text, "pixel_size", 1.0);
  manifest.z_step_size = extract_double(text, "z_step_size", 0.0);
  manifest.morphology_filepath = extract_string(text, "morphology_filepath");
  manifest.morphology_focus_filepath = extract_string(text, "morphology_focus_filepath");

  const auto transcripts_csv = extract_string(text, "transcripts_csv_filepath");
  if (!transcripts_csv.empty()) {
    manifest.transcripts_csv_path = transcripts_csv;
  }
  const auto transcripts_parquet = extract_string(text, "transcripts_parquet_filepath");
  if (!transcripts_parquet.empty()) {
    manifest.transcripts_parquet_path = transcripts_parquet;
  }

  return manifest;
}

// Load one Xenium bundle into TranscriptTable and CellTable objects.
XeniumBundleData load_xenium_bundle(
    const std::string& manifest_path_or_dir,
    const XeniumLoadOptions& options) {
  auto stage_start = std::chrono::steady_clock::now();
  XeniumBundleData data;
  const std::string manifest_path = resolve_manifest_path(manifest_path_or_dir);
  data.manifest = read_xenium_manifest(manifest_path);
  emit_load_progress(options, "Read Xenium manifest: " + manifest_path, stage_start);
  const std::string dataset_dir = dirname_of(manifest_path);

  const std::string parquet_path =
      join_path(dataset_dir, data.manifest.transcripts_parquet_path.empty()
                                 ? "transcripts.parquet"
                                 : data.manifest.transcripts_parquet_path);
  const std::string csv_path =
      join_path(dataset_dir, data.manifest.transcripts_csv_path.empty()
                                 ? "transcripts.csv.gz"
                                 : data.manifest.transcripts_csv_path);

  const bool have_parquet = file_exists(parquet_path);
  const bool have_csv = file_exists(csv_path);
  if (options.prefer_parquet && have_parquet) {
    emit_load_progress(options, "Loading transcripts from parquet", stage_start);
    load_xenium_transcripts_parquet(parquet_path, options, data);
    data.used_parquet = true;
  } else if (have_csv) {
    emit_load_progress(options, "Loading transcripts from CSV", stage_start);
    load_xenium_transcripts_csv(csv_path, options, data);
  } else if (have_parquet) {
    emit_load_progress(options, "Loading transcripts from parquet", stage_start);
    load_xenium_transcripts_parquet(parquet_path, options, data);
    data.used_parquet = true;
  } else {
    throw std::runtime_error(
        "Could not locate Xenium transcripts.parquet or transcripts.csv.gz next to manifest: " +
        manifest_path);
  }

  stage_start = std::chrono::steady_clock::now();
  data.transcripts.finalize();
  emit_load_progress(
      options,
      "Finalized transcript table: " + std::to_string(data.transcripts.size()) +
          " transcripts, " + std::to_string(data.transcripts.num_cells()) +
          " cells, " + std::to_string(data.transcripts.num_genes()) + " genes",
      stage_start);

  if (options.read_cells) {
    const std::string cells_parquet_path = join_path(dataset_dir, "cells.parquet");
    const std::string cells_csv_path = join_path(dataset_dir, "cells.csv.gz");
    const bool have_cells_parquet = file_exists(cells_parquet_path);
    const bool have_cells_csv = file_exists(cells_csv_path);
    if (options.prefer_parquet && have_cells_parquet) {
      emit_load_progress(options, "Loading cells from parquet", stage_start);
      load_xenium_cells_parquet(cells_parquet_path, options, data);
    } else if (have_cells_csv) {
      emit_load_progress(options, "Loading cells from CSV", stage_start);
      load_xenium_cells_csv(cells_csv_path, options, data);
    } else if (have_cells_parquet) {
      emit_load_progress(options, "Loading cells from parquet", stage_start);
      load_xenium_cells_parquet(cells_parquet_path, options, data);
    }
    stage_start = std::chrono::steady_clock::now();
    keep_loaded_cells_only(data);
    emit_load_progress(
        options,
        "Filtered cell metadata to loaded transcripts: " +
            std::to_string(data.cells.size()) + " cells",
        stage_start);
  }

  return data;
}

}  // namespace celladmix
