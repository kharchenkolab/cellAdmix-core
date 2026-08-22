#include <Rcpp.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <arrow/result.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/file_reader.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <zlib.h>

#include "celladmix/clustering.hpp"
#include "celladmix/coherence.hpp"
#include "celladmix/generative.hpp"
#include "celladmix/domains.hpp"
#include "celladmix/input_store.hpp"
#include "celladmix/membrane.hpp"
#include "celladmix/ncv.hpp"
#include "celladmix/nmf_euclidean.hpp"
#include "celladmix/nmf_kl.hpp"
#include "celladmix/pipeline.hpp"
#include "celladmix/pipeline_store.hpp"
#include "celladmix/report_embedding.hpp"
#include "celladmix/run_store.hpp"
#include "celladmix/simulation.hpp"
#include "celladmix/spatial.hpp"
#include "celladmix/tabular.hpp"
#include "celladmix/workflow.hpp"
#include "celladmix/xenium.hpp"

using namespace Rcpp;

namespace {

// Rcpp bridge between the R workflow API and the C++ core implementation.

NumericMatrix matrix_to_r(const celladmix::DenseMatrix& matrix);

// Read an optional single string from an SEXP argument.
std::optional<std::string> optional_string_sexp(SEXP value) {
  if (Rf_isNull(value)) {
    return std::nullopt;
  }
  CharacterVector chr(value);
  if (chr.size() == 0 || CharacterVector::is_na(chr[0])) {
    return std::nullopt;
  }
  return as<std::string>(chr[0]);
}

// Read an optional character vector, dropping NA/empty values and duplicates.
std::vector<std::string> optional_string_vector_sexp(SEXP value) {
  std::vector<std::string> out;
  if (Rf_isNull(value)) {
    return out;
  }
  CharacterVector chr(value);
  out.reserve(static_cast<std::size_t>(chr.size()));
  for (R_xlen_t i = 0; i < chr.size(); ++i) {
    if (CharacterVector::is_na(chr[i])) {
      continue;
    }
    const std::string entry = as<std::string>(chr[i]);
    if (entry.empty()) {
      continue;
    }
    if (std::find(out.begin(), out.end(), entry) == out.end()) {
      out.push_back(entry);
    }
  }
  return out;
}

// Replace cell-type labels from a named character vector where names are cell ids.
void apply_cell_type_overrides(celladmix::CellTable& cells, SEXP cell_types_sexp) {
  if (Rf_isNull(cell_types_sexp)) {
    return;
  }
  CharacterVector cell_types(cell_types_sexp);
  CharacterVector names = cell_types.names();
  if (names.size() != cell_types.size()) {
    stop("cell_type override must be a named character vector with cell ids as names");
  }
  std::unordered_map<std::string, std::string> by_cell;
  by_cell.reserve(static_cast<std::size_t>(cell_types.size()));
  for (R_xlen_t i = 0; i < cell_types.size(); ++i) {
    if (CharacterVector::is_na(names[i]) || CharacterVector::is_na(cell_types[i])) {
      continue;
    }
    const std::string cell_id = as<std::string>(names[i]);
    const std::string cell_type = as<std::string>(cell_types[i]);
    if (!cell_id.empty() && !cell_type.empty()) {
      by_cell[cell_id] = cell_type;
    }
  }
  if (cells.cell_types.size() != cells.cell_ids.size()) {
    cells.cell_types.assign(cells.cell_ids.size(), "");
  }
  for (std::size_t i = 0; i < cells.cell_ids.size(); ++i) {
    const auto it = by_cell.find(cells.cell_ids[i]);
    if (it != by_cell.end()) {
      cells.cell_types[i] = it->second;
    }
  }
}

// Read an optional integer from an SEXP argument.
int optional_int_sexp(SEXP value, int default_value = -1) {
  const int parsed = as<int>(value);
  return IntegerVector::is_na(parsed) ? default_value : parsed;
}

// Read an optional seed from an SEXP argument.
unsigned int optional_seed_sexp(SEXP value, unsigned int default_value = 1U) {
  const int parsed = as<int>(value);
  return IntegerVector::is_na(parsed) ? default_value : static_cast<unsigned int>(parsed);
}

// Unwrap an Arrow result and attach context to any error.
template <class T>
T arrow_unwrap(arrow::Result<T>&& result, const char* context) {
  if (!result.ok()) {
    throw std::runtime_error(std::string(context) + ": " + result.status().ToString());
  }
  return std::move(*result);
}

// Raise a contextual error when an Arrow status is not OK.
void arrow_check(const arrow::Status& status, const char* context) {
  if (!status.ok()) {
    throw std::runtime_error(std::string(context) + ": " + status.ToString());
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

std::shared_ptr<arrow::Table> read_parquet_table(const std::string& path) {
  auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path), "Open parquet file");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open parquet reader");
  auto reader = arrow_unwrap(builder.Build(), "Build parquet reader");
  std::shared_ptr<arrow::Table> table;
  arrow_check(reader->ReadTable(&table), "Read parquet table");
  return table;
}

std::vector<std::string> extract_string_column(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& name) {
  const int idx = table->schema()->GetFieldIndex(name);
  if (idx < 0) stop("Missing parquet column: %s", name);
  const auto chunked = table->column(idx);
  std::vector<std::string> out(static_cast<std::size_t>(chunked->length()));
  int64_t offset = 0;
  for (int c = 0; c < chunked->num_chunks(); ++c) {
    const auto chunk = chunked->chunk(c);
    if (chunk->type_id() == arrow::Type::STRING) {
      const auto arr = std::static_pointer_cast<arrow::StringArray>(chunk);
      for (int64_t i = 0; i < chunk->length(); ++i) {
        out[static_cast<std::size_t>(offset + i)] = arr->GetString(i);
      }
    } else if (chunk->type_id() == arrow::Type::LARGE_STRING) {
      const auto arr = std::static_pointer_cast<arrow::LargeStringArray>(chunk);
      for (int64_t i = 0; i < chunk->length(); ++i) {
        out[static_cast<std::size_t>(offset + i)] = arr->GetString(i);
      }
    } else {
      stop("Unsupported string column type for %s", name);
    }
    offset += chunk->length();
  }
  return out;
}

std::vector<double> extract_double_column(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& name) {
  const int idx = table->schema()->GetFieldIndex(name);
  if (idx < 0) stop("Missing parquet column: %s", name);
  const auto chunked = table->column(idx);
  std::vector<double> out(static_cast<std::size_t>(chunked->length()));
  int64_t offset = 0;
  for (int c = 0; c < chunked->num_chunks(); ++c) {
    const auto chunk = chunked->chunk(c);
    if (chunk->type_id() == arrow::Type::DOUBLE) {
      const auto arr = std::static_pointer_cast<arrow::DoubleArray>(chunk);
      for (int64_t i = 0; i < chunk->length(); ++i) {
        out[static_cast<std::size_t>(offset + i)] = arr->Value(i);
      }
    } else if (chunk->type_id() == arrow::Type::FLOAT) {
      const auto arr = std::static_pointer_cast<arrow::FloatArray>(chunk);
      for (int64_t i = 0; i < chunk->length(); ++i) {
        out[static_cast<std::size_t>(offset + i)] = static_cast<double>(arr->Value(i));
      }
    } else {
      stop("Unsupported numeric column type for %s", name);
    }
    offset += chunk->length();
  }
  return out;
}

std::vector<int> extract_int_column(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& name) {
  const int idx = table->schema()->GetFieldIndex(name);
  if (idx < 0) stop("Missing parquet column: %s", name);
  const auto chunked = table->column(idx);
  std::vector<int> out(static_cast<std::size_t>(chunked->length()));
  int64_t offset = 0;
  for (int c = 0; c < chunked->num_chunks(); ++c) {
    const auto chunk = chunked->chunk(c);
    if (chunk->type_id() == arrow::Type::INT32) {
      const auto arr = std::static_pointer_cast<arrow::Int32Array>(chunk);
      for (int64_t i = 0; i < chunk->length(); ++i) {
        out[static_cast<std::size_t>(offset + i)] = arr->Value(i);
      }
    } else if (chunk->type_id() == arrow::Type::INT64) {
      const auto arr = std::static_pointer_cast<arrow::Int64Array>(chunk);
      for (int64_t i = 0; i < chunk->length(); ++i) {
        out[static_cast<std::size_t>(offset + i)] = static_cast<int>(arr->Value(i));
      }
    } else {
      stop("Unsupported integer column type for %s", name);
    }
    offset += chunk->length();
  }
  return out;
}

DataFrame parquet_table_to_df(const std::string& path) {
  const auto table = read_parquet_table(path);
  List out;
  const auto schema = table->schema();
  for (int col = 0; col < schema->num_fields(); ++col) {
    const auto& field = schema->field(col);
    const auto& name = field->name();
    switch (field->type()->id()) {
      case arrow::Type::STRING:
      case arrow::Type::LARGE_STRING:
        out[name] = wrap(extract_string_column(table, name));
        break;
      case arrow::Type::DOUBLE:
      case arrow::Type::FLOAT:
        out[name] = wrap(extract_double_column(table, name));
        break;
      case arrow::Type::INT32:
      case arrow::Type::INT64:
        out[name] = wrap(extract_int_column(table, name));
        break;
      default:
        stop("Unsupported parquet column type in %s: %s", path, name);
    }
  }
  return DataFrame(out);
}

bool has_suffix(const std::string& value, const std::string& suffix) {
  return value.size() >= suffix.size() &&
      value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string lower_string(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
    return static_cast<char>(std::tolower(ch));
  });
  return value;
}

std::string logical_extension(const std::string& path) {
  std::string lower = lower_string(path);
  if (has_suffix(lower, ".gz")) {
    lower.resize(lower.size() - 3);
  }
  return std::filesystem::path(lower).extension().string();
}

std::vector<std::string> split_csv_header(const std::string& line, char delimiter = ',') {
  std::vector<std::string> out;
  std::string current;
  bool in_quotes = false;
  for (const char ch : line) {
    if (ch == '"') {
      in_quotes = !in_quotes;
      continue;
    }
    if (ch == delimiter && !in_quotes) {
      out.push_back(current);
      current.clear();
      continue;
    }
    if (ch != '\r' && ch != '\n') {
      current.push_back(ch);
    }
  }
  out.push_back(current);
  return out;
}

// Read a small CSV header line while handling gzip-compressed files.
std::vector<std::string> read_csv_header_columns(const std::string& path) {
  std::string line;
  if (path.size() >= 3 && path.substr(path.size() - 3) == ".gz") {
    gzFile file = gzopen(path.c_str(), "rb");
    if (file == nullptr) {
      throw std::runtime_error("failed to open gzip csv: " + path);
    }
    char buffer[16384];
    const char* read = gzgets(file, buffer, sizeof(buffer));
    if (read != nullptr) {
      line = buffer;
    }
    gzclose(file);
  } else {
    std::ifstream input(path);
    std::getline(input, line);
  }
  if (line.empty()) {
    return {};
  }
  const char delimiter = logical_extension(path) == ".tsv" ? '\t' : ',';
  return split_csv_header(line, delimiter);
}

// Read parquet column names without materializing the full table.
std::vector<std::string> read_parquet_column_names(const std::string& path) {
  auto reader = parquet::ParquetFileReader::OpenFile(path, false);
  const auto metadata = reader->metadata();
  const auto schema = metadata->schema();
  std::vector<std::string> out;
  out.reserve(static_cast<std::size_t>(schema->num_columns()));
  for (int i = 0; i < schema->num_columns(); ++i) {
    out.push_back(schema->Column(i)->name());
  }
  return out;
}

// Package one probed input file for the R-side manifest/preflight summary.
List file_probe_to_r(
    const std::string& role,
    const std::string& path,
    const std::vector<std::string>& columns = {}) {
  return List::create(
      _["role"] = role,
      _["path"] = path,
      _["exists"] = std::filesystem::exists(std::filesystem::path(path)),
      _["columns"] = wrap(columns));
}

// Load per-cell strata labels from a cluster or annotation parquet sidecar.
std::vector<std::string> load_cell_label_strata_for_table(
    const std::string& path,
    const celladmix::TranscriptTable& table) {
  const auto label_table = read_parquet_table(path);
  const auto cell_id = extract_string_column(label_table, "cell_id");
  std::vector<std::string> labels;
  const int cluster_idx = label_table->schema()->GetFieldIndex("cluster");
  const int label_idx = label_table->schema()->GetFieldIndex("label");
  if (cluster_idx >= 0) {
    const auto cluster = extract_int_column(label_table, "cluster");
    labels.reserve(cluster.size());
    for (const int value : cluster) {
      labels.push_back("cluster_" + std::to_string(value));
    }
  } else if (label_idx >= 0) {
    labels = extract_string_column(label_table, "label");
  } else {
    throw std::runtime_error("cluster annotation parquet must contain cluster or label column");
  }

  std::unordered_map<std::string, std::string> lookup;
  lookup.reserve(cell_id.size());
  for (std::size_t i = 0; i < cell_id.size(); ++i) {
    lookup[cell_id[i]] = labels[i];
  }

  std::vector<std::string> out(table.num_cells(), "");
  for (std::size_t i = 0; i < table.cells.size(); ++i) {
    const auto it = lookup.find(table.cells[i]);
    if (it != lookup.end()) {
      out[i] = it->second;
    }
  }
  return out;
}

std::vector<std::string> load_cell_label_strata_for_cells(
    const std::string& path,
    const std::vector<std::string>& cells) {
  const auto label_table = read_parquet_table(path);
  const auto cell_id = extract_string_column(label_table, "cell_id");
  std::vector<std::string> labels;
  const int cluster_idx = label_table->schema()->GetFieldIndex("cluster");
  const int label_idx = label_table->schema()->GetFieldIndex("label");
  if (cluster_idx >= 0) {
    const auto cluster = extract_int_column(label_table, "cluster");
    labels.reserve(cluster.size());
    for (const int value : cluster) {
      labels.push_back("cluster_" + std::to_string(value));
    }
  } else if (label_idx >= 0) {
    labels = extract_string_column(label_table, "label");
  } else {
    throw std::runtime_error("cluster annotation parquet must contain cluster or label column");
  }

  std::unordered_map<std::string, std::string> lookup;
  lookup.reserve(cell_id.size());
  for (std::size_t i = 0; i < cell_id.size(); ++i) {
    lookup[cell_id[i]] = labels[i];
  }

  std::vector<std::string> out(cells.size(), "");
  for (std::size_t i = 0; i < cells.size(); ++i) {
    const auto it = lookup.find(cells[i]);
    if (it != lookup.end()) {
      out[i] = it->second;
    }
  }
  return out;
}

// Check whether a data frame exposes one named column.
bool data_frame_has_column(const DataFrame& df, const std::string& name) {
  const CharacterVector names = df.names();
  for (const auto& value : names) {
    if (as<std::string>(value) == name) {
      return true;
    }
  }
  return false;
}

// Extract a required string column from an R data frame.
std::vector<std::string> data_frame_string_column(
    const DataFrame& df,
    const std::string& name) {
  if (!data_frame_has_column(df, name)) {
    stop("Missing required data frame column: %s", name);
  }
  return as<std::vector<std::string>>(df[name]);
}

// Extract a required numeric column from an R data frame.
std::vector<double> data_frame_numeric_column(
    const DataFrame& df,
    const std::string& name) {
  if (!data_frame_has_column(df, name)) {
    stop("Missing required data frame column: %s", name);
  }
  return as<std::vector<double>>(df[name]);
}

// Convert a transcript data frame into the core TranscriptTable plus stable row ids.
std::pair<celladmix::TranscriptTable, std::vector<std::string>> transcript_table_from_df(
    const DataFrame& df) {
  celladmix::TranscriptTable table;
  table.gene_key = data_frame_string_column(df, "gene");
  table.orig_cell_id = data_frame_string_column(df, "cell");
  table.x = data_frame_numeric_column(df, "x");
  table.y = data_frame_numeric_column(df, "y");
  if (data_frame_has_column(df, "z")) {
    table.z = data_frame_numeric_column(df, "z");
  }
  if (data_frame_has_column(df, "transcript_id")) {
    table.transcript_ids = data_frame_string_column(df, "transcript_id");
  }

  std::vector<std::string> row_ids;
  if (data_frame_has_column(df, "mol_id")) {
    row_ids = data_frame_string_column(df, "mol_id");
  } else if (!table.transcript_ids.empty()) {
    row_ids = table.transcript_ids;
  } else {
    row_ids.resize(table.gene_key.size());
    for (std::size_t i = 0; i < row_ids.size(); ++i) {
      row_ids[i] = std::to_string(i + 1);
    }
  }

  table.finalize();
  return std::make_pair(std::move(table), std::move(row_ids));
}

// Build NCVs with the original reference neighborhood convention used in parity checks.
List build_reference_ncv_matrix(
    const DataFrame& df,
    int k,
    SEXP query_ids_sexp,
    int row_sample_n,
    unsigned int seed) {
  if (k <= 0) {
    stop("k must be positive");
  }

  auto parsed = transcript_table_from_df(df);
  celladmix::TranscriptTable table = std::move(parsed.first);
  std::vector<std::string> row_ids = std::move(parsed.second);

  const auto by_cell = table.transcripts_by_cell();
  std::vector<bool> keep_mask(table.size(), false);
  for (std::size_t cell = 0; cell < by_cell.size(); ++cell) {
    if (static_cast<int>(by_cell[cell].size()) <= k) {
      continue;
    }
    for (const int idx : by_cell[cell]) {
      keep_mask[static_cast<std::size_t>(idx)] = true;
    }
  }

  celladmix::TranscriptTable filtered = celladmix::subset_transcripts(
      table,
      std::nullopt,
      nullptr,
      &keep_mask);
  std::vector<std::string> filtered_row_ids;
  filtered_row_ids.reserve(filtered.size());
  for (std::size_t i = 0; i < row_ids.size(); ++i) {
    if (keep_mask[i]) {
      filtered_row_ids.push_back(row_ids[i]);
    }
  }

  std::vector<std::string> genes = filtered.gene_key;
  std::sort(genes.begin(), genes.end());
  genes.erase(std::unique(genes.begin(), genes.end()), genes.end());
  std::unordered_map<std::string, int> gene_to_col;
  gene_to_col.reserve(genes.size());
  for (std::size_t j = 0; j < genes.size(); ++j) {
    gene_to_col.emplace(genes[j], static_cast<int>(j));
  }

  std::vector<int> query_indices;
  std::vector<std::string> query_row_ids;
  if (Rf_isNull(query_ids_sexp)) {
    query_indices.resize(filtered.size());
    query_row_ids = filtered_row_ids;
    for (std::size_t i = 0; i < filtered.size(); ++i) {
      query_indices[i] = static_cast<int>(i);
    }
  } else {
    const std::vector<std::string> requested = as<std::vector<std::string>>(query_ids_sexp);
    std::unordered_map<std::string, int> row_lookup;
    row_lookup.reserve(filtered_row_ids.size());
    for (std::size_t i = 0; i < filtered_row_ids.size(); ++i) {
      row_lookup.emplace(filtered_row_ids[i], static_cast<int>(i));
    }
    for (const auto& id : requested) {
      const auto it = row_lookup.find(id);
      if (it == row_lookup.end()) {
        continue;
      }
      query_indices.push_back(it->second);
      query_row_ids.push_back(id);
    }
  }

  celladmix::DenseMatrix x(static_cast<int>(query_indices.size()), static_cast<int>(genes.size()), 0.0);
  const auto filtered_by_cell = filtered.transcripts_by_cell();
  std::vector<std::vector<int>> query_positions_by_cell(filtered.num_cells());
  for (std::size_t qi = 0; qi < query_indices.size(); ++qi) {
    const int query = query_indices[qi];
    const int cell = filtered.cell_index[static_cast<std::size_t>(query)];
    query_positions_by_cell[static_cast<std::size_t>(cell)].push_back(static_cast<int>(qi));
  }

  for (std::size_t cell = 0; cell < filtered_by_cell.size(); ++cell) {
    const auto& members = filtered_by_cell[cell];
    const auto& query_positions = query_positions_by_cell[cell];
    if (members.empty() || query_positions.empty()) {
      continue;
    }
    const celladmix::SpatialKnnIndex cell_index(filtered, members);
    for (const int qi : query_positions) {
      const int query = query_indices[static_cast<std::size_t>(qi)];
      const auto hits = cell_index.query(query, k + 1, true);
      for (const auto& hit : hits) {
        const auto it = gene_to_col.find(filtered.gene_key[static_cast<std::size_t>(hit.index)]);
        if (it != gene_to_col.end()) {
          x(qi, it->second) += 1.0;
        }
      }
    }
  }

  if (row_sample_n > 0 && row_sample_n < x.rows()) {
    std::mt19937 rng(seed);
    std::vector<int> keep_rows(static_cast<std::size_t>(x.rows()));
    std::iota(keep_rows.begin(), keep_rows.end(), 0);
    std::shuffle(keep_rows.begin(), keep_rows.end(), rng);
    keep_rows.resize(static_cast<std::size_t>(row_sample_n));

    celladmix::DenseMatrix sampled(static_cast<int>(keep_rows.size()), x.cols(), 0.0);
    std::vector<std::string> sampled_row_ids;
    sampled_row_ids.reserve(keep_rows.size());
    for (std::size_t out_i = 0; out_i < keep_rows.size(); ++out_i) {
      const int in_i = keep_rows[out_i];
      for (int j = 0; j < x.cols(); ++j) {
        sampled(static_cast<int>(out_i), j) = x(in_i, j);
      }
      sampled_row_ids.push_back(query_row_ids[static_cast<std::size_t>(in_i)]);
    }
    x = std::move(sampled);
    query_row_ids = std::move(sampled_row_ids);
  }

  std::vector<int> keep_row_idx;
  keep_row_idx.reserve(static_cast<std::size_t>(x.rows()));
  for (int i = 0; i < x.rows(); ++i) {
    double sum = 0.0;
    for (int j = 0; j < x.cols(); ++j) {
      sum += x(i, j);
    }
    if (sum > 0.0) {
      keep_row_idx.push_back(i);
    }
  }
  std::vector<int> keep_col_idx;
  keep_col_idx.reserve(static_cast<std::size_t>(x.cols()));
  for (int j = 0; j < x.cols(); ++j) {
    double sum = 0.0;
    for (int i = 0; i < x.rows(); ++i) {
      sum += x(i, j);
    }
    if (sum > 0.0) {
      keep_col_idx.push_back(j);
    }
  }

  celladmix::DenseMatrix out(static_cast<int>(keep_row_idx.size()), static_cast<int>(keep_col_idx.size()), 0.0);
  CharacterVector out_row_names(static_cast<R_xlen_t>(keep_row_idx.size()));
  CharacterVector out_col_names(static_cast<R_xlen_t>(keep_col_idx.size()));
  for (std::size_t jj = 0; jj < keep_col_idx.size(); ++jj) {
    out_col_names[static_cast<R_xlen_t>(jj)] = genes[static_cast<std::size_t>(keep_col_idx[jj])];
  }
  for (std::size_t ii = 0; ii < keep_row_idx.size(); ++ii) {
    const int src_i = keep_row_idx[ii];
    out_row_names[static_cast<R_xlen_t>(ii)] = query_row_ids[static_cast<std::size_t>(src_i)];
    for (std::size_t jj = 0; jj < keep_col_idx.size(); ++jj) {
      out(static_cast<int>(ii), static_cast<int>(jj)) = x(src_i, keep_col_idx[jj]);
    }
  }

  NumericMatrix out_matrix = matrix_to_r(out);
  out_matrix.attr("dimnames") = List::create(out_row_names, out_col_names);
  return List::create(
      _["x"] = out_matrix,
      _["row_ids"] = out_row_names,
      _["genes"] = out_col_names);
}

// Build NCVs with the current pipeline neighborhood convention used in the fit path.
List build_pipeline_ncv_matrix(
    const DataFrame& df,
    int k,
    SEXP query_ids_sexp) {
  if (k <= 0) {
    stop("k must be positive");
  }

  auto parsed = transcript_table_from_df(df);
  celladmix::TranscriptTable table = std::move(parsed.first);
  std::vector<std::string> row_ids = std::move(parsed.second);

  std::vector<int> query_indices;
  std::vector<std::string> query_row_ids;
  if (Rf_isNull(query_ids_sexp)) {
    query_indices.resize(table.size());
    query_row_ids = row_ids;
    for (std::size_t i = 0; i < table.size(); ++i) {
      query_indices[i] = static_cast<int>(i);
    }
  } else {
    const std::vector<std::string> requested = as<std::vector<std::string>>(query_ids_sexp);
    std::unordered_map<std::string, int> row_lookup;
    row_lookup.reserve(row_ids.size());
    for (std::size_t i = 0; i < row_ids.size(); ++i) {
      row_lookup.emplace(row_ids[i], static_cast<int>(i));
    }
    for (const auto& id : requested) {
      const auto it = row_lookup.find(id);
      if (it == row_lookup.end()) {
        continue;
      }
      query_indices.push_back(it->second);
      query_row_ids.push_back(id);
    }
  }

  celladmix::NcvOptions options;
  options.k = k + 1;
  options.include_self = true;
  options.within_cell = true;
  options.query_indices = query_indices;

  const celladmix::DenseMatrix x = celladmix::build_ncv_matrix(table, options);
  NumericMatrix out_matrix = matrix_to_r(x);
  out_matrix.attr("dimnames") = List::create(
      wrap(query_row_ids),
      wrap(table.genes));
  return List::create(
      _["x"] = out_matrix,
      _["row_ids"] = wrap(query_row_ids),
      _["genes"] = wrap(table.genes));
}

CharacterVector load_run_training_transcript_ids(const std::string& path_or_dir) {
  const auto training_obs = celladmix::load_run_training_obs_ids(path_or_dir);
  if (training_obs.empty()) {
    return CharacterVector(0);
  }

  std::unordered_map<std::int64_t, std::string> obs_to_tx;
  obs_to_tx.reserve(training_obs.size());
  std::unordered_set<std::int64_t> keep_obs(training_obs.begin(), training_obs.end());

  const auto manifest = celladmix::read_run_manifest(path_or_dir);
  const auto table = read_parquet_table(manifest.paths.molecules_parquet);
  const int obs_idx = table->schema()->GetFieldIndex("obs_id");
  const int tx_idx = table->schema()->GetFieldIndex("transcript_id");
  if (obs_idx < 0) {
    stop("molecules.parquet is missing obs_id");
  }
  const auto obs_chunks = table->column(obs_idx);
  std::shared_ptr<arrow::ChunkedArray> tx_chunks;
  if (tx_idx >= 0) {
    tx_chunks = table->column(tx_idx);
  }

  for (int c = 0; c < obs_chunks->num_chunks(); ++c) {
    const auto obs_chunk = obs_chunks->chunk(c);
    std::shared_ptr<arrow::Array> tx_chunk;
    if (tx_idx >= 0) {
      tx_chunk = tx_chunks->chunk(c);
    }
    for (int64_t i = 0; i < obs_chunk->length(); ++i) {
      std::int64_t obs = 0;
      if (obs_chunk->type_id() == arrow::Type::INT64) {
        obs = std::static_pointer_cast<arrow::Int64Array>(obs_chunk)->Value(i);
      } else {
        obs = static_cast<std::int64_t>(
            std::static_pointer_cast<arrow::Int32Array>(obs_chunk)->Value(i));
      }
      if (keep_obs.find(obs) == keep_obs.end()) {
        continue;
      }
      if (tx_idx >= 0) {
        obs_to_tx.emplace(
            obs,
            std::static_pointer_cast<arrow::StringArray>(tx_chunk)->GetString(i));
      } else {
        obs_to_tx.emplace(obs, std::to_string(obs));
      }
    }
    if (obs_to_tx.size() == keep_obs.size()) {
      break;
    }
  }

  CharacterVector out(static_cast<R_xlen_t>(training_obs.size()));
  for (std::size_t i = 0; i < training_obs.size(); ++i) {
    const auto it = obs_to_tx.find(training_obs[i]);
    if (it == obs_to_tx.end()) {
      out[static_cast<R_xlen_t>(i)] = NA_STRING;
    } else {
      out[static_cast<R_xlen_t>(i)] = it->second;
    }
  }
  return out;
}

List export_training_ncv_matrix(const std::string& path_or_dir) {
  const auto run_data = celladmix::load_run_training_cell_data(path_or_dir);
  if (run_data.training_query_indices.empty()) {
    stop("Run does not contain sampled training rows");
  }

  celladmix::NcvOptions options;
  options.k = run_data.manifest.pipeline_options.ncv_k + 1;
  options.include_self = true;
  options.within_cell = true;
  options.query_indices = run_data.training_query_indices;

  const celladmix::DenseMatrix x = celladmix::build_ncv_matrix(run_data.transcripts, options);
  NumericMatrix out_matrix = matrix_to_r(x);

  CharacterVector row_ids(static_cast<R_xlen_t>(run_data.training_obs_ids.size()));
  CharacterVector cell_ids(static_cast<R_xlen_t>(run_data.training_query_indices.size()));
  CharacterVector cell_types(static_cast<R_xlen_t>(run_data.training_query_indices.size()));
  for (std::size_t i = 0; i < run_data.training_query_indices.size(); ++i) {
    row_ids[static_cast<R_xlen_t>(i)] = std::to_string(run_data.training_obs_ids[i]);
    const int query = run_data.training_query_indices[i];
    const int cell = run_data.transcripts.cell_index[static_cast<std::size_t>(query)];
    cell_ids[static_cast<R_xlen_t>(i)] = run_data.cells.cell_ids[static_cast<std::size_t>(cell)];
    if (!run_data.cells.cell_types.empty()) {
      cell_types[static_cast<R_xlen_t>(i)] = run_data.cells.cell_types[static_cast<std::size_t>(cell)];
    } else {
      cell_types[static_cast<R_xlen_t>(i)] = "";
    }
  }
  CharacterVector genes(static_cast<R_xlen_t>(run_data.transcripts.genes.size()));
  for (std::size_t j = 0; j < run_data.transcripts.genes.size(); ++j) {
    genes[static_cast<R_xlen_t>(j)] = run_data.transcripts.genes[j];
  }
  out_matrix.attr("dimnames") = List::create(row_ids, genes);

  return List::create(
      _["x"] = out_matrix,
      _["row_ids"] = row_ids,
      _["genes"] = genes,
      _["cell_id"] = cell_ids,
      _["cell_type"] = cell_types,
      _["ncv_k"] = run_data.manifest.pipeline_options.ncv_k,
      _["run_path"] = run_data.manifest.paths.root_dir);
}

void write_cell_clusters_parquet(
    const std::string& path,
    const celladmix::CellClusteringResult& result,
    int64_t row_group_size = 65536) {
  const auto table_out = arrow::Table::Make(
      arrow::schema({
          arrow::field("cell_id", arrow::utf8()),
          arrow::field("cell_type", arrow::utf8()),
          arrow::field("cluster", arrow::int32()),
          arrow::field("transcript_count", arrow::int32()),
          arrow::field("detected_genes", arrow::int32()),
          arrow::field("x", arrow::float64()),
          arrow::field("y", arrow::float64()),
          arrow::field("z", arrow::float64()),
          arrow::field("analysis_crop", arrow::utf8()),
      }),
      {
          build_string_array(result.cells.cell_ids),
          build_string_array(
              result.cells.cell_types.empty()
                  ? std::vector<std::string>(result.cells.size(), "")
                  : result.cells.cell_types),
          build_int32_array(result.clusters),
          build_int32_array(result.transcript_counts),
          build_int32_array(result.detected_genes),
          build_double_array(result.cells.centroid_x),
          build_double_array(result.cells.centroid_y),
          build_double_array(result.cells.centroid_z),
          build_string_array(result.crop_ids),
      });
  auto sink = arrow_unwrap(arrow::io::FileOutputStream::Open(path), "Open cluster parquet output");
  auto writer = arrow_unwrap(
      parquet::arrow::FileWriter::Open(*table_out->schema(), arrow::default_memory_pool(), sink),
      "Open cluster parquet writer");
  arrow_check(
      writer->WriteTable(
          *table_out,
          std::max<int64_t>(1, std::min<int64_t>(table_out->num_rows(), row_group_size))),
      "Write cluster parquet");
  arrow_check(writer->Close(), "Close cluster parquet writer");
  arrow_check(sink->Close(), "Close cluster parquet sink");
}

void write_cell_labels_parquet(
    const std::string& path,
    const std::vector<std::string>& cell_id,
    const std::vector<std::string>& label,
    int64_t row_group_size = 65536) {
  if (cell_id.size() != label.size()) {
    throw std::runtime_error("cell_id and label vectors must have the same length");
  }
  const auto output_path = std::filesystem::path(path);
  if (output_path.has_parent_path()) {
    std::filesystem::create_directories(output_path.parent_path());
  }
  const auto table_out = arrow::Table::Make(
      arrow::schema({
          arrow::field("cell_id", arrow::utf8()),
          arrow::field("label", arrow::utf8()),
      }),
      {
          build_string_array(cell_id),
          build_string_array(label),
      });
  auto sink = arrow_unwrap(arrow::io::FileOutputStream::Open(path), "Open cell label parquet output");
  auto writer = arrow_unwrap(
      parquet::arrow::FileWriter::Open(*table_out->schema(), arrow::default_memory_pool(), sink),
      "Open cell label parquet writer");
  arrow_check(
      writer->WriteTable(
          *table_out,
          std::max<int64_t>(1, std::min<int64_t>(table_out->num_rows(), row_group_size))),
      "Write cell label parquet");
  arrow_check(writer->Close(), "Close cell label parquet writer");
  arrow_check(sink->Close(), "Close cell label parquet sink");
}

void write_cell_embedding_parquet(
    const std::string& path,
    const celladmix::CellClusteringResult& result,
    int64_t row_group_size = 65536) {
  if (!result.has_umap) {
    throw std::runtime_error("cannot write cell embedding parquet when UMAP was not computed");
  }
  std::vector<double> umap_1(result.cells.size(), 0.0);
  std::vector<double> umap_2(result.cells.size(), 0.0);
  for (std::size_t i = 0; i < result.cells.size(); ++i) {
    if (result.umap.cols() > 0) {
      umap_1[i] = result.umap(static_cast<int>(i), 0);
    }
    if (result.umap.cols() > 1) {
      umap_2[i] = result.umap(static_cast<int>(i), 1);
    }
  }

  const auto table_out = arrow::Table::Make(
      arrow::schema({
          arrow::field("cell_id", arrow::utf8()),
          arrow::field("cluster", arrow::int32()),
          arrow::field("umap_1", arrow::float64()),
          arrow::field("umap_2", arrow::float64()),
          arrow::field("analysis_crop", arrow::utf8()),
      }),
      {
          build_string_array(result.cells.cell_ids),
          build_int32_array(result.clusters),
          build_double_array(umap_1),
          build_double_array(umap_2),
          build_string_array(result.crop_ids),
      });
  auto sink = arrow_unwrap(arrow::io::FileOutputStream::Open(path), "Open embedding parquet output");
  auto writer = arrow_unwrap(
      parquet::arrow::FileWriter::Open(*table_out->schema(), arrow::default_memory_pool(), sink),
      "Open embedding parquet writer");
  arrow_check(
      writer->WriteTable(
          *table_out,
          std::max<int64_t>(1, std::min<int64_t>(table_out->num_rows(), row_group_size))),
      "Write embedding parquet");
  arrow_check(writer->Close(), "Close embedding parquet writer");
  arrow_check(sink->Close(), "Close embedding parquet sink");
}

celladmix::DenseMatrix load_factor_matrix_dense(const celladmix::RunManifest& manifest) {
  const auto table = read_parquet_table(manifest.paths.factors_parquet);
  const auto factor_id = extract_int_column(table, "factor_id");
  const auto gene_idx = extract_int_column(table, "gene_idx");
  const auto loading = extract_double_column(table, "loading");

  celladmix::DenseMatrix out(
      static_cast<int>(manifest.n_factors),
      static_cast<int>(manifest.genes.size()),
      0.0);
  for (std::size_t i = 0; i < loading.size(); ++i) {
    const int factor = factor_id[i] - 1;
    const int gene = gene_idx[i];
    if (factor >= 0 &&
        factor < static_cast<int>(manifest.n_factors) &&
        gene >= 0 &&
        gene < static_cast<int>(manifest.genes.size())) {
      out(factor, gene) = loading[i];
    }
  }
  return out;
}

NumericMatrix load_factor_matrix_from_run(const celladmix::RunManifest& manifest) {
  return matrix_to_r(load_factor_matrix_dense(manifest));
}

void write_bridge_scores_parquet(
    const std::string& path,
    const celladmix::TranscriptTable& table,
    const std::vector<celladmix::BridgeEvidence>& scores) {
  std::vector<std::string> cell_a;
  std::vector<std::string> cell_b;
  std::vector<int> factor;
  std::vector<double> fraction_a;
  std::vector<double> fraction_b;
  std::vector<double> crossing_fraction;
  std::vector<double> score;
  cell_a.reserve(scores.size());
  cell_b.reserve(scores.size());
  factor.reserve(scores.size());
  fraction_a.reserve(scores.size());
  fraction_b.reserve(scores.size());
  crossing_fraction.reserve(scores.size());
  score.reserve(scores.size());
  for (const auto& row : scores) {
    cell_a.push_back(table.cells[static_cast<std::size_t>(row.pair.cell_a)]);
    cell_b.push_back(table.cells[static_cast<std::size_t>(row.pair.cell_b)]);
    factor.push_back(row.factor + 1);
    fraction_a.push_back(row.fraction_a);
    fraction_b.push_back(row.fraction_b);
    crossing_fraction.push_back(row.crossing_fraction);
    score.push_back(row.score);
  }

  const auto table_out = arrow::Table::Make(
      arrow::schema({
          arrow::field("cell_a", arrow::utf8()),
          arrow::field("cell_b", arrow::utf8()),
          arrow::field("factor", arrow::int32()),
          arrow::field("fraction_a", arrow::float64()),
          arrow::field("fraction_b", arrow::float64()),
          arrow::field("crossing_fraction", arrow::float64()),
          arrow::field("score", arrow::float64()),
      }),
      {
          build_string_array(cell_a),
          build_string_array(cell_b),
          build_int32_array(factor),
          build_double_array(fraction_a),
          build_double_array(fraction_b),
          build_double_array(crossing_fraction),
          build_double_array(score),
      });
  auto sink = arrow_unwrap(arrow::io::FileOutputStream::Open(path), "Open bridge score output");
  auto writer = arrow_unwrap(
      parquet::arrow::FileWriter::Open(*table_out->schema(), arrow::default_memory_pool(), sink),
      "Open bridge score writer");
  arrow_check(writer->WriteTable(*table_out, std::max<int64_t>(1, std::min<int64_t>(table_out->num_rows(), 65536))),
              "Write bridge score table");
  arrow_check(writer->Close(), "Close bridge score writer");
  arrow_check(sink->Close(), "Close bridge score sink");
}

void write_bridge_pair_scores_parquet(
    const std::string& path,
    const celladmix::TranscriptTable& table,
    const celladmix::BridgeTestResult& result) {
  std::vector<std::string> target_cell;
  std::vector<std::string> source_cell;
  std::vector<std::string> target_cell_type;
  std::vector<std::string> source_cell_type;
  std::vector<int> factor;
  std::vector<int> factor_count;
  std::vector<int> total_crossing_count;
  std::vector<int> used_in_summary;
  std::vector<double> target_fraction;
  std::vector<double> source_fraction;
  std::vector<double> crossing_fraction;
  std::vector<double> score;
  const auto& scores = result.pair_scores;
  target_cell.reserve(scores.size());
  source_cell.reserve(scores.size());
  target_cell_type.reserve(scores.size());
  source_cell_type.reserve(scores.size());
  factor.reserve(scores.size());
  factor_count.reserve(scores.size());
  total_crossing_count.reserve(scores.size());
  used_in_summary.reserve(scores.size());
  target_fraction.reserve(scores.size());
  source_fraction.reserve(scores.size());
  crossing_fraction.reserve(scores.size());
  score.reserve(scores.size());

  for (const auto& row : scores) {
    target_cell.push_back(table.cells[static_cast<std::size_t>(row.target_cell)]);
    source_cell.push_back(table.cells[static_cast<std::size_t>(row.source_cell)]);
    target_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.target_type)]);
    source_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.source_type)]);
    factor.push_back(row.factor + 1);
    factor_count.push_back(row.factor_count);
    total_crossing_count.push_back(row.total_crossing_count);
    used_in_summary.push_back(row.used_in_summary ? 1 : 0);
    target_fraction.push_back(row.target_fraction);
    source_fraction.push_back(row.source_fraction);
    crossing_fraction.push_back(row.crossing_fraction);
    score.push_back(row.score);
  }

  const auto table_out = arrow::Table::Make(
      arrow::schema({
          arrow::field("target_cell", arrow::utf8()),
          arrow::field("source_cell", arrow::utf8()),
          arrow::field("target_cell_type", arrow::utf8()),
          arrow::field("source_cell_type", arrow::utf8()),
          arrow::field("factor", arrow::int32()),
          arrow::field("factor_count", arrow::int32()),
          arrow::field("target_fraction", arrow::float64()),
          arrow::field("source_fraction", arrow::float64()),
          arrow::field("crossing_fraction", arrow::float64()),
          arrow::field("score", arrow::float64()),
          arrow::field("total_crossing_count", arrow::int32()),
          arrow::field("used_in_summary", arrow::int32()),
      }),
      {
          build_string_array(target_cell),
          build_string_array(source_cell),
          build_string_array(target_cell_type),
          build_string_array(source_cell_type),
          build_int32_array(factor),
          build_int32_array(factor_count),
          build_double_array(target_fraction),
          build_double_array(source_fraction),
          build_double_array(crossing_fraction),
          build_double_array(score),
          build_int32_array(total_crossing_count),
          build_int32_array(used_in_summary),
      });
  auto sink = arrow_unwrap(arrow::io::FileOutputStream::Open(path), "Open bridge pair score output");
  auto writer = arrow_unwrap(
      parquet::arrow::FileWriter::Open(*table_out->schema(), arrow::default_memory_pool(), sink),
      "Open bridge pair score writer");
  arrow_check(writer->WriteTable(*table_out, std::max<int64_t>(1, std::min<int64_t>(table_out->num_rows(), 65536))),
              "Write bridge pair score table");
  arrow_check(writer->Close(), "Close bridge pair score writer");
  arrow_check(sink->Close(), "Close bridge pair score sink");
}

void write_bridge_summary_parquet(
    const std::string& path,
    const celladmix::BridgeTestResult& result) {
  std::vector<std::string> target_cell_type;
  std::vector<std::string> source_cell_type;
  std::vector<int> factor;
  std::vector<int> n_pairs;
  std::vector<double> mean_score;
  std::vector<double> mean_null_score;
  std::vector<double> q75_score;
  std::vector<double> p_value;
  std::vector<double> neg_log10_p;
  const auto& summaries = result.summaries;
  target_cell_type.reserve(summaries.size());
  source_cell_type.reserve(summaries.size());
  factor.reserve(summaries.size());
  n_pairs.reserve(summaries.size());
  mean_score.reserve(summaries.size());
  mean_null_score.reserve(summaries.size());
  q75_score.reserve(summaries.size());
  p_value.reserve(summaries.size());
  neg_log10_p.reserve(summaries.size());

  for (const auto& row : summaries) {
    target_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.target_type)]);
    source_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.source_type)]);
    factor.push_back(row.factor + 1);
    n_pairs.push_back(row.n_pairs);
    mean_score.push_back(row.mean_score);
    mean_null_score.push_back(row.mean_null_score);
    q75_score.push_back(row.q75_score);
    p_value.push_back(row.p_value);
    neg_log10_p.push_back(row.neg_log10_p);
  }

  const auto table_out = arrow::Table::Make(
      arrow::schema({
          arrow::field("target_cell_type", arrow::utf8()),
          arrow::field("source_cell_type", arrow::utf8()),
          arrow::field("factor", arrow::int32()),
          arrow::field("n_pairs", arrow::int32()),
          arrow::field("mean_score", arrow::float64()),
          arrow::field("mean_null_score", arrow::float64()),
          arrow::field("q75_score", arrow::float64()),
          arrow::field("p_value", arrow::float64()),
          arrow::field("neg_log10_p", arrow::float64()),
      }),
      {
          build_string_array(target_cell_type),
          build_string_array(source_cell_type),
          build_int32_array(factor),
          build_int32_array(n_pairs),
          build_double_array(mean_score),
          build_double_array(mean_null_score),
          build_double_array(q75_score),
          build_double_array(p_value),
          build_double_array(neg_log10_p),
      });
  auto sink = arrow_unwrap(arrow::io::FileOutputStream::Open(path), "Open bridge summary output");
  auto writer = arrow_unwrap(
      parquet::arrow::FileWriter::Open(*table_out->schema(), arrow::default_memory_pool(), sink),
      "Open bridge summary writer");
  arrow_check(writer->WriteTable(*table_out, std::max<int64_t>(1, std::min<int64_t>(table_out->num_rows(), 65536))),
              "Write bridge summary table");
  arrow_check(writer->Close(), "Close bridge summary writer");
  arrow_check(sink->Close(), "Close bridge summary sink");
}

void write_membrane_pair_scores_parquet(
    const std::string& path,
    const celladmix::TranscriptTable& table,
    const celladmix::MembraneTestResult& result) {
  std::vector<std::string> target_cell;
  std::vector<std::string> source_cell;
  std::vector<std::string> target_cell_type;
  std::vector<std::string> source_cell_type;
  std::vector<int> factor;
  std::vector<int> factor_count;
  std::vector<int> scored_molecules;
  std::vector<int> used_in_summary;
  std::vector<double> mean_score;
  std::vector<double> fraction_positive;
  std::vector<double> mean_directional_weight;
  const auto& scores = result.pair_scores;
  target_cell.reserve(scores.size());
  source_cell.reserve(scores.size());
  target_cell_type.reserve(scores.size());
  source_cell_type.reserve(scores.size());
  factor.reserve(scores.size());
  factor_count.reserve(scores.size());
  scored_molecules.reserve(scores.size());
  used_in_summary.reserve(scores.size());
  mean_score.reserve(scores.size());
  fraction_positive.reserve(scores.size());
  mean_directional_weight.reserve(scores.size());

  for (const auto& row : scores) {
    target_cell.push_back(table.cells[static_cast<std::size_t>(row.target_cell)]);
    source_cell.push_back(table.cells[static_cast<std::size_t>(row.source_cell)]);
    target_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.target_type)]);
    source_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.source_type)]);
    factor.push_back(row.factor + 1);
    factor_count.push_back(row.factor_count);
    scored_molecules.push_back(row.scored_molecules);
    used_in_summary.push_back(row.used_in_summary ? 1 : 0);
    mean_score.push_back(row.mean_score);
    fraction_positive.push_back(row.fraction_positive);
    mean_directional_weight.push_back(row.mean_directional_weight);
  }

  const auto table_out = arrow::Table::Make(
      arrow::schema({
          arrow::field("target_cell", arrow::utf8()),
          arrow::field("source_cell", arrow::utf8()),
          arrow::field("target_cell_type", arrow::utf8()),
          arrow::field("source_cell_type", arrow::utf8()),
          arrow::field("factor", arrow::int32()),
          arrow::field("factor_count", arrow::int32()),
          arrow::field("scored_molecules", arrow::int32()),
          arrow::field("mean_score", arrow::float64()),
          arrow::field("fraction_positive", arrow::float64()),
          arrow::field("mean_directional_weight", arrow::float64()),
          arrow::field("used_in_summary", arrow::int32()),
      }),
      {
          build_string_array(target_cell),
          build_string_array(source_cell),
          build_string_array(target_cell_type),
          build_string_array(source_cell_type),
          build_int32_array(factor),
          build_int32_array(factor_count),
          build_int32_array(scored_molecules),
          build_double_array(mean_score),
          build_double_array(fraction_positive),
          build_double_array(mean_directional_weight),
          build_int32_array(used_in_summary),
      });
  auto sink = arrow_unwrap(arrow::io::FileOutputStream::Open(path), "Open membrane pair score output");
  auto writer = arrow_unwrap(
      parquet::arrow::FileWriter::Open(*table_out->schema(), arrow::default_memory_pool(), sink),
      "Open membrane pair score writer");
  arrow_check(writer->WriteTable(*table_out, std::max<int64_t>(1, std::min<int64_t>(table_out->num_rows(), 65536))),
              "Write membrane pair score table");
  arrow_check(writer->Close(), "Close membrane pair score writer");
  arrow_check(sink->Close(), "Close membrane pair score sink");
}

void write_membrane_summary_parquet(
    const std::string& path,
    const celladmix::MembraneTestResult& result) {
  std::vector<std::string> target_cell_type;
  std::vector<std::string> source_cell_type;
  std::vector<int> factor;
  std::vector<int> n_pairs;
  std::vector<double> mean_score;
  std::vector<double> q75_score;
  std::vector<double> fraction_positive_pairs;
  std::vector<double> p_value;
  std::vector<double> neg_log10_p;
  const auto& summaries = result.summaries;
  target_cell_type.reserve(summaries.size());
  source_cell_type.reserve(summaries.size());
  factor.reserve(summaries.size());
  n_pairs.reserve(summaries.size());
  mean_score.reserve(summaries.size());
  q75_score.reserve(summaries.size());
  fraction_positive_pairs.reserve(summaries.size());
  p_value.reserve(summaries.size());
  neg_log10_p.reserve(summaries.size());

  for (const auto& row : summaries) {
    target_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.target_type)]);
    source_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.source_type)]);
    factor.push_back(row.factor + 1);
    n_pairs.push_back(row.n_pairs);
    mean_score.push_back(row.mean_score);
    q75_score.push_back(row.q75_score);
    fraction_positive_pairs.push_back(row.fraction_positive_pairs);
    p_value.push_back(row.p_value);
    neg_log10_p.push_back(row.neg_log10_p);
  }

  const auto table_out = arrow::Table::Make(
      arrow::schema({
          arrow::field("target_cell_type", arrow::utf8()),
          arrow::field("source_cell_type", arrow::utf8()),
          arrow::field("factor", arrow::int32()),
          arrow::field("n_pairs", arrow::int32()),
          arrow::field("mean_score", arrow::float64()),
          arrow::field("q75_score", arrow::float64()),
          arrow::field("fraction_positive_pairs", arrow::float64()),
          arrow::field("p_value", arrow::float64()),
          arrow::field("neg_log10_p", arrow::float64()),
      }),
      {
          build_string_array(target_cell_type),
          build_string_array(source_cell_type),
          build_int32_array(factor),
          build_int32_array(n_pairs),
          build_double_array(mean_score),
          build_double_array(q75_score),
          build_double_array(fraction_positive_pairs),
          build_double_array(p_value),
          build_double_array(neg_log10_p),
      });
  auto sink = arrow_unwrap(arrow::io::FileOutputStream::Open(path), "Open membrane summary output");
  auto writer = arrow_unwrap(
      parquet::arrow::FileWriter::Open(*table_out->schema(), arrow::default_memory_pool(), sink),
      "Open membrane summary writer");
  arrow_check(writer->WriteTable(*table_out, std::max<int64_t>(1, std::min<int64_t>(table_out->num_rows(), 65536))),
              "Write membrane summary table");
  arrow_check(writer->Close(), "Close membrane summary writer");
  arrow_check(sink->Close(), "Close membrane summary sink");
}

void write_coherence_cell_scores_parquet(
    const std::string& path,
    const celladmix::TranscriptTable& table,
    const celladmix::CoherenceTestResult& result) {
  std::vector<std::string> target_cell;
  std::vector<std::string> target_cell_type;
  std::vector<std::string> source_cell_type;
  std::vector<int> factor;
  std::vector<int> factor_count;
  std::vector<int> active_count;
  std::vector<int> used_in_summary;
  std::vector<double> mean_raw_score;
  std::vector<double> mean_coherence_score;
  std::vector<double> mean_score;
  std::vector<double> mean_null_score;
  std::vector<double> score_delta;
  std::vector<double> q75_score;
  std::vector<double> active_fraction;
  std::vector<double> mean_edge_weight;
  std::vector<int> largest_patch_count;
  std::vector<double> largest_patch_fraction;
  std::vector<double> patch_score;
  std::vector<double> mean_null_patch_score;
  std::vector<double> patch_score_delta;
  std::vector<double> source_log_enrichment;
  std::vector<double> source_probability;
  const auto& scores = result.cell_scores;
  target_cell.reserve(scores.size());
  target_cell_type.reserve(scores.size());
  source_cell_type.reserve(scores.size());
  factor.reserve(scores.size());
  factor_count.reserve(scores.size());
  active_count.reserve(scores.size());
  used_in_summary.reserve(scores.size());
  mean_raw_score.reserve(scores.size());
  mean_coherence_score.reserve(scores.size());
  mean_score.reserve(scores.size());
  mean_null_score.reserve(scores.size());
  score_delta.reserve(scores.size());
  q75_score.reserve(scores.size());
  active_fraction.reserve(scores.size());
  mean_edge_weight.reserve(scores.size());
  largest_patch_count.reserve(scores.size());
  largest_patch_fraction.reserve(scores.size());
  patch_score.reserve(scores.size());
  mean_null_patch_score.reserve(scores.size());
  patch_score_delta.reserve(scores.size());
  source_log_enrichment.reserve(scores.size());
  source_probability.reserve(scores.size());

  for (const auto& row : scores) {
    target_cell.push_back(table.cells[static_cast<std::size_t>(row.target_cell)]);
    target_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.target_type)]);
    source_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.source_type)]);
    factor.push_back(row.factor + 1);
    factor_count.push_back(row.factor_count);
    active_count.push_back(row.active_count);
    used_in_summary.push_back(row.used_in_summary ? 1 : 0);
    mean_raw_score.push_back(row.mean_raw_score);
    mean_coherence_score.push_back(row.mean_coherence_score);
    mean_score.push_back(row.mean_score);
    mean_null_score.push_back(row.mean_null_score);
    score_delta.push_back(row.score_delta);
    q75_score.push_back(row.q75_score);
    active_fraction.push_back(row.active_fraction);
    mean_edge_weight.push_back(row.mean_edge_weight);
    largest_patch_count.push_back(row.largest_patch_count);
    largest_patch_fraction.push_back(row.largest_patch_fraction);
    patch_score.push_back(row.patch_score);
    mean_null_patch_score.push_back(row.mean_null_patch_score);
    patch_score_delta.push_back(row.patch_score_delta);
    source_log_enrichment.push_back(row.source_log_enrichment);
    source_probability.push_back(row.source_probability);
  }

  const auto table_out = arrow::Table::Make(
      arrow::schema({
          arrow::field("target_cell", arrow::utf8()),
          arrow::field("target_cell_type", arrow::utf8()),
          arrow::field("source_cell_type", arrow::utf8()),
          arrow::field("factor", arrow::int32()),
          arrow::field("factor_count", arrow::int32()),
          arrow::field("active_count", arrow::int32()),
          arrow::field("mean_raw_score", arrow::float64()),
          arrow::field("mean_coherence_score", arrow::float64()),
          arrow::field("mean_score", arrow::float64()),
          arrow::field("mean_null_score", arrow::float64()),
          arrow::field("score_delta", arrow::float64()),
          arrow::field("q75_score", arrow::float64()),
          arrow::field("active_fraction", arrow::float64()),
          arrow::field("mean_edge_weight", arrow::float64()),
          arrow::field("largest_patch_count", arrow::int32()),
          arrow::field("largest_patch_fraction", arrow::float64()),
          arrow::field("patch_score", arrow::float64()),
          arrow::field("mean_null_patch_score", arrow::float64()),
          arrow::field("patch_score_delta", arrow::float64()),
          arrow::field("source_log_enrichment", arrow::float64()),
          arrow::field("source_probability", arrow::float64()),
          arrow::field("used_in_summary", arrow::int32()),
      }),
      {
          build_string_array(target_cell),
          build_string_array(target_cell_type),
          build_string_array(source_cell_type),
          build_int32_array(factor),
          build_int32_array(factor_count),
          build_int32_array(active_count),
          build_double_array(mean_raw_score),
          build_double_array(mean_coherence_score),
          build_double_array(mean_score),
          build_double_array(mean_null_score),
          build_double_array(score_delta),
          build_double_array(q75_score),
          build_double_array(active_fraction),
          build_double_array(mean_edge_weight),
          build_int32_array(largest_patch_count),
          build_double_array(largest_patch_fraction),
          build_double_array(patch_score),
          build_double_array(mean_null_patch_score),
          build_double_array(patch_score_delta),
          build_double_array(source_log_enrichment),
          build_double_array(source_probability),
          build_int32_array(used_in_summary),
      });
  auto sink = arrow_unwrap(arrow::io::FileOutputStream::Open(path), "Open coherence cell score output");
  auto writer = arrow_unwrap(
      parquet::arrow::FileWriter::Open(*table_out->schema(), arrow::default_memory_pool(), sink),
      "Open coherence cell score writer");
  arrow_check(writer->WriteTable(*table_out, std::max<int64_t>(1, std::min<int64_t>(table_out->num_rows(), 65536))),
              "Write coherence cell score table");
  arrow_check(writer->Close(), "Close coherence cell score writer");
  arrow_check(sink->Close(), "Close coherence cell score sink");
}

void write_coherence_summary_parquet(
    const std::string& path,
    const celladmix::CoherenceTestResult& result) {
  std::vector<std::string> target_cell_type;
  std::vector<std::string> source_cell_type;
  std::vector<int> factor;
  std::vector<int> n_cells;
  std::vector<int> n_molecules;
  std::vector<int> n_active_molecules;
  std::vector<double> active_fraction;
  std::vector<double> mean_score;
  std::vector<double> mean_null_score;
  std::vector<double> mean_delta_score;
  std::vector<double> q75_score;
  std::vector<double> mean_patch_score;
  std::vector<double> mean_null_patch_score;
  std::vector<double> mean_delta_patch_score;
  std::vector<double> mean_largest_patch_fraction;
  std::vector<double> patch_p_value;
  std::vector<double> patch_neg_log10_p;
  std::vector<double> p_value;
  std::vector<double> neg_log10_p;
  std::vector<double> source_log_enrichment;
  std::vector<double> source_probability;
  const auto& summaries = result.summaries;
  target_cell_type.reserve(summaries.size());
  source_cell_type.reserve(summaries.size());
  factor.reserve(summaries.size());
  n_cells.reserve(summaries.size());
  n_molecules.reserve(summaries.size());
  n_active_molecules.reserve(summaries.size());
  active_fraction.reserve(summaries.size());
  mean_score.reserve(summaries.size());
  mean_null_score.reserve(summaries.size());
  mean_delta_score.reserve(summaries.size());
  q75_score.reserve(summaries.size());
  mean_patch_score.reserve(summaries.size());
  mean_null_patch_score.reserve(summaries.size());
  mean_delta_patch_score.reserve(summaries.size());
  mean_largest_patch_fraction.reserve(summaries.size());
  patch_p_value.reserve(summaries.size());
  patch_neg_log10_p.reserve(summaries.size());
  p_value.reserve(summaries.size());
  neg_log10_p.reserve(summaries.size());
  source_log_enrichment.reserve(summaries.size());
  source_probability.reserve(summaries.size());

  for (const auto& row : summaries) {
    target_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.target_type)]);
    source_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.source_type)]);
    factor.push_back(row.factor + 1);
    n_cells.push_back(row.n_cells);
    n_molecules.push_back(row.n_molecules);
    n_active_molecules.push_back(row.n_active_molecules);
    active_fraction.push_back(row.active_fraction);
    mean_score.push_back(row.mean_score);
    mean_null_score.push_back(row.mean_null_score);
    mean_delta_score.push_back(row.mean_delta_score);
    q75_score.push_back(row.q75_score);
    mean_patch_score.push_back(row.mean_patch_score);
    mean_null_patch_score.push_back(row.mean_null_patch_score);
    mean_delta_patch_score.push_back(row.mean_delta_patch_score);
    mean_largest_patch_fraction.push_back(row.mean_largest_patch_fraction);
    patch_p_value.push_back(row.patch_p_value);
    patch_neg_log10_p.push_back(row.patch_neg_log10_p);
    p_value.push_back(row.p_value);
    neg_log10_p.push_back(row.neg_log10_p);
    source_log_enrichment.push_back(row.source_log_enrichment);
    source_probability.push_back(row.source_probability);
  }

  const auto table_out = arrow::Table::Make(
      arrow::schema({
          arrow::field("target_cell_type", arrow::utf8()),
          arrow::field("source_cell_type", arrow::utf8()),
          arrow::field("factor", arrow::int32()),
          arrow::field("n_cells", arrow::int32()),
          arrow::field("n_molecules", arrow::int32()),
          arrow::field("n_active_molecules", arrow::int32()),
          arrow::field("active_fraction", arrow::float64()),
          arrow::field("mean_score", arrow::float64()),
          arrow::field("mean_null_score", arrow::float64()),
          arrow::field("mean_delta_score", arrow::float64()),
          arrow::field("q75_score", arrow::float64()),
          arrow::field("mean_patch_score", arrow::float64()),
          arrow::field("mean_null_patch_score", arrow::float64()),
          arrow::field("mean_delta_patch_score", arrow::float64()),
          arrow::field("mean_largest_patch_fraction", arrow::float64()),
          arrow::field("patch_p_value", arrow::float64()),
          arrow::field("patch_neg_log10_p", arrow::float64()),
          arrow::field("p_value", arrow::float64()),
          arrow::field("neg_log10_p", arrow::float64()),
          arrow::field("source_log_enrichment", arrow::float64()),
          arrow::field("source_probability", arrow::float64()),
      }),
      {
          build_string_array(target_cell_type),
          build_string_array(source_cell_type),
          build_int32_array(factor),
          build_int32_array(n_cells),
          build_int32_array(n_molecules),
          build_int32_array(n_active_molecules),
          build_double_array(active_fraction),
          build_double_array(mean_score),
          build_double_array(mean_null_score),
          build_double_array(mean_delta_score),
          build_double_array(q75_score),
          build_double_array(mean_patch_score),
          build_double_array(mean_null_patch_score),
          build_double_array(mean_delta_patch_score),
          build_double_array(mean_largest_patch_fraction),
          build_double_array(patch_p_value),
          build_double_array(patch_neg_log10_p),
          build_double_array(p_value),
          build_double_array(neg_log10_p),
          build_double_array(source_log_enrichment),
          build_double_array(source_probability),
      });
  auto sink = arrow_unwrap(arrow::io::FileOutputStream::Open(path), "Open coherence summary output");
  auto writer = arrow_unwrap(
      parquet::arrow::FileWriter::Open(*table_out->schema(), arrow::default_memory_pool(), sink),
      "Open coherence summary writer");
  arrow_check(writer->WriteTable(*table_out, std::max<int64_t>(1, std::min<int64_t>(table_out->num_rows(), 65536))),
              "Write coherence summary table");
  arrow_check(writer->Close(), "Close coherence summary writer");
  arrow_check(sink->Close(), "Close coherence summary sink");
}

std::optional<celladmix::CropBox> optional_bbox_sexp(SEXP value) {
  if (Rf_isNull(value)) {
    return std::nullopt;
  }
  const NumericVector bbox(value);
  if (bbox.size() != 4 && bbox.size() != 6) {
    stop("bbox must contain 4 or 6 numeric values");
  }
  celladmix::CropBox crop;
  crop.crop_id = "region";
  crop.xmin = bbox[0];
  crop.xmax = bbox[1];
  crop.ymin = bbox[2];
  crop.ymax = bbox[3];
  crop.has_z = bbox.size() == 6;
  if (crop.has_z) {
    crop.zmin = bbox[4];
    crop.zmax = bbox[5];
  }
  return crop;
}

std::vector<celladmix::CropBox> crop_boxes_from_df(SEXP crops_sexp) {
  if (Rf_isNull(crops_sexp)) {
    return {};
  }
  const DataFrame df = as<DataFrame>(crops_sexp);
  if (!df.containsElementNamed("crop_id") ||
      !df.containsElementNamed("xmin") ||
      !df.containsElementNamed("xmax") ||
      !df.containsElementNamed("ymin") ||
      !df.containsElementNamed("ymax")) {
    stop("crops must contain crop_id/xmin/xmax/ymin/ymax columns");
  }
  const auto crop_id = as<std::vector<std::string>>(df["crop_id"]);
  const auto xmin = as<std::vector<double>>(df["xmin"]);
  const auto xmax = as<std::vector<double>>(df["xmax"]);
  const auto ymin = as<std::vector<double>>(df["ymin"]);
  const auto ymax = as<std::vector<double>>(df["ymax"]);
  std::vector<double> zmin(crop_id.size(), 0.0);
  std::vector<double> zmax(crop_id.size(), 0.0);
  std::vector<bool> has_z(crop_id.size(), false);
  if (df.containsElementNamed("zmin") && df.containsElementNamed("zmax")) {
    zmin = as<std::vector<double>>(df["zmin"]);
    zmax = as<std::vector<double>>(df["zmax"]);
    for (std::size_t i = 0; i < crop_id.size(); ++i) {
      has_z[i] = std::isfinite(zmin[i]) || std::isfinite(zmax[i]);
    }
  }
  std::vector<celladmix::CropBox> out;
  out.reserve(crop_id.size());
  for (std::size_t i = 0; i < crop_id.size(); ++i) {
    out.push_back(celladmix::CropBox{
        crop_id[i], xmin[i], xmax[i], ymin[i], ymax[i], zmin[i], zmax[i], has_z[i]});
  }
  return out;
}

std::string compose_training_stratum(
    const std::string& crop_id,
    const std::string& cell_type) {
  if (!crop_id.empty() && !cell_type.empty()) {
    return crop_id + "|" + cell_type;
  }
  if (!crop_id.empty()) {
    return crop_id;
  }
  if (!cell_type.empty()) {
    return cell_type;
  }
  return {};
}

std::vector<std::string> training_cell_strata_from_cells(
    const celladmix::TranscriptTable& table,
    const celladmix::CellTable& cells,
    const std::vector<std::string>* cell_crop_ids = nullptr) {
  const bool has_crop = cell_crop_ids != nullptr && cell_crop_ids->size() == cells.cell_ids.size();
  const bool has_cell_type = !cells.cell_types.empty();
  if (!has_crop && !has_cell_type) {
    return {};
  }

  std::unordered_map<std::string, std::size_t> lookup;
  lookup.reserve(cells.cell_ids.size());
  for (std::size_t i = 0; i < cells.cell_ids.size(); ++i) {
    lookup.emplace(cells.cell_ids[i], i);
  }

  std::vector<std::string> out(table.cells.size(), "");
  for (std::size_t i = 0; i < table.cells.size(); ++i) {
    const auto it = lookup.find(table.cells[i]);
    if (it == lookup.end()) {
      continue;
    }
    const std::size_t idx = it->second;
    out[i] = compose_training_stratum(
        has_crop ? (*cell_crop_ids)[idx] : "",
        has_cell_type ? cells.cell_types[idx] : "");
  }
  return out;
}

std::vector<std::string> training_cell_strata_from_cell_metadata(
    const celladmix::CellTable& cells,
    const std::vector<std::string>& cell_crop_ids = {}) {
  const bool has_crop = cell_crop_ids.size() == cells.cell_ids.size();
  const bool has_cell_type = !cells.cell_types.empty();
  if (!has_crop && !has_cell_type) {
    return {};
  }
  std::vector<std::string> out(cells.cell_ids.size(), "");
  for (std::size_t i = 0; i < cells.cell_ids.size(); ++i) {
    out[i] = compose_training_stratum(
        has_crop ? cell_crop_ids[i] : "",
        has_cell_type ? cells.cell_types[i] : "");
  }
  return out;
}

std::vector<std::string> subset_cell_crop_ids(
    const celladmix::CellTable& subset_cells,
    const celladmix::CellTable& source_cells,
    const std::vector<std::string>& source_crop_ids) {
  if (source_crop_ids.size() != source_cells.cell_ids.size()) {
    return {};
  }

  std::unordered_map<std::string, std::string> lookup;
  lookup.reserve(source_cells.cell_ids.size());
  for (std::size_t i = 0; i < source_cells.cell_ids.size(); ++i) {
    lookup.emplace(source_cells.cell_ids[i], source_crop_ids[i]);
  }

  std::vector<std::string> out(subset_cells.cell_ids.size(), "");
  for (std::size_t i = 0; i < subset_cells.cell_ids.size(); ++i) {
    const auto it = lookup.find(subset_cells.cell_ids[i]);
    if (it != lookup.end()) {
      out[i] = it->second;
    }
  }
  return out;
}

List manifest_to_r(const celladmix::XeniumManifest& manifest) {
  return List::create(
      _["run_name"] = manifest.run_name,
      _["pixel_size"] = manifest.pixel_size,
      _["z_step_size"] = manifest.z_step_size,
      _["morphology_filepath"] = manifest.morphology_filepath,
      _["morphology_focus_filepath"] = manifest.morphology_focus_filepath,
      _["transcripts_csv_path"] = manifest.transcripts_csv_path,
      _["transcripts_parquet_path"] = manifest.transcripts_parquet_path);
}

DataFrame run_source_files_to_r(const std::vector<celladmix::RunSourceFile>& files) {
  CharacterVector role(static_cast<R_xlen_t>(files.size()));
  CharacterVector path(static_cast<R_xlen_t>(files.size()));
  LogicalVector exists(static_cast<R_xlen_t>(files.size()));
  for (std::size_t i = 0; i < files.size(); ++i) {
    role[static_cast<R_xlen_t>(i)] = files[i].role;
    path[static_cast<R_xlen_t>(i)] = files[i].path;
    exists[static_cast<R_xlen_t>(i)] = files[i].exists;
  }
  List out = List::create(
      _["role"] = role,
      _["path"] = path,
      _["exists"] = exists);
  out.attr("class") = "data.frame";
  out.attr("row.names") = IntegerVector::create(NA_INTEGER, -static_cast<int>(files.size()));
  return static_cast<DataFrame>(out);
}

List run_manifest_to_r(const celladmix::RunManifest& manifest) {
  List source = List::create(
      _["type"] = manifest.source.type,
      _["path"] = manifest.source.path,
      _["used_parquet"] = manifest.source.used_parquet,
      _["files"] = run_source_files_to_r(manifest.source.files));
  List paths = List::create(
      _["root_dir"] = manifest.paths.root_dir,
      _["run_json"] = manifest.paths.run_json,
      _["factors_parquet"] = manifest.paths.factors_parquet,
      _["cells_parquet"] = manifest.paths.cells_parquet,
      _["molecules_parquet"] = manifest.paths.molecules_parquet,
      _["training_rows_parquet"] = manifest.paths.training_rows_parquet,
      _["correction_summary_parquet"] =
          celladmix::correction_summary_parquet_path(manifest.paths.root_dir),
      _["report_dir"] = std::filesystem::path(manifest.paths.root_dir).append("report").string(),
      _["training_molecule_umap_parquet"] =
          celladmix::training_molecule_umap_parquet_path(manifest.paths.root_dir),
      _["scores_dir"] = manifest.paths.scores_dir,
      _["corrected_dir"] = manifest.paths.corrected_dir);
  List pipeline_options = List::create(
      _["ncv_k"] = manifest.pipeline_options.ncv_k,
      _["rank"] = manifest.pipeline_options.rank,
      _["graph_k"] = manifest.pipeline_options.graph_k,
      _["same_label_ratio"] = manifest.pipeline_options.same_label_ratio,
      _["nmf_iterations"] = manifest.pipeline_options.nmf_iterations,
      _["nmf_init"] = manifest.pipeline_options.nmf_init,
      _["nmf_variant"] = manifest.pipeline_options.nmf_variant,
      _["molecule_scoring"] = manifest.pipeline_options.molecule_scoring,
      _["nmf_n_runs"] = manifest.pipeline_options.nmf_n_runs,
      _["nmf_train_max_rows"] = manifest.pipeline_options.nmf_train_max_rows,
      _["nmf_min_molecules"] = manifest.pipeline_options.nmf_min_molecules,
      _["num_threads"] = manifest.pipeline_options.num_threads,
      _["return_ncv"] = manifest.pipeline_options.return_ncv,
      _["seed"] = manifest.pipeline_options.seed,
      _["training_scope_cell_types"] = wrap(manifest.pipeline_options.training_scope_cell_types));
  List storage_options = List::create(
      _["tile_size"] = manifest.storage_options.tile_size,
      _["parquet_row_group_size"] = manifest.storage_options.parquet_row_group_size);
  List nmf_diagnostics = List::create(
      _["final_objective"] = manifest.nmf_diagnostics.final_objective,
      _["selected_seed"] = static_cast<double>(manifest.nmf_diagnostics.selected_seed),
      _["selected_run"] = static_cast<double>(manifest.nmf_diagnostics.selected_run + 1),
      _["candidate_final_objectives"] = wrap(manifest.nmf_diagnostics.candidate_final_objectives),
      _["candidate_best_match_correlations"] =
          wrap(manifest.nmf_diagnostics.candidate_best_match_correlations),
      _["selected_factor_stability"] = wrap(manifest.nmf_diagnostics.selected_factor_stability),
      _["candidate_final_objective_mean"] =
          manifest.nmf_diagnostics.candidate_final_objective_mean,
      _["candidate_final_objective_sd"] =
          manifest.nmf_diagnostics.candidate_final_objective_sd,
      _["candidate_best_match_correlation_mean"] =
          manifest.nmf_diagnostics.candidate_best_match_correlation_mean,
      _["stability_metric"] = manifest.nmf_diagnostics.stability_metric,
      _["stability_comparison_runs"] = manifest.nmf_diagnostics.stability_comparison_runs,
      _["stable_factor_count"] = manifest.nmf_diagnostics.stable_factor_count,
      _["stability_threshold"] = manifest.nmf_diagnostics.stability_threshold);

  return List::create(
      _["path"] = manifest.paths.root_dir,
      _["run_type"] = manifest.run_type,
      _["package_version"] = manifest.package_version,
      _["annotation_hash"] = manifest.annotation_hash,
      _["source"] = source,
      _["paths"] = paths,
      _["pipeline_options"] = pipeline_options,
      _["storage_options"] = storage_options,
      _["nmf_diagnostics"] = nmf_diagnostics,
      _["nmf_final_objective"] = manifest.nmf_diagnostics.final_objective,
      _["nmf_selected_seed"] = static_cast<double>(manifest.nmf_diagnostics.selected_seed),
      _["nmf_selected_run"] = static_cast<double>(manifest.nmf_diagnostics.selected_run + 1),
      _["nmf_candidate_final_objectives"] =
          wrap(manifest.nmf_diagnostics.candidate_final_objectives),
      _["nmf_candidate_best_match_correlations"] =
          wrap(manifest.nmf_diagnostics.candidate_best_match_correlations),
      _["nmf_factor_stability"] = wrap(manifest.nmf_diagnostics.selected_factor_stability),
      _["nmf_candidate_final_objective_mean"] =
          manifest.nmf_diagnostics.candidate_final_objective_mean,
      _["nmf_candidate_final_objective_sd"] =
          manifest.nmf_diagnostics.candidate_final_objective_sd,
      _["nmf_candidate_best_match_correlation_mean"] =
          manifest.nmf_diagnostics.candidate_best_match_correlation_mean,
      _["nmf_stability_metric"] = manifest.nmf_diagnostics.stability_metric,
      _["nmf_stability_comparison_runs"] = manifest.nmf_diagnostics.stability_comparison_runs,
      _["nmf_stable_factor_count"] = manifest.nmf_diagnostics.stable_factor_count,
      _["nmf_stability_threshold"] = manifest.nmf_diagnostics.stability_threshold,
      _["analysis_crop"] = manifest.analysis_crop.has_value() ? wrap(*manifest.analysis_crop) : R_NilValue,
      _["parent_run"] = manifest.parent_run.has_value() ? wrap(*manifest.parent_run) : R_NilValue,
      _["genes"] = manifest.genes,
      _["nmf_gene_weights"] = manifest.nmf_gene_weights,
      _["nmf_transform_target_row_sum"] = manifest.nmf_transform_target_row_sum,
      _["crop_ids"] = manifest.crop_ids,
      _["n_transcripts"] = static_cast<double>(manifest.n_transcripts),
      _["n_cells"] = static_cast<double>(manifest.n_cells),
      _["n_factors"] = static_cast<double>(manifest.n_factors),
      _["n_training_rows"] = static_cast<double>(manifest.n_training_rows),
      _["has_z"] = manifest.has_z,
      _["has_qv"] = manifest.has_qv,
      _["has_transcript_id"] = manifest.has_transcript_id,
      _["has_cell_type"] = manifest.has_cell_type,
      _["has_sample_id"] = manifest.has_sample_id,
      _["has_fov_id"] = manifest.has_fov_id,
      _["has_nucleus_id"] = manifest.has_nucleus_id,
      _["has_overlaps_nucleus"] = manifest.has_overlaps_nucleus,
      _["has_nucleus_distance"] = manifest.has_nucleus_distance);
}

List input_store_manifest_to_r(const celladmix::InputStoreManifest& manifest) {
  const auto paths = celladmix::make_input_store_paths(manifest.root_dir);
  List path_list = List::create(
      _["root_dir"] = paths.root_dir,
      _["manifest_json"] = paths.manifest_json,
      _["genes_parquet"] = paths.genes_parquet,
      _["cells_parquet"] = paths.cells_parquet,
      _["counts_parquet"] = paths.counts_parquet,
      _["molecules_parquet"] = paths.molecules_parquet,
      _["cell_offsets_parquet"] = paths.cell_offsets_parquet);
  List capabilities = List::create(
      _["has_molecules"] = manifest.has_molecules,
      _["has_molecule_rows"] = manifest.has_molecule_rows,
      _["has_cell_gene_counts"] = manifest.has_cell_gene_counts,
      _["has_cell_offsets"] = manifest.has_cell_offsets,
      _["has_spatial_tiles"] = manifest.has_spatial_tiles,
      _["has_labels"] = manifest.has_labels,
      _["has_factor_scores"] = manifest.has_factor_scores);
  return List::create(
      _["path"] = manifest.root_dir,
      _["format_version"] = manifest.format_version,
      _["source_type"] = manifest.source_type,
      _["store_mode"] = manifest.store_mode,
      _["source_path"] = manifest.source_path,
      _["source_fingerprint"] = manifest.source_fingerprint,
      _["filter_signature"] = manifest.filter_signature,
      _["paths"] = path_list,
      _["capabilities"] = capabilities,
      _["n_molecules"] = static_cast<double>(manifest.n_molecules),
      _["n_cells"] = static_cast<double>(manifest.n_cells),
      _["n_genes"] = static_cast<double>(manifest.n_genes),
      _["has_z"] = manifest.has_z,
      _["has_qv"] = manifest.has_qv,
      _["has_transcript_id"] = manifest.has_transcript_id,
      _["has_cell_type"] = manifest.has_cell_type,
      _["has_sample_id"] = manifest.has_sample_id,
      _["has_fov_id"] = manifest.has_fov_id,
      _["has_overlaps_nucleus"] = manifest.has_overlaps_nucleus,
      _["has_nucleus_distance"] = manifest.has_nucleus_distance);
}

std::string resolve_xenium_manifest_path_local(const std::string& path_or_dir) {
  const auto path = std::filesystem::path(path_or_dir);
  if (std::filesystem::is_directory(path)) {
    return (path / "experiment.xenium").string();
  }
  return path.string();
}

std::string resolve_relative_to_parent(
    const std::filesystem::path& parent,
    const std::string& child) {
  if (child.empty()) {
    return {};
  }
  const auto child_path = std::filesystem::path(child);
  if (child_path.is_absolute()) {
    return child_path.string();
  }
  return (parent / child_path).string();
}

std::vector<celladmix::RunSourceFile> build_xenium_source_files(
    const std::string& source_path,
    const celladmix::XeniumBundleData& bundle) {
  const auto manifest_path = std::filesystem::absolute(
      std::filesystem::path(resolve_xenium_manifest_path_local(source_path)));
  const auto dataset_dir = manifest_path.parent_path();

  auto make_file = [](const std::string& role, const std::string& path) {
    celladmix::RunSourceFile file;
    file.role = role;
    file.path = path;
    file.exists = !path.empty() && std::filesystem::exists(std::filesystem::path(path));
    return file;
  };

  std::vector<celladmix::RunSourceFile> files;
  files.push_back(make_file("experiment_manifest", manifest_path.string()));

  const auto transcripts_parquet = resolve_relative_to_parent(
      dataset_dir,
      bundle.manifest.transcripts_parquet_path);
  const auto transcripts_csv = resolve_relative_to_parent(
      dataset_dir,
      bundle.manifest.transcripts_csv_path);
  const auto morphology = resolve_relative_to_parent(
      dataset_dir,
      bundle.manifest.morphology_filepath);
  const auto morphology_focus = resolve_relative_to_parent(
      dataset_dir,
      bundle.manifest.morphology_focus_filepath);
  const auto cells_parquet = (dataset_dir / "cells.parquet").string();
  const auto cells_csv = (dataset_dir / "cells.csv.gz").string();

  files.push_back(make_file("transcripts_parquet", transcripts_parquet));
  files.push_back(make_file("transcripts_csv", transcripts_csv));
  if (!bundle.manifest.morphology_filepath.empty()) {
    files.push_back(make_file("morphology_image", morphology));
  }
  if (!bundle.manifest.morphology_focus_filepath.empty()) {
    files.push_back(make_file("morphology_focus_image", morphology_focus));
  }
  if (std::filesystem::exists(std::filesystem::path(cells_parquet))) {
    files.push_back(make_file("cells_parquet", cells_parquet));
  }
  if (std::filesystem::exists(std::filesystem::path(cells_csv))) {
    files.push_back(make_file("cells_csv", cells_csv));
  }
  return files;
}

std::vector<celladmix::RunSourceFile> build_input_store_source_files(
    const std::string& store_dir,
    const celladmix::InputStoreManifest& manifest) {
  const auto paths = celladmix::make_input_store_paths(store_dir);
  auto make_file = [](const std::string& role, const std::string& path) {
    celladmix::RunSourceFile file;
    file.role = role;
    file.path = path;
    file.exists = !path.empty() && std::filesystem::exists(std::filesystem::path(path));
    return file;
  };

  std::vector<celladmix::RunSourceFile> files;
  files.push_back(make_file("input_store_manifest", paths.manifest_json));
  files.push_back(make_file("genes", paths.genes_parquet));
  files.push_back(make_file("cells", paths.cells_parquet));
  files.push_back(make_file("cell_gene_counts", paths.counts_parquet));
  if (manifest.has_molecules) {
    files.push_back(make_file("molecules", paths.molecules_parquet));
  }
  if (manifest.has_cell_offsets) {
    files.push_back(make_file("cell_offsets", paths.cell_offsets_parquet));
  }
  return files;
}

std::vector<std::string> cell_crop_ids_from_transcripts(
    const celladmix::TranscriptTable& table,
    const std::vector<std::string>& transcript_crop_ids) {
  std::vector<std::string> out(table.num_cells(), "");
  if (transcript_crop_ids.size() != table.size()) {
    return out;
  }
  for (std::size_t row = 0; row < table.size(); ++row) {
    const std::size_t cell = static_cast<std::size_t>(table.cell_index[row]);
    if (out[cell].empty()) {
      out[cell] = transcript_crop_ids[row];
    }
  }
  return out;
}

std::vector<celladmix::RunSourceFile> build_tabular_source_files(
    const celladmix::TabularSourceSpec& source,
    bool used_parquet) {
  auto make_file = [](const std::string& role, const std::string& path) {
    celladmix::RunSourceFile file;
    file.role = role;
    file.path = path;
    file.exists = !path.empty() && std::filesystem::exists(std::filesystem::path(path));
    return file;
  };

  std::vector<celladmix::RunSourceFile> files;
  files.push_back(make_file(
      used_parquet ? "molecules_parquet" : "molecules_csv",
      std::filesystem::absolute(std::filesystem::path(source.molecules_path)).string()));
  if (!source.segmentation_mask_path.empty()) {
    files.push_back(make_file(
        "segmentation_mask",
        std::filesystem::absolute(std::filesystem::path(source.segmentation_mask_path)).string()));
  }
  return files;
}

List probe_xenium_source_to_r(const std::string& source_path) {
  const std::string manifest_path = resolve_xenium_manifest_path_local(source_path);
  const auto manifest = celladmix::read_xenium_manifest(manifest_path);
  const auto manifest_abs = std::filesystem::absolute(std::filesystem::path(manifest_path));
  const auto dataset_dir = manifest_abs.parent_path();

  const std::string transcripts_parquet =
      resolve_relative_to_parent(dataset_dir, manifest.transcripts_parquet_path);
  const std::string transcripts_csv =
      resolve_relative_to_parent(dataset_dir, manifest.transcripts_csv_path);
  const std::string cells_parquet = (dataset_dir / "cells.parquet").string();
  const std::string cells_csv_gz = (dataset_dir / "cells.csv.gz").string();
  const std::string cells_csv = (dataset_dir / "cells.csv").string();
  const std::string morphology =
      resolve_relative_to_parent(dataset_dir, manifest.morphology_filepath);
  const std::string morphology_focus =
      resolve_relative_to_parent(dataset_dir, manifest.morphology_focus_filepath);

  const bool has_transcripts_parquet = !transcripts_parquet.empty() &&
      std::filesystem::exists(std::filesystem::path(transcripts_parquet));
  const bool has_transcripts_csv = !transcripts_csv.empty() &&
      std::filesystem::exists(std::filesystem::path(transcripts_csv));
  const bool has_cells_parquet = std::filesystem::exists(std::filesystem::path(cells_parquet));
  const bool has_cells_csv_gz = std::filesystem::exists(std::filesystem::path(cells_csv_gz));
  const bool has_cells_csv = std::filesystem::exists(std::filesystem::path(cells_csv));

  const std::string preferred_transcripts =
      has_transcripts_parquet ? transcripts_parquet : transcripts_csv;
  const std::string preferred_cells = has_cells_parquet
      ? cells_parquet
      : (has_cells_csv_gz ? cells_csv_gz : cells_csv);

  std::vector<std::string> transcript_columns;
  if (has_transcripts_parquet) {
    transcript_columns = read_parquet_column_names(transcripts_parquet);
  } else if (has_transcripts_csv) {
    transcript_columns = read_csv_header_columns(transcripts_csv);
  }

  std::vector<std::string> cell_columns;
  if (has_cells_parquet) {
    cell_columns = read_parquet_column_names(cells_parquet);
  } else if (has_cells_csv_gz) {
    cell_columns = read_csv_header_columns(cells_csv_gz);
  } else if (has_cells_csv) {
    cell_columns = read_csv_header_columns(cells_csv);
  }

  return List::create(
      _["manifest"] = manifest_to_r(manifest),
      _["bundle_dir"] = dataset_dir.string(),
      _["manifest_path"] = manifest_abs.string(),
      _["transcripts"] = List::create(
          _["preferred_path"] = preferred_transcripts,
          _["has_parquet"] = has_transcripts_parquet,
          _["has_csv"] = has_transcripts_csv,
          _["columns"] = wrap(transcript_columns)),
      _["cells"] = List::create(
          _["preferred_path"] = preferred_cells,
          _["has_parquet"] = has_cells_parquet,
          _["has_csv"] = has_cells_csv_gz || has_cells_csv,
          _["columns"] = wrap(cell_columns)),
      _["morphology"] = file_probe_to_r("morphology_image", morphology),
      _["morphology_focus"] = file_probe_to_r("morphology_focus_image", morphology_focus));
}

List probe_tabular_source_to_r(
    const std::string& molecules_path,
    const std::string& segmentation_mask_path) {
  const auto abs_path = std::filesystem::absolute(std::filesystem::path(molecules_path));
  const auto ext = logical_extension(abs_path.string());
  std::vector<std::string> columns;
  if (ext == ".parquet" || ext == ".pq") {
    columns = read_parquet_column_names(abs_path.string());
  } else {
    columns = read_csv_header_columns(abs_path.string());
  }

  SEXP segmentation_mask = R_NilValue;
  if (!segmentation_mask_path.empty()) {
    segmentation_mask = file_probe_to_r(
        "segmentation_mask",
        std::filesystem::absolute(std::filesystem::path(segmentation_mask_path)).string());
  }

  return List::create(
      _["molecules"] = List::create(
          _["path"] = abs_path.string(),
          _["exists"] = std::filesystem::exists(abs_path),
          _["columns"] = wrap(columns),
          _["is_parquet"] = (ext == ".parquet" || ext == ".pq")),
      _["segmentation_mask"] = segmentation_mask);
}

void annotate_transcripts_from_cells(
    celladmix::TranscriptTable& table,
    const celladmix::CellTable& cells) {
  if (cells.cell_ids.empty()) {
    return;
  }
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
    if (it == cell_lookup.end()) {
      continue;
    }
    const std::size_t cell_idx = it->second;
    if (add_cell_types) table.cell_types[i] = cells.cell_types[cell_idx];
    if (add_sample_ids) table.sample_ids[i] = cells.sample_ids[cell_idx];
    if (add_fov_ids) table.fov_ids[i] = cells.fov_ids[cell_idx];
  }
}

IntegerVector labels_to_r(const std::vector<int>& labels) {
  IntegerVector out(static_cast<R_xlen_t>(labels.size()));
  for (std::size_t i = 0; i < labels.size(); ++i) {
    out[static_cast<R_xlen_t>(i)] = labels[i] + 1;
  }
  return out;
}

IntegerVector indices_to_r1(const std::vector<int>& indices) {
  IntegerVector out(static_cast<R_xlen_t>(indices.size()));
  for (std::size_t i = 0; i < indices.size(); ++i) {
    out[static_cast<R_xlen_t>(i)] = indices[i] + 1;
  }
  return out;
}

NumericMatrix matrix_to_r(const celladmix::DenseMatrix& matrix) {
  NumericMatrix out(matrix.rows(), matrix.cols());
  for (int i = 0; i < matrix.rows(); ++i) {
    for (int j = 0; j < matrix.cols(); ++j) {
      out(i, j) = matrix(i, j);
    }
  }
  return out;
}

celladmix::DenseMatrix dense_matrix_from_r(SEXP matrix_sexp) {
  if (!Rf_isMatrix(matrix_sexp) || !Rf_isNumeric(matrix_sexp)) {
    stop("x must be a numeric matrix");
  }
  NumericMatrix matrix(matrix_sexp);
  celladmix::DenseMatrix out(matrix.nrow(), matrix.ncol(), 0.0);
  for (int i = 0; i < matrix.nrow(); ++i) {
    for (int j = 0; j < matrix.ncol(); ++j) {
      out(i, j) = matrix(i, j);
    }
  }
  return out;
}

celladmix::SparseRowMatrix sparse_row_matrix_from_dgc(SEXP matrix_sexp) {
  if (!Rf_isS4(matrix_sexp)) {
    stop("x must be a dgCMatrix");
  }
  S4 matrix(matrix_sexp);
  CharacterVector classes = matrix.attr("class");
  bool is_dgc = false;
  for (R_xlen_t i = 0; i < classes.size(); ++i) {
    if (as<std::string>(classes[i]) == "dgCMatrix") {
      is_dgc = true;
      break;
    }
  }
  if (!is_dgc) {
    stop("x must be a dgCMatrix");
  }

  IntegerVector dims = matrix.slot("Dim");
  IntegerVector p = matrix.slot("p");
  IntegerVector i = matrix.slot("i");
  NumericVector x = matrix.slot("x");

  const int nrow = dims[0];
  const int ncol = dims[1];
  std::vector<int> row_counts(static_cast<std::size_t>(nrow), 0);
  for (R_xlen_t idx = 0; idx < i.size(); ++idx) {
    row_counts[static_cast<std::size_t>(i[idx])] += 1;
  }

  std::vector<int> indptr(static_cast<std::size_t>(nrow + 1), 0);
  for (int row = 0; row < nrow; ++row) {
    indptr[static_cast<std::size_t>(row + 1)] =
        indptr[static_cast<std::size_t>(row)] + row_counts[static_cast<std::size_t>(row)];
  }

  std::vector<int> indices(static_cast<std::size_t>(x.size()), 0);
  std::vector<double> values(static_cast<std::size_t>(x.size()), 0.0);
  std::vector<int> next = indptr;

  for (int col = 0; col < ncol; ++col) {
    for (int ptr = p[col]; ptr < p[col + 1]; ++ptr) {
      const int row = i[ptr];
      const int offset = next[static_cast<std::size_t>(row)]++;
      indices[static_cast<std::size_t>(offset)] = col;
      values[static_cast<std::size_t>(offset)] = x[ptr];
    }
  }

  return celladmix::SparseRowMatrix(
      nrow,
      ncol,
      std::move(indptr),
      std::move(indices),
      std::move(values));
}

NumericMatrix counts_to_r(
    const celladmix::DenseMatrix& counts_by_cell_gene,
    const std::vector<std::string>& cell_names,
    const std::vector<std::string>& gene_names) {
  NumericMatrix out(static_cast<R_xlen_t>(gene_names.size()), static_cast<R_xlen_t>(cell_names.size()));
  for (int cell = 0; cell < counts_by_cell_gene.rows(); ++cell) {
    for (int gene = 0; gene < counts_by_cell_gene.cols(); ++gene) {
      out(gene, cell) = counts_by_cell_gene(cell, gene);
    }
  }
  out.attr("dimnames") = List::create(wrap(gene_names), wrap(cell_names));
  return out;
}

List sparse_counts_to_r_list(const celladmix::CellCountMatrix& counts) {
  return List::create(
      _["p"] = wrap(counts.indptr),
      _["i"] = wrap(counts.indices),
      _["x"] = wrap(counts.values),
      _["genes"] = wrap(counts.genes),
      _["cells"] = wrap(counts.cells.cell_ids),
      _["cell_x"] = wrap(counts.cells.centroid_x),
      _["cell_y"] = wrap(counts.cells.centroid_y),
      _["cell_type"] = counts.cells.cell_types.empty()
          ? wrap(std::vector<std::string>(counts.cells.cell_ids.size(), ""))
          : wrap(counts.cells.cell_types));
}

std::vector<std::size_t> select_transcript_rows(
    std::size_t total_size,
    const std::optional<std::string>& crop_id = std::nullopt,
    const std::vector<std::string>* crop_ids = nullptr,
    const std::vector<bool>* keep_mask = nullptr,
    int sample_n = -1,
    unsigned int seed = 1U) {
  std::vector<std::size_t> selected;
  selected.reserve(total_size);
  for (std::size_t i = 0; i < total_size; ++i) {
    if (crop_id.has_value()) {
      if (crop_ids == nullptr || crop_ids->size() != total_size) {
        throw std::runtime_error("crop_id selection requires crop ids matching transcript rows");
      }
      if ((*crop_ids)[i] != *crop_id) {
        continue;
      }
    }
    if (keep_mask != nullptr) {
      if (keep_mask->size() != total_size) {
        throw std::runtime_error("keep_mask must match transcript rows");
      }
      if (!(*keep_mask)[i]) {
        continue;
      }
    }
    selected.push_back(i);
  }
  if (sample_n > 0 && static_cast<std::size_t>(sample_n) < selected.size()) {
    std::mt19937 rng(seed);
    std::shuffle(selected.begin(), selected.end(), rng);
    selected.resize(static_cast<std::size_t>(sample_n));
    std::sort(selected.begin(), selected.end());
  }
  return selected;
}

DataFrame transcripts_to_df_selected(
    const celladmix::TranscriptTable& table,
    const std::vector<std::size_t>& indices,
    const std::optional<std::string>& analysis_crop = std::nullopt,
    const std::vector<int>* labels = nullptr,
    const celladmix::DenseMatrix* factor_scores = nullptr,
    const std::vector<std::string>* analysis_crop_values = nullptr,
    const std::vector<double>* factor_margin = nullptr) {
  if (labels != nullptr && labels->size() != table.size()) {
    stop("labels must match transcript rows");
  }
  if (factor_scores != nullptr && factor_scores->rows() != static_cast<int>(table.size())) {
    stop("factor_scores must have one row per transcript");
  }
  if (analysis_crop_values != nullptr && analysis_crop_values->size() != table.size()) {
    stop("analysis_crop_values must match transcript rows");
  }
  if (factor_margin != nullptr && factor_margin->size() != table.size()) {
    stop("factor_margin must match transcript rows");
  }

  std::vector<std::string> gene;
  std::vector<std::string> cell;
  std::vector<std::string> celltype;
  std::vector<std::string> transcript_id;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<std::string> sample_id;
  std::vector<std::string> fov_id;
  std::vector<double> qv;
  std::vector<std::string> nucleus_id;
  std::vector<std::string> crop_value;
  std::vector<double> factor_margin_value;
  IntegerVector factor_label(static_cast<R_xlen_t>(indices.size()));

  const bool has_cell_type = !table.cell_types.empty();
  const bool has_tx_id = !table.transcript_ids.empty();
  const bool has_sample = !table.sample_ids.empty();
  const bool has_fov = !table.fov_ids.empty();
  const bool has_qv = !table.qv.empty();
  const bool has_nucleus = !table.orig_nucleus_ids.empty();
  const bool has_labels = labels != nullptr;
  const bool has_scores = factor_scores != nullptr;

  gene.reserve(indices.size());
  cell.reserve(indices.size());
  celltype.reserve(indices.size());
  transcript_id.reserve(indices.size());
  x.reserve(indices.size());
  y.reserve(indices.size());
  z.reserve(indices.size());
  sample_id.reserve(indices.size());
  fov_id.reserve(indices.size());
  qv.reserve(indices.size());
  nucleus_id.reserve(indices.size());
  crop_value.reserve(indices.size());
  factor_margin_value.reserve(indices.size());

  for (std::size_t pos = 0; pos < indices.size(); ++pos) {
    const std::size_t i = indices[pos];
    gene.push_back(table.gene_key[i]);
    cell.push_back(table.orig_cell_id[i]);
    x.push_back(table.x[i]);
    y.push_back(table.y[i]);
    z.push_back(table.z[i]);
    if (has_cell_type) celltype.push_back(table.cell_types[i]);
    if (has_tx_id) transcript_id.push_back(table.transcript_ids[i]);
    if (has_sample) sample_id.push_back(table.sample_ids[i]);
    if (has_fov) fov_id.push_back(table.fov_ids[i]);
    if (has_qv) qv.push_back(table.qv[i]);
    if (has_nucleus) nucleus_id.push_back(table.orig_nucleus_ids[i]);
    if (analysis_crop_values != nullptr) {
      crop_value.push_back((*analysis_crop_values)[i]);
    } else if (analysis_crop.has_value()) {
      crop_value.push_back(*analysis_crop);
    }
    if (has_labels) factor_label[static_cast<R_xlen_t>(pos)] = (*labels)[i] + 1;
    if (factor_margin != nullptr) factor_margin_value.push_back((*factor_margin)[i]);
  }

  List out = List::create(
      _["gene"] = gene,
      _["cell"] = cell,
      _["celltype"] = has_cell_type ? wrap(celltype) : wrap(std::vector<std::string>(indices.size(), "")),
      _["transcript_id"] = has_tx_id ? wrap(transcript_id) : wrap(std::vector<std::string>(indices.size(), "")),
      _["x"] = x,
      _["y"] = y,
      _["z"] = z);
  if (has_sample) out["sample_id"] = sample_id;
  if (has_fov) out["fov_id"] = fov_id;
  if (has_qv) out["qv"] = qv;
  if (has_nucleus) out["orig_nucleus_id"] = nucleus_id;
  if (analysis_crop_values != nullptr || analysis_crop.has_value()) out["analysis_crop"] = crop_value;
  if (has_labels) out["factor_label"] = factor_label;
  if (factor_margin != nullptr) out["factor_margin"] = factor_margin_value;
  if (has_scores) {
    for (int factor = 0; factor < factor_scores->cols(); ++factor) {
      NumericVector score(static_cast<R_xlen_t>(indices.size()));
      for (std::size_t pos = 0; pos < indices.size(); ++pos) {
        score[static_cast<R_xlen_t>(pos)] = (*factor_scores)(static_cast<int>(indices[pos]), factor);
      }
      out[std::string("factor_") + std::to_string(factor + 1) + "_score"] = score;
    }
  }
  return DataFrame(out);
}

DataFrame cells_to_df(
    const celladmix::CellTable& cells,
    const std::vector<std::string>& analysis_crop = {}) {
  const std::vector<std::string> cell_types = cells.cell_types.empty()
      ? std::vector<std::string>(cells.size(), "")
      : cells.cell_types;
  List out = List::create(
      _["cell_id"] = cells.cell_ids,
      _["cell_type"] = cell_types,
      _["x"] = cells.centroid_x,
      _["y"] = cells.centroid_y,
      _["z"] = cells.centroid_z);
  if (!cells.sample_ids.empty()) out["sample_id"] = cells.sample_ids;
  if (!cells.fov_ids.empty()) out["fov_id"] = cells.fov_ids;
  if (!analysis_crop.empty()) out["analysis_crop"] = analysis_crop;
  return DataFrame(out);
}

DataFrame cell_clusters_to_df(const celladmix::CellClusteringResult& result) {
  std::vector<std::string> cell_type = result.cells.cell_types.empty()
      ? std::vector<std::string>(result.cells.size(), "")
      : result.cells.cell_types;
  return DataFrame::create(
      _["cell_id"] = result.cells.cell_ids,
      _["cell_type"] = cell_type,
      _["cluster"] = wrap(result.clusters),
      _["transcript_count"] = wrap(result.transcript_counts),
      _["detected_genes"] = wrap(result.detected_genes),
      _["x"] = wrap(result.cells.centroid_x),
      _["y"] = wrap(result.cells.centroid_y),
      _["z"] = wrap(result.cells.centroid_z),
      _["analysis_crop"] = wrap(result.crop_ids));
}

SEXP cell_embedding_to_df(const celladmix::CellClusteringResult& result) {
  if (!result.has_umap) {
    return R_NilValue;
  }
  NumericVector umap_1(static_cast<R_xlen_t>(result.cells.size()));
  NumericVector umap_2(static_cast<R_xlen_t>(result.cells.size()));
  for (std::size_t i = 0; i < result.cells.size(); ++i) {
    if (result.umap.cols() > 0) {
      umap_1[static_cast<R_xlen_t>(i)] = result.umap(static_cast<int>(i), 0);
    }
    if (result.umap.cols() > 1) {
      umap_2[static_cast<R_xlen_t>(i)] = result.umap(static_cast<int>(i), 1);
    }
  }
  return DataFrame::create(
      _["cell_id"] = result.cells.cell_ids,
      _["cluster"] = wrap(result.clusters),
      _["umap_1"] = umap_1,
      _["umap_2"] = umap_2,
      _["analysis_crop"] = wrap(result.crop_ids));
}

DataFrame bridge_scores_to_df(
    const celladmix::TranscriptTable& table,
    const std::vector<celladmix::BridgeEvidence>& scores) {
  std::vector<std::string> cell_a_id;
  std::vector<std::string> cell_b_id;
  IntegerVector factor(static_cast<R_xlen_t>(scores.size()));
  NumericVector fraction_a(static_cast<R_xlen_t>(scores.size()));
  NumericVector fraction_b(static_cast<R_xlen_t>(scores.size()));
  NumericVector crossing_fraction(static_cast<R_xlen_t>(scores.size()));
  NumericVector score(static_cast<R_xlen_t>(scores.size()));

  cell_a_id.reserve(scores.size());
  cell_b_id.reserve(scores.size());

  for (std::size_t i = 0; i < scores.size(); ++i) {
    const auto& evidence = scores[i];
    cell_a_id.push_back(table.cells[static_cast<std::size_t>(evidence.pair.cell_a)]);
    cell_b_id.push_back(table.cells[static_cast<std::size_t>(evidence.pair.cell_b)]);
    factor[static_cast<R_xlen_t>(i)] = evidence.factor + 1;
    fraction_a[static_cast<R_xlen_t>(i)] = evidence.fraction_a;
    fraction_b[static_cast<R_xlen_t>(i)] = evidence.fraction_b;
    crossing_fraction[static_cast<R_xlen_t>(i)] = evidence.crossing_fraction;
    score[static_cast<R_xlen_t>(i)] = evidence.score;
  }

  return DataFrame::create(
      _["cell_a"] = cell_a_id,
      _["cell_b"] = cell_b_id,
      _["factor"] = factor,
      _["fraction_a"] = fraction_a,
      _["fraction_b"] = fraction_b,
      _["crossing_fraction"] = crossing_fraction,
      _["score"] = score);
}

DataFrame bridge_pair_scores_to_df(
    const celladmix::TranscriptTable& table,
    const celladmix::BridgeTestResult& result) {
  const auto& scores = result.pair_scores;
  std::vector<std::string> target_cell;
  std::vector<std::string> source_cell;
  std::vector<std::string> target_cell_type;
  std::vector<std::string> source_cell_type;
  IntegerVector factor(static_cast<R_xlen_t>(scores.size()));
  IntegerVector factor_count(static_cast<R_xlen_t>(scores.size()));
  NumericVector target_fraction(static_cast<R_xlen_t>(scores.size()));
  NumericVector source_fraction(static_cast<R_xlen_t>(scores.size()));
  NumericVector crossing_fraction(static_cast<R_xlen_t>(scores.size()));
  NumericVector score(static_cast<R_xlen_t>(scores.size()));
  IntegerVector total_crossing_count(static_cast<R_xlen_t>(scores.size()));
  LogicalVector used_in_summary(static_cast<R_xlen_t>(scores.size()));

  target_cell.reserve(scores.size());
  source_cell.reserve(scores.size());
  target_cell_type.reserve(scores.size());
  source_cell_type.reserve(scores.size());

  for (std::size_t i = 0; i < scores.size(); ++i) {
    const auto& row = scores[i];
    target_cell.push_back(table.cells[static_cast<std::size_t>(row.target_cell)]);
    source_cell.push_back(table.cells[static_cast<std::size_t>(row.source_cell)]);
    target_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.target_type)]);
    source_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.source_type)]);
    factor[static_cast<R_xlen_t>(i)] = row.factor + 1;
    factor_count[static_cast<R_xlen_t>(i)] = row.factor_count;
    target_fraction[static_cast<R_xlen_t>(i)] = row.target_fraction;
    source_fraction[static_cast<R_xlen_t>(i)] = row.source_fraction;
    crossing_fraction[static_cast<R_xlen_t>(i)] = row.crossing_fraction;
    score[static_cast<R_xlen_t>(i)] = row.score;
    total_crossing_count[static_cast<R_xlen_t>(i)] = row.total_crossing_count;
    used_in_summary[static_cast<R_xlen_t>(i)] = row.used_in_summary;
  }

  return DataFrame::create(
      _["target_cell"] = target_cell,
      _["source_cell"] = source_cell,
      _["target_cell_type"] = target_cell_type,
      _["source_cell_type"] = source_cell_type,
      _["factor"] = factor,
      _["factor_count"] = factor_count,
      _["target_fraction"] = target_fraction,
      _["source_fraction"] = source_fraction,
      _["crossing_fraction"] = crossing_fraction,
      _["score"] = score,
      _["total_crossing_count"] = total_crossing_count,
      _["used_in_summary"] = used_in_summary);
}

DataFrame bridge_summary_to_df(const celladmix::BridgeTestResult& result) {
  const auto& summaries = result.summaries;
  std::vector<std::string> target_cell_type;
  std::vector<std::string> source_cell_type;
  IntegerVector factor(static_cast<R_xlen_t>(summaries.size()));
  IntegerVector n_pairs(static_cast<R_xlen_t>(summaries.size()));
  NumericVector mean_score(static_cast<R_xlen_t>(summaries.size()));
  NumericVector mean_null_score(static_cast<R_xlen_t>(summaries.size()));
  NumericVector q75_score(static_cast<R_xlen_t>(summaries.size()));
  NumericVector p_value(static_cast<R_xlen_t>(summaries.size()));
  NumericVector neg_log10_p(static_cast<R_xlen_t>(summaries.size()));
  target_cell_type.reserve(summaries.size());
  source_cell_type.reserve(summaries.size());

  for (std::size_t i = 0; i < summaries.size(); ++i) {
    const auto& row = summaries[i];
    target_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.target_type)]);
    source_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.source_type)]);
    factor[static_cast<R_xlen_t>(i)] = row.factor + 1;
    n_pairs[static_cast<R_xlen_t>(i)] = row.n_pairs;
    mean_score[static_cast<R_xlen_t>(i)] = row.mean_score;
    mean_null_score[static_cast<R_xlen_t>(i)] = row.mean_null_score;
    q75_score[static_cast<R_xlen_t>(i)] = row.q75_score;
    p_value[static_cast<R_xlen_t>(i)] = row.p_value;
    neg_log10_p[static_cast<R_xlen_t>(i)] = row.neg_log10_p;
  }

  return DataFrame::create(
      _["target_cell_type"] = target_cell_type,
      _["source_cell_type"] = source_cell_type,
      _["factor"] = factor,
      _["n_pairs"] = n_pairs,
      _["mean_score"] = mean_score,
      _["mean_null_score"] = mean_null_score,
      _["q75_score"] = q75_score,
      _["p_value"] = p_value,
      _["neg_log10_p"] = neg_log10_p);
}

DataFrame membrane_pair_scores_to_df(
    const celladmix::TranscriptTable& table,
    const celladmix::MembraneTestResult& result) {
  const auto& scores = result.pair_scores;
  std::vector<std::string> target_cell;
  std::vector<std::string> source_cell;
  std::vector<std::string> target_cell_type;
  std::vector<std::string> source_cell_type;
  IntegerVector factor(static_cast<R_xlen_t>(scores.size()));
  IntegerVector factor_count(static_cast<R_xlen_t>(scores.size()));
  IntegerVector scored_molecules(static_cast<R_xlen_t>(scores.size()));
  NumericVector mean_score(static_cast<R_xlen_t>(scores.size()));
  NumericVector fraction_positive(static_cast<R_xlen_t>(scores.size()));
  NumericVector mean_directional_weight(static_cast<R_xlen_t>(scores.size()));
  LogicalVector used_in_summary(static_cast<R_xlen_t>(scores.size()));

  target_cell.reserve(scores.size());
  source_cell.reserve(scores.size());
  target_cell_type.reserve(scores.size());
  source_cell_type.reserve(scores.size());

  for (std::size_t i = 0; i < scores.size(); ++i) {
    const auto& row = scores[i];
    target_cell.push_back(table.cells[static_cast<std::size_t>(row.target_cell)]);
    source_cell.push_back(table.cells[static_cast<std::size_t>(row.source_cell)]);
    target_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.target_type)]);
    source_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.source_type)]);
    factor[static_cast<R_xlen_t>(i)] = row.factor + 1;
    factor_count[static_cast<R_xlen_t>(i)] = row.factor_count;
    scored_molecules[static_cast<R_xlen_t>(i)] = row.scored_molecules;
    mean_score[static_cast<R_xlen_t>(i)] = row.mean_score;
    fraction_positive[static_cast<R_xlen_t>(i)] = row.fraction_positive;
    mean_directional_weight[static_cast<R_xlen_t>(i)] = row.mean_directional_weight;
    used_in_summary[static_cast<R_xlen_t>(i)] = row.used_in_summary;
  }

  return DataFrame::create(
      _["target_cell"] = target_cell,
      _["source_cell"] = source_cell,
      _["target_cell_type"] = target_cell_type,
      _["source_cell_type"] = source_cell_type,
      _["factor"] = factor,
      _["factor_count"] = factor_count,
      _["scored_molecules"] = scored_molecules,
      _["mean_score"] = mean_score,
      _["fraction_positive"] = fraction_positive,
      _["mean_directional_weight"] = mean_directional_weight,
      _["used_in_summary"] = used_in_summary);
}

DataFrame membrane_summary_to_df(const celladmix::MembraneTestResult& result) {
  const auto& summaries = result.summaries;
  std::vector<std::string> target_cell_type;
  std::vector<std::string> source_cell_type;
  IntegerVector factor(static_cast<R_xlen_t>(summaries.size()));
  IntegerVector n_pairs(static_cast<R_xlen_t>(summaries.size()));
  NumericVector mean_score(static_cast<R_xlen_t>(summaries.size()));
  NumericVector q75_score(static_cast<R_xlen_t>(summaries.size()));
  NumericVector fraction_positive_pairs(static_cast<R_xlen_t>(summaries.size()));
  NumericVector p_value(static_cast<R_xlen_t>(summaries.size()));
  NumericVector neg_log10_p(static_cast<R_xlen_t>(summaries.size()));

  target_cell_type.reserve(summaries.size());
  source_cell_type.reserve(summaries.size());
  for (std::size_t i = 0; i < summaries.size(); ++i) {
    const auto& row = summaries[i];
    target_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.target_type)]);
    source_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.source_type)]);
    factor[static_cast<R_xlen_t>(i)] = row.factor + 1;
    n_pairs[static_cast<R_xlen_t>(i)] = row.n_pairs;
    mean_score[static_cast<R_xlen_t>(i)] = row.mean_score;
    q75_score[static_cast<R_xlen_t>(i)] = row.q75_score;
    fraction_positive_pairs[static_cast<R_xlen_t>(i)] = row.fraction_positive_pairs;
    p_value[static_cast<R_xlen_t>(i)] = row.p_value;
    neg_log10_p[static_cast<R_xlen_t>(i)] = row.neg_log10_p;
  }

  return DataFrame::create(
      _["target_cell_type"] = target_cell_type,
      _["source_cell_type"] = source_cell_type,
      _["factor"] = factor,
      _["n_pairs"] = n_pairs,
      _["mean_score"] = mean_score,
      _["q75_score"] = q75_score,
      _["fraction_positive_pairs"] = fraction_positive_pairs,
      _["p_value"] = p_value,
      _["neg_log10_p"] = neg_log10_p);
}

DataFrame coherence_cell_scores_to_df(
    const celladmix::TranscriptTable& table,
    const celladmix::CoherenceTestResult& result) {
  const auto& scores = result.cell_scores;
  std::vector<std::string> target_cell;
  std::vector<std::string> target_cell_type;
  std::vector<std::string> source_cell_type;
  IntegerVector factor(static_cast<R_xlen_t>(scores.size()));
  IntegerVector factor_count(static_cast<R_xlen_t>(scores.size()));
  IntegerVector active_count(static_cast<R_xlen_t>(scores.size()));
  NumericVector mean_raw_score(static_cast<R_xlen_t>(scores.size()));
  NumericVector mean_coherence_score(static_cast<R_xlen_t>(scores.size()));
  NumericVector mean_score(static_cast<R_xlen_t>(scores.size()));
  NumericVector mean_null_score(static_cast<R_xlen_t>(scores.size()));
  NumericVector score_delta(static_cast<R_xlen_t>(scores.size()));
  NumericVector q75_score(static_cast<R_xlen_t>(scores.size()));
  NumericVector active_fraction(static_cast<R_xlen_t>(scores.size()));
  NumericVector mean_edge_weight(static_cast<R_xlen_t>(scores.size()));
  IntegerVector largest_patch_count(static_cast<R_xlen_t>(scores.size()));
  NumericVector largest_patch_fraction(static_cast<R_xlen_t>(scores.size()));
  NumericVector patch_score(static_cast<R_xlen_t>(scores.size()));
  NumericVector mean_null_patch_score(static_cast<R_xlen_t>(scores.size()));
  NumericVector patch_score_delta(static_cast<R_xlen_t>(scores.size()));
  NumericVector source_log_enrichment(static_cast<R_xlen_t>(scores.size()));
  NumericVector source_probability(static_cast<R_xlen_t>(scores.size()));
  LogicalVector used_in_summary(static_cast<R_xlen_t>(scores.size()));

  target_cell.reserve(scores.size());
  target_cell_type.reserve(scores.size());
  source_cell_type.reserve(scores.size());
  for (std::size_t i = 0; i < scores.size(); ++i) {
    const auto& row = scores[i];
    target_cell.push_back(table.cells[static_cast<std::size_t>(row.target_cell)]);
    target_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.target_type)]);
    source_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.source_type)]);
    factor[static_cast<R_xlen_t>(i)] = row.factor + 1;
    factor_count[static_cast<R_xlen_t>(i)] = row.factor_count;
    active_count[static_cast<R_xlen_t>(i)] = row.active_count;
    mean_raw_score[static_cast<R_xlen_t>(i)] = row.mean_raw_score;
    mean_coherence_score[static_cast<R_xlen_t>(i)] = row.mean_coherence_score;
    mean_score[static_cast<R_xlen_t>(i)] = row.mean_score;
    mean_null_score[static_cast<R_xlen_t>(i)] = row.mean_null_score;
    score_delta[static_cast<R_xlen_t>(i)] = row.score_delta;
    q75_score[static_cast<R_xlen_t>(i)] = row.q75_score;
    active_fraction[static_cast<R_xlen_t>(i)] = row.active_fraction;
    mean_edge_weight[static_cast<R_xlen_t>(i)] = row.mean_edge_weight;
    largest_patch_count[static_cast<R_xlen_t>(i)] = row.largest_patch_count;
    largest_patch_fraction[static_cast<R_xlen_t>(i)] = row.largest_patch_fraction;
    patch_score[static_cast<R_xlen_t>(i)] = row.patch_score;
    mean_null_patch_score[static_cast<R_xlen_t>(i)] = row.mean_null_patch_score;
    patch_score_delta[static_cast<R_xlen_t>(i)] = row.patch_score_delta;
    source_log_enrichment[static_cast<R_xlen_t>(i)] = row.source_log_enrichment;
    source_probability[static_cast<R_xlen_t>(i)] = row.source_probability;
    used_in_summary[static_cast<R_xlen_t>(i)] = row.used_in_summary;
  }

  return DataFrame::create(
      _["target_cell"] = target_cell,
      _["target_cell_type"] = target_cell_type,
      _["source_cell_type"] = source_cell_type,
      _["factor"] = factor,
      _["factor_count"] = factor_count,
      _["active_count"] = active_count,
      _["mean_raw_score"] = mean_raw_score,
      _["mean_coherence_score"] = mean_coherence_score,
      _["mean_score"] = mean_score,
      _["mean_null_score"] = mean_null_score,
      _["score_delta"] = score_delta,
      _["q75_score"] = q75_score,
      _["active_fraction"] = active_fraction,
      _["mean_edge_weight"] = mean_edge_weight,
      _["largest_patch_count"] = largest_patch_count,
      _["largest_patch_fraction"] = largest_patch_fraction,
      _["patch_score"] = patch_score,
      _["mean_null_patch_score"] = mean_null_patch_score,
      _["patch_score_delta"] = patch_score_delta,
      _["source_log_enrichment"] = source_log_enrichment,
      _["source_probability"] = source_probability,
      _["used_in_summary"] = used_in_summary);
}

DataFrame coherence_summary_to_df(const celladmix::CoherenceTestResult& result) {
  const auto& summaries = result.summaries;
  std::vector<std::string> target_cell_type;
  std::vector<std::string> source_cell_type;
  IntegerVector factor(static_cast<R_xlen_t>(summaries.size()));
  IntegerVector n_cells(static_cast<R_xlen_t>(summaries.size()));
  IntegerVector n_molecules(static_cast<R_xlen_t>(summaries.size()));
  IntegerVector n_active_molecules(static_cast<R_xlen_t>(summaries.size()));
  NumericVector active_fraction(static_cast<R_xlen_t>(summaries.size()));
  NumericVector mean_score(static_cast<R_xlen_t>(summaries.size()));
  NumericVector mean_null_score(static_cast<R_xlen_t>(summaries.size()));
  NumericVector mean_delta_score(static_cast<R_xlen_t>(summaries.size()));
  NumericVector q75_score(static_cast<R_xlen_t>(summaries.size()));
  NumericVector mean_patch_score(static_cast<R_xlen_t>(summaries.size()));
  NumericVector mean_null_patch_score(static_cast<R_xlen_t>(summaries.size()));
  NumericVector mean_delta_patch_score(static_cast<R_xlen_t>(summaries.size()));
  NumericVector mean_largest_patch_fraction(static_cast<R_xlen_t>(summaries.size()));
  NumericVector patch_p_value(static_cast<R_xlen_t>(summaries.size()));
  NumericVector patch_neg_log10_p(static_cast<R_xlen_t>(summaries.size()));
  NumericVector p_value(static_cast<R_xlen_t>(summaries.size()));
  NumericVector neg_log10_p(static_cast<R_xlen_t>(summaries.size()));
  NumericVector source_log_enrichment(static_cast<R_xlen_t>(summaries.size()));
  NumericVector source_probability(static_cast<R_xlen_t>(summaries.size()));

  target_cell_type.reserve(summaries.size());
  source_cell_type.reserve(summaries.size());
  for (std::size_t i = 0; i < summaries.size(); ++i) {
    const auto& row = summaries[i];
    target_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.target_type)]);
    source_cell_type.push_back(result.cell_types[static_cast<std::size_t>(row.source_type)]);
    factor[static_cast<R_xlen_t>(i)] = row.factor + 1;
    n_cells[static_cast<R_xlen_t>(i)] = row.n_cells;
    n_molecules[static_cast<R_xlen_t>(i)] = row.n_molecules;
    n_active_molecules[static_cast<R_xlen_t>(i)] = row.n_active_molecules;
    active_fraction[static_cast<R_xlen_t>(i)] = row.active_fraction;
    mean_score[static_cast<R_xlen_t>(i)] = row.mean_score;
    mean_null_score[static_cast<R_xlen_t>(i)] = row.mean_null_score;
    mean_delta_score[static_cast<R_xlen_t>(i)] = row.mean_delta_score;
    q75_score[static_cast<R_xlen_t>(i)] = row.q75_score;
    mean_patch_score[static_cast<R_xlen_t>(i)] = row.mean_patch_score;
    mean_null_patch_score[static_cast<R_xlen_t>(i)] = row.mean_null_patch_score;
    mean_delta_patch_score[static_cast<R_xlen_t>(i)] = row.mean_delta_patch_score;
    mean_largest_patch_fraction[static_cast<R_xlen_t>(i)] = row.mean_largest_patch_fraction;
    patch_p_value[static_cast<R_xlen_t>(i)] = row.patch_p_value;
    patch_neg_log10_p[static_cast<R_xlen_t>(i)] = row.patch_neg_log10_p;
    p_value[static_cast<R_xlen_t>(i)] = row.p_value;
    neg_log10_p[static_cast<R_xlen_t>(i)] = row.neg_log10_p;
    source_log_enrichment[static_cast<R_xlen_t>(i)] = row.source_log_enrichment;
    source_probability[static_cast<R_xlen_t>(i)] = row.source_probability;
  }

  return DataFrame::create(
      _["target_cell_type"] = target_cell_type,
      _["source_cell_type"] = source_cell_type,
      _["factor"] = factor,
      _["n_cells"] = n_cells,
      _["n_molecules"] = n_molecules,
      _["n_active_molecules"] = n_active_molecules,
      _["active_fraction"] = active_fraction,
      _["mean_score"] = mean_score,
      _["mean_null_score"] = mean_null_score,
      _["mean_delta_score"] = mean_delta_score,
      _["q75_score"] = q75_score,
      _["mean_patch_score"] = mean_patch_score,
      _["mean_null_patch_score"] = mean_null_patch_score,
      _["mean_delta_patch_score"] = mean_delta_patch_score,
      _["mean_largest_patch_fraction"] = mean_largest_patch_fraction,
      _["patch_p_value"] = patch_p_value,
      _["patch_neg_log10_p"] = patch_neg_log10_p,
      _["p_value"] = p_value,
      _["neg_log10_p"] = neg_log10_p,
      _["source_log_enrichment"] = source_log_enrichment,
      _["source_probability"] = source_probability);
}

}  // namespace

extern "C" SEXP _cellAdmixCore_celladmix_core_version() {
  return wrap(std::string(celladmix::kCelladmixVersion));
}

extern "C" SEXP _cellAdmixCore_celladmix_weighted_nmf_matrix(
    SEXP matrix_sexp,
    SEXP column_weights_sexp,
    SEXP rank_sexp,
    SEXP max_iterations_sexp,
    SEXP init_mode_sexp,
    SEXP row_groups_sexp,
    SEXP n_runs_sexp,
    SEXP num_threads_sexp,
    SEXP tolerance_sexp,
    SEXP update_epsilon_sexp,
    SEXP renormalize_each_iteration_sexp,
    SEXP lee_style_epsilon_sexp,
    SEXP random_init_sexp,
    SEXP seed_sexp) {
  try {
    const auto x = dense_matrix_from_r(matrix_sexp);
    std::vector<double> column_weights;
    if (Rf_isNull(column_weights_sexp)) {
      column_weights = celladmix::default_column_weights(x);
    } else {
      column_weights = as<std::vector<double>>(column_weights_sexp);
    }

    celladmix::WeightedNmfOptions options;
    options.rank = as<int>(rank_sexp);
    options.max_iterations = as<int>(max_iterations_sexp);
    const std::string init_mode = as<std::string>(init_mode_sexp);
    if (init_mode != "random" && init_mode != "cluster") {
      stop("init_mode must be 'random' or 'cluster'");
    }
    if (init_mode == "cluster") {
      if (Rf_isNull(row_groups_sexp)) {
        stop("row_groups must be provided when init_mode='cluster'");
      }
      options.init_groups = as<std::vector<int>>(row_groups_sexp);
    }
    options.n_runs = as<int>(n_runs_sexp);
    options.num_threads = as<int>(num_threads_sexp);
    options.tolerance = as<double>(tolerance_sexp);
    options.update_epsilon = as<double>(update_epsilon_sexp);
    options.renormalize_each_iteration = as<bool>(renormalize_each_iteration_sexp);
    options.lee_style_epsilon = as<bool>(lee_style_epsilon_sexp);
    options.random_init = as<std::string>(random_init_sexp);
    if (options.random_init != "legacy" && options.random_init != "rnmf") {
      stop("random_init must be 'legacy' or 'rnmf'");
    }
    options.seed = static_cast<unsigned int>(as<int>(seed_sexp));

    const auto fit = celladmix::weighted_nmf(x, column_weights, options);
    const double final_loss = fit.losses.empty() ? NA_REAL : fit.losses.back();
    return List::create(
        _["w"] = matrix_to_r(fit.w),
        _["h"] = matrix_to_r(fit.h),
        _["losses"] = wrap(fit.losses),
        _["final_loss"] = final_loss,
        _["selected_seed"] = static_cast<double>(fit.selected_seed),
        _["selected_run"] = static_cast<double>(fit.selected_run + 1),
        _["candidate_final_losses"] = wrap(fit.candidate_final_losses),
        _["selected_factor_stability"] = wrap(fit.selected_factor_stability),
        _["column_weights"] = wrap(column_weights));
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_weighted_nmf_matrix: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_weighted_nmf_sparse_matrix(
    SEXP matrix_sexp,
    SEXP column_weights_sexp,
    SEXP rank_sexp,
    SEXP max_iterations_sexp,
    SEXP init_mode_sexp,
    SEXP row_groups_sexp,
    SEXP n_runs_sexp,
    SEXP num_threads_sexp,
    SEXP tolerance_sexp,
    SEXP update_epsilon_sexp,
    SEXP renormalize_each_iteration_sexp,
    SEXP lee_style_epsilon_sexp,
    SEXP random_init_sexp,
    SEXP seed_sexp) {
  try {
    const auto x = sparse_row_matrix_from_dgc(matrix_sexp);
    std::vector<double> column_weights;
    if (Rf_isNull(column_weights_sexp)) {
      column_weights = celladmix::default_column_weights(x);
    } else {
      column_weights = as<std::vector<double>>(column_weights_sexp);
    }

    celladmix::WeightedNmfOptions options;
    options.rank = as<int>(rank_sexp);
    options.max_iterations = as<int>(max_iterations_sexp);
    const std::string init_mode = as<std::string>(init_mode_sexp);
    if (init_mode != "random" && init_mode != "cluster") {
      stop("init_mode must be 'random' or 'cluster'");
    }
    if (init_mode == "cluster") {
      if (Rf_isNull(row_groups_sexp)) {
        stop("row_groups must be provided when init_mode='cluster'");
      }
      options.init_groups = as<std::vector<int>>(row_groups_sexp);
    }
    options.n_runs = as<int>(n_runs_sexp);
    options.num_threads = as<int>(num_threads_sexp);
    options.tolerance = as<double>(tolerance_sexp);
    options.update_epsilon = as<double>(update_epsilon_sexp);
    options.renormalize_each_iteration = as<bool>(renormalize_each_iteration_sexp);
    options.lee_style_epsilon = as<bool>(lee_style_epsilon_sexp);
    options.random_init = as<std::string>(random_init_sexp);
    if (options.random_init != "legacy" && options.random_init != "rnmf") {
      stop("random_init must be 'legacy' or 'rnmf'");
    }
    options.seed = static_cast<unsigned int>(as<int>(seed_sexp));

    const auto fit = celladmix::weighted_nmf(x, column_weights, options);
    const double final_loss = fit.losses.empty() ? NA_REAL : fit.losses.back();
    return List::create(
        _["w"] = matrix_to_r(fit.w),
        _["h"] = matrix_to_r(fit.h),
        _["losses"] = wrap(fit.losses),
        _["final_loss"] = final_loss,
        _["selected_seed"] = static_cast<double>(fit.selected_seed),
        _["selected_run"] = static_cast<double>(fit.selected_run + 1),
        _["candidate_final_losses"] = wrap(fit.candidate_final_losses),
        _["selected_factor_stability"] = wrap(fit.selected_factor_stability),
        _["column_weights"] = wrap(column_weights));
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_weighted_nmf_sparse_matrix: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_sparse_nmf_matrix(
    SEXP matrix_sexp,
    SEXP column_weights_sexp,
    SEXP rank_sexp,
    SEXP max_iterations_sexp,
    SEXP init_mode_sexp,
    SEXP row_groups_sexp,
    SEXP n_runs_sexp,
    SEXP num_threads_sexp,
    SEXP tolerance_sexp,
    SEXP update_epsilon_sexp,
    SEXP random_init_sexp,
    SEXP loss_mode_sexp,
    SEXP h_l1_penalty_sexp,
    SEXP h_diversity_penalty_sexp,
    SEXP seed_sexp) {
  try {
    const auto x = sparse_row_matrix_from_dgc(matrix_sexp);
    std::vector<double> column_weights;
    if (Rf_isNull(column_weights_sexp)) {
      column_weights = celladmix::default_column_weights(x);
    } else {
      column_weights = as<std::vector<double>>(column_weights_sexp);
    }

    celladmix::SparseNmfOptions options;
    options.rank = as<int>(rank_sexp);
    options.max_iterations = as<int>(max_iterations_sexp);
    options.init_mode = as<std::string>(init_mode_sexp);
    if (options.init_mode != "random" && options.init_mode != "cluster") {
      stop("init_mode must be 'random' or 'cluster'");
    }
    if (options.init_mode == "cluster") {
      if (Rf_isNull(row_groups_sexp)) {
        stop("row_groups must be provided when init_mode='cluster'");
      }
      options.init_groups = as<std::vector<int>>(row_groups_sexp);
    }
    options.n_runs = as<int>(n_runs_sexp);
    options.num_threads = as<int>(num_threads_sexp);
    options.tolerance = as<double>(tolerance_sexp);
    options.update_epsilon = as<double>(update_epsilon_sexp);
    options.random_init = as<std::string>(random_init_sexp);
    options.loss_mode = as<std::string>(loss_mode_sexp);
    options.h_l1_penalty = as<double>(h_l1_penalty_sexp);
    options.h_diversity_penalty = as<double>(h_diversity_penalty_sexp);
    options.seed = static_cast<unsigned int>(as<int>(seed_sexp));

    const auto fit = celladmix::sparse_nmf(x, column_weights, options);
    return List::create(
        _["w"] = matrix_to_r(fit.w),
        _["h"] = matrix_to_r(fit.h),
        _["losses"] = wrap(fit.losses),
        _["final_objective"] = fit.final_objective,
        _["selected_seed"] = static_cast<double>(fit.selected_seed),
        _["selected_run"] = static_cast<double>(fit.selected_run + 1),
        _["candidate_final_objectives"] = wrap(fit.candidate_final_objectives),
        _["candidate_best_match_correlations"] = wrap(fit.candidate_best_match_correlations),
        _["selected_factor_stability"] = wrap(fit.selected_factor_stability),
        _["candidate_final_objective_mean"] = fit.candidate_final_objective_mean,
        _["candidate_final_objective_sd"] = fit.candidate_final_objective_sd,
        _["candidate_best_match_correlation_mean"] = fit.candidate_best_match_correlation_mean,
        _["column_weights"] = wrap(column_weights),
        _["loss_mode"] = options.loss_mode,
        _["h_l1_penalty"] = options.h_l1_penalty,
        _["h_diversity_penalty"] = options.h_diversity_penalty);
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_sparse_nmf_matrix: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_build_reference_ncv_matrix(
    SEXP transcripts_sexp,
    SEXP k_sexp,
    SEXP query_ids_sexp,
    SEXP row_sample_n_sexp,
    SEXP seed_sexp) {
  try {
    return build_reference_ncv_matrix(
        as<DataFrame>(transcripts_sexp),
        as<int>(k_sexp),
        query_ids_sexp,
        optional_int_sexp(row_sample_n_sexp, -1),
        optional_seed_sexp(seed_sexp, 1U));
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_build_reference_ncv_matrix: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_build_pipeline_ncv_matrix(
    SEXP transcripts_sexp,
    SEXP k_sexp,
    SEXP query_ids_sexp) {
  try {
    return build_pipeline_ncv_matrix(
        as<DataFrame>(transcripts_sexp),
        as<int>(k_sexp),
        query_ids_sexp);
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_build_pipeline_ncv_matrix: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_load_run_training_transcript_ids(
    SEXP path_sexp) {
  try {
    return load_run_training_transcript_ids(as<std::string>(path_sexp));
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_load_run_training_transcript_ids: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_export_training_ncv_matrix(
    SEXP path_sexp) {
  try {
    return export_training_ncv_matrix(as<std::string>(path_sexp));
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_export_training_ncv_matrix: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_read_xenium_manifest(SEXP path_sexp) {
  try {
    return manifest_to_r(celladmix::read_xenium_manifest(as<std::string>(path_sexp)));
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_read_xenium_manifest: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_probe_xenium(SEXP path_sexp) {
  try {
    return probe_xenium_source_to_r(as<std::string>(path_sexp));
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_probe_xenium: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_probe_tabular(
    SEXP molecules_path_sexp,
    SEXP segmentation_mask_path_sexp) {
  try {
    return probe_tabular_source_to_r(
        as<std::string>(molecules_path_sexp),
        optional_string_sexp(segmentation_mask_path_sexp).value_or(""));
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_probe_tabular: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_read_input_store(SEXP store_dir_sexp) {
  try {
    const auto manifest = celladmix::read_input_store_manifest(as<std::string>(store_dir_sexp));
    auto out = input_store_manifest_to_r(manifest);
    out.attr("class") = CharacterVector::create("celladmix_store");
    return out;
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_read_input_store: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_collect_input_store_counts(SEXP store_dir_sexp) {
  try {
    const auto counts = celladmix::load_input_store_counts(as<std::string>(store_dir_sexp));
    return sparse_counts_to_r_list(counts);
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_collect_input_store_counts: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_identify_domains_store(
    SEXP store_dir_sexp,
    SEXP cell_types_sexp,
    SEXP scales_sexp,
    SEXP transform_sexp,
    SEXP include_self_sexp,
    SEXP self_weight_sexp,
    SEXP include_center_type_sexp,
    SEXP center_type_weight_sexp,
    SEXP n_domains_sexp,
    SEXP smooth_k_sexp,
    SEXP smooth_lambda_sexp,
    SEXP max_edge_distance_sexp,
    SEXP min_region_size_sexp,
    SEXP seed_sexp,
    SEXP num_threads_sexp) {
  try {
    auto all_cells = celladmix::load_input_store_cells(as<std::string>(store_dir_sexp));
    apply_cell_type_overrides(all_cells, cell_types_sexp);
    const std::size_t original_n = all_cells.cell_ids.size();
    if (all_cells.cell_types.size() != original_n) {
      stop("Domain detection requires cell annotations; provide annotations/annotation_col or prepare a store with cell_type");
    }

    celladmix::CellTable cells;
    cells.cell_ids.reserve(original_n);
    cells.cell_types.reserve(original_n);
    cells.centroid_x.reserve(original_n);
    cells.centroid_y.reserve(original_n);
    cells.centroid_z.reserve(original_n);
    if (all_cells.sample_ids.size() == original_n) cells.sample_ids.reserve(original_n);
    if (all_cells.fov_ids.size() == original_n) cells.fov_ids.reserve(original_n);
    for (std::size_t i = 0; i < original_n; ++i) {
      const auto& label = all_cells.cell_types[i];
      if (label.empty()) {
        continue;
      }
      cells.cell_ids.push_back(all_cells.cell_ids[i]);
      cells.cell_types.push_back(label);
      cells.centroid_x.push_back(all_cells.centroid_x[i]);
      cells.centroid_y.push_back(all_cells.centroid_y[i]);
      if (all_cells.centroid_z.size() == original_n) {
        cells.centroid_z.push_back(all_cells.centroid_z[i]);
      }
      if (all_cells.sample_ids.size() == original_n) {
        cells.sample_ids.push_back(all_cells.sample_ids[i]);
      }
      if (all_cells.fov_ids.size() == original_n) {
        cells.fov_ids.push_back(all_cells.fov_ids[i]);
      }
    }
    const std::size_t n = cells.cell_ids.size();
    if (n == 0) {
      stop("Domain detection found no annotated cells");
    }

    std::unordered_map<std::string, int> type_index;
    std::vector<std::string> type_names;
    std::vector<int> annotation_ids(n, -1);
    type_index.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
      const auto& label = cells.cell_types[i];
      if (label.empty()) {
        stop("Domain detection encountered an empty annotation for cell '%s'", cells.cell_ids[i].c_str());
      }
      auto it = type_index.find(label);
      if (it == type_index.end()) {
        const int idx = static_cast<int>(type_names.size());
        type_index.emplace(label, idx);
        type_names.push_back(label);
        annotation_ids[i] = idx;
      } else {
        annotation_ids[i] = it->second;
      }
    }

    celladmix::DomainOptions options;
    options.scales = as<std::vector<int>>(scales_sexp);
    options.transform = as<std::string>(transform_sexp);
    options.include_self = as<bool>(include_self_sexp);
    options.self_weight = as<double>(self_weight_sexp);
    options.include_center_type = as<bool>(include_center_type_sexp);
    options.center_type_weight = as<double>(center_type_weight_sexp);
    options.n_domains = as<int>(n_domains_sexp);
    options.smooth_k = as<int>(smooth_k_sexp);
    options.smooth_lambda = as<double>(smooth_lambda_sexp);
    const double max_edge_distance = as<double>(max_edge_distance_sexp);
    options.max_edge_distance = NumericVector::is_na(max_edge_distance) ? -1.0 : max_edge_distance;
    options.min_region_size = as<int>(min_region_size_sexp);
    options.seed = optional_seed_sexp(seed_sexp);
    options.num_threads = as<int>(num_threads_sexp);

    const auto result = celladmix::identify_domains(
        cells,
        annotation_ids,
        static_cast<int>(type_names.size()),
        options);

    CharacterVector cell_id(static_cast<R_xlen_t>(n));
    CharacterVector sample_id(static_cast<R_xlen_t>(n));
    CharacterVector fov_id(static_cast<R_xlen_t>(n));
    CharacterVector annotation(static_cast<R_xlen_t>(n));
    NumericVector x(static_cast<R_xlen_t>(n));
    NumericVector y(static_cast<R_xlen_t>(n));
    NumericVector z(static_cast<R_xlen_t>(n));
    IntegerVector domain(static_cast<R_xlen_t>(n));
    IntegerVector region(static_cast<R_xlen_t>(n));
    for (std::size_t i = 0; i < n; ++i) {
      cell_id[static_cast<R_xlen_t>(i)] = cells.cell_ids[i];
      if (cells.sample_ids.size() == n) {
        sample_id[static_cast<R_xlen_t>(i)] = cells.sample_ids[i];
      } else {
        sample_id[static_cast<R_xlen_t>(i)] = NA_STRING;
      }
      if (cells.fov_ids.size() == n) {
        fov_id[static_cast<R_xlen_t>(i)] = cells.fov_ids[i];
      } else {
        fov_id[static_cast<R_xlen_t>(i)] = NA_STRING;
      }
      annotation[static_cast<R_xlen_t>(i)] = cells.cell_types[i];
      x[static_cast<R_xlen_t>(i)] = cells.centroid_x[i];
      y[static_cast<R_xlen_t>(i)] = cells.centroid_y[i];
      z[static_cast<R_xlen_t>(i)] =
          cells.centroid_z.size() == n ? cells.centroid_z[i] : 0.0;
      domain[static_cast<R_xlen_t>(i)] = result.labels[i] + 1;
      region[static_cast<R_xlen_t>(i)] = result.region_ids[i] + 1;
    }

    NumericMatrix composition = matrix_to_r(result.domain_composition);
    CharacterVector domain_names(result.domain_composition.rows());
    for (int i = 0; i < result.domain_composition.rows(); ++i) {
      domain_names[i] = "D" + std::to_string(i + 1);
    }
    composition.attr("dimnames") = List::create(domain_names, type_names);

    IntegerVector domain_sizes(result.domain_sizes.begin(), result.domain_sizes.end());
    IntegerVector region_sizes(result.region_sizes.begin(), result.region_sizes.end());
    domain_sizes.attr("names") = domain_names;
    CharacterVector region_names(result.region_sizes.size());
    for (R_xlen_t i = 0; i < region_names.size(); ++i) {
      region_names[i] = "R" + std::to_string(i + 1);
    }
    region_sizes.attr("names") = region_names;

    return List::create(
        _["cell_id"] = cell_id,
        _["x"] = x,
        _["y"] = y,
        _["z"] = z,
        _["sample_id"] = sample_id,
        _["fov_id"] = fov_id,
        _["annotation"] = annotation,
        _["domain"] = domain,
        _["region"] = region,
        _["domain_composition"] = composition,
        _["domain_sizes"] = domain_sizes,
        _["region_sizes"] = region_sizes,
        _["cell_type_levels"] = type_names,
      _["spatial_coherence"] = result.spatial_coherence,
      _["boundary_fraction"] = result.boundary_fraction,
      _["kmeans_iterations"] = result.kmeans_iterations,
      _["smooth_iterations"] = result.smooth_iterations,
      _["n_store_cells"] = static_cast<double>(original_n),
      _["n_dropped_unannotated"] = static_cast<double>(original_n - n));
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_identify_domains_store: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_build_xenium_store(
    SEXP path_sexp,
    SEXP crops_sexp,
    SEXP cell_filter_sexp,
    SEXP gene_filter_sexp,
    SEXP min_qv_sexp,
    SEXP keep_unassigned_sexp,
    SEXP keep_non_gene_sexp,
    SEXP read_cells_sexp,
    SEXP prefer_parquet_sexp,
    SEXP store_dir_sexp,
    SEXP materialize_molecules_sexp,
    SEXP force_sexp,
    SEXP num_threads_sexp,
    SEXP row_group_size_sexp,
    SEXP verbose_sexp) {
  try {
    const auto start = std::chrono::steady_clock::now();
    auto stage_start = start;
    const bool verbose = as<bool>(verbose_sexp);
    auto log_info = [&](const std::string& message, double stage_sec = -1.0) {
      if (!verbose) return;
      const auto now = std::chrono::system_clock::now();
      const auto tt = std::chrono::system_clock::to_time_t(now);
      std::tm tm{};
#ifdef _WIN32
      localtime_s(&tm, &tt);
#else
      localtime_r(&tt, &tm);
#endif
      std::ostringstream line;
      line << "[INFO " << std::put_time(&tm, "%H:%M:%S") << " +"
           << std::fixed << std::setprecision(3)
           << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
           << "s] " << message;
      if (stage_sec >= 0.0) line << " (" << std::fixed << std::setprecision(3) << stage_sec << "s)";
      Rcpp::Rcout << line.str() << std::endl;
    };

    celladmix::XeniumLoadOptions load_options;
    load_options.crops = crop_boxes_from_df(crops_sexp);
    load_options.cell_filter = optional_string_vector_sexp(cell_filter_sexp);
    load_options.gene_filter = optional_string_vector_sexp(gene_filter_sexp);
    const double min_qv = as<double>(min_qv_sexp);
    if (!NumericVector::is_na(min_qv)) load_options.min_qv = min_qv;
    load_options.keep_unassigned = as<bool>(keep_unassigned_sexp);
    load_options.keep_non_gene = as<bool>(keep_non_gene_sexp);
    load_options.read_cells = as<bool>(read_cells_sexp);
    load_options.prefer_parquet = as<bool>(prefer_parquet_sexp);
    load_options.progress = log_info;

    celladmix::InputStoreBuildOptions store_options;
    store_options.store_dir = as<std::string>(store_dir_sexp);
    store_options.materialize_molecules = as<bool>(materialize_molecules_sexp);
    store_options.force = as<bool>(force_sexp);
    store_options.num_threads = as<int>(num_threads_sexp);
    store_options.parquet_row_group_size = as<int>(row_group_size_sexp);

    stage_start = std::chrono::steady_clock::now();
    const auto manifest = celladmix::build_xenium_input_store(
        as<std::string>(path_sexp),
        load_options,
        store_options);
    log_info(
        "Built Xenium input store: " + std::to_string(manifest.n_molecules) +
            " molecules, " + std::to_string(manifest.n_cells) + " cells",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());
    auto out = input_store_manifest_to_r(manifest);
    out.attr("class") = CharacterVector::create("celladmix_store");
    return out;
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_build_xenium_store: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_build_tabular_store(
    SEXP molecules_path_sexp,
    SEXP crops_sexp,
    SEXP min_qv_sexp,
    SEXP keep_unassigned_sexp,
    SEXP x_col_sexp,
    SEXP y_col_sexp,
    SEXP z_col_sexp,
    SEXP gene_col_sexp,
    SEXP qv_col_sexp,
    SEXP cell_id_col_sexp,
    SEXP segmentation_mask_path_sexp,
    SEXP cell_type_col_sexp,
    SEXP cell_metadata_path_sexp,
    SEXP cell_metadata_cell_id_col_sexp,
    SEXP cell_metadata_cell_type_col_sexp,
    SEXP sample_id_col_sexp,
    SEXP fov_id_col_sexp,
    SEXP sample_id_sexp,
    SEXP fov_id_sexp,
    SEXP store_dir_sexp,
    SEXP materialize_molecules_sexp,
    SEXP force_sexp,
    SEXP num_threads_sexp,
    SEXP row_group_size_sexp,
    SEXP verbose_sexp) {
  try {
    const auto start = std::chrono::steady_clock::now();
    const bool verbose = as<bool>(verbose_sexp);
    auto log_info = [&](const std::string& message, double stage_sec = -1.0) {
      if (!verbose) return;
      const auto now = std::chrono::system_clock::now();
      const auto tt = std::chrono::system_clock::to_time_t(now);
      std::tm tm{};
#ifdef _WIN32
      localtime_s(&tm, &tt);
#else
      localtime_r(&tt, &tm);
#endif
      std::ostringstream line;
      line << "[INFO " << std::put_time(&tm, "%H:%M:%S") << " +"
           << std::fixed << std::setprecision(3)
           << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
           << "s] " << message;
      if (stage_sec >= 0.0) line << " (" << std::fixed << std::setprecision(3) << stage_sec << "s)";
      Rcpp::Rcout << line.str() << std::endl;
    };

    celladmix::TabularSourceSpec source;
    source.molecules_path = as<std::string>(molecules_path_sexp);
    source.x_col = as<std::string>(x_col_sexp);
    source.y_col = as<std::string>(y_col_sexp);
    source.z_col = optional_string_sexp(z_col_sexp).value_or("");
    source.gene_col = as<std::string>(gene_col_sexp);
    source.qv_col = optional_string_sexp(qv_col_sexp).value_or("");
    source.cell_id_col = optional_string_sexp(cell_id_col_sexp).value_or("");
    source.segmentation_mask_path = optional_string_sexp(segmentation_mask_path_sexp).value_or("");
    source.cell_type_col = optional_string_sexp(cell_type_col_sexp).value_or("");
    source.cell_metadata_path = optional_string_sexp(cell_metadata_path_sexp).value_or("");
    source.cell_metadata_cell_id_col =
        optional_string_sexp(cell_metadata_cell_id_col_sexp).value_or("cell");
    source.cell_metadata_cell_type_col =
        optional_string_sexp(cell_metadata_cell_type_col_sexp).value_or("");
    source.sample_id_col = optional_string_sexp(sample_id_col_sexp).value_or("");
    source.fov_id_col = optional_string_sexp(fov_id_col_sexp).value_or("");
    source.sample_id = optional_string_sexp(sample_id_sexp).value_or("");
    source.fov_id = optional_string_sexp(fov_id_sexp).value_or("");

    celladmix::TabularLoadOptions load_options;
    load_options.crops = crop_boxes_from_df(crops_sexp);
    const double min_qv = as<double>(min_qv_sexp);
    if (!NumericVector::is_na(min_qv)) load_options.min_qv = min_qv;
    load_options.keep_unassigned = as<bool>(keep_unassigned_sexp);

    celladmix::InputStoreBuildOptions store_options;
    store_options.store_dir = as<std::string>(store_dir_sexp);
    store_options.materialize_molecules = as<bool>(materialize_molecules_sexp);
    store_options.force = as<bool>(force_sexp);
    store_options.num_threads = as<int>(num_threads_sexp);
    store_options.parquet_row_group_size = as<int>(row_group_size_sexp);

    const auto stage_start = std::chrono::steady_clock::now();
    const auto manifest = celladmix::build_tabular_input_store(source, load_options, store_options);
    log_info(
        "Built tabular input store: " + std::to_string(manifest.n_molecules) +
            " molecules, " + std::to_string(manifest.n_cells) + " cells",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());
    auto out = input_store_manifest_to_r(manifest);
    out.attr("class") = CharacterVector::create("celladmix_store");
    return out;
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_build_tabular_store: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_cluster_store(
    SEXP store_dir_sexp,
    SEXP min_molecules_sexp,
    SEXP min_genes_sexp,
    SEXP cells_max_sexp,
    SEXP n_variable_genes_sexp,
    SEXP pca_dims_sexp,
    SEXP graph_k_sexp,
    SEXP cluster_resolution_sexp,
    SEXP compute_umap_sexp,
    SEXP umap_neighbors_sexp,
    SEXP umap_epochs_sexp,
    SEXP num_threads_sexp,
    SEXP umap_parallel_optimization_sexp,
    SEXP normalization_scale_sexp,
    SEXP seed_sexp,
    SEXP clusters_out_sexp,
    SEXP embedding_out_sexp,
    SEXP verbose_sexp) {
  try {
    const auto cluster_start = std::chrono::steady_clock::now();
    const bool verbose = as<bool>(verbose_sexp);
    auto log_info = [&](const std::string& message, double stage_sec = -1.0) {
      if (!verbose) return;
      const auto now = std::chrono::system_clock::now();
      const auto tt = std::chrono::system_clock::to_time_t(now);
      std::tm tm{};
#ifdef _WIN32
      localtime_s(&tm, &tt);
#else
      localtime_r(&tt, &tm);
#endif
      std::ostringstream line;
      line << "[INFO " << std::put_time(&tm, "%H:%M:%S") << " +"
           << std::fixed << std::setprecision(3)
           << std::chrono::duration<double>(std::chrono::steady_clock::now() - cluster_start).count()
           << "s] " << message;
      if (stage_sec >= 0.0) line << " (" << std::fixed << std::setprecision(3) << stage_sec << "s)";
      Rcpp::Rcout << line.str() << std::endl;
    };
    auto stage_start = std::chrono::steady_clock::now();
    const std::string store_dir = as<std::string>(store_dir_sexp);
    const auto store_manifest = celladmix::read_input_store_manifest(store_dir);
    if (!store_manifest.has_cell_gene_counts) {
      throw std::runtime_error("cell clustering requires an input store with cell-gene counts");
    }
    const auto counts = celladmix::load_input_store_counts(store_dir);
    log_info(
        "Loaded input-store cell-gene counts: " + std::to_string(counts.transcript_counts.size()) +
            " cells",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    celladmix::CellClusteringOptions options;
    options.min_molecules = as<int>(min_molecules_sexp);
    options.min_genes = as<int>(min_genes_sexp);
    options.cells_max = optional_int_sexp(cells_max_sexp, -1);
    options.n_variable_genes = as<int>(n_variable_genes_sexp);
    options.pca_dims = as<int>(pca_dims_sexp);
    options.graph_k = as<int>(graph_k_sexp);
    options.cluster_resolution = as<double>(cluster_resolution_sexp);
    options.compute_umap = as<bool>(compute_umap_sexp);
    options.umap_neighbors = as<int>(umap_neighbors_sexp);
    options.umap_epochs = as<int>(umap_epochs_sexp);
    options.num_threads = as<int>(num_threads_sexp);
    options.umap_parallel_optimization = as<bool>(umap_parallel_optimization_sexp);
    options.normalization_scale = as<double>(normalization_scale_sexp);
    options.seed = optional_seed_sexp(seed_sexp, 1U);
    options.progress = log_info;

    stage_start = std::chrono::steady_clock::now();
    const auto result = celladmix::cluster_cell_counts(counts, options);
    log_info(
        "Finished store-backed cell clustering",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    stage_start = std::chrono::steady_clock::now();
    const std::string clusters_out = as<std::string>(clusters_out_sexp);
    const std::string embedding_out = as<std::string>(embedding_out_sexp);
    std::filesystem::create_directories(std::filesystem::path(clusters_out).parent_path());
    write_cell_clusters_parquet(clusters_out, result);
    if (result.has_umap) {
      std::filesystem::create_directories(std::filesystem::path(embedding_out).parent_path());
      write_cell_embedding_parquet(embedding_out, result);
    }
    log_info("Wrote cell clustering outputs", std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    return List::create(
        _["clusters"] = cell_clusters_to_df(result),
        _["embedding"] = cell_embedding_to_df(result),
        _["clusters_path"] = clusters_out,
        _["embedding_path"] = result.has_umap ? wrap(embedding_out) : R_NilValue,
        _["n_cells"] = static_cast<double>(result.cells.size()),
        _["n_clusters"] = static_cast<double>(
            result.clusters.empty() ? 0 : *std::max_element(result.clusters.begin(), result.clusters.end())),
        _["n_variable_genes"] = static_cast<double>(result.variable_genes.size()),
        _["variable_genes"] = wrap(result.variable_genes),
        _["pca_variance_explained"] = wrap(result.pca_variance_explained));
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_cluster_store: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_cluster_counts_matrix(
    SEXP p_sexp,
    SEXP i_sexp,
    SEXP x_sexp,
    SEXP genes_sexp,
    SEXP cells_sexp,
    SEXP min_molecules_sexp,
    SEXP min_genes_sexp,
    SEXP cells_max_sexp,
    SEXP n_variable_genes_sexp,
    SEXP pca_dims_sexp,
    SEXP graph_k_sexp,
    SEXP cluster_resolution_sexp,
    SEXP compute_umap_sexp,
    SEXP umap_neighbors_sexp,
    SEXP umap_epochs_sexp,
    SEXP num_threads_sexp,
    SEXP umap_parallel_optimization_sexp,
    SEXP normalization_scale_sexp,
    SEXP seed_sexp) {
  try {
    celladmix::CellCountMatrix counts;
    counts.indptr = as<std::vector<int>>(p_sexp);
    counts.indices = as<std::vector<int>>(i_sexp);
    counts.values = as<std::vector<double>>(x_sexp);
    counts.genes = as<std::vector<std::string>>(genes_sexp);
    counts.cells.cell_ids = as<std::vector<std::string>>(cells_sexp);

    counts.transcript_counts.resize(counts.cells.cell_ids.size(), 0);
    counts.detected_genes.resize(counts.cells.cell_ids.size(), 0);
    for (std::size_t c = 0; c + 1 < counts.indptr.size(); ++c) {
      double total = 0.0;
      for (int p = counts.indptr[c]; p < counts.indptr[c + 1]; ++p) {
        total += counts.values[static_cast<std::size_t>(p)];
      }
      counts.transcript_counts[c] = static_cast<int>(total);
      counts.detected_genes[c] = counts.indptr[c + 1] - counts.indptr[c];
    }
    counts.crop_ids.assign(counts.cells.cell_ids.size(), "");
    counts.cells.centroid_x.assign(counts.cells.cell_ids.size(), 0.0);
    counts.cells.centroid_y.assign(counts.cells.cell_ids.size(), 0.0);
    counts.cells.centroid_z.assign(counts.cells.cell_ids.size(), 0.0);

    celladmix::CellClusteringOptions options;
    options.min_molecules = as<int>(min_molecules_sexp);
    options.min_genes = as<int>(min_genes_sexp);
    options.cells_max = optional_int_sexp(cells_max_sexp, -1);
    options.n_variable_genes = as<int>(n_variable_genes_sexp);
    options.pca_dims = as<int>(pca_dims_sexp);
    options.graph_k = as<int>(graph_k_sexp);
    options.cluster_resolution = as<double>(cluster_resolution_sexp);
    options.compute_umap = as<bool>(compute_umap_sexp);
    options.umap_neighbors = as<int>(umap_neighbors_sexp);
    options.umap_epochs = as<int>(umap_epochs_sexp);
    options.num_threads = as<int>(num_threads_sexp);
    options.umap_parallel_optimization = as<bool>(umap_parallel_optimization_sexp);
    options.normalization_scale = as<double>(normalization_scale_sexp);
    options.seed = optional_seed_sexp(seed_sexp, 1U);

    const auto result = celladmix::cluster_cell_counts(counts, options);
    return List::create(
        _["clusters"] = cell_clusters_to_df(result),
        _["embedding"] = cell_embedding_to_df(result),
        _["n_cells"] = static_cast<double>(result.cells.size()),
        _["n_clusters"] = static_cast<double>(
            result.clusters.empty() ? 0 : *std::max_element(result.clusters.begin(), result.clusters.end())),
        _["n_variable_genes"] = static_cast<double>(result.variable_genes.size()),
        _["variable_genes"] = wrap(result.variable_genes),
        _["pca_variance_explained"] = wrap(result.pca_variance_explained));
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_cluster_counts_matrix: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_cluster_run_counts(
    SEXP run_dir_sexp,
    SEXP min_molecules_sexp,
    SEXP min_genes_sexp,
    SEXP cells_max_sexp,
    SEXP n_variable_genes_sexp,
    SEXP pca_dims_sexp,
    SEXP graph_k_sexp,
    SEXP cluster_resolution_sexp,
    SEXP compute_umap_sexp,
    SEXP umap_neighbors_sexp,
    SEXP umap_epochs_sexp,
    SEXP num_threads_sexp,
    SEXP umap_parallel_optimization_sexp,
    SEXP normalization_scale_sexp,
    SEXP seed_sexp,
    SEXP verbose_sexp) {
  try {
    const auto cluster_start = std::chrono::steady_clock::now();
    const bool verbose = as<bool>(verbose_sexp);
    auto log_info = [&](const std::string& message, double stage_sec = -1.0) {
      if (!verbose) return;
      const auto now = std::chrono::system_clock::now();
      const auto tt = std::chrono::system_clock::to_time_t(now);
      std::tm tm{};
#ifdef _WIN32
      localtime_s(&tm, &tt);
#else
      localtime_r(&tt, &tm);
#endif
      std::ostringstream line;
      line << "[INFO " << std::put_time(&tm, "%H:%M:%S") << " +"
           << std::fixed << std::setprecision(3)
           << std::chrono::duration<double>(std::chrono::steady_clock::now() - cluster_start).count()
           << "s] " << message;
      if (stage_sec >= 0.0) line << " (" << std::fixed << std::setprecision(3) << stage_sec << "s)";
      Rcpp::Rcout << line.str() << std::endl;
    };

    auto stage_start = std::chrono::steady_clock::now();
    const std::string run_dir = as<std::string>(run_dir_sexp);
    const auto counts = celladmix::collect_run_counts(run_dir);
    log_info(
        "Loaded run cell-gene counts: " + std::to_string(counts.transcript_counts.size()) +
            " cells",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    celladmix::CellClusteringOptions options;
    options.min_molecules = as<int>(min_molecules_sexp);
    options.min_genes = as<int>(min_genes_sexp);
    options.cells_max = optional_int_sexp(cells_max_sexp, -1);
    options.n_variable_genes = as<int>(n_variable_genes_sexp);
    options.pca_dims = as<int>(pca_dims_sexp);
    options.graph_k = as<int>(graph_k_sexp);
    options.cluster_resolution = as<double>(cluster_resolution_sexp);
    options.compute_umap = as<bool>(compute_umap_sexp);
    options.umap_neighbors = as<int>(umap_neighbors_sexp);
    options.umap_epochs = as<int>(umap_epochs_sexp);
    options.num_threads = as<int>(num_threads_sexp);
    options.umap_parallel_optimization = as<bool>(umap_parallel_optimization_sexp);
    options.normalization_scale = as<double>(normalization_scale_sexp);
    options.seed = optional_seed_sexp(seed_sexp, 1U);
    options.progress = log_info;

    stage_start = std::chrono::steady_clock::now();
    const auto result = celladmix::cluster_cell_counts(counts, options);
    log_info(
        "Finished run-count cell-state embedding",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    return List::create(
        _["clusters"] = cell_clusters_to_df(result),
        _["embedding"] = cell_embedding_to_df(result),
        _["n_cells"] = static_cast<double>(result.cells.size()),
        _["n_clusters"] = static_cast<double>(
            result.clusters.empty() ? 0 : *std::max_element(result.clusters.begin(), result.clusters.end())),
        _["n_variable_genes"] = static_cast<double>(result.variable_genes.size()),
        _["variable_genes"] = wrap(result.variable_genes),
        _["pca_variance_explained"] = wrap(result.pca_variance_explained));
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_cluster_run_counts: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_fit_store_run(
    SEXP store_dir_sexp,
    SEXP analysis_crop_sexp,
    SEXP ncv_k_sexp,
    SEXP rank_sexp,
    SEXP graph_k_sexp,
    SEXP same_label_ratio_sexp,
    SEXP nmf_iterations_sexp,
    SEXP nmf_init_sexp,
    SEXP nmf_variant_sexp,
    SEXP molecule_scoring_sexp,
    SEXP nmf_n_runs_sexp,
    SEXP nmf_train_max_rows_sexp,
    SEXP nmf_min_molecules_sexp,
    SEXP num_threads_sexp,
    SEXP training_labels_path_sexp,
    SEXP use_cell_type_training_sexp,
    SEXP training_scope_cell_types_sexp,
    SEXP annotation_hash_sexp,
    SEXP seed_sexp,
    SEXP out_dir_sexp,
    SEXP tile_size_sexp,
    SEXP row_group_size_sexp,
    SEXP report_ncv_umap_sexp,
    SEXP verbose_sexp,
    SEXP nmf_fixed_h_sexp,
    SEXP nmf_fixed_h_genes_sexp) {
  try {
    const auto fit_start = std::chrono::steady_clock::now();
    auto stage_start = fit_start;
    const bool verbose = as<bool>(verbose_sexp);
    auto log_info = [&](const std::string& message, double stage_sec = -1.0) {
      if (!verbose) {
        return;
      }
      const auto now = std::chrono::system_clock::now();
      const auto tt = std::chrono::system_clock::to_time_t(now);
      std::tm tm{};
#ifdef _WIN32
      localtime_s(&tm, &tt);
#else
      localtime_r(&tt, &tm);
#endif
      std::ostringstream line;
      line << "[INFO " << std::put_time(&tm, "%H:%M:%S") << " +"
           << std::fixed << std::setprecision(3)
           << std::chrono::duration<double>(std::chrono::steady_clock::now() - fit_start).count()
           << "s] " << message;
      if (stage_sec >= 0.0) {
        line << " (" << std::fixed << std::setprecision(3) << stage_sec << "s)";
      }
      Rcpp::Rcout << line.str() << std::endl;
    };

    log_info("Starting store-backed fit");
    const std::string store_dir = as<std::string>(store_dir_sexp);
    const auto analysis_crop = optional_string_sexp(analysis_crop_sexp);
    const auto training_labels_path = optional_string_sexp(training_labels_path_sexp);

    celladmix::BasicPipelineOptions pipeline_options;
    pipeline_options.ncv_k = as<int>(ncv_k_sexp);
    pipeline_options.rank = as<int>(rank_sexp);
    pipeline_options.graph_k = as<int>(graph_k_sexp);
    pipeline_options.same_label_ratio = as<double>(same_label_ratio_sexp);
    pipeline_options.nmf_iterations = as<int>(nmf_iterations_sexp);
    pipeline_options.nmf_init = as<std::string>(nmf_init_sexp);
    pipeline_options.nmf_variant = as<std::string>(nmf_variant_sexp);
    pipeline_options.molecule_scoring = as<std::string>(molecule_scoring_sexp);
    pipeline_options.nmf_n_runs = as<int>(nmf_n_runs_sexp);
    pipeline_options.nmf_train_max_rows = optional_int_sexp(nmf_train_max_rows_sexp, -1);
    pipeline_options.nmf_min_molecules = as<int>(nmf_min_molecules_sexp);
    pipeline_options.num_threads = as<int>(num_threads_sexp);
    pipeline_options.return_ncv = false;
    pipeline_options.seed = optional_seed_sexp(seed_sexp, 1U);
    pipeline_options.training_scope_cell_types =
        optional_string_vector_sexp(training_scope_cell_types_sexp);
    pipeline_options.annotation_hash =
        optional_string_sexp(annotation_hash_sexp).value_or("");

    celladmix::RunStorageOptions storage_options;
    storage_options.tile_size = as<double>(tile_size_sexp);
    storage_options.parquet_row_group_size = as<int>(row_group_size_sexp);

    const auto store_manifest = celladmix::read_input_store_manifest(store_dir);
    if (!store_manifest.has_molecule_rows || !store_manifest.has_cell_offsets) {
      throw std::runtime_error("store-backed fit requires an input store with molecule rows and cell offsets");
    }
    const auto counts = celladmix::load_input_store_counts(store_dir);
    if (nmf_fixed_h_sexp != R_NilValue) {
      const Rcpp::NumericMatrix h_in(nmf_fixed_h_sexp);
      const auto h_genes = optional_string_vector_sexp(nmf_fixed_h_genes_sexp);
      if (static_cast<int>(h_genes.size()) != h_in.ncol()) {
        throw std::runtime_error("nmf_fixed_h_genes must name every column of nmf_fixed_h");
      }
      std::unordered_map<std::string, int> gene_index;
      for (std::size_t g = 0; g < counts.genes.size(); ++g) {
        gene_index[counts.genes[g]] = static_cast<int>(g);
      }
      const int rank_in = h_in.nrow();
      pipeline_options.nmf_fixed_h.assign(
          static_cast<std::size_t>(rank_in) * counts.genes.size(), 0.0);
      int matched = 0;
      for (int c = 0; c < h_in.ncol(); ++c) {
        const auto hit = gene_index.find(h_genes[static_cast<std::size_t>(c)]);
        if (hit == gene_index.end()) {
          continue;
        }
        ++matched;
        for (int f = 0; f < rank_in; ++f) {
          pipeline_options.nmf_fixed_h[
              static_cast<std::size_t>(f) * counts.genes.size() +
              static_cast<std::size_t>(hit->second)] = h_in(f, c);
        }
      }
      if (matched < 2) {
        throw std::runtime_error("nmf_fixed_h genes do not match the input store gene vocabulary");
      }
    }
    stage_start = std::chrono::steady_clock::now();
    if (training_labels_path.has_value()) {
      pipeline_options.training_cell_strata =
          load_cell_label_strata_for_cells(*training_labels_path, counts.cells.cell_ids);
    } else if (as<bool>(use_cell_type_training_sexp)) {
      pipeline_options.training_cell_strata =
          training_cell_strata_from_cell_metadata(counts.cells, counts.crop_ids);
    } else {
      pipeline_options.training_cell_strata.assign(counts.cells.size(), "__all_cells__");
    }
    log_info(
        "Resolved store training strata for " + std::to_string(counts.cells.size()) + " cells",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    celladmix::RunSourceInfo source;
    source.type = store_manifest.source_type.empty() ? "input_store" : store_manifest.source_type;
    source.path = store_manifest.source_path.empty() ? store_dir : store_manifest.source_path;
    source.used_parquet = true;
    source.files = build_input_store_source_files(store_dir, store_manifest);

    const auto fit = celladmix::run_basic_pipeline_store(
        store_dir,
        source,
        pipeline_options,
        storage_options,
        as<std::string>(out_dir_sexp),
        analysis_crop,
        as<bool>(report_ncv_umap_sexp),
        &counts,
        verbose);

    auto out = run_manifest_to_r(fit.manifest);
    out["h"] = matrix_to_r(fit.nmf.h);
    out["genes"] = fit.manifest.genes;
    out["nmf_final_objective"] = fit.nmf.final_objective;
    out["nmf_selected_seed"] = static_cast<double>(fit.nmf.selected_seed);
    out["nmf_selected_run"] = static_cast<double>(fit.nmf.selected_run + 1);
    out["nmf_candidate_final_objectives"] = wrap(fit.nmf.candidate_final_objectives);
    out["nmf_candidate_best_match_correlations"] = wrap(fit.nmf.candidate_best_match_correlations);
    out["nmf_factor_stability"] = wrap(fit.nmf.selected_factor_stability);
    out["nmf_candidate_final_objective_mean"] = fit.nmf.candidate_final_objective_mean;
    out["nmf_candidate_final_objective_sd"] = fit.nmf.candidate_final_objective_sd;
    out["nmf_candidate_best_match_correlation_mean"] = fit.nmf.candidate_best_match_correlation_mean;
    out["nmf_stability_comparison_runs"] = fit.nmf.stability_comparison_runs;
    out["nmf_stable_factor_count"] = fit.nmf.stable_factor_count;
    out["nmf_stability_threshold"] = fit.nmf.stability_threshold;
    return out;
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_fit_store_run: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_read_run(SEXP path_sexp) {
  try {
    const auto manifest = celladmix::read_run_manifest(as<std::string>(path_sexp));
    auto out = run_manifest_to_r(manifest);
    out["h"] = load_factor_matrix_from_run(manifest);
    out["genes"] = manifest.genes;
    return out;
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_read_run: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_collect_run_transcripts(
    SEXP path_sexp,
    SEXP analysis_crop_sexp,
    SEXP bbox_sexp,
    SEXP sample_n_sexp,
    SEXP seed_sexp) {
  try {
    celladmix::RunLoadOptions options;
    options.crop_id = optional_string_sexp(analysis_crop_sexp);
    options.region = optional_bbox_sexp(bbox_sexp);
    options.sample_n = optional_int_sexp(sample_n_sexp, -1);
    options.seed = optional_seed_sexp(seed_sexp, 1U);
    const auto run_data = celladmix::load_run_data(as<std::string>(path_sexp), options);
    return transcripts_to_df_selected(
        run_data.transcripts,
        select_transcript_rows(run_data.transcripts.size()),
        std::nullopt,
        &run_data.labels,
        nullptr,
        &run_data.transcript_crop_ids,
        &run_data.factor_margin);
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_collect_run_transcripts: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_collect_run_cells(SEXP path_sexp) {
  try {
    const auto manifest = celladmix::read_run_manifest(as<std::string>(path_sexp));
    return parquet_table_to_df(manifest.paths.cells_parquet);
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_collect_run_cells: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_collect_correction_summary(SEXP path_sexp) {
  try {
    const auto summary_path =
        celladmix::correction_summary_parquet_path(as<std::string>(path_sexp));
    if (!std::filesystem::exists(summary_path)) {
      stop("No correction summary found for run: %s", as<std::string>(path_sexp));
    }
    return parquet_table_to_df(summary_path);
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_collect_correction_summary: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_build_report_data(
    SEXP path_sexp,
    SEXP build_ncv_umap_sexp,
    SEXP force_sexp,
    SEXP umap_neighbors_sexp,
    SEXP umap_epochs_sexp,
    SEXP normalization_scale_sexp,
    SEXP pca_dims_sexp,
    SEXP num_threads_sexp,
    SEXP seed_sexp,
    SEXP verbose_sexp) {
  try {
    const auto report_start = std::chrono::steady_clock::now();
    auto stage_start = report_start;
    const bool verbose = as<bool>(verbose_sexp);
    const int report_threads = std::max(1, as<int>(num_threads_sexp));
    auto log_info = [&](const std::string& message, double stage_sec = -1.0) {
      if (!verbose) {
        return;
      }
      const auto now = std::chrono::system_clock::now();
      const auto tt = std::chrono::system_clock::to_time_t(now);
      std::tm tm{};
#ifdef _WIN32
      localtime_s(&tm, &tt);
#else
      localtime_r(&tt, &tm);
#endif
      std::ostringstream line;
      line << "[INFO " << std::put_time(&tm, "%H:%M:%S") << " +"
           << std::fixed << std::setprecision(3)
           << std::chrono::duration<double>(std::chrono::steady_clock::now() - report_start).count()
           << "s] " << message;
      if (stage_sec >= 0.0) {
        line << " (" << std::fixed << std::setprecision(3) << stage_sec << "s)";
      }
      Rcpp::Rcout << line.str() << std::endl;
    };
    const std::string run_path = as<std::string>(path_sexp);
    List out;
    if (!as<bool>(build_ncv_umap_sexp)) {
      log_info("Skipped report data build");
      return out;
    }

    log_info("Starting report data build");
    stage_start = std::chrono::steady_clock::now();
    const auto manifest = celladmix::read_run_manifest(run_path);
    const std::string parquet_path =
        celladmix::training_molecule_umap_parquet_path(manifest.paths.root_dir);
    log_info(
        "Loaded run manifest",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    if (!as<bool>(force_sexp) &&
        std::filesystem::exists(std::filesystem::path(parquet_path))) {
      out["ncv_umap_parquet"] = parquet_path;
      log_info("Reused existing NCV UMAP report sidecar");
      return out;
    }

    stage_start = std::chrono::steady_clock::now();
    const auto run_data = celladmix::load_run_training_cell_data(run_path);
    log_info(
        "Loaded training-cell run data: " + std::to_string(run_data.transcripts.size()) +
            " transcripts across " + std::to_string(run_data.cells.size()) + " cells",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    celladmix::BasicPipelineResult fit;
    fit.training_query_indices = run_data.training_query_indices;
    fit.labels = run_data.labels;
    fit.nmf.h = load_factor_matrix_dense(manifest);

    stage_start = std::chrono::steady_clock::now();
    celladmix::NcvOptions ncv_options;
    ncv_options.k = manifest.pipeline_options.ncv_k + 1;
    ncv_options.include_self = true;
    ncv_options.within_cell = true;
    ncv_options.query_indices = run_data.training_query_indices;
    celladmix::NcvFeatureTransform transform;
    transform.mode = manifest.pipeline_options.nmf_variant;
    transform.gene_weights = manifest.nmf_gene_weights;
    transform.target_row_sum = manifest.nmf_transform_target_row_sum;
    std::string molecule_scoring = manifest.pipeline_options.molecule_scoring;
    if (molecule_scoring == "auto") {
      molecule_scoring = manifest.pipeline_options.nmf_variant == "ls_nmf"
          ? "gene_loadings"
          : "ncv_projection";
    }
    if (molecule_scoring == "gene_loadings") {
      fit.factor_scores = celladmix::project_gene_loadings_to_factors(
          run_data.transcripts,
          fit.nmf.h,
          ncv_options,
          report_threads);
    } else if (manifest.pipeline_options.nmf_variant == "ls_nmf") {
      fit.factor_scores = celladmix::project_ncv_to_factors(
          run_data.transcripts,
          fit.nmf.h,
          ncv_options,
          report_threads);
    } else {
      fit.factor_scores = celladmix::project_ncv_to_factors_kl(
          run_data.transcripts,
          fit.nmf.h,
          ncv_options,
          report_threads,
          10,
          1e-4,
          1e-10,
          transform);
    }
    log_info(
        "Projected NCVs for report data",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    stage_start = std::chrono::steady_clock::now();
    fit.nmf.w = fit.factor_scores;
    log_info(
        "Extracted sampled NCV factor scores",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    stage_start = std::chrono::steady_clock::now();
    celladmix::write_training_molecule_umap_report(
        manifest.paths.root_dir,
        run_data.transcripts,
        fit,
        manifest.pipeline_options.ncv_k,
        as<int>(umap_neighbors_sexp),
        as<int>(umap_epochs_sexp),
        optional_seed_sexp(seed_sexp, manifest.pipeline_options.seed),
        as<double>(normalization_scale_sexp),
        as<int>(pca_dims_sexp),
        report_threads,
        &run_data.training_obs_ids);
    log_info(
        "Built NCV UMAP report sidecar",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    out["ncv_umap_parquet"] = parquet_path;
    log_info("Finished report data build");
    return out;
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_build_report_data: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_collect_training_molecules(SEXP path_sexp) {
  try {
    const auto report = celladmix::load_training_molecule_umap_report(as<std::string>(path_sexp));
    List out = List::create(
        _["training_rank"] = wrap(report.training_rank),
        _["obs_id"] = wrap(report.obs_id),
        _["cell_id"] = wrap(report.cell_id),
        _["gene"] = wrap(report.gene),
        _["x"] = wrap(report.x),
        _["y"] = wrap(report.y),
        _["z"] = wrap(report.z),
        _["factor_label"] = wrap(report.factor_label),
        _["factor_margin"] = wrap(report.factor_margin),
        _["dominant_component"] = wrap(report.dominant_component),
        _["dominant_component_margin"] = wrap(report.dominant_component_margin),
        _["projected_component"] = wrap(report.projected_component),
        _["projected_component_margin"] = wrap(report.projected_component_margin),
        _["umap_x"] = wrap(report.umap_x),
        _["umap_y"] = wrap(report.umap_y));
    for (int component = 0; component < report.component_weights.cols(); ++component) {
      NumericVector values(report.component_weights.rows());
      for (int row = 0; row < report.component_weights.rows(); ++row) {
        values[row] = report.component_weights(row, component);
      }
      out.push_back(values, "component_" + std::to_string(component + 1));
    }
    for (int component = 0; component < report.projected_component_weights.cols(); ++component) {
      NumericVector values(report.projected_component_weights.rows());
      for (int row = 0; row < report.projected_component_weights.rows(); ++row) {
        values[row] = report.projected_component_weights(row, component);
      }
      out.push_back(values, "projected_component_" + std::to_string(component + 1));
    }
    out.attr("class") = "data.frame";
    out.attr("row.names") = IntegerVector::create(NA_INTEGER, -static_cast<int>(report.obs_id.size()));
    return out;
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_collect_training_molecules: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_collect_run_counts(
    SEXP path_sexp,
    SEXP analysis_crop_sexp,
    SEXP bbox_sexp) {
  try {
    const std::string path = as<std::string>(path_sexp);
    celladmix::RunLoadOptions options;
    options.crop_id = optional_string_sexp(analysis_crop_sexp);
    options.region = optional_bbox_sexp(bbox_sexp);
    const auto run_data = celladmix::load_run_data(path, options);
    const auto full_cells = celladmix::load_run_cells(path);
    const auto subset_counts = celladmix::build_cell_gene_counts(run_data.transcripts);
    celladmix::DenseMatrix counts(
      static_cast<int>(full_cells.size()),
      static_cast<int>(run_data.transcripts.genes.size()),
      0.0
    );
    std::unordered_map<std::string, int> full_cell_index;
    full_cell_index.reserve(full_cells.cell_ids.size());
    for (std::size_t i = 0; i < full_cells.cell_ids.size(); ++i) {
      full_cell_index.emplace(full_cells.cell_ids[i], static_cast<int>(i));
    }
    for (int cell = 0; cell < subset_counts.rows(); ++cell) {
      const auto it = full_cell_index.find(run_data.transcripts.cells[static_cast<std::size_t>(cell)]);
      if (it == full_cell_index.end()) {
        continue;
      }
      for (int gene = 0; gene < subset_counts.cols(); ++gene) {
        counts(it->second, gene) = subset_counts(cell, gene);
      }
    }
    return counts_to_r(counts, full_cells.cell_ids, run_data.transcripts.genes);
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_collect_run_counts: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_collect_run_counts_sparse(
    SEXP path_sexp,
    SEXP analysis_crop_sexp,
    SEXP bbox_sexp) {
  try {
    celladmix::RunLoadOptions options;
    options.crop_id = optional_string_sexp(analysis_crop_sexp);
    options.region = optional_bbox_sexp(bbox_sexp);
    const auto counts = celladmix::collect_run_counts(as<std::string>(path_sexp), options);
    return sparse_counts_to_r_list(counts);
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_collect_run_counts_sparse: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_bridge_scores_run(
    SEXP path_sexp,
    SEXP analysis_crop_sexp,
    SEXP bbox_sexp,
    SEXP cell_types_sexp,
    SEXP candidate_mode_sexp,
    SEXP candidate_k_sexp,
    SEXP cell_candidate_k_sexp,
    SEXP candidate_pairs_per_type_pair_sexp,
    SEXP crossing_k_sexp,
    SEXP min_type_pair_contacts_sexp,
    SEXP min_factor_molecules_sexp,
    SEXP min_pairs_sexp,
    SEXP max_cells_per_type_pair_sexp,
    SEXP null_iterations_sexp,
    SEXP null_max_iterations_sexp,
    SEXP cell_candidate_halo_sexp,
    SEXP null_step_sexp,
    SEXP null_pool_fraction_sexp,
    SEXP num_threads_sexp,
    SEXP seed_sexp,
    SEXP compute_null_sexp,
    SEXP verbose_sexp,
    SEXP ensemble_member_sexp) {
  try {
    const auto bridge_start = std::chrono::steady_clock::now();
    auto stage_start = bridge_start;
    const bool verbose = as<bool>(verbose_sexp);
    auto log_info = [&](const std::string& message, double stage_sec = -1.0) {
      if (!verbose) {
        return;
      }
      const auto now = std::chrono::system_clock::now();
      const auto tt = std::chrono::system_clock::to_time_t(now);
      std::tm tm{};
#ifdef _WIN32
      localtime_s(&tm, &tt);
#else
      localtime_r(&tt, &tm);
#endif
      std::ostringstream line;
      line << "[INFO " << std::put_time(&tm, "%H:%M:%S") << " +"
           << std::fixed << std::setprecision(3)
           << std::chrono::duration<double>(std::chrono::steady_clock::now() - bridge_start).count()
           << "s] " << message;
      if (stage_sec >= 0.0) {
        line << " (" << std::fixed << std::setprecision(3) << stage_sec << "s)";
      }
      Rcpp::Rcout << line.str() << std::endl;
    };

    celladmix::RunLoadOptions options;
    options.crop_id = optional_string_sexp(analysis_crop_sexp);
    options.region = optional_bbox_sexp(bbox_sexp);
    options.ensemble_member = as<int>(ensemble_member_sexp);
    auto run_data = celladmix::load_run_bridge_data(as<std::string>(path_sexp), options);
    apply_cell_type_overrides(run_data.cells, cell_types_sexp);
    log_info(
        "Loaded bridge run data: molecules=" + std::to_string(run_data.transcripts.size()) +
            ", cells=" + std::to_string(run_data.transcripts.num_cells()),
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());
    celladmix::BridgeTestOptions bridge_options;
    bridge_options.candidate_mode = as<std::string>(candidate_mode_sexp);
    bridge_options.candidate_k = as<int>(candidate_k_sexp);
    bridge_options.cell_candidate_k = as<int>(cell_candidate_k_sexp);
    bridge_options.candidate_pairs_per_type_pair = as<int>(candidate_pairs_per_type_pair_sexp);
    bridge_options.crossing_k = as<int>(crossing_k_sexp);
    bridge_options.min_type_pair_contacts = as<int>(min_type_pair_contacts_sexp);
    bridge_options.min_factor_molecules = as<int>(min_factor_molecules_sexp);
    bridge_options.min_pairs = as<int>(min_pairs_sexp);
    bridge_options.max_cells_per_type_pair = as<int>(max_cells_per_type_pair_sexp);
    bridge_options.null_iterations = as<int>(null_iterations_sexp);
    bridge_options.null_max_iterations = as<int>(null_max_iterations_sexp);
    bridge_options.cell_candidate_halo = as<double>(cell_candidate_halo_sexp);
    bridge_options.null_step = as<double>(null_step_sexp);
    bridge_options.null_pool_fraction = as<double>(null_pool_fraction_sexp);
    bridge_options.num_threads = as<int>(num_threads_sexp);
    bridge_options.seed = static_cast<unsigned int>(as<int>(seed_sexp));
    bridge_options.compute_null = as<bool>(compute_null_sexp);
    bridge_options.progress = log_info;

    stage_start = std::chrono::steady_clock::now();
    const auto result = celladmix::run_bridge_test(
        run_data.transcripts,
        run_data.cells,
        run_data.labels,
        bridge_options);
    log_info(
        "Ran bridge test: score_rows=" + std::to_string(result.pair_scores.size()) +
            ", summaries=" + std::to_string(result.summaries.size()),
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());
    auto suffix = options.crop_id.has_value() ? "_" + *options.crop_id : std::string();
    if (options.ensemble_member >= 0) {
      suffix += "_m" + std::to_string(options.ensemble_member);
    }
    const auto pair_score_path =
        std::filesystem::path(run_data.manifest.paths.scores_dir) /
        (std::string("bridge_pair_scores") + suffix + ".parquet");
    const auto summary_path =
        std::filesystem::path(run_data.manifest.paths.scores_dir) /
        (std::string("bridge_summary") + suffix + ".parquet");
    stage_start = std::chrono::steady_clock::now();
    write_bridge_pair_scores_parquet(pair_score_path.string(), run_data.transcripts, result);
    write_bridge_summary_parquet(summary_path.string(), result);
    log_info(
        "Wrote bridge outputs",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    stage_start = std::chrono::steady_clock::now();
    auto pair_scores = bridge_pair_scores_to_df(run_data.transcripts, result);
    auto summary = bridge_summary_to_df(result);
    log_info(
        "Materialized bridge results for R",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());
    pair_scores.attr("path") = pair_score_path.string();
    summary.attr("path") = summary_path.string();
    auto out = List::create(
        _["scores"] = pair_scores,
        _["summary"] = summary,
        _["cell_types"] = result.cell_types,
        _["paths"] = List::create(
            _["scores"] = pair_score_path.string(),
            _["summary"] = summary_path.string()));
    out.attr("class") = CharacterVector::create("celladmix_bridge_result", "list");
    return out;
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_bridge_scores_run: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_membrane_scores_run(
    SEXP path_sexp,
    SEXP analysis_crop_sexp,
    SEXP bbox_sexp,
    SEXP image_path_sexp,
    SEXP pixel_size_sexp,
    SEXP x_offset_sexp,
    SEXP y_offset_sexp,
    SEXP epsilon_sexp,
    SEXP cell_types_sexp,
    SEXP cell_candidate_k_sexp,
    SEXP candidate_pairs_per_type_pair_sexp,
    SEXP cell_candidate_halo_sexp,
    SEXP min_factor_molecules_sexp,
    SEXP min_pairs_sexp,
    SEXP max_cells_per_type_pair_sexp,
    SEXP control_distance_fraction_sexp,
    SEXP line_samples_sexp,
    SEXP num_threads_sexp,
    SEXP seed_sexp,
    SEXP verbose_sexp,
    SEXP ensemble_member_sexp) {
  try {
    const auto membrane_start = std::chrono::steady_clock::now();
    auto stage_start = membrane_start;
    const bool verbose = as<bool>(verbose_sexp);
    auto log_info = [&](const std::string& message, double stage_sec = -1.0) {
      if (!verbose) {
        return;
      }
      const auto now = std::chrono::system_clock::now();
      const auto tt = std::chrono::system_clock::to_time_t(now);
      std::tm tm{};
#ifdef _WIN32
      localtime_s(&tm, &tt);
#else
      localtime_r(&tt, &tm);
#endif
      std::ostringstream line;
      line << "[INFO " << std::put_time(&tm, "%H:%M:%S") << " +"
           << std::fixed << std::setprecision(3)
           << std::chrono::duration<double>(std::chrono::steady_clock::now() - membrane_start).count()
           << "s] " << message;
      if (stage_sec >= 0.0) {
        line << " (" << std::fixed << std::setprecision(3) << stage_sec << "s)";
      }
      Rcpp::Rcout << line.str() << std::endl;
    };

    celladmix::RunLoadOptions options;
    options.crop_id = optional_string_sexp(analysis_crop_sexp);
    options.region = optional_bbox_sexp(bbox_sexp);
    options.ensemble_member = as<int>(ensemble_member_sexp);
    auto run_data = celladmix::load_run_bridge_data(as<std::string>(path_sexp), options);
    apply_cell_type_overrides(run_data.cells, cell_types_sexp);
    log_info(
        "Loaded membrane run data: molecules=" + std::to_string(run_data.transcripts.size()) +
            ", cells=" + std::to_string(run_data.transcripts.num_cells()),
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    celladmix::MembraneImageOptions image_options;
    image_options.image_path = as<std::string>(image_path_sexp);
    image_options.pixel_size = as<double>(pixel_size_sexp);
    image_options.x_offset = as<double>(x_offset_sexp);
    image_options.y_offset = as<double>(y_offset_sexp);
    image_options.epsilon = as<double>(epsilon_sexp);

    celladmix::MembraneTestOptions membrane_options;
    membrane_options.cell_candidate_k = as<int>(cell_candidate_k_sexp);
    membrane_options.candidate_pairs_per_type_pair = as<int>(candidate_pairs_per_type_pair_sexp);
    membrane_options.cell_candidate_halo = as<double>(cell_candidate_halo_sexp);
    membrane_options.min_factor_molecules = as<int>(min_factor_molecules_sexp);
    membrane_options.min_pairs = as<int>(min_pairs_sexp);
    membrane_options.max_cells_per_type_pair = as<int>(max_cells_per_type_pair_sexp);
    membrane_options.control_distance_fraction = as<double>(control_distance_fraction_sexp);
    membrane_options.line_samples = as<int>(line_samples_sexp);
    membrane_options.num_threads = as<int>(num_threads_sexp);
    membrane_options.seed = static_cast<unsigned int>(as<int>(seed_sexp));
    membrane_options.progress = log_info;

    stage_start = std::chrono::steady_clock::now();
    const auto result = celladmix::run_membrane_test(
        run_data.transcripts,
        run_data.cells,
        run_data.labels,
        image_options,
        membrane_options);
    log_info(
        "Ran membrane test: score_rows=" + std::to_string(result.pair_scores.size()) +
            ", summaries=" + std::to_string(result.summaries.size()),
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    std::string suffix;
    if (options.crop_id.has_value()) {
      suffix = "_" + *options.crop_id;
    } else if (options.region.has_value()) {
      std::ostringstream key;
      key << std::fixed << std::setprecision(6)
          << options.region->xmin << "_" << options.region->xmax << "_"
          << options.region->ymin << "_" << options.region->ymax;
      if (options.region->has_z) {
        key << "_" << options.region->zmin << "_" << options.region->zmax;
      }
      std::ostringstream hashed;
      hashed << std::hex << std::hash<std::string>{}(key.str());
      suffix = "_bbox_" + hashed.str();
    }
    if (options.ensemble_member >= 0) {
      suffix += "_m" + std::to_string(options.ensemble_member);
    }
    const auto pair_score_path =
        std::filesystem::path(run_data.manifest.paths.scores_dir) /
        (std::string("membrane_pair_scores") + suffix + ".parquet");
    const auto summary_path =
        std::filesystem::path(run_data.manifest.paths.scores_dir) /
        (std::string("membrane_summary") + suffix + ".parquet");
    stage_start = std::chrono::steady_clock::now();
    write_membrane_pair_scores_parquet(pair_score_path.string(), run_data.transcripts, result);
    write_membrane_summary_parquet(summary_path.string(), result);
    log_info(
        "Wrote membrane outputs",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    stage_start = std::chrono::steady_clock::now();
    auto pair_scores = membrane_pair_scores_to_df(run_data.transcripts, result);
    auto summary = membrane_summary_to_df(result);
    log_info(
        "Materialized membrane results for R",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());
    pair_scores.attr("path") = pair_score_path.string();
    summary.attr("path") = summary_path.string();
    auto out = List::create(
        _["scores"] = pair_scores,
        _["summary"] = summary,
        _["cell_types"] = result.cell_types,
        _["image"] = List::create(
            _["path"] = image_options.image_path,
            _["pixel_size"] = image_options.pixel_size,
            _["x_offset"] = image_options.x_offset,
            _["y_offset"] = image_options.y_offset),
        _["paths"] = List::create(
            _["scores"] = pair_score_path.string(),
            _["summary"] = summary_path.string()));
    out.attr("class") = CharacterVector::create("celladmix_membrane_result", "list");
    return out;
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_membrane_scores_run: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_coherence_scores_run(
    SEXP path_sexp,
    SEXP analysis_crop_sexp,
    SEXP bbox_sexp,
    SEXP cell_types_sexp,
    SEXP stain_mode_sexp,
    SEXP image_path_sexp,
    SEXP pixel_size_sexp,
    SEXP x_offset_sexp,
    SEXP y_offset_sexp,
    SEXP stain_max_pixels_sexp,
    SEXP k_neighbors_sexp,
    SEXP min_factor_molecules_sexp,
    SEXP min_cells_sexp,
    SEXP max_neighbor_distance_sexp,
    SEXP distance_sigma_sexp,
    SEXP lambda_coherence_sexp,
    SEXP beta_margin_sexp,
    SEXP score_threshold_sexp,
    SEXP source_pseudocount_sexp,
    SEXP min_source_log_enrichment_sexp,
    SEXP include_self_source_sexp,
    SEXP compute_null_sexp,
    SEXP null_method_sexp,
    SEXP null_iterations_sexp,
    SEXP null_exclude_factor_sexp,
    SEXP null_match_nucleus_sexp,
    SEXP null_match_density_sexp,
    SEXP null_nucleus_distance_bins_sexp,
    SEXP null_density_bins_sexp,
    SEXP seed_sexp,
    SEXP normalize_membrane_sexp,
    SEXP membrane_low_quantile_sexp,
    SEXP membrane_high_quantile_sexp,
    SEXP membrane_alpha_sexp,
    SEXP line_samples_sexp,
    SEXP patch_edge_weight_min_sexp,
    SEXP num_threads_sexp,
    SEXP verbose_sexp,
    SEXP ensemble_member_sexp) {
  try {
    const auto coherence_start = std::chrono::steady_clock::now();
    auto stage_start = coherence_start;
    const bool verbose = as<bool>(verbose_sexp);
    auto log_info = [&](const std::string& message, double stage_sec = -1.0) {
      if (!verbose) {
        return;
      }
      const auto now = std::chrono::system_clock::now();
      const auto tt = std::chrono::system_clock::to_time_t(now);
      std::tm tm{};
#ifdef _WIN32
      localtime_s(&tm, &tt);
#else
      localtime_r(&tt, &tm);
#endif
      std::ostringstream line;
      line << "[INFO " << std::put_time(&tm, "%H:%M:%S") << " +"
           << std::fixed << std::setprecision(3)
           << std::chrono::duration<double>(std::chrono::steady_clock::now() - coherence_start).count()
           << "s] " << message;
      if (stage_sec >= 0.0) {
        line << " (" << std::fixed << std::setprecision(3) << stage_sec << "s)";
      }
      Rcpp::Rcout << line.str() << std::endl;
    };

    celladmix::RunLoadOptions load_options;
    load_options.crop_id = optional_string_sexp(analysis_crop_sexp);
    load_options.region = optional_bbox_sexp(bbox_sexp);
    load_options.ensemble_member = as<int>(ensemble_member_sexp);
    auto run_data = celladmix::load_run_data(as<std::string>(path_sexp), load_options);
    apply_cell_type_overrides(run_data.cells, cell_types_sexp);
    log_info(
        "Loaded coherence run data: molecules=" + std::to_string(run_data.transcripts.size()) +
            ", cells=" + std::to_string(run_data.transcripts.num_cells()),
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    const std::string stain_mode = as<std::string>(stain_mode_sexp);
    std::optional<celladmix::Image2D> membrane_image;
    if (stain_mode == "membrane_barrier") {
      const std::string image_path = as<std::string>(image_path_sexp);
      if (image_path.empty()) {
        stop("image_path is required when stain_mode = 'membrane_barrier'");
      }
      if (run_data.transcripts.size() == 0) {
        stop("No molecules available for image-bounded coherence scoring");
      }
      double xmin = *std::min_element(run_data.transcripts.x.begin(), run_data.transcripts.x.end());
      double xmax = *std::max_element(run_data.transcripts.x.begin(), run_data.transcripts.x.end());
      double ymin = *std::min_element(run_data.transcripts.y.begin(), run_data.transcripts.y.end());
      double ymax = *std::max_element(run_data.transcripts.y.begin(), run_data.transcripts.y.end());
      const double pad = std::max(1.0, as<double>(pixel_size_sexp));
      xmin -= pad;
      xmax += pad;
      ymin -= pad;
      ymax += pad;
      celladmix::MembraneImageOptions image_options;
      image_options.image_path = image_path;
      image_options.pixel_size = as<double>(pixel_size_sexp);
      image_options.x_offset = as<double>(x_offset_sexp);
      image_options.y_offset = as<double>(y_offset_sexp);
      stage_start = std::chrono::steady_clock::now();
      membrane_image = celladmix::read_stain_image_crop(
          image_options,
          xmin,
          xmax,
          ymin,
          ymax,
          as<int>(stain_max_pixels_sexp));
      log_info(
          "Read membrane barrier image crop: width=" + std::to_string(membrane_image->width) +
              ", height=" + std::to_string(membrane_image->height),
          std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());
    } else if (stain_mode != "none") {
      stop("stain_mode must be 'none' or 'membrane_barrier'");
    }

    celladmix::CoherenceTestOptions coherence_options;
    coherence_options.k_neighbors = as<int>(k_neighbors_sexp);
    coherence_options.min_factor_molecules = as<int>(min_factor_molecules_sexp);
    coherence_options.min_cells = as<int>(min_cells_sexp);
    coherence_options.max_neighbor_distance = as<double>(max_neighbor_distance_sexp);
    coherence_options.distance_sigma = as<double>(distance_sigma_sexp);
    coherence_options.lambda_coherence = as<double>(lambda_coherence_sexp);
    coherence_options.beta_margin = as<double>(beta_margin_sexp);
    coherence_options.score_threshold = as<double>(score_threshold_sexp);
    coherence_options.source_pseudocount = as<double>(source_pseudocount_sexp);
    coherence_options.min_source_log_enrichment = as<double>(min_source_log_enrichment_sexp);
    coherence_options.include_self_source = as<bool>(include_self_source_sexp);
    coherence_options.compute_null = as<bool>(compute_null_sexp);
    coherence_options.null_method = as<std::string>(null_method_sexp);
    coherence_options.null_iterations = as<int>(null_iterations_sexp);
    coherence_options.null_exclude_factor = as<bool>(null_exclude_factor_sexp);
    coherence_options.null_match_nucleus = as<bool>(null_match_nucleus_sexp);
    coherence_options.null_match_density = as<bool>(null_match_density_sexp);
    coherence_options.null_nucleus_distance_bins = as<int>(null_nucleus_distance_bins_sexp);
    coherence_options.null_density_bins = as<int>(null_density_bins_sexp);
    coherence_options.seed = static_cast<unsigned int>(as<int>(seed_sexp));
    coherence_options.use_membrane_barrier = stain_mode == "membrane_barrier";
    coherence_options.normalize_membrane = as<bool>(normalize_membrane_sexp);
    coherence_options.membrane_low_quantile = as<double>(membrane_low_quantile_sexp);
    coherence_options.membrane_high_quantile = as<double>(membrane_high_quantile_sexp);
    coherence_options.membrane_alpha = as<double>(membrane_alpha_sexp);
    coherence_options.line_samples = as<int>(line_samples_sexp);
    coherence_options.patch_edge_weight_min = as<double>(patch_edge_weight_min_sexp);
    coherence_options.num_threads = as<int>(num_threads_sexp);
    coherence_options.progress = log_info;

    stage_start = std::chrono::steady_clock::now();
    const auto result = celladmix::run_coherence_test(
        run_data.transcripts,
        run_data.cells,
        run_data.labels,
        run_data.factor_margin,
        membrane_image.has_value() ? &(*membrane_image) : nullptr,
        coherence_options);
    log_info(
        "Ran coherence test: score_rows=" + std::to_string(result.cell_scores.size()) +
            ", summaries=" + std::to_string(result.summaries.size()),
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    std::string suffix;
    if (stain_mode != "none") {
      suffix += "_" + stain_mode;
    }
    if (load_options.crop_id.has_value()) {
      suffix += "_" + *load_options.crop_id;
    } else if (load_options.region.has_value()) {
      std::ostringstream key;
      key << std::fixed << std::setprecision(6)
          << load_options.region->xmin << "_" << load_options.region->xmax << "_"
          << load_options.region->ymin << "_" << load_options.region->ymax;
      if (load_options.region->has_z) {
        key << "_" << load_options.region->zmin << "_" << load_options.region->zmax;
      }
      std::ostringstream hashed;
      hashed << std::hex << std::hash<std::string>{}(key.str());
      suffix += "_bbox_" + hashed.str();
    }
    if (load_options.ensemble_member >= 0) {
      suffix += "_m" + std::to_string(load_options.ensemble_member);
    }
    const auto score_path =
        std::filesystem::path(run_data.manifest.paths.scores_dir) /
        (std::string("coherence_cell_scores") + suffix + ".parquet");
    const auto summary_path =
        std::filesystem::path(run_data.manifest.paths.scores_dir) /
        (std::string("coherence_summary") + suffix + ".parquet");
    stage_start = std::chrono::steady_clock::now();
    write_coherence_cell_scores_parquet(score_path.string(), run_data.transcripts, result);
    write_coherence_summary_parquet(summary_path.string(), result);
    log_info(
        "Wrote coherence outputs",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());

    stage_start = std::chrono::steady_clock::now();
    auto scores = coherence_cell_scores_to_df(run_data.transcripts, result);
    auto summary = coherence_summary_to_df(result);
    log_info(
        "Materialized coherence results for R",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());
    scores.attr("path") = score_path.string();
    summary.attr("path") = summary_path.string();
    auto out = List::create(
        _["scores"] = scores,
        _["summary"] = summary,
        _["cell_types"] = result.cell_types,
        _["stain_mode"] = stain_mode,
        _["membrane_normalized"] = result.membrane_normalized,
        _["membrane_scale_low"] = result.membrane_scale_low,
        _["membrane_scale_high"] = result.membrane_scale_high,
        _["null_matching"] = List::create(
            _["method"] = result.null_method,
            _["nucleus_overlap"] = result.null_match_nucleus_overlap_used,
            _["nucleus_distance"] = result.null_match_nucleus_distance_used,
            _["density"] = result.null_match_density_used,
            _["nucleus_distance_bins"] = result.null_nucleus_distance_bins,
            _["density_bins"] = result.null_density_bins),
        _["paths"] = List::create(
            _["scores"] = score_path.string(),
            _["summary"] = summary_path.string()));
    out.attr("class") = CharacterVector::create("celladmix_coherence_result", "list");
    return out;
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_coherence_scores_run: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_read_stain_image_crop(
    SEXP image_path_sexp,
    SEXP pixel_size_sexp,
    SEXP x_offset_sexp,
    SEXP y_offset_sexp,
    SEXP bbox_sexp,
    SEXP max_pixels_sexp) {
  try {
    auto bbox = optional_bbox_sexp(bbox_sexp);
    if (!bbox.has_value()) {
      stop("bbox is required");
    }
    celladmix::MembraneImageOptions image_options;
    image_options.image_path = as<std::string>(image_path_sexp);
    image_options.pixel_size = as<double>(pixel_size_sexp);
    image_options.x_offset = as<double>(x_offset_sexp);
    image_options.y_offset = as<double>(y_offset_sexp);
    const auto image = celladmix::read_stain_image_crop(
        image_options,
        bbox->xmin,
        bbox->xmax,
        bbox->ymin,
        bbox->ymax,
        as<int>(max_pixels_sexp));
    NumericMatrix values(image.height, image.width);
    for (int row = 0; row < image.height; ++row) {
      for (int col = 0; col < image.width; ++col) {
        values(row, col) =
            image.values[static_cast<std::size_t>(row) * static_cast<std::size_t>(image.width) +
                         static_cast<std::size_t>(col)];
      }
    }
    return List::create(
        _["values"] = values,
        _["xlim"] = NumericVector::create(bbox->xmin, bbox->xmax),
        _["ylim"] = NumericVector::create(bbox->ymin, bbox->ymax),
        _["pixel_size"] = image.pixel_size,
        _["width"] = image.width,
        _["height"] = image.height,
        _["path"] = image_options.image_path);
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_read_stain_image_crop: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_cell_neighbor_type_counts(
    SEXP x_sexp,
    SEXP y_sexp,
    SEXP type_codes_sexp,
    SEXP n_types_sexp,
    SEXP k_sexp) {
  try {
    const std::vector<double> x = as<std::vector<double>>(x_sexp);
    const std::vector<double> y = as<std::vector<double>>(y_sexp);
    const std::vector<int> type_codes = as<std::vector<int>>(type_codes_sexp);
    const auto counts = celladmix::cell_neighbor_type_counts(
        x, y, type_codes, as<int>(n_types_sexp), as<int>(k_sexp));
    return matrix_to_r(counts);
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_cell_neighbor_type_counts: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_fit_generative(
    SEXP counts_indptr_sexp,
    SEXP counts_indices_sexp,
    SEXP counts_values_sexp,
    SEXP n_genes_sexp,
    SEXP cell_ids_sexp,
    SEXP x_sexp,
    SEXP y_sexp,
    SEXP type_codes_sexp,
    SEXP n_types_sexp,
    SEXP molecules_parquet_sexp,
    SEXP cells_parquet_sexp,
    SEXP pairs_sexp,
    SEXP factor_to_type_sexp,
    SEXP programs_flat_sexp,
    SEXP program_type_sexp,
    SEXP options_sexp) {
  try {
    std::vector<celladmix::GenerativePairSpec> pairs;
    const List pair_list(pairs_sexp);
    for (R_xlen_t i = 0; i < pair_list.size(); ++i) {
      const List d(pair_list[i]);
      celladmix::GenerativePairSpec ps;
      ps.source_type = as<int>(d["source_type"]);
      ps.target_type = as<int>(d["target_type"]);
      ps.pool = as<std::vector<int>>(d["pool"]);
      ps.strict = as<std::vector<int>>(d["strict"]);
      if (d.containsElementNamed("guide")) {
        ps.guide = as<std::vector<int>>(d["guide"]);
      }
      if (d.containsElementNamed("exposure")) {
        ps.exposure = as<std::vector<double>>(d["exposure"]);
      }
      pairs.push_back(std::move(ps));
    }
    celladmix::GenerativeOptions opt;
    const List options(options_sexp);
    const auto set_double = [&](const char* key, double& field) {
      if (options.containsElementNamed(key)) field = as<double>(options[key]);
    };
    const auto set_int = [&](const char* key, int& field) {
      if (options.containsElementNamed(key)) field = as<int>(options[key]);
    };
    const auto set_bool = [&](const char* key, bool& field) {
      if (options.containsElementNamed(key)) field = as<bool>(options[key]);
    };
    set_double("alpha_prior_strength", opt.alpha_prior_strength);
    set_double("alpha_cap", opt.alpha_cap);
    set_double("lambda_max", opt.lambda_max);
    set_double("ambient_prior_strength", opt.ambient_prior_strength);
    set_double("ambient_cap", opt.ambient_cap);
    set_double("floor_total", opt.floor_total);
    set_double("induced_z", opt.induced_z);
    set_double("induced_min_excess", opt.induced_min_excess);
    set_double("profile_cv", opt.profile_cv);
    set_double("rho_shape", opt.rho_shape);
    set_double("rho_cap", opt.rho_cap);
    set_double("near_um", opt.near_um);
    set_double("near_min_molecules", opt.near_min_molecules);
    set_double("dose_weight", opt.dose_weight);
    set_double("far_um", opt.far_um);
    set_double("far_min_molecules", opt.far_min_molecules);
    set_int("em_iterations", opt.em_iterations);
    set_int("topup_em_iterations", opt.topup_em_iterations);
    set_int("topup_passes", opt.topup_passes);
    set_int("outer_rounds", opt.outer_rounds);
    set_int("neighbor_k", opt.neighbor_k);
    set_bool("use_ambient", opt.use_ambient);
    set_bool("use_induced", opt.use_induced);
    set_int("num_threads", opt.num_threads);
    set_int("n_programs", opt.n_programs);
    if (options.containsElementNamed("init_mode")) {
      opt.init_mode = as<std::string>(options["init_mode"]);
    }
    if (options.containsElementNamed("seed")) {
      opt.seed = static_cast<unsigned int>(as<int>(options["seed"]));
    }

    const auto res = celladmix::fit_generative(
        as<std::vector<int>>(counts_indptr_sexp),
        as<std::vector<int>>(counts_indices_sexp),
        as<std::vector<double>>(counts_values_sexp),
        as<int>(n_genes_sexp),
        as<std::vector<std::string>>(cell_ids_sexp),
        as<std::vector<double>>(x_sexp),
        as<std::vector<double>>(y_sexp),
        as<std::vector<int>>(type_codes_sexp),
        as<int>(n_types_sexp),
        as<std::string>(molecules_parquet_sexp),
        as<std::string>(cells_parquet_sexp),
        pairs,
        as<std::vector<int>>(factor_to_type_sexp),
        as<std::vector<double>>(programs_flat_sexp),
        as<std::vector<int>>(program_type_sexp),
        opt);

    List pair_cells(res.pair_cells.size());
    List pair_dose(res.pair_dose.size());
    List pair_alpha(res.pair_alpha.size());
    List pair_rho(res.pair_rho.size());
    List pair_induced(res.pair_induced.size());
    for (std::size_t j = 0; j < res.pair_cells.size(); ++j) {
      pair_cells[j] = wrap(res.pair_cells[j]);
      pair_dose[j] = wrap(res.pair_dose[j]);
      pair_alpha[j] = wrap(res.pair_alpha[j]);
      pair_rho[j] = wrap(res.pair_rho[j]);
      pair_induced[j] = wrap(res.pair_induced[j]);
    }
    std::vector<int> ind_pair;
    std::vector<int> ind_gene;
    std::vector<double> ind_excess;
    std::vector<double> ind_expected;
    std::vector<double> ind_z;
    for (const auto& row : res.induced) {
      ind_pair.push_back(row.pair);
      ind_gene.push_back(row.gene);
      ind_excess.push_back(row.excess);
      ind_expected.push_back(row.expected);
      ind_z.push_back(row.z);
    }
    std::vector<int> sum_pair;
    std::vector<double> sum_prior;
    std::vector<double> sum_post;
    std::vector<double> sum_ind;
    std::vector<double> sum_dose;
    for (const auto& row : res.pairs) {
      sum_pair.push_back(row.pair);
      sum_prior.push_back(row.prior_molecules);
      sum_post.push_back(row.posterior_molecules);
      sum_ind.push_back(row.induced_molecules);
      sum_dose.push_back(row.mean_dose_exposed);
    }
    return List::create(
        Named("removed") = wrap(res.removed),
        Named("removed_without_retention") = wrap(res.removed_without_retention),
        Named("n_programs_used") = wrap(res.n_programs_used),
        Named("pair_cells") = pair_cells,
        Named("pair_dose") = pair_dose,
        Named("pair_alpha") = pair_alpha,
        Named("pair_rho") = pair_rho,
        Named("pair_induced") = pair_induced,
        Named("ambient_scale") = wrap(res.ambient_scale),
        Named("factor_alignment") = wrap(res.factor_alignment),
        Named("whole_cell_profiles") = res.whole_cell_profiles,
        Named("induced") = List::create(
            Named("pair") = wrap(ind_pair),
            Named("gene") = wrap(ind_gene),
            Named("excess") = wrap(ind_excess),
            Named("expected") = wrap(ind_expected),
            Named("z") = wrap(ind_z)),
        Named("pairs") = List::create(
            Named("pair") = wrap(sum_pair),
            Named("prior_molecules") = wrap(sum_prior),
            Named("posterior_molecules") = wrap(sum_post),
            Named("induced_molecules") = wrap(sum_ind),
            Named("mean_dose_exposed") = wrap(sum_dose)));
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_fit_generative: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_cell_nearest_type_distance(
    SEXP x_sexp,
    SEXP y_sexp,
    SEXP type_codes_sexp,
    SEXP n_types_sexp) {
  try {
    const std::vector<double> x = as<std::vector<double>>(x_sexp);
    const std::vector<double> y = as<std::vector<double>>(y_sexp);
    const std::vector<int> type_codes = as<std::vector<int>>(type_codes_sexp);
    const auto dist = celladmix::cell_nearest_type_distance(
        x, y, type_codes, as<int>(n_types_sexp));
    return matrix_to_r(dist);
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_cell_nearest_type_distance: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_write_cell_labels(
    SEXP path_sexp,
    SEXP labels_sexp,
    SEXP row_group_size_sexp) {
  try {
    CharacterVector labels(labels_sexp);
    CharacterVector names(labels.attr("names"));
    if (names.size() != labels.size()) {
      stop("labels must be a named character vector with cell IDs as names");
    }
    std::vector<std::string> cell_id;
    std::vector<std::string> label;
    cell_id.reserve(static_cast<std::size_t>(labels.size()));
    label.reserve(static_cast<std::size_t>(labels.size()));
    for (R_xlen_t i = 0; i < labels.size(); ++i) {
      if (CharacterVector::is_na(names[i]) || CharacterVector::is_na(labels[i])) {
        continue;
      }
      const std::string cell = as<std::string>(names[i]);
      const std::string value = as<std::string>(labels[i]);
      if (!cell.empty() && !value.empty()) {
        cell_id.push_back(cell);
        label.push_back(value);
      }
    }
    if (cell_id.empty()) {
      stop("labels did not contain any usable named entries");
    }
    write_cell_labels_parquet(
        as<std::string>(path_sexp),
        cell_id,
        label,
        as<int>(row_group_size_sexp));
    return List::create(
        _["path"] = as<std::string>(path_sexp),
        _["n_labels"] = static_cast<double>(label.size()));
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_write_cell_labels: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_correct_run(
    SEXP path_sexp,
    SEXP rules_sexp,
    SEXP out_dir_sexp,
    SEXP cell_types_sexp,
    SEXP rule_member_sexp,
    SEXP min_votes_sexp) {
  try {
    auto run_data = celladmix::load_run_data(as<std::string>(path_sexp));
    const DataFrame rules = as<DataFrame>(rules_sexp);
    if (!rules.containsElementNamed("factor") || !rules.containsElementNamed("target_cell_type")) {
      stop("rules must contain factor and target_cell_type columns");
    }
    const auto factor = as<std::vector<int>>(rules["factor"]);
    const auto target_cell_type = as<std::vector<std::string>>(rules["target_cell_type"]);

    // Ensemble voting: each rule belongs to one member labeling (-1 = the
    // run's own labels); a molecule is removed when at least min_votes
    // member labelings remove it.
    std::vector<int> rule_member(factor.size(), -1);
    if (!Rf_isNull(rule_member_sexp)) {
      const IntegerVector member_vec(rule_member_sexp);
      if (static_cast<std::size_t>(member_vec.size()) == factor.size()) {
        for (R_xlen_t i = 0; i < member_vec.size(); ++i) {
          rule_member[static_cast<std::size_t>(i)] =
              member_vec[i] == NA_INTEGER ? -1 : member_vec[i];
        }
      } else if (member_vec.size() != 0) {
        stop("rule member vector must match the rule count");
      }
    }
    const int min_votes = Rf_isNull(min_votes_sexp)
        ? 1
        : std::max(1, as<int>(min_votes_sexp));

    std::unordered_map<std::string, std::string> cell_type_override;
    if (!Rf_isNull(cell_types_sexp)) {
      const CharacterVector labels(cell_types_sexp);
      const SEXP names_sexp = Rf_getAttrib(cell_types_sexp, R_NamesSymbol);
      if (Rf_isNull(names_sexp)) {
        stop("cell_types override must be a named character vector");
      }
      const CharacterVector names(names_sexp);
      if (names.size() != labels.size()) {
        stop("cell_types override names must match labels");
      }
      cell_type_override.reserve(static_cast<std::size_t>(labels.size()));
      for (R_xlen_t i = 0; i < labels.size(); ++i) {
        if (CharacterVector::is_na(names[i]) || CharacterVector::is_na(labels[i])) {
          continue;
        }
        const std::string cell = as<std::string>(names[i]);
        const std::string label = as<std::string>(labels[i]);
        if (!cell.empty() && !label.empty()) {
          cell_type_override.emplace(cell, label);
        }
      }
    }

    std::vector<std::string> transcript_cell_types;
    if (!cell_type_override.empty()) {
      transcript_cell_types.assign(run_data.transcripts.size(), "");
      for (std::size_t i = 0; i < run_data.transcripts.size(); ++i) {
        std::string cell = run_data.transcripts.orig_cell_id[i];
        if (!run_data.transcripts.cell_index.empty() &&
            run_data.transcripts.cell_index[i] >= 0 &&
            static_cast<std::size_t>(run_data.transcripts.cell_index[i]) <
                run_data.transcripts.cells.size()) {
          cell = run_data.transcripts.cells[
              static_cast<std::size_t>(run_data.transcripts.cell_index[i])];
        }
        const auto found = cell_type_override.find(cell);
        if (found != cell_type_override.end()) {
          transcript_cell_types[i] = found->second;
        }
      }

      if (!run_data.cells.cell_ids.empty()) {
        run_data.cells.cell_types.assign(run_data.cells.cell_ids.size(), "");
        for (std::size_t i = 0; i < run_data.cells.cell_ids.size(); ++i) {
          const auto found = cell_type_override.find(run_data.cells.cell_ids[i]);
          if (found != cell_type_override.end()) {
            run_data.cells.cell_types[i] = found->second;
          }
        }
      }
    }

    if (run_data.labels.size() != run_data.transcripts.size()) {
      stop("labels length must match TranscriptTable size");
    }
    std::map<int, std::vector<std::size_t>> rules_by_member;
    for (std::size_t rule = 0; rule < factor.size(); ++rule) {
      rules_by_member[rule_member[rule]].push_back(rule);
    }

    auto mark_rule_removed = [&](
        const std::vector<int>& labels,
        std::size_t rule,
        std::vector<bool>& removed) {
      if (!transcript_cell_types.empty()) {
        for (std::size_t i = 0; i < run_data.transcripts.size(); ++i) {
          if (labels[i] == factor[rule] - 1 &&
              transcript_cell_types[i] == target_cell_type[rule]) {
            removed[i] = true;
          }
        }
      } else {
        const auto rule_keep = celladmix::apply_removal_rule(
            run_data.transcripts,
            labels,
            factor[rule] - 1,
            target_cell_type[rule]);
        for (std::size_t i = 0; i < rule_keep.size(); ++i) {
          if (!rule_keep[i]) {
            removed[i] = true;
          }
        }
      }
    };

    std::vector<int> votes(run_data.transcripts.size(), 0);
    for (const auto& [member, rule_ids] : rules_by_member) {
      std::vector<int> member_label_vec;
      const std::vector<int>* labels_ptr = &run_data.labels;
      if (member >= 0) {
        const auto member_labels = celladmix::load_ensemble_member_labels(
            as<std::string>(path_sexp), member);
        member_label_vec.resize(run_data.transcripts.size());
        for (std::size_t i = 0; i < run_data.transcripts.size(); ++i) {
          member_label_vec[i] = member_labels.label_for(run_data.obs_ids[i]);
        }
        labels_ptr = &member_label_vec;
      }
      std::vector<bool> removed(run_data.transcripts.size(), false);
      for (const auto rule : rule_ids) {
        mark_rule_removed(*labels_ptr, rule, removed);
      }
      for (std::size_t i = 0; i < removed.size(); ++i) {
        if (removed[i]) {
          votes[i] += 1;
        }
      }
    }

    // min_votes is interpreted against the caller's intended member count:
    // members that contributed no rules simply never vote, so a threshold
    // above the voting-member count removes nothing.
    const int n_members = static_cast<int>(rules_by_member.size());
    std::vector<bool> keep_mask(run_data.transcripts.size(), true);
    for (std::size_t i = 0; i < keep_mask.size(); ++i) {
      keep_mask[i] = votes[i] < min_votes;
    }

    const auto manifest = celladmix::write_corrected_run(
        as<std::string>(out_dir_sexp),
        as<std::string>(path_sexp),
        run_data,
        keep_mask);
    int kept = 0;
    for (bool keep : keep_mask) {
      kept += keep ? 1 : 0;
    }
    IntegerVector vote_histogram(n_members);
    for (std::size_t i = 0; i < votes.size(); ++i) {
      if (votes[i] > 0 && votes[i] <= n_members) {
        vote_histogram[votes[i] - 1] += 1;
      }
    }
    auto out = run_manifest_to_r(manifest);
    out["h"] = load_factor_matrix_from_run(manifest);
    out["genes"] = manifest.genes;
    out["n_removed"] = static_cast<double>(keep_mask.size() - kept);
    out["n_members"] = n_members;
    out["min_votes"] = min_votes;
    out["vote_histogram"] = vote_histogram;
    out["correction_summary"] =
        parquet_table_to_df(celladmix::correction_summary_parquet_path(manifest.paths.root_dir));
    return out;
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_correct_run: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_ensemble_prepare(
    SEXP path_sexp,
    SEXP num_threads_sexp) {
  try {
    const int count = celladmix::ensure_ensemble_labels(
        as<std::string>(path_sexp),
        as<int>(num_threads_sexp));
    return Rcpp::wrap(count);
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_ensemble_prepare: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_ensemble_member_fractions(
    SEXP path_sexp,
    SEXP member_sexp) {
  try {
    const auto fractions = celladmix::ensemble_member_cell_fractions(
        as<std::string>(path_sexp),
        as<int>(member_sexp));
    const auto cells = celladmix::load_run_cells(as<std::string>(path_sexp));
    return List::create(
        _["fractions"] = matrix_to_r(fractions),
        _["cell_id"] = cells.cell_ids);
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_ensemble_member_fractions: unknown C++ exception");
  }
  return R_NilValue;
}

extern "C" SEXP _cellAdmixCore_celladmix_simulate_nsclc(
    SEXP transcripts_per_cell_sexp,
    SEXP admixture_per_target_cell_sexp,
    SEXP seed_sexp) {
  try {
    celladmix::SimulationParams params;
    params.transcripts_per_cell = as<int>(transcripts_per_cell_sexp);
    params.admixture_per_target_cell = as<int>(admixture_per_target_cell_sexp);
    params.seed = static_cast<unsigned int>(as<int>(seed_sexp));
    const auto sim = celladmix::simulate_nsclc_admixture(params);
    return List::create(
        _["transcripts"] = transcripts_to_df_selected(
            sim.transcripts, select_transcript_rows(sim.transcripts.size())),
        _["cells"] = cells_to_df(sim.cells),
        _["malignant_markers"] = sim.malignant_markers,
        _["fibro_markers"] = sim.fibro_markers,
        _["true_admixture"] = sim.true_admixture);
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("celladmix_simulate_nsclc: unknown C++ exception");
  }
  return R_NilValue;
}
