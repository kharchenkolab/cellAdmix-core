#include "celladmix/run_store.hpp"

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <arrow/result.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/file_reader.h>
#include <parquet/statistics.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "celladmix/graph.hpp"
#include "celladmix/ncv.hpp"
#include "celladmix/nmf_kl.hpp"
#include "celladmix/pipeline_common.hpp"
#include "celladmix/workflow.hpp"

namespace celladmix {

NmfRunDiagnostics nmf_diagnostics_from_fit(const SparseNmfResult& fit) {
  NmfRunDiagnostics out;
  out.final_objective = fit.final_objective;
  out.selected_seed = fit.selected_seed;
  out.selected_run = fit.selected_run;
  out.candidate_final_objectives = fit.candidate_final_objectives;
  out.candidate_best_match_correlations = fit.candidate_best_match_correlations;
  out.selected_factor_stability = fit.selected_factor_stability;
  out.candidate_final_objective_mean = fit.candidate_final_objective_mean;
  out.candidate_final_objective_sd = fit.candidate_final_objective_sd;
  out.candidate_best_match_correlation_mean = fit.candidate_best_match_correlation_mean;
  out.stability_comparison_runs = fit.stability_comparison_runs;
  out.stable_factor_count = fit.stable_factor_count;
  out.stability_threshold = fit.stability_threshold;
  return out;
}

// Persisted run I/O for molecules, cells, factors, manifests, and cropped reloads.

namespace {

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

// Accept either a run directory or a run.json path and return the run directory.
std::filesystem::path resolve_run_dir(const std::string& path_or_dir) {
  const auto path = std::filesystem::path(path_or_dir);
  if (std::filesystem::is_directory(path)) {
    return path;
  }
  return path.parent_path();
}

// Resolve a possibly relative manifest path against the run root.
std::filesystem::path resolve_path(const std::filesystem::path& root, const std::string& path) {
  const auto candidate = std::filesystem::path(path);
  if (candidate.is_absolute()) {
    return candidate;
  }
  return root / candidate;
}

// Create an output directory if it does not already exist.
void ensure_directory(const std::filesystem::path& path) {
  std::error_code ec;
  std::filesystem::create_directories(path, ec);
  if (ec) {
    throw std::runtime_error(
        "Could not create directory '" + path.string() + "': " + ec.message());
  }
}

// Materialize a string Arrow array from a standard vector.
std::shared_ptr<arrow::Array> build_string_array(const std::vector<std::string>& values) {
  arrow::StringBuilder builder;
  arrow_check(builder.Reserve(static_cast<int64_t>(values.size())), "Reserve string builder");
  for (const auto& value : values) {
    arrow_check(builder.Append(value), "Append string");
  }
  return arrow_unwrap(builder.Finish(), "Finish string array");
}

// Materialize a double Arrow array from a standard vector.
std::shared_ptr<arrow::Array> build_double_array(const std::vector<double>& values) {
  arrow::DoubleBuilder builder;
  arrow_check(builder.Reserve(static_cast<int64_t>(values.size())), "Reserve double builder");
  for (double value : values) {
    arrow_check(builder.Append(value), "Append double");
  }
  return arrow_unwrap(builder.Finish(), "Finish double array");
}

// Materialize an int32 Arrow array from a standard vector.
std::shared_ptr<arrow::Array> build_int32_array(const std::vector<int>& values) {
  arrow::Int32Builder builder;
  arrow_check(builder.Reserve(static_cast<int64_t>(values.size())), "Reserve int32 builder");
  for (int value : values) {
    arrow_check(builder.Append(value), "Append int32");
  }
  return arrow_unwrap(builder.Finish(), "Finish int32 array");
}

// Materialize an int64 Arrow array from a standard vector.
std::shared_ptr<arrow::Array> build_int64_array(const std::vector<std::int64_t>& values) {
  arrow::Int64Builder builder;
  arrow_check(builder.Reserve(static_cast<int64_t>(values.size())), "Reserve int64 builder");
  for (std::int64_t value : values) {
    arrow_check(builder.Append(value), "Append int64");
  }
  return arrow_unwrap(builder.Finish(), "Finish int64 array");
}

// Write one Arrow table to parquet with the requested row-group size.
void write_parquet_table(
    const std::shared_ptr<arrow::Table>& table,
    const std::filesystem::path& path,
    int row_group_size) {
  auto sink = arrow_unwrap(
      arrow::io::FileOutputStream::Open(path.string()),
      "Open parquet output");
  auto writer = arrow_unwrap(
      parquet::arrow::FileWriter::Open(*table->schema(), arrow::default_memory_pool(), sink),
      "Open parquet writer");
  if (table->schema()->metadata()) {
    arrow_check(writer->AddKeyValueMetadata(table->schema()->metadata()), "Attach parquet metadata");
  }
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

// Compute the top-two factor margin per transcript.
std::vector<double> compute_factor_margin(const DenseMatrix& factor_scores) {
  std::vector<double> out(static_cast<std::size_t>(factor_scores.rows()), 0.0);
  for (int row = 0; row < factor_scores.rows(); ++row) {
    double best = -std::numeric_limits<double>::infinity();
    double second = -std::numeric_limits<double>::infinity();
    for (int col = 0; col < factor_scores.cols(); ++col) {
      const double value = factor_scores(row, col);
      if (value >= best) {
        second = best;
        best = value;
      } else if (value > second) {
        second = value;
      }
    }
    if (!std::isfinite(second)) {
      second = 0.0;
    }
    out[static_cast<std::size_t>(row)] = best - second;
  }
  return out;
}

// Convert one coordinate to the integer tile index used for parquet partition columns.
int tile_index(double value, double tile_size) {
  return static_cast<int>(std::floor(value / tile_size));
}

// Expand per-row crop metadata for molecule parquet output.
std::vector<std::string> crop_names_for_rows(
    std::size_t n_rows,
    const std::optional<std::string>& analysis_crop,
    const std::vector<std::string>* transcript_crop_ids) {
  if (analysis_crop.has_value()) {
    return std::vector<std::string>(n_rows, *analysis_crop);
  }
  if (transcript_crop_ids != nullptr) {
    if (transcript_crop_ids->size() != n_rows) {
      throw std::runtime_error("transcript_crop_ids must match transcript rows");
    }
    return *transcript_crop_ids;
  }
  return std::vector<std::string>(n_rows, "");
}

// Collect unique non-empty strings while preserving first-seen order.
std::vector<std::string> unique_nonempty(const std::vector<std::string>& values) {
  std::vector<std::string> out;
  out.reserve(values.size());
  for (const auto& value : values) {
    if (value.empty()) {
      continue;
    }
    if (std::find(out.begin(), out.end(), value) == out.end()) {
      out.push_back(value);
    }
  }
  return out;
}

struct MoleculeWriteColumns {
  std::vector<std::int64_t> obs_id;
  std::vector<int> crop_idx;
  std::vector<int> tile_x;
  std::vector<int> tile_y;
  std::vector<int> tile_z;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<int> gene_idx;
  std::vector<int> cell_idx;
  std::vector<int> factor_label;
  std::vector<double> factor_margin;
  std::vector<double> qv;
  std::vector<int> overlaps_nucleus;
  std::vector<double> nucleus_distance;
  std::vector<std::string> transcript_id;
};

// Assemble sorted molecule parquet columns, optionally applying a keep mask.
MoleculeWriteColumns build_molecule_columns(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    const std::vector<double>& factor_margin,
    const std::vector<std::string>& crop_names,
    double tile_size,
    const std::vector<bool>* keep_mask = nullptr) {
  if (labels.size() != table.size()) {
    throw std::runtime_error("labels must match TranscriptTable size");
  }
  if (!factor_margin.empty() && factor_margin.size() != table.size()) {
    throw std::runtime_error("factor_margin must match TranscriptTable size");
  }
  if (crop_names.size() != table.size()) {
    throw std::runtime_error("crop_names must match TranscriptTable size");
  }
  if (tile_size <= 0.0) {
    throw std::runtime_error("tile_size must be positive");
  }

  const auto crop_ids = unique_nonempty(crop_names);
  std::unordered_map<std::string, int> crop_lookup;
  crop_lookup.reserve(crop_ids.size());
  for (std::size_t i = 0; i < crop_ids.size(); ++i) {
    crop_lookup.emplace(crop_ids[i], static_cast<int>(i));
  }

  std::vector<std::size_t> order;
  order.reserve(table.size());
  for (std::size_t i = 0; i < table.size(); ++i) {
    if (keep_mask != nullptr) {
      if (keep_mask->size() != table.size()) {
        throw std::runtime_error("keep_mask must match TranscriptTable size");
      }
      if (!(*keep_mask)[i]) {
        continue;
      }
    }
    order.push_back(i);
  }

  auto crop_idx_for = [&](std::size_t i) {
    const auto it = crop_lookup.find(crop_names[i]);
    return it == crop_lookup.end() ? -1 : it->second;
  };

  std::stable_sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
    const int crop_l = crop_idx_for(lhs);
    const int crop_r = crop_idx_for(rhs);
    if (crop_l != crop_r) return crop_l < crop_r;
    const int tile_x_l = tile_index(table.x[lhs], tile_size);
    const int tile_x_r = tile_index(table.x[rhs], tile_size);
    if (tile_x_l != tile_x_r) return tile_x_l < tile_x_r;
    const int tile_y_l = tile_index(table.y[lhs], tile_size);
    const int tile_y_r = tile_index(table.y[rhs], tile_size);
    if (tile_y_l != tile_y_r) return tile_y_l < tile_y_r;
    const int tile_z_l = tile_index(table.z[lhs], tile_size);
    const int tile_z_r = tile_index(table.z[rhs], tile_size);
    if (tile_z_l != tile_z_r) return tile_z_l < tile_z_r;
    return lhs < rhs;
  });

  MoleculeWriteColumns out;
  out.obs_id.reserve(order.size());
  out.crop_idx.reserve(order.size());
  out.tile_x.reserve(order.size());
  out.tile_y.reserve(order.size());
  out.tile_z.reserve(order.size());
  out.x.reserve(order.size());
  out.y.reserve(order.size());
  out.z.reserve(order.size());
  out.gene_idx.reserve(order.size());
  out.cell_idx.reserve(order.size());
  out.factor_label.reserve(order.size());
  out.factor_margin.reserve(order.size());
  if (!table.qv.empty()) out.qv.reserve(order.size());
  if (!table.overlaps_nucleus.empty()) out.overlaps_nucleus.reserve(order.size());
  if (!table.nucleus_distance.empty()) out.nucleus_distance.reserve(order.size());
  if (!table.transcript_ids.empty()) out.transcript_id.reserve(order.size());

  for (const std::size_t i : order) {
    out.obs_id.push_back(static_cast<std::int64_t>(i));
    out.crop_idx.push_back(crop_idx_for(i));
    out.tile_x.push_back(tile_index(table.x[i], tile_size));
    out.tile_y.push_back(tile_index(table.y[i], tile_size));
    out.tile_z.push_back(tile_index(table.z[i], tile_size));
    out.x.push_back(table.x[i]);
    out.y.push_back(table.y[i]);
    out.z.push_back(table.z[i]);
    out.gene_idx.push_back(table.gene_index[i]);
    out.cell_idx.push_back(table.cell_index[i]);
    out.factor_label.push_back(labels[i]);
    out.factor_margin.push_back(
        factor_margin.empty() ? 0.0 : factor_margin[i]);
    if (!table.qv.empty()) out.qv.push_back(table.qv[i]);
    if (!table.overlaps_nucleus.empty()) out.overlaps_nucleus.push_back(table.overlaps_nucleus[i]);
    if (!table.nucleus_distance.empty()) out.nucleus_distance.push_back(table.nucleus_distance[i]);
    if (!table.transcript_ids.empty()) out.transcript_id.push_back(table.transcript_ids[i]);
  }
  return out;
}

// Write molecule-level results, coordinates, and factor assignments to parquet.
void write_molecules_parquet(
    const std::filesystem::path& path,
    const TranscriptTable& table,
    const std::vector<int>& labels,
    const std::vector<double>& factor_margin,
    const std::vector<std::string>& crop_names,
    const RunStorageOptions& storage_options,
    const std::vector<bool>* keep_mask = nullptr) {
  const auto columns = build_molecule_columns(
      table,
      labels,
      factor_margin,
      crop_names,
      storage_options.tile_size,
      keep_mask);

  std::vector<std::shared_ptr<arrow::Field>> fields = {
      arrow::field("obs_id", arrow::int64()),
      arrow::field("crop_idx", arrow::int32()),
      arrow::field("tile_x", arrow::int32()),
      arrow::field("tile_y", arrow::int32()),
      arrow::field("tile_z", arrow::int32()),
      arrow::field("x", arrow::float64()),
      arrow::field("y", arrow::float64()),
      arrow::field("z", arrow::float64()),
      arrow::field("gene_idx", arrow::int32()),
      arrow::field("cell_idx", arrow::int32()),
      arrow::field("factor_label", arrow::int32()),
      arrow::field("factor_margin", arrow::float64()),
  };

  std::vector<std::shared_ptr<arrow::Array>> arrays = {
      build_int64_array(columns.obs_id),
      build_int32_array(columns.crop_idx),
      build_int32_array(columns.tile_x),
      build_int32_array(columns.tile_y),
      build_int32_array(columns.tile_z),
      build_double_array(columns.x),
      build_double_array(columns.y),
      build_double_array(columns.z),
      build_int32_array(columns.gene_idx),
      build_int32_array(columns.cell_idx),
      build_int32_array(columns.factor_label),
      build_double_array(columns.factor_margin),
  };

  if (!table.qv.empty()) {
    fields.push_back(arrow::field("qv", arrow::float64()));
    arrays.push_back(build_double_array(columns.qv));
  }
  if (!table.overlaps_nucleus.empty()) {
    fields.push_back(arrow::field("overlaps_nucleus", arrow::int32()));
    arrays.push_back(build_int32_array(columns.overlaps_nucleus));
  }
  if (!table.nucleus_distance.empty()) {
    fields.push_back(arrow::field("nucleus_distance", arrow::float64()));
    arrays.push_back(build_double_array(columns.nucleus_distance));
  }
  if (!table.transcript_ids.empty()) {
    fields.push_back(arrow::field("transcript_id", arrow::utf8()));
    arrays.push_back(build_string_array(columns.transcript_id));
  }

  const auto table_out = arrow::Table::Make(std::make_shared<arrow::Schema>(fields), arrays);
  write_parquet_table(table_out, path, storage_options.parquet_row_group_size);
}

// Persist the training-row identity map used to trace sampled NCVs back to molecules.
void write_training_rows_parquet(
    const std::filesystem::path& path,
    const std::vector<int>& training_query_indices,
    int row_group_size) {
  std::vector<std::int64_t> obs_id;
  std::vector<int> training_rank;
  obs_id.reserve(training_query_indices.size());
  training_rank.reserve(training_query_indices.size());

  for (std::size_t i = 0; i < training_query_indices.size(); ++i) {
    obs_id.push_back(static_cast<std::int64_t>(training_query_indices[i]));
    training_rank.push_back(static_cast<int>(i));
  }

  const auto table_out = arrow::Table::Make(
      arrow::schema({
          arrow::field("training_rank", arrow::int32()),
          arrow::field("obs_id", arrow::int64()),
      }),
      {
          build_int32_array(training_rank),
          build_int64_array(obs_id),
      });
  write_parquet_table(table_out, path, row_group_size);
}

// Persist per-cell molecule removal counts for a corrected run.
void write_correction_summary_parquet(
    const std::filesystem::path& path,
    const RunData& run_data,
    const std::vector<bool>& keep_mask,
    int row_group_size) {
  std::vector<std::string> cell_id = run_data.cells.cell_ids;
  std::vector<std::string> cell_type = run_data.cells.cell_types;
  if (cell_id.empty()) {
    cell_id = run_data.transcripts.cells;
    cell_type.assign(cell_id.size(), "");
  }
  if (cell_type.size() != cell_id.size()) {
    cell_type.assign(cell_id.size(), "");
  }

  std::unordered_map<std::string, int> cell_lookup;
  cell_lookup.reserve(cell_id.size());
  for (std::size_t i = 0; i < cell_id.size(); ++i) {
    cell_lookup.emplace(cell_id[i], static_cast<int>(i));
  }

  std::vector<int> before(cell_id.size(), 0);
  std::vector<int> after(cell_id.size(), 0);
  std::vector<int> removed(cell_id.size(), 0);
  for (std::size_t i = 0; i < run_data.transcripts.size(); ++i) {
    std::string id = run_data.transcripts.orig_cell_id[i];
    if (!run_data.transcripts.cell_index.empty() &&
        run_data.transcripts.cell_index[i] >= 0 &&
        static_cast<std::size_t>(run_data.transcripts.cell_index[i]) < run_data.transcripts.cells.size()) {
      id = run_data.transcripts.cells[static_cast<std::size_t>(run_data.transcripts.cell_index[i])];
    }
    const auto it = cell_lookup.find(id);
    if (it == cell_lookup.end()) {
      continue;
    }
    const auto idx = static_cast<std::size_t>(it->second);
    before[idx] += 1;
    if (keep_mask[i]) {
      after[idx] += 1;
    } else {
      removed[idx] += 1;
    }
    if (cell_type[idx].empty() && run_data.transcripts.has_cell_types()) {
      cell_type[idx] = run_data.transcripts.cell_types[i];
    }
  }

  std::vector<double> fraction_removed(cell_id.size(), 0.0);
  std::vector<int> modified(cell_id.size(), 0);
  for (std::size_t i = 0; i < cell_id.size(); ++i) {
    if (before[i] > 0) {
      fraction_removed[i] = static_cast<double>(removed[i]) / static_cast<double>(before[i]);
    }
    modified[i] = removed[i] > 0 ? 1 : 0;
  }

  const auto table_out = arrow::Table::Make(
      arrow::schema({
          arrow::field("cell_id", arrow::utf8()),
          arrow::field("cell_type", arrow::utf8()),
          arrow::field("n_molecules_before", arrow::int32()),
          arrow::field("n_molecules_after", arrow::int32()),
          arrow::field("n_removed", arrow::int32()),
          arrow::field("fraction_removed", arrow::float64()),
          arrow::field("modified", arrow::int32()),
      }),
      {
          build_string_array(cell_id),
          build_string_array(cell_type),
          build_int32_array(before),
          build_int32_array(after),
          build_int32_array(removed),
          build_double_array(fraction_removed),
          build_int32_array(modified),
      });
  write_parquet_table(table_out, path, row_group_size);
}

// Persist factor loadings as a long factor-by-gene parquet table.
void write_factors_parquet(
    const std::filesystem::path& path,
    const TranscriptTable& table,
    const BasicPipelineResult& fit,
    int row_group_size) {
  std::vector<int> factor_id;
  std::vector<int> gene_idx;
  std::vector<std::string> gene;
  std::vector<double> loading;
  const std::size_t total =
      static_cast<std::size_t>(fit.nmf.h.rows()) * static_cast<std::size_t>(fit.nmf.h.cols());
  factor_id.reserve(total);
  gene_idx.reserve(total);
  gene.reserve(total);
  loading.reserve(total);

  for (int factor = 0; factor < fit.nmf.h.rows(); ++factor) {
    for (int gene_i = 0; gene_i < fit.nmf.h.cols(); ++gene_i) {
      factor_id.push_back(factor + 1);
      gene_idx.push_back(gene_i);
      gene.push_back(table.genes[static_cast<std::size_t>(gene_i)]);
      loading.push_back(fit.nmf.h(factor, gene_i));
    }
  }

  const auto table_out = arrow::Table::Make(
      arrow::schema({
          arrow::field("factor_id", arrow::int32()),
          arrow::field("gene_idx", arrow::int32()),
          arrow::field("gene", arrow::utf8()),
          arrow::field("loading", arrow::float64()),
      }),
      {
          build_int32_array(factor_id),
          build_int32_array(gene_idx),
          build_string_array(gene),
          build_double_array(loading),
      });
  write_parquet_table(table_out, path, row_group_size);
}

// Build a cell table aligned to the write order expected by cells.parquet.
CellTable build_cells_for_write(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    const CellTable* source_cells,
    bool include_all_source_cells) {
  const auto summaries = summarize_cells(table, &labels);
  std::unordered_map<std::string, CellFactorSummary> summary_by_id;
  summary_by_id.reserve(summaries.size());
  for (const auto& summary : summaries) {
    summary_by_id.emplace(summary.cell_id, summary);
  }

  CellTable out;
  std::vector<std::string> ordered_cell_ids;
  if (source_cells != nullptr && include_all_source_cells) {
    ordered_cell_ids = source_cells->cell_ids;
  } else {
    ordered_cell_ids = table.cells;
  }

  const int n_factors = summaries.empty()
      ? 0
      : static_cast<int>(summaries.front().factor_fractions.size());

  out.cell_ids.reserve(ordered_cell_ids.size());
  out.centroid_x.reserve(ordered_cell_ids.size());
  out.centroid_y.reserve(ordered_cell_ids.size());
  out.centroid_z.reserve(ordered_cell_ids.size());
  if (source_cells != nullptr && !source_cells->cell_types.empty()) {
    out.cell_types.reserve(ordered_cell_ids.size());
  }
  if (source_cells != nullptr && !source_cells->sample_ids.empty()) {
    out.sample_ids.reserve(ordered_cell_ids.size());
  }
  if (source_cells != nullptr && !source_cells->fov_ids.empty()) {
    out.fov_ids.reserve(ordered_cell_ids.size());
  }

  // Store factor fractions in parallel vectors using fov_ids/sample_ids? no.
  (void)n_factors;
  for (const auto& cell_id : ordered_cell_ids) {
    out.cell_ids.push_back(cell_id);
    const auto it = summary_by_id.find(cell_id);
    if (source_cells != nullptr) {
      const auto src_it = std::find(source_cells->cell_ids.begin(), source_cells->cell_ids.end(), cell_id);
      if (src_it != source_cells->cell_ids.end()) {
        const std::size_t idx = static_cast<std::size_t>(std::distance(source_cells->cell_ids.begin(), src_it));
        out.centroid_x.push_back(source_cells->centroid_x[idx]);
        out.centroid_y.push_back(source_cells->centroid_y[idx]);
        out.centroid_z.push_back(source_cells->centroid_z[idx]);
        if (!source_cells->cell_types.empty()) out.cell_types.push_back(source_cells->cell_types[idx]);
        if (!source_cells->sample_ids.empty()) out.sample_ids.push_back(source_cells->sample_ids[idx]);
        if (!source_cells->fov_ids.empty()) out.fov_ids.push_back(source_cells->fov_ids[idx]);
        continue;
      }
    }
    if (it != summary_by_id.end()) {
      out.centroid_x.push_back(it->second.centroid_x);
      out.centroid_y.push_back(it->second.centroid_y);
      out.centroid_z.push_back(it->second.centroid_z);
      if (source_cells != nullptr && !source_cells->cell_types.empty()) {
        out.cell_types.push_back(it->second.cell_type);
      }
      if (source_cells != nullptr && !source_cells->sample_ids.empty()) {
        out.sample_ids.push_back("");
      }
      if (source_cells != nullptr && !source_cells->fov_ids.empty()) {
        out.fov_ids.push_back("");
      }
    } else {
      out.centroid_x.push_back(0.0);
      out.centroid_y.push_back(0.0);
      out.centroid_z.push_back(0.0);
      if (source_cells != nullptr && !source_cells->cell_types.empty()) out.cell_types.push_back("");
      if (source_cells != nullptr && !source_cells->sample_ids.empty()) out.sample_ids.push_back("");
      if (source_cells != nullptr && !source_cells->fov_ids.empty()) out.fov_ids.push_back("");
    }
  }

  return out;
}

// Write per-cell consensus factors and summaries to parquet.
void write_cells_parquet(
    const std::filesystem::path& path,
    const TranscriptTable& table,
    const std::vector<int>& labels,
    const CellTable* source_cells,
    bool include_all_source_cells,
    const std::vector<std::string>* ordered_cell_ids_override,
    int row_group_size) {
  const auto summaries = summarize_cells(table, &labels);
  std::unordered_map<std::string, CellFactorSummary> summary_by_id;
  summary_by_id.reserve(summaries.size());
  int n_factors = 0;
  for (const auto& summary : summaries) {
    summary_by_id.emplace(summary.cell_id, summary);
    n_factors = std::max<int>(n_factors, static_cast<int>(summary.factor_fractions.size()));
  }

  std::vector<std::string> ordered_cell_ids;
  if (ordered_cell_ids_override != nullptr) {
    ordered_cell_ids = *ordered_cell_ids_override;
  } else if (source_cells != nullptr && include_all_source_cells) {
    ordered_cell_ids = source_cells->cell_ids;
  } else {
    ordered_cell_ids = table.cells;
  }

  std::unordered_map<std::string, std::size_t> source_lookup;
  if (source_cells != nullptr) {
    source_lookup.reserve(source_cells->cell_ids.size());
    for (std::size_t i = 0; i < source_cells->cell_ids.size(); ++i) {
      source_lookup.emplace(source_cells->cell_ids[i], i);
    }
  }

  std::vector<int> cell_idx;
  std::vector<std::string> cell_id;
  std::vector<std::string> cell_type;
  std::vector<std::string> sample_id;
  std::vector<std::string> fov_id;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<int> transcript_count;
  std::vector<int> dominant_factor;
  std::vector<double> dominant_fraction;
  std::vector<std::vector<double>> factor_fractions(static_cast<std::size_t>(n_factors));

  const bool write_cell_type = source_cells != nullptr && !source_cells->cell_types.empty();
  const bool write_sample = source_cells != nullptr && !source_cells->sample_ids.empty();
  const bool write_fov = source_cells != nullptr && !source_cells->fov_ids.empty();

  for (std::size_t i = 0; i < ordered_cell_ids.size(); ++i) {
    const auto& id = ordered_cell_ids[i];
    const auto summary_it = summary_by_id.find(id);
    cell_idx.push_back(static_cast<int>(i));
    cell_id.push_back(id);

    std::size_t source_index = std::numeric_limits<std::size_t>::max();
    const auto src_it = source_lookup.find(id);
    if (src_it != source_lookup.end()) {
      source_index = src_it->second;
    }

    if (source_index != std::numeric_limits<std::size_t>::max()) {
      x.push_back(source_cells->centroid_x[source_index]);
      y.push_back(source_cells->centroid_y[source_index]);
      z.push_back(source_cells->centroid_z[source_index]);
      if (write_cell_type) cell_type.push_back(source_cells->cell_types[source_index]);
      if (write_sample) sample_id.push_back(source_cells->sample_ids[source_index]);
      if (write_fov) fov_id.push_back(source_cells->fov_ids[source_index]);
    } else if (summary_it != summary_by_id.end()) {
      x.push_back(summary_it->second.centroid_x);
      y.push_back(summary_it->second.centroid_y);
      z.push_back(summary_it->second.centroid_z);
      if (write_cell_type) cell_type.push_back(summary_it->second.cell_type);
      if (write_sample) sample_id.push_back("");
      if (write_fov) fov_id.push_back("");
    } else {
      x.push_back(0.0);
      y.push_back(0.0);
      z.push_back(0.0);
      if (write_cell_type) cell_type.push_back("");
      if (write_sample) sample_id.push_back("");
      if (write_fov) fov_id.push_back("");
    }

    if (summary_it != summary_by_id.end()) {
      transcript_count.push_back(summary_it->second.transcript_count);
      dominant_factor.push_back(summary_it->second.dominant_factor + 1);
      dominant_fraction.push_back(summary_it->second.dominant_fraction);
      for (int factor = 0; factor < n_factors; ++factor) {
        factor_fractions[static_cast<std::size_t>(factor)].push_back(
            factor < static_cast<int>(summary_it->second.factor_fractions.size())
            ? summary_it->second.factor_fractions[static_cast<std::size_t>(factor)]
            : 0.0);
      }
    } else {
      transcript_count.push_back(0);
      dominant_factor.push_back(-1);
      dominant_fraction.push_back(0.0);
      for (int factor = 0; factor < n_factors; ++factor) {
        factor_fractions[static_cast<std::size_t>(factor)].push_back(0.0);
      }
    }
  }

  std::vector<std::shared_ptr<arrow::Field>> fields = {
      arrow::field("cell_idx", arrow::int32()),
      arrow::field("cell_id", arrow::utf8()),
      arrow::field("x", arrow::float64()),
      arrow::field("y", arrow::float64()),
      arrow::field("z", arrow::float64()),
      arrow::field("transcript_count", arrow::int32()),
      arrow::field("dominant_factor", arrow::int32()),
      arrow::field("dominant_fraction", arrow::float64()),
  };
  std::vector<std::shared_ptr<arrow::Array>> arrays = {
      build_int32_array(cell_idx),
      build_string_array(cell_id),
      build_double_array(x),
      build_double_array(y),
      build_double_array(z),
      build_int32_array(transcript_count),
      build_int32_array(dominant_factor),
      build_double_array(dominant_fraction),
  };
  if (write_cell_type) {
    fields.push_back(arrow::field("cell_type", arrow::utf8()));
    arrays.push_back(build_string_array(cell_type));
  }
  if (write_sample) {
    fields.push_back(arrow::field("sample_id", arrow::utf8()));
    arrays.push_back(build_string_array(sample_id));
  }
  if (write_fov) {
    fields.push_back(arrow::field("fov_id", arrow::utf8()));
    arrays.push_back(build_string_array(fov_id));
  }
  for (int factor = 0; factor < n_factors; ++factor) {
    fields.push_back(
        arrow::field("factor_" + std::to_string(factor + 1) + "_fraction", arrow::float64()));
    arrays.push_back(build_double_array(factor_fractions[static_cast<std::size_t>(factor)]));
  }

  const auto table_out = arrow::Table::Make(std::make_shared<arrow::Schema>(fields), arrays);
  write_parquet_table(table_out, path, row_group_size);
}

// Serialize the run manifest to JSON.
nlohmann::json manifest_to_json(const RunManifest& manifest) {
  nlohmann::json source_files = nlohmann::json::array();
  for (const auto& file : manifest.source.files) {
    source_files.push_back(
        {
            {"role", file.role},
            {"path", file.path},
            {"exists", file.exists},
        });
  }
  return nlohmann::json{
      {"format_version", manifest.format_version},
      {"package_version", manifest.package_version},
      {"annotation_hash", manifest.annotation_hash},
      {"run_type", manifest.run_type},
      {"source",
       {
           {"type", manifest.source.type},
           {"path", manifest.source.path},
           {"used_parquet", manifest.source.used_parquet},
           {"files", source_files},
       }},
      {"paths",
       {
           {"run_json", std::filesystem::path(manifest.paths.run_json).filename().string()},
           {"factors_parquet", std::filesystem::path(manifest.paths.factors_parquet).filename().string()},
           {"cells_parquet", std::filesystem::path(manifest.paths.cells_parquet).filename().string()},
           {"molecules_parquet", std::filesystem::path(manifest.paths.molecules_parquet).filename().string()},
           {"training_rows_parquet", std::filesystem::path(manifest.paths.training_rows_parquet).filename().string()},
           {"scores_dir", std::filesystem::path(manifest.paths.scores_dir).filename().string()},
           {"corrected_dir", std::filesystem::path(manifest.paths.corrected_dir).filename().string()},
       }},
      {"pipeline_options",
       {
           {"ncv_k", manifest.pipeline_options.ncv_k},
           {"rank", manifest.pipeline_options.rank},
           {"graph_k", manifest.pipeline_options.graph_k},
           {"same_label_ratio", manifest.pipeline_options.same_label_ratio},
           {"nmf_iterations", manifest.pipeline_options.nmf_iterations},
           {"nmf_init", manifest.pipeline_options.nmf_init},
           {"nmf_variant", manifest.pipeline_options.nmf_variant},
           {"molecule_scoring", manifest.pipeline_options.molecule_scoring},
           {"nmf_n_runs", manifest.pipeline_options.nmf_n_runs},
           {"nmf_train_max_rows", manifest.pipeline_options.nmf_train_max_rows},
           {"nmf_min_molecules", manifest.pipeline_options.nmf_min_molecules},
           {"num_threads", manifest.pipeline_options.num_threads},
           {"return_ncv", manifest.pipeline_options.return_ncv},
           {"seed", manifest.pipeline_options.seed},
           {"training_scope_cell_types", manifest.pipeline_options.training_scope_cell_types},
       }},
      {"storage_options",
       {
           {"tile_size", manifest.storage_options.tile_size},
           {"parquet_row_group_size", manifest.storage_options.parquet_row_group_size},
       }},
      {"analysis_crop", manifest.analysis_crop.has_value() ? nlohmann::json(*manifest.analysis_crop) : nlohmann::json(nullptr)},
      {"parent_run", manifest.parent_run.has_value() ? nlohmann::json(*manifest.parent_run) : nlohmann::json(nullptr)},
      {"genes", manifest.genes},
      {"nmf_gene_weights", manifest.nmf_gene_weights},
      {"nmf_diagnostics",
       {
           {"final_objective", manifest.nmf_diagnostics.final_objective},
           {"selected_seed", manifest.nmf_diagnostics.selected_seed},
           {"selected_run", manifest.nmf_diagnostics.selected_run},
           {"candidate_final_objectives", manifest.nmf_diagnostics.candidate_final_objectives},
           {"candidate_best_match_correlations", manifest.nmf_diagnostics.candidate_best_match_correlations},
           {"selected_factor_stability", manifest.nmf_diagnostics.selected_factor_stability},
           {"candidate_final_objective_mean", manifest.nmf_diagnostics.candidate_final_objective_mean},
           {"candidate_final_objective_sd", manifest.nmf_diagnostics.candidate_final_objective_sd},
           {"candidate_best_match_correlation_mean", manifest.nmf_diagnostics.candidate_best_match_correlation_mean},
           {"stability_metric", manifest.nmf_diagnostics.stability_metric},
           {"stability_comparison_runs", manifest.nmf_diagnostics.stability_comparison_runs},
           {"stable_factor_count", manifest.nmf_diagnostics.stable_factor_count},
           {"stability_threshold", manifest.nmf_diagnostics.stability_threshold},
       }},
      {"nmf_transform_target_row_sum", manifest.nmf_transform_target_row_sum},
      {"crop_ids", manifest.crop_ids},
      {"n_transcripts", manifest.n_transcripts},
      {"n_cells", manifest.n_cells},
      {"n_factors", manifest.n_factors},
      {"n_training_rows", manifest.n_training_rows},
      {"has_z", manifest.has_z},
      {"has_qv", manifest.has_qv},
      {"has_transcript_id", manifest.has_transcript_id},
      {"has_cell_type", manifest.has_cell_type},
      {"has_sample_id", manifest.has_sample_id},
      {"has_fov_id", manifest.has_fov_id},
      {"has_nucleus_id", manifest.has_nucleus_id},
      {"has_overlaps_nucleus", manifest.has_overlaps_nucleus},
      {"has_nucleus_distance", manifest.has_nucleus_distance},
  };
}

// Parse a run manifest JSON file back into strongly typed paths and options.
RunManifest manifest_from_json(const nlohmann::json& json, const std::filesystem::path& root_dir) {
  RunManifest manifest;
  manifest.format_version = json.at("format_version").get<std::string>();
  manifest.package_version = json.value("package_version", std::string{});
  manifest.annotation_hash = json.value("annotation_hash", std::string{});
  manifest.run_type = json.at("run_type").get<std::string>();
  const auto& source = json.at("source");
  manifest.source.type = source.at("type").get<std::string>();
  manifest.source.path = source.at("path").get<std::string>();
  manifest.source.used_parquet = source.at("used_parquet").get<bool>();
  if (source.contains("files")) {
    for (const auto& file_json : source.at("files")) {
      RunSourceFile file;
      file.role = file_json.at("role").get<std::string>();
      file.path = file_json.at("path").get<std::string>();
      file.exists = file_json.value("exists", false);
      manifest.source.files.push_back(std::move(file));
    }
  }
  const auto& paths = json.at("paths");
  manifest.paths.root_dir = root_dir.string();
  manifest.paths.run_json = resolve_path(root_dir, paths.at("run_json").get<std::string>()).string();
  manifest.paths.factors_parquet =
      resolve_path(root_dir, paths.at("factors_parquet").get<std::string>()).string();
  manifest.paths.cells_parquet =
      resolve_path(root_dir, paths.at("cells_parquet").get<std::string>()).string();
  manifest.paths.molecules_parquet =
      resolve_path(root_dir, paths.at("molecules_parquet").get<std::string>()).string();
  if (paths.contains("training_rows_parquet")) {
    manifest.paths.training_rows_parquet =
        resolve_path(root_dir, paths.at("training_rows_parquet").get<std::string>()).string();
  }
  manifest.paths.scores_dir =
      resolve_path(root_dir, paths.at("scores_dir").get<std::string>()).string();
  manifest.paths.corrected_dir =
      resolve_path(root_dir, paths.at("corrected_dir").get<std::string>()).string();
  const auto& pipeline_options = json.at("pipeline_options");
  manifest.pipeline_options.ncv_k = pipeline_options.at("ncv_k").get<int>();
  manifest.pipeline_options.rank = pipeline_options.at("rank").get<int>();
  manifest.pipeline_options.graph_k = pipeline_options.at("graph_k").get<int>();
  manifest.pipeline_options.same_label_ratio = pipeline_options.at("same_label_ratio").get<double>();
  manifest.pipeline_options.nmf_iterations = pipeline_options.at("nmf_iterations").get<int>();
  manifest.pipeline_options.nmf_init = pipeline_options.value("nmf_init", std::string("auto"));
  manifest.pipeline_options.nmf_variant = pipeline_options.value("nmf_variant", std::string("kl"));
  manifest.pipeline_options.molecule_scoring =
      pipeline_options.value("molecule_scoring", std::string("auto"));
  manifest.pipeline_options.nmf_n_runs = pipeline_options.value("nmf_n_runs", 1);
  manifest.pipeline_options.nmf_train_max_rows = pipeline_options.at("nmf_train_max_rows").get<int>();
  manifest.pipeline_options.nmf_min_molecules = pipeline_options.value("nmf_min_molecules", 10);
  manifest.pipeline_options.num_threads = pipeline_options.value("num_threads", 1);
  manifest.pipeline_options.return_ncv = pipeline_options.at("return_ncv").get<bool>();
  manifest.pipeline_options.seed = pipeline_options.at("seed").get<unsigned int>();
  manifest.pipeline_options.training_scope_cell_types =
      pipeline_options.value("training_scope_cell_types", std::vector<std::string>{});
  const auto& storage_options = json.at("storage_options");
  manifest.storage_options.tile_size = storage_options.at("tile_size").get<double>();
  manifest.storage_options.parquet_row_group_size =
      storage_options.at("parquet_row_group_size").get<int>();
  if (!json.at("analysis_crop").is_null()) {
    manifest.analysis_crop = json.at("analysis_crop").get<std::string>();
  }
  if (!json.at("parent_run").is_null()) {
    manifest.parent_run = json.at("parent_run").get<std::string>();
  }
  manifest.genes = json.at("genes").get<std::vector<std::string>>();
  manifest.nmf_gene_weights = json.value("nmf_gene_weights", std::vector<double>{});
  if (json.contains("nmf_diagnostics") && json.at("nmf_diagnostics").is_object()) {
    const auto& diag = json.at("nmf_diagnostics");
    manifest.nmf_diagnostics.final_objective = diag.value("final_objective", 0.0);
    manifest.nmf_diagnostics.selected_seed = diag.value("selected_seed", 1U);
    manifest.nmf_diagnostics.selected_run = diag.value("selected_run", 0);
    manifest.nmf_diagnostics.candidate_final_objectives =
        diag.value("candidate_final_objectives", std::vector<double>{});
    manifest.nmf_diagnostics.candidate_best_match_correlations =
        diag.value("candidate_best_match_correlations", std::vector<double>{});
    manifest.nmf_diagnostics.selected_factor_stability =
        diag.value("selected_factor_stability", std::vector<double>{});
    manifest.nmf_diagnostics.candidate_final_objective_mean =
        diag.value("candidate_final_objective_mean", 0.0);
    manifest.nmf_diagnostics.candidate_final_objective_sd =
        diag.value("candidate_final_objective_sd", 0.0);
    manifest.nmf_diagnostics.candidate_best_match_correlation_mean =
        diag.value("candidate_best_match_correlation_mean", 0.0);
    manifest.nmf_diagnostics.stability_metric =
        diag.value("stability_metric", std::string{"best_match_legacy"});
    manifest.nmf_diagnostics.stability_comparison_runs =
        diag.value("stability_comparison_runs", 0);
    manifest.nmf_diagnostics.stable_factor_count =
        diag.value("stable_factor_count", 0);
    manifest.nmf_diagnostics.stability_threshold =
        diag.value("stability_threshold", 0.0);
  }
  manifest.nmf_transform_target_row_sum = json.value("nmf_transform_target_row_sum", 0.0);
  manifest.crop_ids = json.at("crop_ids").get<std::vector<std::string>>();
  manifest.n_transcripts = json.at("n_transcripts").get<std::size_t>();
  manifest.n_cells = json.at("n_cells").get<std::size_t>();
  manifest.n_factors = json.at("n_factors").get<std::size_t>();
  manifest.n_training_rows = json.value("n_training_rows", static_cast<std::size_t>(0));
  manifest.has_z = json.at("has_z").get<bool>();
  manifest.has_qv = json.at("has_qv").get<bool>();
  manifest.has_transcript_id = json.at("has_transcript_id").get<bool>();
  manifest.has_cell_type = json.at("has_cell_type").get<bool>();
  manifest.has_sample_id = json.at("has_sample_id").get<bool>();
  manifest.has_fov_id = json.at("has_fov_id").get<bool>();
  manifest.has_nucleus_id = json.value("has_nucleus_id", false);
  manifest.has_overlaps_nucleus = json.value("has_overlaps_nucleus", false);
  manifest.has_nucleus_distance = json.value("has_nucleus_distance", false);
  return manifest;
}

// Persist the manifest JSON into the run directory.
void write_manifest_json_impl(const RunManifest& manifest) {
  std::ofstream out(manifest.paths.run_json);
  if (!out) {
    throw std::runtime_error("Could not open run.json for writing: " + manifest.paths.run_json);
  }
  out << manifest_to_json(manifest).dump(2) << "\n";
}

// Read an arbitrary parquet file into one Arrow table.
std::shared_ptr<arrow::Table> read_parquet_table(const std::string& path) {
  auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path), "Open parquet file");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open parquet reader");
  auto reader = arrow_unwrap(builder.Build(), "Build parquet reader");
  std::shared_ptr<arrow::Table> table;
  arrow_check(reader->ReadTable(&table), "Read parquet table");
  return table;
}

// Look up one required column index in an Arrow schema.
int find_column_index(const std::shared_ptr<arrow::Schema>& schema, const std::string& name) {
  const int idx = schema->GetFieldIndex(name);
  if (idx < 0) {
    throw std::runtime_error("Missing required column: " + name);
  }
  return idx;
}

// Check whether an Arrow schema exposes one named column.
bool has_column(const std::shared_ptr<arrow::Schema>& schema, const std::string& name) {
  return schema->GetFieldIndex(name) >= 0;
}

// Extract one UTF-8 column into a standard string vector.
std::vector<std::string> extract_string_column(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& name) {
  const int idx = find_column_index(table->schema(), name);
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
    } else {
      const auto arr = std::static_pointer_cast<arrow::LargeStringArray>(chunk);
      for (int64_t i = 0; i < chunk->length(); ++i) {
        out[static_cast<std::size_t>(offset + i)] = arr->GetString(i);
      }
    }
    offset += chunk->length();
  }
  return out;
}

// Extract one numeric column into doubles, accepting float or integer storage.
std::vector<double> extract_double_column(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& name) {
  const int idx = find_column_index(table->schema(), name);
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
    } else {
      const auto arr = std::static_pointer_cast<arrow::FloatArray>(chunk);
      for (int64_t i = 0; i < chunk->length(); ++i) {
        out[static_cast<std::size_t>(offset + i)] = static_cast<double>(arr->Value(i));
      }
    }
    offset += chunk->length();
  }
  return out;
}

// Extract one integer column into int32 values, accepting int64 storage as well.
std::vector<int> extract_int32_column(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& name) {
  const int idx = find_column_index(table->schema(), name);
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
    } else {
      const auto arr = std::static_pointer_cast<arrow::Int64Array>(chunk);
      for (int64_t i = 0; i < chunk->length(); ++i) {
        out[static_cast<std::size_t>(offset + i)] = static_cast<int>(arr->Value(i));
      }
    }
    offset += chunk->length();
  }
  return out;
}

// Lightweight typed view over one numeric Arrow array chunk.
struct NumericArrayView {
  const arrow::DoubleArray* double_arr = nullptr;
  const arrow::FloatArray* float_arr = nullptr;
  const arrow::Int64Array* int64_arr = nullptr;
  const arrow::Int32Array* int32_arr = nullptr;

  explicit NumericArrayView(const std::shared_ptr<arrow::Array>& array) {
    if (!array) return;
    switch (array->type_id()) {
      case arrow::Type::DOUBLE:
        double_arr = static_cast<const arrow::DoubleArray*>(array.get());
        break;
      case arrow::Type::FLOAT:
        float_arr = static_cast<const arrow::FloatArray*>(array.get());
        break;
      case arrow::Type::INT64:
        int64_arr = static_cast<const arrow::Int64Array*>(array.get());
        break;
      case arrow::Type::INT32:
        int32_arr = static_cast<const arrow::Int32Array*>(array.get());
        break;
      default:
        throw std::runtime_error("Unsupported numeric column type in run store");
    }
  }

  double value(int64_t i) const {
    if (double_arr) return double_arr->Value(i);
    if (float_arr) return static_cast<double>(float_arr->Value(i));
    if (int64_arr) return static_cast<double>(int64_arr->Value(i));
    return static_cast<double>(int32_arr->Value(i));
  }

  int int_value(int64_t i) const {
    if (int64_arr) return static_cast<int>(int64_arr->Value(i));
    if (int32_arr) return int32_arr->Value(i);
    return static_cast<int>(std::llround(value(i)));
  }
};

// Lightweight typed view over one Arrow string array chunk.
struct StringArrayView {
  const arrow::StringArray* string_arr = nullptr;
  const arrow::LargeStringArray* large_string_arr = nullptr;

  explicit StringArrayView(const std::shared_ptr<arrow::Array>& array) {
    if (!array) return;
    if (array->type_id() == arrow::Type::STRING) {
      string_arr = static_cast<const arrow::StringArray*>(array.get());
    } else if (array->type_id() == arrow::Type::LARGE_STRING) {
      large_string_arr = static_cast<const arrow::LargeStringArray*>(array.get());
    } else {
      throw std::runtime_error("Unsupported string column type in run store");
    }
  }

  std::string value(int64_t i) const {
    if (string_arr) return string_arr->GetString(i);
    if (large_string_arr) return large_string_arr->GetString(i);
    return {};
  }
};

// Read row-group min/max statistics for numeric pruning.
std::optional<std::pair<double, double>> parquet_numeric_minmax(
    const parquet::ColumnChunkMetaData& column) {
  const auto stats = column.statistics();
  if (stats == nullptr || !stats->HasMinMax()) {
    return std::nullopt;
  }
  switch (column.type()) {
    case parquet::Type::DOUBLE: {
      const auto typed = std::static_pointer_cast<parquet::DoubleStatistics>(stats);
      return std::make_pair(typed->min(), typed->max());
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
      return std::nullopt;
  }
}

// Test whether one row group overlaps a requested numeric range.
bool row_group_overlaps_numeric(
    const parquet::RowGroupMetaData& row_group,
    int column_idx,
    double min_value,
    double max_value) {
  if (column_idx < 0) {
    return true;
  }
  const auto column = row_group.ColumnChunk(column_idx);
  if (column == nullptr) {
    return true;
  }
  const auto minmax = parquet_numeric_minmax(*column);
  if (!minmax.has_value()) {
    return true;
  }
  return !(minmax->second < min_value || minmax->first > max_value);
}

// Select molecule parquet row groups that overlap the requested crop and region.
std::vector<int> select_molecule_row_groups(
    const RunManifest& manifest,
    const std::optional<int>& crop_idx_filter,
    const std::optional<CropBox>& region) {
  auto input = arrow_unwrap(
      arrow::io::ReadableFile::Open(manifest.paths.molecules_parquet),
      "Open molecules parquet");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open molecules reader");
  auto reader = arrow_unwrap(builder.Build(), "Build molecules reader");

  std::shared_ptr<arrow::Schema> schema;
  arrow_check(reader->GetSchema(&schema), "Get molecules schema");

  const int crop_idx_col = schema->GetFieldIndex("crop_idx");
  const int tile_x_col = schema->GetFieldIndex("tile_x");
  const int tile_y_col = schema->GetFieldIndex("tile_y");
  const int tile_z_col = schema->GetFieldIndex("tile_z");

  const int min_tile_x = region.has_value()
      ? tile_index(region->xmin, manifest.storage_options.tile_size)
      : std::numeric_limits<int>::min();
  const int max_tile_x = region.has_value()
      ? tile_index(region->xmax, manifest.storage_options.tile_size)
      : std::numeric_limits<int>::max();
  const int min_tile_y = region.has_value()
      ? tile_index(region->ymin, manifest.storage_options.tile_size)
      : std::numeric_limits<int>::min();
  const int max_tile_y = region.has_value()
      ? tile_index(region->ymax, manifest.storage_options.tile_size)
      : std::numeric_limits<int>::max();
  const int min_tile_z = region.has_value() && region->has_z
      ? tile_index(region->zmin, manifest.storage_options.tile_size)
      : std::numeric_limits<int>::min();
  const int max_tile_z = region.has_value() && region->has_z
      ? tile_index(region->zmax, manifest.storage_options.tile_size)
      : std::numeric_limits<int>::max();

  auto parquet_reader = parquet::ParquetFileReader::OpenFile(manifest.paths.molecules_parquet, false);
  const auto metadata = parquet_reader->metadata();
  std::vector<int> row_groups;
  row_groups.reserve(static_cast<std::size_t>(metadata->num_row_groups()));
  for (int rg = 0; rg < metadata->num_row_groups(); ++rg) {
    const auto row_group = metadata->RowGroup(rg);
    if (crop_idx_filter.has_value() &&
        !row_group_overlaps_numeric(*row_group, crop_idx_col, *crop_idx_filter, *crop_idx_filter)) {
      continue;
    }
    if (region.has_value()) {
      if (!row_group_overlaps_numeric(*row_group, tile_x_col, min_tile_x, max_tile_x)) continue;
      if (!row_group_overlaps_numeric(*row_group, tile_y_col, min_tile_y, max_tile_y)) continue;
      if (region->has_z && !row_group_overlaps_numeric(*row_group, tile_z_col, min_tile_z, max_tile_z)) continue;
    }
    row_groups.push_back(rg);
  }
  return row_groups;
}

// Resolve a crop name into the integer crop index stored in molecules.parquet.
std::optional<int> crop_index_from_name(
    const RunManifest& manifest,
    const std::optional<std::string>& crop_id) {
  if (!crop_id.has_value()) {
    return std::nullopt;
  }
  const auto it = std::find(manifest.crop_ids.begin(), manifest.crop_ids.end(), *crop_id);
  if (it == manifest.crop_ids.end()) {
    return -999999;
  }
  return static_cast<int>(std::distance(manifest.crop_ids.begin(), it));
}

}  // namespace

// Public wrapper that writes the run manifest JSON.
void write_manifest_json(const RunManifest& manifest) {
  write_manifest_json_impl(manifest);
}

// Build the standard set of output paths for one run directory.
RunPaths make_run_paths(const std::string& root_dir) {
  const auto root = std::filesystem::absolute(std::filesystem::path(root_dir));
  return RunPaths{
      root.string(),
      (root / "run.json").string(),
      (root / "factors.parquet").string(),
      (root / "cells.parquet").string(),
      (root / "molecules.parquet").string(),
      (root / "training_rows.parquet").string(),
      (root / "scores").string(),
      (root / "corrected").string(),
  };
}

std::string correction_summary_parquet_path(const std::string& run_path_or_dir) {
  const auto root = resolve_run_dir(run_path_or_dir);
  return (root / "correction_summary.parquet").string();
}

// Persist the core results of one fit into a new run directory.
RunManifest write_basic_run(
    const std::string& root_dir,
    const RunSourceInfo& source,
    const TranscriptTable& table,
    const BasicPipelineResult& fit,
    const BasicPipelineOptions& pipeline_options,
    const RunStorageOptions& storage_options,
    const std::optional<std::string>& analysis_crop,
    const std::vector<std::string>* transcript_crop_ids,
    const CellTable* source_cells) {
  const RunPaths paths = make_run_paths(root_dir);
  ensure_directory(paths.root_dir);
  ensure_directory(paths.scores_dir);
  ensure_directory(paths.corrected_dir);

  const auto crop_names = crop_names_for_rows(table.size(), analysis_crop, transcript_crop_ids);
  const auto factor_margin = compute_factor_margin(fit.factor_scores);

  write_factors_parquet(
      paths.factors_parquet,
      table,
      fit,
      storage_options.parquet_row_group_size);
  if (!fit.nmf.candidate_h.empty()) {
    write_ensemble_h_parquet(
        ensemble_h_parquet_path(paths.root_dir),
        fit.nmf.candidate_h,
        storage_options.parquet_row_group_size);
  }
  write_cells_parquet(
      paths.cells_parquet,
      table,
      fit.labels,
      source_cells,
      /*include_all_source_cells=*/true,
      &table.cells,
      storage_options.parquet_row_group_size);
  write_molecules_parquet(
      paths.molecules_parquet,
      table,
      fit.labels,
      factor_margin,
      crop_names,
      storage_options);
  write_training_rows_parquet(
      paths.training_rows_parquet,
      fit.training_query_indices,
      storage_options.parquet_row_group_size);

  RunManifest manifest;
  manifest.source = source;
  manifest.paths = paths;
  manifest.pipeline_options = pipeline_options;
  manifest.storage_options = storage_options;
  manifest.analysis_crop = analysis_crop;
  manifest.genes = table.genes;
  manifest.nmf_gene_weights = fit.nmf_transform.gene_weights;
  manifest.package_version = kCelladmixVersion;
  manifest.annotation_hash = pipeline_options.annotation_hash;
  manifest.nmf_diagnostics = nmf_diagnostics_from_fit(fit.nmf);
  manifest.nmf_transform_target_row_sum = fit.nmf_transform.target_row_sum;
  manifest.crop_ids = unique_nonempty(crop_names);
  manifest.n_transcripts = table.size();
  manifest.n_cells = source_cells != nullptr ? source_cells->size() : table.num_cells();
  manifest.n_factors = static_cast<std::size_t>(fit.nmf.h.rows());
  manifest.n_training_rows = fit.training_query_indices.size();
  manifest.has_z = true;
  manifest.has_qv = !table.qv.empty();
  manifest.has_transcript_id = !table.transcript_ids.empty();
  manifest.has_cell_type = source_cells != nullptr
      ? !source_cells->cell_types.empty()
      : !table.cell_types.empty();
  manifest.has_sample_id = source_cells != nullptr
      ? !source_cells->sample_ids.empty()
      : !table.sample_ids.empty();
  manifest.has_fov_id = source_cells != nullptr
      ? !source_cells->fov_ids.empty()
      : !table.fov_ids.empty();
  manifest.has_nucleus_id = !table.orig_nucleus_ids.empty();
  manifest.has_overlaps_nucleus = !table.overlaps_nucleus.empty();
  manifest.has_nucleus_distance = !table.nucleus_distance.empty();
  write_manifest_json(manifest);
  return manifest;
}

// Read a run manifest from either run.json or a run directory.
RunManifest read_run_manifest(const std::string& path_or_dir) {
  const auto root = resolve_run_dir(path_or_dir);
  const auto run_json = std::filesystem::is_directory(std::filesystem::path(path_or_dir))
      ? root / "run.json"
      : std::filesystem::path(path_or_dir);
  std::ifstream in(run_json);
  if (!in) {
    throw std::runtime_error("Could not open run manifest: " + run_json.string());
  }
  nlohmann::json json;
  in >> json;
  return manifest_from_json(json, root);
}

// Load only the persisted cell table from a run.
CellTable load_run_cells(const std::string& path_or_dir) {
  const auto manifest = read_run_manifest(path_or_dir);
  const auto table = read_parquet_table(manifest.paths.cells_parquet);
  CellTable cells;
  cells.cell_ids = extract_string_column(table, "cell_id");
  cells.centroid_x = extract_double_column(table, "x");
  cells.centroid_y = extract_double_column(table, "y");
  cells.centroid_z = extract_double_column(table, "z");
  if (has_column(table->schema(), "cell_type")) {
    cells.cell_types = extract_string_column(table, "cell_type");
  }
  if (has_column(table->schema(), "sample_id")) {
    cells.sample_ids = extract_string_column(table, "sample_id");
  }
  if (has_column(table->schema(), "fov_id")) {
    cells.fov_ids = extract_string_column(table, "fov_id");
  }
  return cells;
}

// Load the sampled training molecule ids stored for one run.
std::vector<std::int64_t> load_run_training_obs_ids(const std::string& path_or_dir) {
  const auto manifest = read_run_manifest(path_or_dir);
  if (manifest.paths.training_rows_parquet.empty() ||
      !std::filesystem::exists(manifest.paths.training_rows_parquet)) {
    return {};
  }
  const auto table = read_parquet_table(manifest.paths.training_rows_parquet);
  const int obs_idx = table->schema()->GetFieldIndex("obs_id");
  if (obs_idx < 0) {
    return {};
  }
  const auto chunked = table->column(obs_idx);
  std::vector<std::int64_t> out(static_cast<std::size_t>(chunked->length()));
  int64_t offset = 0;
  for (int c = 0; c < chunked->num_chunks(); ++c) {
    const auto chunk = chunked->chunk(c);
    if (chunk->type_id() == arrow::Type::INT64) {
      const auto arr = std::static_pointer_cast<arrow::Int64Array>(chunk);
      for (int64_t i = 0; i < chunk->length(); ++i) {
        out[static_cast<std::size_t>(offset + i)] = arr->Value(i);
      }
    } else {
      const auto arr = std::static_pointer_cast<arrow::Int32Array>(chunk);
      for (int64_t i = 0; i < chunk->length(); ++i) {
        out[static_cast<std::size_t>(offset + i)] = static_cast<std::int64_t>(arr->Value(i));
      }
    }
    offset += chunk->length();
  }
  return out;
}

// Load only cells that contain sampled training molecules for report-side NCV UMAP rebuilds.
TrainingRunData load_run_training_cell_data(const std::string& path_or_dir) {
  TrainingRunData out;
  out.manifest = read_run_manifest(path_or_dir);
  const CellTable all_cells = load_run_cells(path_or_dir);
  out.training_obs_ids = load_run_training_obs_ids(path_or_dir);
  if (out.training_obs_ids.empty()) {
    return out;
  }

  std::unordered_map<std::int64_t, int> training_rank_by_obs;
  training_rank_by_obs.reserve(out.training_obs_ids.size());
  for (std::size_t rank = 0; rank < out.training_obs_ids.size(); ++rank) {
    training_rank_by_obs.emplace(out.training_obs_ids[rank], static_cast<int>(rank));
  }

  auto input = arrow_unwrap(
      arrow::io::ReadableFile::Open(out.manifest.paths.molecules_parquet),
      "Open molecules parquet");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open molecules reader");
  auto reader = arrow_unwrap(builder.Build(), "Build molecules reader");
  reader->set_use_threads(true);
  reader->set_batch_size(65536);

  std::shared_ptr<arrow::Schema> schema;
  arrow_check(reader->GetSchema(&schema), "Get molecules schema");
  const int obs_col = find_column_index(schema, "obs_id");
  const int cell_col = find_column_index(schema, "cell_idx");
  const int x_col = find_column_index(schema, "x");
  const int y_col = find_column_index(schema, "y");
  const int z_col = find_column_index(schema, "z");
  const int gene_col = find_column_index(schema, "gene_idx");
  const int label_col = find_column_index(schema, "factor_label");
  const int margin_col = find_column_index(schema, "factor_margin");
  const int qv_col = schema->GetFieldIndex("qv");
  const int tx_col = schema->GetFieldIndex("transcript_id");

  auto parquet_reader = parquet::ParquetFileReader::OpenFile(out.manifest.paths.molecules_parquet, false);
  const auto metadata = parquet_reader->metadata();
  std::vector<int> row_groups(static_cast<std::size_t>(metadata->num_row_groups()));
  std::iota(row_groups.begin(), row_groups.end(), 0);

  std::vector<char> keep_cell_mask(all_cells.cell_ids.size(), 0);
  std::vector<int> training_cell_by_rank(out.training_obs_ids.size(), -1);
  auto first_pass = arrow_unwrap(
      reader->GetRecordBatchReader(row_groups, std::vector<int>{obs_col, cell_col}),
      "Create training-cell discovery reader");
  std::shared_ptr<arrow::RecordBatch> batch;
  while (true) {
    arrow_check(first_pass->ReadNext(&batch), "Read training-cell discovery batch");
    if (!batch) break;
    NumericArrayView obs_view(batch->column(0));
    NumericArrayView cell_view(batch->column(1));
    for (int64_t row = 0; row < batch->num_rows(); ++row) {
      const auto obs_id = static_cast<std::int64_t>(obs_view.value(row));
      const auto it = training_rank_by_obs.find(obs_id);
      if (it == training_rank_by_obs.end()) {
        continue;
      }
      const int cell = cell_view.int_value(row);
      if (cell < 0 || cell >= static_cast<int>(all_cells.cell_ids.size())) {
        throw std::runtime_error("cell_idx out of range in molecules.parquet");
      }
      keep_cell_mask[static_cast<std::size_t>(cell)] = 1;
      training_cell_by_rank[static_cast<std::size_t>(it->second)] = cell;
    }
  }
  for (std::size_t rank = 0; rank < training_cell_by_rank.size(); ++rank) {
    if (training_cell_by_rank[rank] < 0) {
      throw std::runtime_error("training obs_id was not found in persisted molecule rows");
    }
  }

  std::vector<int> local_cell_by_global(all_cells.cell_ids.size(), -1);
  for (std::size_t cell = 0; cell < all_cells.cell_ids.size(); ++cell) {
    if (!keep_cell_mask[cell]) {
      continue;
    }
    local_cell_by_global[cell] = static_cast<int>(out.cells.cell_ids.size());
    out.cells.cell_ids.push_back(all_cells.cell_ids[cell]);
    out.cells.centroid_x.push_back(all_cells.centroid_x[cell]);
    out.cells.centroid_y.push_back(all_cells.centroid_y[cell]);
    out.cells.centroid_z.push_back(all_cells.centroid_z[cell]);
    if (!all_cells.cell_types.empty()) out.cells.cell_types.push_back(all_cells.cell_types[cell]);
    if (!all_cells.sample_ids.empty()) out.cells.sample_ids.push_back(all_cells.sample_ids[cell]);
    if (!all_cells.fov_ids.empty()) out.cells.fov_ids.push_back(all_cells.fov_ids[cell]);
  }

  std::vector<int> projected = {obs_col, x_col, y_col, z_col, gene_col, cell_col, label_col, margin_col};
  const bool has_qv = qv_col >= 0;
  const bool has_tx_id = tx_col >= 0;
  if (has_qv) projected.push_back(qv_col);
  if (has_tx_id) projected.push_back(tx_col);
  out.training_query_indices.assign(out.training_obs_ids.size(), -1);

  auto second_pass = arrow_unwrap(
      reader->GetRecordBatchReader(row_groups, projected),
      "Create training-cell molecule reader");
  while (true) {
    arrow_check(second_pass->ReadNext(&batch), "Read training-cell molecule batch");
    if (!batch) break;
    int col = 0;
    NumericArrayView obs_view(batch->column(col++));
    NumericArrayView x_view(batch->column(col++));
    NumericArrayView y_view(batch->column(col++));
    NumericArrayView z_view(batch->column(col++));
    NumericArrayView gene_view(batch->column(col++));
    NumericArrayView cell_view(batch->column(col++));
    NumericArrayView label_view(batch->column(col++));
    NumericArrayView margin_view(batch->column(col++));
    std::optional<NumericArrayView> qv_view;
    if (has_qv) qv_view.emplace(batch->column(col++));
    std::optional<StringArrayView> tx_view;
    if (has_tx_id) tx_view.emplace(batch->column(col++));

    for (int64_t row = 0; row < batch->num_rows(); ++row) {
      const int global_cell = cell_view.int_value(row);
      if (global_cell < 0 || global_cell >= static_cast<int>(local_cell_by_global.size())) {
        throw std::runtime_error("cell_idx out of range in molecules.parquet");
      }
      const int local_cell = local_cell_by_global[static_cast<std::size_t>(global_cell)];
      if (local_cell < 0) {
        continue;
      }
      const int local_row = static_cast<int>(out.transcripts.x.size());
      const auto obs_id = static_cast<std::int64_t>(obs_view.value(row));
      out.obs_ids.push_back(obs_id);
      out.transcripts.x.push_back(x_view.value(row));
      out.transcripts.y.push_back(y_view.value(row));
      out.transcripts.z.push_back(z_view.value(row));
      out.transcripts.gene_index.push_back(gene_view.int_value(row));
      out.transcripts.cell_index.push_back(local_cell);
      if (qv_view.has_value()) out.transcripts.qv.push_back(qv_view->value(row));
      if (tx_view.has_value()) out.transcripts.transcript_ids.push_back(tx_view->value(row));
      out.labels.push_back(label_view.int_value(row));
      out.factor_margin.push_back(margin_view.value(row));
      const auto rank_it = training_rank_by_obs.find(obs_id);
      if (rank_it != training_rank_by_obs.end()) {
        out.training_query_indices[static_cast<std::size_t>(rank_it->second)] = local_row;
      }
    }
  }

  for (std::size_t rank = 0; rank < out.training_query_indices.size(); ++rank) {
    if (out.training_query_indices[rank] < 0) {
      throw std::runtime_error("training obs_id was not found in loaded training-cell molecules");
    }
  }
  out.transcripts.genes = out.manifest.genes;
  out.transcripts.cells = out.cells.cell_ids;
  return out;
}

// Load a run back into memory with optional crop, region, and subsampling filters.
BridgeRunData load_run_bridge_data(const std::string& path_or_dir, const RunLoadOptions& options) {
  BridgeRunData out;
  out.manifest = read_run_manifest(path_or_dir);
  out.cells = load_run_cells(path_or_dir);

  const auto crop_idx_filter = crop_index_from_name(out.manifest, options.crop_id);
  if (crop_idx_filter.has_value() && *crop_idx_filter < -1) {
    out.transcripts.cells = out.cells.cell_ids;
    return out;
  }

  auto input = arrow_unwrap(
      arrow::io::ReadableFile::Open(out.manifest.paths.molecules_parquet),
      "Open molecules parquet");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open molecules reader");
  auto reader = arrow_unwrap(builder.Build(), "Build molecules reader");
  std::shared_ptr<arrow::Schema> schema;
  arrow_check(reader->GetSchema(&schema), "Get molecules schema");

  const bool member_labels_active = options.ensemble_member >= 0;
  EnsembleMemberLabels member_labels;
  if (member_labels_active) {
    member_labels = load_ensemble_member_labels(path_or_dir, options.ensemble_member);
  }

  std::vector<int> projected = {
      find_column_index(schema, "crop_idx"),
      find_column_index(schema, "x"),
      find_column_index(schema, "y"),
      find_column_index(schema, "z"),
      find_column_index(schema, "cell_idx"),
      find_column_index(schema, "factor_label"),
  };
  if (member_labels_active) {
    projected.push_back(find_column_index(schema, "obs_id"));
  }

  const auto row_groups = select_molecule_row_groups(out.manifest, crop_idx_filter, options.region);
  auto batch_reader = arrow_unwrap(
      reader->GetRecordBatchReader(row_groups, projected),
      "Get bridge molecule record batch reader");

  std::shared_ptr<arrow::RecordBatch> batch;
  while (true) {
    arrow_check(batch_reader->ReadNext(&batch), "Read bridge molecule record batch");
    if (!batch) {
      break;
    }
    int col = 0;
    NumericArrayView crop_view(batch->column(col++));
    NumericArrayView x_view(batch->column(col++));
    NumericArrayView y_view(batch->column(col++));
    NumericArrayView z_view(batch->column(col++));
    NumericArrayView cell_view(batch->column(col++));
    NumericArrayView label_view(batch->column(col++));
    std::shared_ptr<arrow::Array> obs_array;
    if (member_labels_active) {
      obs_array = batch->column(col++);
    }
    NumericArrayView obs_view(obs_array);

    for (int64_t i = 0; i < batch->num_rows(); ++i) {
      const int crop_idx = crop_view.int_value(i);
      if (crop_idx_filter.has_value() && crop_idx != *crop_idx_filter) {
        continue;
      }
      const double x = x_view.value(i);
      const double y = y_view.value(i);
      const double z = z_view.value(i);
      if (options.region.has_value()) {
        if (x < options.region->xmin || x > options.region->xmax ||
            y < options.region->ymin || y > options.region->ymax) {
          continue;
        }
        if (options.region->has_z &&
            (z < options.region->zmin || z > options.region->zmax)) {
          continue;
        }
      }

      const int cell_idx = cell_view.int_value(i);
      if (cell_idx < 0 || cell_idx >= static_cast<int>(out.cells.cell_ids.size())) {
        throw std::runtime_error("cell_idx out of range in molecules.parquet");
      }
      out.transcripts.x.push_back(x);
      out.transcripts.y.push_back(y);
      out.transcripts.z.push_back(z);
      out.transcripts.cell_index.push_back(cell_idx);
      out.labels.push_back(member_labels_active
          ? member_labels.label_for(static_cast<std::int64_t>(obs_view.value(i)))
          : label_view.int_value(i));
    }
  }

  out.transcripts.cells = out.cells.cell_ids;
  return out;
}

CellCountMatrix collect_run_counts(const std::string& path_or_dir, const RunLoadOptions& options) {
  const auto manifest = read_run_manifest(path_or_dir);
  CellCountMatrix out;
  out.cells = load_run_cells(path_or_dir);
  out.genes = manifest.genes;
  out.transcript_counts.assign(out.cells.size(), 0);
  out.crop_ids.assign(out.cells.size(), "");

  const auto crop_idx_filter = crop_index_from_name(manifest, options.crop_id);
  out.indptr.reserve(out.cells.size() + 1U);
  out.indptr.push_back(0);
  if (crop_idx_filter.has_value() && *crop_idx_filter < -1) {
    out.indptr.resize(out.cells.size() + 1U, 0);
    return out;
  }

  auto input = arrow_unwrap(
      arrow::io::ReadableFile::Open(manifest.paths.molecules_parquet),
      "Open molecules parquet");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open molecules reader");
  auto reader = arrow_unwrap(builder.Build(), "Build molecules reader");
  reader->set_use_threads(true);
  reader->set_batch_size(65536);

  std::shared_ptr<arrow::Schema> schema;
  arrow_check(reader->GetSchema(&schema), "Get molecules schema");

  std::vector<int> projected;
  auto add_column = [&](const std::string& name) {
    const int pos = static_cast<int>(projected.size());
    projected.push_back(find_column_index(schema, name));
    return pos;
  };
  const int crop_pos = crop_idx_filter.has_value() ? add_column("crop_idx") : -1;
  const int x_pos = options.region.has_value() ? add_column("x") : -1;
  const int y_pos = options.region.has_value() ? add_column("y") : -1;
  const int z_pos = options.region.has_value() ? add_column("z") : -1;
  const int gene_pos = add_column("gene_idx");
  const int cell_pos = add_column("cell_idx");

  const auto row_groups = select_molecule_row_groups(manifest, crop_idx_filter, options.region);
  auto batch_reader = arrow_unwrap(
      reader->GetRecordBatchReader(row_groups, projected),
      "Get sparse count molecule record batch reader");

  std::vector<std::unordered_map<int, int>> counts_by_cell(out.cells.size());
  std::shared_ptr<arrow::RecordBatch> batch;
  while (true) {
    arrow_check(batch_reader->ReadNext(&batch), "Read sparse count molecule batch");
    if (!batch) {
      break;
    }

    std::optional<NumericArrayView> crop_view;
    std::optional<NumericArrayView> x_view;
    std::optional<NumericArrayView> y_view;
    std::optional<NumericArrayView> z_view;
    if (crop_pos >= 0) crop_view.emplace(batch->column(crop_pos));
    if (x_pos >= 0) x_view.emplace(batch->column(x_pos));
    if (y_pos >= 0) y_view.emplace(batch->column(y_pos));
    if (z_pos >= 0) z_view.emplace(batch->column(z_pos));
    NumericArrayView gene_view(batch->column(gene_pos));
    NumericArrayView cell_view(batch->column(cell_pos));

    for (int64_t row = 0; row < batch->num_rows(); ++row) {
      if (crop_view.has_value() && crop_view->int_value(row) != *crop_idx_filter) {
        continue;
      }
      if (options.region.has_value()) {
        const double x = x_view->value(row);
        const double y = y_view->value(row);
        const double z = z_view->value(row);
        if (x < options.region->xmin || x > options.region->xmax ||
            y < options.region->ymin || y > options.region->ymax) {
          continue;
        }
        if (options.region->has_z &&
            (z < options.region->zmin || z > options.region->zmax)) {
          continue;
        }
      }

      const int gene = gene_view.int_value(row);
      const int cell = cell_view.int_value(row);
      if (gene < 0 || gene >= static_cast<int>(out.genes.size())) {
        throw std::runtime_error("gene_idx out of range in molecules.parquet");
      }
      if (cell < 0 || cell >= static_cast<int>(out.cells.size())) {
        throw std::runtime_error("cell_idx out of range in molecules.parquet");
      }
      counts_by_cell[static_cast<std::size_t>(cell)][gene] += 1;
      out.transcript_counts[static_cast<std::size_t>(cell)] += 1;
    }
  }

  std::vector<std::pair<int, int>> entries;
  for (std::size_t cell = 0; cell < counts_by_cell.size(); ++cell) {
    entries.clear();
    entries.reserve(counts_by_cell[cell].size());
    for (const auto& entry : counts_by_cell[cell]) {
      entries.emplace_back(entry.first, entry.second);
    }
    std::sort(entries.begin(), entries.end(), [](const auto& left, const auto& right) {
      return left.first < right.first;
    });
    for (const auto& entry : entries) {
      out.indices.push_back(entry.first);
      out.values.push_back(static_cast<double>(entry.second));
    }
    if (out.indices.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::runtime_error("Sparse run count matrix exceeds R dgCMatrix index capacity");
    }
    out.indptr.push_back(static_cast<int>(out.indices.size()));
  }
  return out;
}

RunData load_run_data(const std::string& path_or_dir, const RunLoadOptions& options) {
  RunData out;
  out.manifest = read_run_manifest(path_or_dir);
  out.cells = load_run_cells(path_or_dir);

  const auto crop_idx_filter = crop_index_from_name(out.manifest, options.crop_id);
  if (crop_idx_filter.has_value() && *crop_idx_filter < -1) {
    return out;
  }

  auto input = arrow_unwrap(
      arrow::io::ReadableFile::Open(out.manifest.paths.molecules_parquet),
      "Open molecules parquet");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open molecules reader");
  auto reader = arrow_unwrap(builder.Build(), "Build molecules reader");
  std::shared_ptr<arrow::Schema> schema;
  arrow_check(reader->GetSchema(&schema), "Get molecules schema");

  std::vector<int> projected = {
      find_column_index(schema, "obs_id"),
      find_column_index(schema, "crop_idx"),
      find_column_index(schema, "x"),
      find_column_index(schema, "y"),
      find_column_index(schema, "z"),
      find_column_index(schema, "gene_idx"),
      find_column_index(schema, "cell_idx"),
      find_column_index(schema, "factor_label"),
      find_column_index(schema, "factor_margin"),
  };
  const bool has_qv = has_column(schema, "qv");
  const bool has_tx_id = has_column(schema, "transcript_id");
  const bool has_overlaps_nucleus = has_column(schema, "overlaps_nucleus");
  const bool has_nucleus_distance = has_column(schema, "nucleus_distance");
  if (has_qv) projected.push_back(find_column_index(schema, "qv"));
  if (has_overlaps_nucleus) projected.push_back(find_column_index(schema, "overlaps_nucleus"));
  if (has_nucleus_distance) projected.push_back(find_column_index(schema, "nucleus_distance"));
  if (has_tx_id) projected.push_back(find_column_index(schema, "transcript_id"));

  const auto row_groups = select_molecule_row_groups(out.manifest, crop_idx_filter, options.region);
  auto batch_reader = arrow_unwrap(
      reader->GetRecordBatchReader(row_groups, projected),
      "Get molecules record batch reader");

  std::vector<std::string> sample_ids;
  std::vector<std::string> fov_ids;
  std::vector<int> source_cell_indices;
  std::vector<int> source_gene_indices;
  std::vector<std::string> crop_names;
  std::vector<std::int64_t> obs_ids;
  std::shared_ptr<arrow::RecordBatch> batch;

  while (true) {
    arrow_check(batch_reader->ReadNext(&batch), "Read molecules record batch");
    if (!batch) {
      break;
    }
    int col = 0;
    NumericArrayView obs_view(batch->column(col++));
    NumericArrayView crop_view(batch->column(col++));
    NumericArrayView x_view(batch->column(col++));
    NumericArrayView y_view(batch->column(col++));
    NumericArrayView z_view(batch->column(col++));
    NumericArrayView gene_view(batch->column(col++));
    NumericArrayView cell_view(batch->column(col++));
    NumericArrayView label_view(batch->column(col++));
    NumericArrayView margin_view(batch->column(col++));
    std::shared_ptr<arrow::Array> qv_array;
    std::shared_ptr<arrow::Array> overlaps_nucleus_array;
    std::shared_ptr<arrow::Array> nucleus_distance_array;
    std::shared_ptr<arrow::Array> tx_array;
    if (has_qv) qv_array = batch->column(col++);
    if (has_overlaps_nucleus) overlaps_nucleus_array = batch->column(col++);
    if (has_nucleus_distance) nucleus_distance_array = batch->column(col++);
    if (has_tx_id) tx_array = batch->column(col++);
    NumericArrayView qv_view(qv_array);
    NumericArrayView overlaps_nucleus_view(overlaps_nucleus_array);
    NumericArrayView nucleus_distance_view(nucleus_distance_array);
    StringArrayView tx_view(tx_array);

    for (int64_t i = 0; i < batch->num_rows(); ++i) {
      const int crop_idx = crop_view.int_value(i);
      if (crop_idx_filter.has_value() && crop_idx != *crop_idx_filter) {
        continue;
      }
      const double x = x_view.value(i);
      const double y = y_view.value(i);
      const double z = z_view.value(i);
      if (options.region.has_value()) {
        if (x < options.region->xmin || x > options.region->xmax ||
            y < options.region->ymin || y > options.region->ymax) {
          continue;
        }
        if (options.region->has_z &&
            (z < options.region->zmin || z > options.region->zmax)) {
          continue;
        }
      }

      out.transcripts.x.push_back(x);
      out.transcripts.y.push_back(y);
      out.transcripts.z.push_back(z);
      const int gene_idx = gene_view.int_value(i);
      const int cell_idx = cell_view.int_value(i);
      if (gene_idx < 0 || gene_idx >= static_cast<int>(out.manifest.genes.size())) {
        throw std::runtime_error("gene_idx out of range in molecules.parquet");
      }
      if (cell_idx < 0 || cell_idx >= static_cast<int>(out.cells.cell_ids.size())) {
        throw std::runtime_error("cell_idx out of range in molecules.parquet");
      }
      out.transcripts.gene_key.push_back(out.manifest.genes[static_cast<std::size_t>(gene_idx)]);
      out.transcripts.orig_cell_id.push_back(out.cells.cell_ids[static_cast<std::size_t>(cell_idx)]);
      if (!out.cells.cell_types.empty()) {
        out.transcripts.cell_types.push_back(out.cells.cell_types[static_cast<std::size_t>(cell_idx)]);
      }
      if (!out.cells.sample_ids.empty()) {
        out.transcripts.sample_ids.push_back(out.cells.sample_ids[static_cast<std::size_t>(cell_idx)]);
      }
      if (!out.cells.fov_ids.empty()) {
        out.transcripts.fov_ids.push_back(out.cells.fov_ids[static_cast<std::size_t>(cell_idx)]);
      }
      if (has_qv) {
        out.transcripts.qv.push_back(qv_view.value(i));
      }
      if (has_overlaps_nucleus) {
        out.transcripts.overlaps_nucleus.push_back(overlaps_nucleus_view.int_value(i));
      }
      if (has_nucleus_distance) {
        out.transcripts.nucleus_distance.push_back(nucleus_distance_view.value(i));
      }
      if (has_tx_id) {
        out.transcripts.transcript_ids.push_back(tx_view.value(i));
      }
      out.labels.push_back(label_view.int_value(i));
      out.factor_margin.push_back(margin_view.value(i));
      out.transcript_crop_ids.push_back(
          crop_idx >= 0 && crop_idx < static_cast<int>(out.manifest.crop_ids.size())
          ? out.manifest.crop_ids[static_cast<std::size_t>(crop_idx)]
          : "");
      const auto obs_id = static_cast<std::int64_t>(obs_view.value(i));
      obs_ids.push_back(obs_id);
      out.obs_ids.push_back(obs_id);
      source_gene_indices.push_back(gene_idx);
      source_cell_indices.push_back(cell_idx);
    }
  }

  if (options.ensemble_member >= 0) {
    const auto member_labels =
        load_ensemble_member_labels(path_or_dir, options.ensemble_member);
    for (std::size_t i = 0; i < out.labels.size(); ++i) {
      out.labels[i] = member_labels.label_for(out.obs_ids[i]);
    }
  }

  if (options.sample_n > 0 && static_cast<std::size_t>(options.sample_n) < out.transcripts.size()) {
    std::vector<std::size_t> order(out.transcripts.size());
    std::iota(order.begin(), order.end(), 0U);
    std::mt19937 rng(options.seed);
    std::shuffle(order.begin(), order.end(), rng);
    order.resize(static_cast<std::size_t>(options.sample_n));
    std::sort(order.begin(), order.end());

    TranscriptTable sampled;
    sampled.genes = out.manifest.genes;
    std::vector<int> sampled_labels;
    std::vector<double> sampled_margin;
    std::vector<int> sampled_source_cells;
    std::vector<int> sampled_source_genes;
    std::vector<std::string> sampled_crop_ids;
    std::vector<std::int64_t> sampled_obs_ids;
    for (const auto idx : order) {
      sampled.x.push_back(out.transcripts.x[idx]);
      sampled.y.push_back(out.transcripts.y[idx]);
      sampled.z.push_back(out.transcripts.z[idx]);
      sampled.gene_key.push_back(out.transcripts.gene_key[idx]);
      sampled.orig_cell_id.push_back(out.transcripts.orig_cell_id[idx]);
      if (!out.transcripts.transcript_ids.empty()) sampled.transcript_ids.push_back(out.transcripts.transcript_ids[idx]);
      if (!out.transcripts.sample_ids.empty()) sampled.sample_ids.push_back(out.transcripts.sample_ids[idx]);
      if (!out.transcripts.fov_ids.empty()) sampled.fov_ids.push_back(out.transcripts.fov_ids[idx]);
      if (!out.transcripts.cell_types.empty()) sampled.cell_types.push_back(out.transcripts.cell_types[idx]);
      if (!out.transcripts.qv.empty()) sampled.qv.push_back(out.transcripts.qv[idx]);
      if (!out.transcripts.overlaps_nucleus.empty()) sampled.overlaps_nucleus.push_back(out.transcripts.overlaps_nucleus[idx]);
      if (!out.transcripts.nucleus_distance.empty()) sampled.nucleus_distance.push_back(out.transcripts.nucleus_distance[idx]);
      sampled_labels.push_back(out.labels[idx]);
      sampled_margin.push_back(out.factor_margin[idx]);
      sampled_source_cells.push_back(source_cell_indices[idx]);
      sampled_source_genes.push_back(source_gene_indices[idx]);
      sampled_crop_ids.push_back(out.transcript_crop_ids[idx]);
      sampled_obs_ids.push_back(out.obs_ids[idx]);
    }
    out.transcripts = std::move(sampled);
    out.obs_ids = std::move(sampled_obs_ids);
    out.labels = std::move(sampled_labels);
    out.factor_margin = std::move(sampled_margin);
    source_cell_indices = std::move(sampled_source_cells);
    source_gene_indices = std::move(sampled_source_genes);
    out.transcript_crop_ids = std::move(sampled_crop_ids);
  }

  out.transcripts.genes = out.manifest.genes;
  out.transcripts.cells = out.cells.cell_ids;
  out.transcripts.gene_index = source_gene_indices;
  out.transcripts.cell_index = source_cell_indices;
  return out;
}

// Write a corrected child run after removing a subset of molecules.
RunManifest write_corrected_run(
    const std::string& root_dir,
    const std::string& parent_run_path_or_dir,
    const RunData& run_data,
    const std::vector<bool>& keep_mask) {
  if (keep_mask.size() != run_data.transcripts.size()) {
    throw std::runtime_error("keep_mask must match run transcript count");
  }

  const RunPaths paths = make_run_paths(root_dir);
  ensure_directory(paths.root_dir);
  ensure_directory(paths.scores_dir);
  ensure_directory(paths.corrected_dir);

  std::vector<int> kept_labels;
  kept_labels.reserve(run_data.transcripts.size());
  std::vector<double> kept_margin;
  kept_margin.reserve(run_data.transcripts.size());
  for (std::size_t i = 0; i < keep_mask.size(); ++i) {
    if (!keep_mask[i]) continue;
    kept_labels.push_back(run_data.labels[i]);
    kept_margin.push_back(run_data.factor_margin[i]);
  }

  const auto corrected_table = subset_transcripts(
      run_data.transcripts,
      std::nullopt,
      nullptr,
      &keep_mask);

  write_cells_parquet(
      paths.cells_parquet,
      corrected_table,
      kept_labels,
      &run_data.cells,
      /*include_all_source_cells=*/true,
      &run_data.transcripts.cells,
      run_data.manifest.storage_options.parquet_row_group_size);
  write_molecules_parquet(
      paths.molecules_parquet,
      run_data.transcripts,
      run_data.labels,
      run_data.factor_margin,
      run_data.transcript_crop_ids,
      run_data.manifest.storage_options,
      &keep_mask);
  write_correction_summary_parquet(
      correction_summary_parquet_path(paths.root_dir),
      run_data,
      keep_mask,
      run_data.manifest.storage_options.parquet_row_group_size);
  std::filesystem::copy_file(
      run_data.manifest.paths.factors_parquet,
      paths.factors_parquet,
      std::filesystem::copy_options::overwrite_existing);
  if (!run_data.manifest.paths.training_rows_parquet.empty() &&
      std::filesystem::exists(run_data.manifest.paths.training_rows_parquet)) {
    std::filesystem::copy_file(
        run_data.manifest.paths.training_rows_parquet,
        paths.training_rows_parquet,
        std::filesystem::copy_options::overwrite_existing);
  }

  RunManifest manifest = run_data.manifest;
  manifest.run_type = "corrected";
  manifest.paths = paths;
  manifest.parent_run = std::filesystem::absolute(resolve_run_dir(parent_run_path_or_dir)).string();
  manifest.n_transcripts = corrected_table.size();
  manifest.n_cells = run_data.cells.size();
  write_manifest_json(manifest);
  return manifest;
}

// --- Ensemble member pool (per-restart factorizations and labelings) ---

namespace {

std::vector<int> read_int32_parquet_column(
    const std::string& path,
    const std::string& name) {
  auto input = arrow_unwrap(
      arrow::io::ReadableFile::Open(path),
      "Open parquet column source");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open parquet column reader");
  auto reader = arrow_unwrap(builder.Build(), "Build parquet column reader");
  std::shared_ptr<arrow::Schema> schema;
  arrow_check(reader->GetSchema(&schema), "Get parquet column schema");
  std::vector<int> row_groups(static_cast<std::size_t>(reader->num_row_groups()));
  std::iota(row_groups.begin(), row_groups.end(), 0);
  auto batch_reader = arrow_unwrap(
      reader->GetRecordBatchReader(row_groups, {find_column_index(schema, name)}),
      "Get parquet column batch reader");
  std::vector<int> out;
  std::shared_ptr<arrow::RecordBatch> batch;
  while (true) {
    arrow_check(batch_reader->ReadNext(&batch), "Read parquet column batch");
    if (!batch) break;
    NumericArrayView view(batch->column(0));
    for (int64_t i = 0; i < batch->num_rows(); ++i) {
      out.push_back(view.int_value(i));
    }
  }
  return out;
}

std::vector<std::int64_t> read_int64_parquet_column(
    const std::string& path,
    const std::string& name) {
  auto input = arrow_unwrap(
      arrow::io::ReadableFile::Open(path),
      "Open parquet column source");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open parquet column reader");
  auto reader = arrow_unwrap(builder.Build(), "Build parquet column reader");
  std::shared_ptr<arrow::Schema> schema;
  arrow_check(reader->GetSchema(&schema), "Get parquet column schema");
  std::vector<int> row_groups(static_cast<std::size_t>(reader->num_row_groups()));
  std::iota(row_groups.begin(), row_groups.end(), 0);
  auto batch_reader = arrow_unwrap(
      reader->GetRecordBatchReader(row_groups, {find_column_index(schema, name)}),
      "Get parquet column batch reader");
  std::vector<std::int64_t> out;
  std::shared_ptr<arrow::RecordBatch> batch;
  while (true) {
    arrow_check(batch_reader->ReadNext(&batch), "Read parquet column batch");
    if (!batch) break;
    NumericArrayView view(batch->column(0));
    for (int64_t i = 0; i < batch->num_rows(); ++i) {
      out.push_back(static_cast<std::int64_t>(view.value(i)));
    }
  }
  return out;
}

}  // namespace

std::string ensemble_h_parquet_path(const std::string& run_path_or_dir) {
  return (resolve_run_dir(run_path_or_dir) / "ensemble_h.parquet").string();
}

std::string ensemble_labels_parquet_path(
    const std::string& run_path_or_dir,
    int member) {
  return (resolve_run_dir(run_path_or_dir) /
      ("ensemble_labels_m" + std::to_string(member) + ".parquet")).string();
}

void write_ensemble_h_parquet(
    const std::string& path,
    const std::vector<DenseMatrix>& candidate_h,
    int row_group_size) {
  std::vector<int> member;
  std::vector<int> factor_id;
  std::vector<int> gene_idx;
  std::vector<double> loading;
  std::size_t total = 0;
  for (const auto& h : candidate_h) {
    total += static_cast<std::size_t>(h.rows()) * static_cast<std::size_t>(h.cols());
  }
  member.reserve(total);
  factor_id.reserve(total);
  gene_idx.reserve(total);
  loading.reserve(total);
  for (std::size_t m = 0; m < candidate_h.size(); ++m) {
    const auto& h = candidate_h[m];
    for (int factor = 0; factor < h.rows(); ++factor) {
      for (int gene_i = 0; gene_i < h.cols(); ++gene_i) {
        member.push_back(static_cast<int>(m));
        factor_id.push_back(factor + 1);
        gene_idx.push_back(gene_i);
        loading.push_back(h(factor, gene_i));
      }
    }
  }
  const auto table_out = arrow::Table::Make(
      arrow::schema({
          arrow::field("member", arrow::int32()),
          arrow::field("factor_id", arrow::int32()),
          arrow::field("gene_idx", arrow::int32()),
          arrow::field("loading", arrow::float64()),
      }),
      {
          build_int32_array(member),
          build_int32_array(factor_id),
          build_int32_array(gene_idx),
          build_double_array(loading),
      });
  write_parquet_table(table_out, std::filesystem::path(path), row_group_size);
}

std::vector<DenseMatrix> load_ensemble_h(const std::string& run_path_or_dir) {
  const auto path = ensemble_h_parquet_path(run_path_or_dir);
  if (!std::filesystem::exists(path)) {
    return {};
  }
  const auto member = read_int32_parquet_column(path, "member");
  const auto factor_id = read_int32_parquet_column(path, "factor_id");
  const auto gene_idx = read_int32_parquet_column(path, "gene_idx");

  auto input = arrow_unwrap(
      arrow::io::ReadableFile::Open(path),
      "Open ensemble H parquet");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open ensemble H reader");
  auto reader = arrow_unwrap(builder.Build(), "Build ensemble H reader");
  std::shared_ptr<arrow::Schema> schema;
  arrow_check(reader->GetSchema(&schema), "Get ensemble H schema");
  std::vector<int> row_groups(static_cast<std::size_t>(reader->num_row_groups()));
  std::iota(row_groups.begin(), row_groups.end(), 0);
  auto batch_reader = arrow_unwrap(
      reader->GetRecordBatchReader(row_groups, {find_column_index(schema, "loading")}),
      "Get ensemble H batch reader");
  std::vector<double> loading;
  loading.reserve(member.size());
  std::shared_ptr<arrow::RecordBatch> batch;
  while (true) {
    arrow_check(batch_reader->ReadNext(&batch), "Read ensemble H batch");
    if (!batch) break;
    NumericArrayView view(batch->column(0));
    for (int64_t i = 0; i < batch->num_rows(); ++i) {
      loading.push_back(view.value(i));
    }
  }

  int n_members = 0;
  int rank = 0;
  int n_genes = 0;
  for (std::size_t i = 0; i < member.size(); ++i) {
    n_members = std::max(n_members, member[i] + 1);
    rank = std::max(rank, factor_id[i]);
    n_genes = std::max(n_genes, gene_idx[i] + 1);
  }
  std::vector<DenseMatrix> out(
      static_cast<std::size_t>(n_members),
      DenseMatrix(rank, n_genes, 0.0));
  for (std::size_t i = 0; i < member.size(); ++i) {
    out[static_cast<std::size_t>(member[i])](factor_id[i] - 1, gene_idx[i]) =
        loading[i];
  }
  return out;
}

int ensemble_member_count(const std::string& run_path_or_dir) {
  int count = 0;
  while (std::filesystem::exists(
      ensemble_labels_parquet_path(run_path_or_dir, count))) {
    ++count;
  }
  return count;
}

int ensure_ensemble_labels(
    const std::string& run_path_or_dir,
    int num_threads) {
  const auto pool = load_ensemble_h(run_path_or_dir);
  if (pool.empty()) {
    return 0;
  }
  const int n_members = static_cast<int>(pool.size());
  bool all_exist = true;
  for (int m = 0; m < n_members; ++m) {
    if (!std::filesystem::exists(ensemble_labels_parquet_path(run_path_or_dir, m))) {
      all_exist = false;
      break;
    }
  }
  if (all_exist) {
    return n_members;
  }

  const auto run = load_run_data(run_path_or_dir);
  const auto& opts = run.manifest.pipeline_options;
  NcvOptions ncv_options;
  ncv_options.k = opts.ncv_k + 1;
  ncv_options.include_self = true;
  ncv_options.within_cell = true;
  const std::string molecule_scoring = resolve_molecule_scoring(opts);
  NcvFeatureTransform transform;
  transform.mode = opts.nmf_variant;
  transform.gene_weights = run.manifest.nmf_gene_weights;
  transform.target_row_sum = run.manifest.nmf_transform_target_row_sum;
  const int threads = std::max(1, num_threads);
  const int selected = run.manifest.nmf_diagnostics.selected_run;

  for (int m = 0; m < n_members; ++m) {
    const auto member_path = ensemble_labels_parquet_path(run_path_or_dir, m);
    if (std::filesystem::exists(member_path)) {
      continue;
    }
    std::vector<int> labels;
    if (m == selected) {
      labels = run.labels;
    } else {
      const DenseMatrix scores = molecule_scoring == "gene_loadings"
          ? project_gene_loadings_to_factors(
                run.transcripts, pool[static_cast<std::size_t>(m)], ncv_options, threads)
          : is_weighted_ls_variant(opts.nmf_variant)
              ? project_ncv_to_factors(
                    run.transcripts, pool[static_cast<std::size_t>(m)], ncv_options, threads)
              : project_ncv_to_factors_kl(
                    run.transcripts, pool[static_cast<std::size_t>(m)], ncv_options, threads,
                    10, 1e-4, 1e-10, transform);
      labels = assign_factors_per_cell(
          run.transcripts, scores, opts.graph_k, opts.same_label_ratio, 20, threads);
    }
    const auto table_out = arrow::Table::Make(
        arrow::schema({arrow::field("label", arrow::int32())}),
        {build_int32_array(labels)});
    write_parquet_table(
        table_out,
        std::filesystem::path(member_path),
        run.manifest.storage_options.parquet_row_group_size);
  }
  return n_members;
}

int EnsembleMemberLabels::label_for(std::int64_t obs_id) const {
  const auto it = std::lower_bound(obs_ids.begin(), obs_ids.end(), obs_id);
  if (it == obs_ids.end() || *it != obs_id) {
    throw std::runtime_error("obs id missing from ensemble member labels");
  }
  return labels[static_cast<std::size_t>(it - obs_ids.begin())];
}

EnsembleMemberLabels load_ensemble_member_labels(
    const std::string& run_path_or_dir,
    int member) {
  const auto manifest = read_run_manifest(run_path_or_dir);
  const auto raw_labels = read_int32_parquet_column(
      ensemble_labels_parquet_path(run_path_or_dir, member), "label");
  const auto raw_obs = read_int64_parquet_column(
      manifest.paths.molecules_parquet, "obs_id");
  if (raw_labels.size() != raw_obs.size()) {
    throw std::runtime_error(
        "ensemble member labels do not match molecules.parquet row count");
  }
  std::vector<std::size_t> order(raw_obs.size());
  std::iota(order.begin(), order.end(), 0U);
  std::sort(order.begin(), order.end(), [&](std::size_t left, std::size_t right) {
    return raw_obs[left] < raw_obs[right];
  });
  EnsembleMemberLabels out;
  out.obs_ids.reserve(raw_obs.size());
  out.labels.reserve(raw_obs.size());
  for (const auto idx : order) {
    out.obs_ids.push_back(raw_obs[idx]);
    out.labels.push_back(raw_labels[idx]);
  }
  return out;
}

DenseMatrix ensemble_member_cell_fractions(
    const std::string& run_path_or_dir,
    int member) {
  const auto manifest = read_run_manifest(run_path_or_dir);
  const auto labels = read_int32_parquet_column(
      ensemble_labels_parquet_path(run_path_or_dir, member), "label");
  const auto cell_idx = read_int32_parquet_column(
      manifest.paths.molecules_parquet, "cell_idx");
  if (labels.size() != cell_idx.size()) {
    throw std::runtime_error(
        "ensemble member labels do not match molecules.parquet row count");
  }
  const auto cells = load_run_cells(run_path_or_dir);
  const int n_cells = static_cast<int>(cells.size());
  const int rank = static_cast<int>(manifest.n_factors);
  DenseMatrix out(n_cells, rank, 0.0);
  std::vector<double> totals(static_cast<std::size_t>(n_cells), 0.0);
  for (std::size_t i = 0; i < labels.size(); ++i) {
    const int cell = cell_idx[i];
    if (cell < 0 || cell >= n_cells) {
      continue;
    }
    totals[static_cast<std::size_t>(cell)] += 1.0;
    const int factor = labels[i];
    if (factor >= 0 && factor < rank) {
      out(cell, factor) += 1.0;
    }
  }
  for (int cell = 0; cell < n_cells; ++cell) {
    const double total = totals[static_cast<std::size_t>(cell)];
    if (total <= 0.0) continue;
    for (int factor = 0; factor < rank; ++factor) {
      out(cell, factor) /= total;
    }
  }
  return out;
}

}  // namespace celladmix
