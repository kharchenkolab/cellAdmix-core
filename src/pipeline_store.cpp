// Cell-blocked fit implementation over the normalized input store.

#include "celladmix/pipeline_store.hpp"

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "celladmix/graph.hpp"
#include "celladmix/input_store.hpp"
#include "celladmix/ncv.hpp"
#include "celladmix/nmf_euclidean.hpp"
#include "celladmix/nmf_kl.hpp"
#include "celladmix/pipeline_common.hpp"
#include "celladmix/report_embedding.hpp"

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
  for (const double value : values) {
    arrow_check(builder.Append(value), "Append double");
  }
  return arrow_unwrap(builder.Finish(), "Finish double array");
}

std::shared_ptr<arrow::Array> build_int32_array(const std::vector<int>& values) {
  arrow::Int32Builder builder;
  arrow_check(builder.Reserve(static_cast<int64_t>(values.size())), "Reserve int32 builder");
  for (const int value : values) {
    arrow_check(builder.Append(value), "Append int32");
  }
  return arrow_unwrap(builder.Finish(), "Finish int32 array");
}

std::shared_ptr<arrow::Array> build_int64_array(const std::vector<std::int64_t>& values) {
  arrow::Int64Builder builder;
  arrow_check(builder.Reserve(static_cast<int64_t>(values.size())), "Reserve int64 builder");
  for (const auto value : values) {
    arrow_check(builder.Append(value), "Append int64");
  }
  return arrow_unwrap(builder.Finish(), "Finish int64 array");
}

void write_parquet_table(
    const std::shared_ptr<arrow::Table>& table,
    const std::filesystem::path& path,
    int row_group_size) {
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

void emit_info(
    const std::chrono::steady_clock::time_point& start,
    const std::string& message,
    double stage_sec) {
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
  if (stage_sec >= 0.0) {
    line << " (" << std::fixed << std::setprecision(3) << stage_sec << "s)";
  }
  std::cout << line.str() << std::endl;
}

struct CompressedTrainingNcv {
  SparseRowMatrix matrix;
  std::vector<int> kept_rows;
  std::vector<int> kept_cols;
};

CompressedTrainingNcv compact_training_ncv(const SparseRowMatrix& x) {
  const auto col_sums = x.col_sums();
  std::vector<int> col_map(static_cast<std::size_t>(x.cols()), -1);
  std::vector<int> kept_cols;
  kept_cols.reserve(static_cast<std::size_t>(x.cols()));
  for (int col = 0; col < x.cols(); ++col) {
    if (col_sums[static_cast<std::size_t>(col)] <= 0.0) {
      continue;
    }
    col_map[static_cast<std::size_t>(col)] = static_cast<int>(kept_cols.size());
    kept_cols.push_back(col);
  }

  std::vector<int> kept_rows;
  std::vector<int> indptr;
  std::vector<int> indices;
  std::vector<double> values;
  indptr.push_back(0);
  for (int row = 0; row < x.rows(); ++row) {
    const int row_begin = static_cast<int>(indices.size());
    for (int p = x.indptr()[static_cast<std::size_t>(row)];
         p < x.indptr()[static_cast<std::size_t>(row + 1)];
         ++p) {
      const int mapped = col_map[static_cast<std::size_t>(x.indices()[static_cast<std::size_t>(p)])];
      if (mapped < 0) {
        continue;
      }
      indices.push_back(mapped);
      values.push_back(x.values()[static_cast<std::size_t>(p)]);
    }
    if (static_cast<int>(indices.size()) == row_begin) {
      continue;
    }
    kept_rows.push_back(row);
    indptr.push_back(static_cast<int>(indices.size()));
  }
  return CompressedTrainingNcv{
      SparseRowMatrix(
          static_cast<int>(kept_rows.size()),
          static_cast<int>(kept_cols.size()),
          std::move(indptr),
          std::move(indices),
          std::move(values)),
      std::move(kept_rows),
      std::move(kept_cols)};
}

std::vector<int> subset_int_vector(const std::vector<int>& values, const std::vector<int>& keep) {
  std::vector<int> out;
  out.reserve(keep.size());
  for (const int index : keep) {
    out.push_back(values[static_cast<std::size_t>(index)]);
  }
  return out;
}

std::vector<int> weighted_capacity_allocation(
    const std::vector<double>& weights,
    const std::vector<int>& capacities,
    int budget,
    std::mt19937& rng) {
  std::vector<int> allocation(weights.size(), 0);
  if (budget <= 0 || weights.empty()) {
    return allocation;
  }
  std::vector<std::size_t> active;
  for (std::size_t i = 0; i < capacities.size(); ++i) {
    if (capacities[i] > 0) active.push_back(i);
  }
  int remaining = budget;
  while (remaining > 0 && !active.empty()) {
    double total_weight = 0.0;
    for (const auto idx : active) total_weight += std::max(weights[idx], 0.0);
    if (total_weight <= 0.0) total_weight = static_cast<double>(active.size());

    std::vector<double> fractional(weights.size(), 0.0);
    int assigned = 0;
    for (const auto idx : active) {
      const int remaining_capacity = capacities[idx] - allocation[idx];
      if (remaining_capacity <= 0) continue;
      const double expected =
          static_cast<double>(remaining) * std::max(weights[idx], 0.0) / total_weight;
      const int add = std::min(remaining_capacity, static_cast<int>(std::floor(expected)));
      allocation[idx] += add;
      assigned += add;
      fractional[idx] = expected - std::floor(expected);
    }

    int leftover = remaining - assigned;
    std::shuffle(active.begin(), active.end(), rng);
    std::stable_sort(active.begin(), active.end(), [&](std::size_t lhs, std::size_t rhs) {
      if (fractional[lhs] == fractional[rhs]) {
        return (capacities[lhs] - allocation[lhs]) > (capacities[rhs] - allocation[rhs]);
      }
      return fractional[lhs] > fractional[rhs];
    });
    for (const auto idx : active) {
      if (leftover == 0) break;
      if (capacities[idx] - allocation[idx] <= 0) continue;
      allocation[idx] += 1;
      assigned += 1;
      leftover -= 1;
    }
    if (assigned == 0) break;
    remaining -= assigned;
    active.erase(
        std::remove_if(
            active.begin(),
            active.end(),
            [&](std::size_t idx) { return capacities[idx] - allocation[idx] <= 0; }),
        active.end());
  }
  return allocation;
}

std::vector<std::string> infer_training_cell_strata(
    const CellTable& cells,
    const std::vector<std::string>& provided) {
  if (!provided.empty()) {
    if (provided.size() != cells.size()) {
      throw std::runtime_error("training_cell_strata must match the number of cells");
    }
    return provided;
  }
  if (cells.cell_types.empty()) {
    return std::vector<std::string>(cells.size(), "__all_cells__");
  }
  std::vector<std::string> out(cells.size(), "__unlabeled__");
  for (std::size_t i = 0; i < cells.size(); ++i) {
    if (!cells.cell_types[i].empty()) out[i] = cells.cell_types[i];
  }
  return out;
}

struct TrainingSelection {
  std::vector<int> global_rows;
  std::vector<int> cell_indices;
};

struct StoreFitScope {
  CellCountMatrix counts;
  std::vector<int> global_cell_by_local;
  std::vector<int> local_cell_by_global;
  std::vector<std::vector<int>> crop_rows_by_cell;
  std::size_t n_molecules = 0;
};

void append_cell_metadata(CellTable& out, const CellTable& source, int cell) {
  const auto idx = static_cast<std::size_t>(cell);
  out.cell_ids.push_back(source.cell_ids[idx]);
  out.centroid_x.push_back(source.centroid_x[idx]);
  out.centroid_y.push_back(source.centroid_y[idx]);
  out.centroid_z.push_back(source.centroid_z[idx]);
  if (!source.cell_types.empty()) out.cell_types.push_back(source.cell_types[idx]);
  if (!source.sample_ids.empty()) out.sample_ids.push_back(source.sample_ids[idx]);
  if (!source.fov_ids.empty()) out.fov_ids.push_back(source.fov_ids[idx]);
}

StoreFitScope build_full_fit_scope(const CellCountMatrix& counts) {
  StoreFitScope scope;
  scope.counts = counts;
  scope.global_cell_by_local.resize(counts.cells.size());
  std::iota(scope.global_cell_by_local.begin(), scope.global_cell_by_local.end(), 0);
  scope.local_cell_by_global.resize(counts.cells.size());
  std::iota(scope.local_cell_by_global.begin(), scope.local_cell_by_global.end(), 0);
  scope.n_molecules = 0;
  for (const int n : counts.transcript_counts) {
    scope.n_molecules += static_cast<std::size_t>(std::max(n, 0));
  }
  return scope;
}

StoreFitScope build_crop_fit_scope(
    const std::string& store_dir,
    const CellCountMatrix& counts,
    const std::vector<InputStoreCellOffset>& offsets,
    const std::string& crop_id) {
  StoreFitScope scope;
  scope.counts.genes = counts.genes;
  scope.local_cell_by_global.assign(counts.cells.size(), -1);
  std::vector<std::unordered_map<int, int>> gene_counts_by_cell;
  constexpr std::int64_t target_block_rows = 512 * 1024;

  auto ensure_local_cell = [&](int global_cell) -> int {
    if (global_cell < 0 || global_cell >= static_cast<int>(counts.cells.size())) {
      throw std::runtime_error("input-store molecule cell_idx is out of range");
    }
    int& local = scope.local_cell_by_global[static_cast<std::size_t>(global_cell)];
    if (local >= 0) {
      return local;
    }
    local = static_cast<int>(scope.counts.cells.size());
    scope.global_cell_by_local.push_back(global_cell);
    append_cell_metadata(scope.counts.cells, counts.cells, global_cell);
    scope.counts.crop_ids.push_back(crop_id);
    scope.counts.transcript_counts.push_back(0);
    scope.crop_rows_by_cell.emplace_back();
    gene_counts_by_cell.emplace_back();
    return local;
  };

  int cell = 0;
  while (cell < static_cast<int>(offsets.size())) {
    while (cell < static_cast<int>(offsets.size()) && offsets[static_cast<std::size_t>(cell)].start < 0) {
      ++cell;
    }
    if (cell >= static_cast<int>(offsets.size())) break;
    std::int64_t range_start = offsets[static_cast<std::size_t>(cell)].start;
    std::int64_t range_end = offsets[static_cast<std::size_t>(cell)].end;
    ++cell;
    while (cell < static_cast<int>(offsets.size())) {
      if (offsets[static_cast<std::size_t>(cell)].start < 0) {
        ++cell;
        continue;
      }
      const std::int64_t candidate_end = offsets[static_cast<std::size_t>(cell)].end;
      if (candidate_end - range_start > target_block_rows && range_end > range_start) {
        break;
      }
      range_end = candidate_end;
      ++cell;
    }

    const auto block = load_input_store_molecule_range(store_dir, range_start, range_end);
    if (block.size() == 0) continue;
    if (block.crop_id.size() != block.size()) {
      throw std::runtime_error("analysis_crop requires crop_id metadata in the input store");
    }
    for (std::size_t row = 0; row < block.size(); ++row) {
      if (block.crop_id[row] != crop_id) {
        continue;
      }
      const int local_cell = ensure_local_cell(block.cell_idx[row]);
      const auto local_idx = static_cast<std::size_t>(local_cell);
      scope.crop_rows_by_cell[local_idx].push_back(static_cast<int>(block.row_index[row]));
      scope.counts.transcript_counts[local_idx] += 1;
      gene_counts_by_cell[local_idx][block.gene_idx[row]] += 1;
      scope.n_molecules += 1;
    }
  }

  scope.counts.indptr.reserve(scope.counts.cells.size() + 1);
  scope.counts.indptr.push_back(0);
  for (auto& gene_counts : gene_counts_by_cell) {
    std::vector<std::pair<int, int>> entries;
    entries.reserve(gene_counts.size());
    for (const auto& kv : gene_counts) entries.emplace_back(kv.first, kv.second);
    std::sort(entries.begin(), entries.end());
    for (const auto& entry : entries) {
      scope.counts.indices.push_back(entry.first);
      scope.counts.values.push_back(static_cast<double>(entry.second));
    }
    scope.counts.indptr.push_back(static_cast<int>(scope.counts.indices.size()));
  }
  return scope;
}

std::vector<std::string> resolve_scope_cell_strata(
    const StoreFitScope& scope,
    const std::vector<std::string>& provided,
    std::size_t full_cell_count) {
  if (provided.empty()) {
    return infer_training_cell_strata(scope.counts.cells, provided);
  }
  if (provided.size() == scope.counts.cells.size()) {
    return infer_training_cell_strata(scope.counts.cells, provided);
  }
  if (provided.size() != full_cell_count) {
    throw std::runtime_error("training_cell_strata must match the full store cells or the fit-scope cells");
  }
  std::vector<std::string> subset;
  subset.reserve(scope.global_cell_by_local.size());
  for (const int global_cell : scope.global_cell_by_local) {
    subset.push_back(provided[static_cast<std::size_t>(global_cell)]);
  }
  return infer_training_cell_strata(scope.counts.cells, subset);
}

std::vector<std::string> unique_nonempty(const std::vector<std::string>& values) {
  std::vector<std::string> out;
  for (const auto& value : values) {
    if (value.empty()) continue;
    if (std::find(out.begin(), out.end(), value) == out.end()) out.push_back(value);
  }
  return out;
}

std::string join_strings(const std::vector<std::string>& values, const char* sep) {
  std::ostringstream out;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i > 0) out << sep;
    out << values[i];
  }
  return out.str();
}

std::size_t filter_training_strata_to_cell_types(
    const CellTable& cells,
    std::vector<std::string>& cell_strata,
    const std::vector<std::string>& requested_cell_types) {
  const auto requested = unique_nonempty(requested_cell_types);
  if (requested.empty()) {
    return cells.size();
  }
  if (cell_strata.size() != cells.size()) {
    throw std::runtime_error("training_cell_strata must match the scoped cells before applying training scope");
  }
  std::unordered_set<std::string> requested_set(requested.begin(), requested.end());
  std::unordered_set<std::string> seen;
  std::size_t selected_cells = 0;
  for (std::size_t i = 0; i < cells.size(); ++i) {
    const auto& stratum = cell_strata[i];
    if (requested_set.find(stratum) == requested_set.end()) {
      cell_strata[i].clear();
      continue;
    }
    seen.insert(stratum);
    ++selected_cells;
  }

  std::vector<std::string> missing;
  for (const auto& cell_type : requested) {
    if (seen.find(cell_type) == seen.end()) missing.push_back(cell_type);
  }
  if (!missing.empty()) {
    throw std::runtime_error(
        "scope values were not found in the training annotation: " + join_strings(missing, ", "));
  }
  if (selected_cells == 0) {
    throw std::runtime_error("scope did not select any cells for NMF training");
  }
  return selected_cells;
}

TrainingSelection sample_training_queries_from_offsets(
    const std::vector<InputStoreCellOffset>& offsets,
    const std::vector<std::string>& cell_strata,
    int max_rows,
    unsigned int seed,
    int min_molecules) {
  std::mt19937 rng(seed);
  const bool has_explicit_subset =
      !cell_strata.empty() &&
      std::any_of(cell_strata.begin(), cell_strata.end(), [](const std::string& value) {
        return value.empty();
      });

  std::vector<int> eligible_cells;
  for (std::size_t cell = 0; cell < offsets.size(); ++cell) {
    const auto& offset = offsets[cell];
    const int size = offset.start < 0 || offset.end <= offset.start
        ? 0
        : static_cast<int>(offset.end - offset.start);
    if (size <= 0) continue;
    if (has_explicit_subset && cell_strata[cell].empty()) continue;
    if (size >= std::max(min_molecules, 1)) eligible_cells.push_back(static_cast<int>(cell));
  }
  if (eligible_cells.empty()) {
    if (has_explicit_subset) {
      throw std::runtime_error("explicit training strata did not label any cells above the training size threshold");
    }
    for (std::size_t cell = 0; cell < offsets.size(); ++cell) {
      if (offsets[cell].start >= 0 && offsets[cell].end > offsets[cell].start) {
        eligible_cells.push_back(static_cast<int>(cell));
      }
    }
  }

  int eligible_rows = 0;
  for (const int cell : eligible_cells) {
    eligible_rows += static_cast<int>(offsets[static_cast<std::size_t>(cell)].end -
        offsets[static_cast<std::size_t>(cell)].start);
  }
  if (eligible_rows == 0) return {};

  std::vector<int> allocation(offsets.size(), 0);
  if (max_rows <= 0 || max_rows >= eligible_rows) {
    for (const int cell : eligible_cells) {
      allocation[static_cast<std::size_t>(cell)] =
          static_cast<int>(offsets[static_cast<std::size_t>(cell)].end -
              offsets[static_cast<std::size_t>(cell)].start);
    }
  } else if (max_rows < static_cast<int>(eligible_cells.size())) {
    auto candidate_cells = eligible_cells;
    std::shuffle(candidate_cells.begin(), candidate_cells.end(), rng);
    candidate_cells.resize(static_cast<std::size_t>(max_rows));
    for (const int cell : candidate_cells) {
      allocation[static_cast<std::size_t>(cell)] = 1;
    }
  } else {
    std::vector<std::string> strata_order;
    std::unordered_map<std::string, int> stratum_lookup;
    std::vector<std::vector<int>> cells_by_stratum;
    for (const int cell : eligible_cells) {
      const std::string key = cell_strata[static_cast<std::size_t>(cell)];
      const auto it = stratum_lookup.find(key);
      if (it == stratum_lookup.end()) {
        const int stratum = static_cast<int>(strata_order.size());
        stratum_lookup.emplace(key, stratum);
        strata_order.push_back(key);
        cells_by_stratum.push_back({cell});
      } else {
        cells_by_stratum[static_cast<std::size_t>(it->second)].push_back(cell);
      }
    }

    int remaining = max_rows;
    for (const int cell : eligible_cells) {
      allocation[static_cast<std::size_t>(cell)] = 1;
      remaining -= 1;
      if (remaining == 0) break;
    }
    if (remaining > 0) {
      std::vector<double> stratum_weights(cells_by_stratum.size(), 0.0);
      std::vector<int> stratum_capacities(cells_by_stratum.size(), 0);
      std::vector<std::vector<int>> residual_by_cell(cells_by_stratum.size());
      for (std::size_t stratum = 0; stratum < cells_by_stratum.size(); ++stratum) {
        int total_residual = 0;
        for (const int cell : cells_by_stratum[stratum]) {
          const int size = static_cast<int>(offsets[static_cast<std::size_t>(cell)].end -
              offsets[static_cast<std::size_t>(cell)].start);
          const int residual = std::max(0, size - allocation[static_cast<std::size_t>(cell)]);
          residual_by_cell[stratum].push_back(residual);
          total_residual += residual;
        }
        stratum_capacities[stratum] = total_residual;
        stratum_weights[stratum] = total_residual > 0 ? std::sqrt(static_cast<double>(total_residual)) : 0.0;
      }
      const auto extra_by_stratum =
          weighted_capacity_allocation(stratum_weights, stratum_capacities, remaining, rng);
      for (std::size_t stratum = 0; stratum < cells_by_stratum.size(); ++stratum) {
        if (extra_by_stratum[stratum] <= 0) continue;
        std::vector<double> cell_weights(cells_by_stratum[stratum].size(), 0.0);
        std::vector<int> cell_capacities(cells_by_stratum[stratum].size(), 0);
        for (std::size_t local = 0; local < cells_by_stratum[stratum].size(); ++local) {
          const int residual = residual_by_cell[stratum][local];
          cell_capacities[local] = residual;
          cell_weights[local] = residual > 0 ? std::sqrt(static_cast<double>(residual)) : 0.0;
        }
        const auto extra_by_cell =
            weighted_capacity_allocation(cell_weights, cell_capacities, extra_by_stratum[stratum], rng);
        for (std::size_t local = 0; local < cells_by_stratum[stratum].size(); ++local) {
          allocation[static_cast<std::size_t>(cells_by_stratum[stratum][local])] += extra_by_cell[local];
        }
      }
    }
  }

  TrainingSelection out;
  for (std::size_t cell = 0; cell < offsets.size(); ++cell) {
    const int take = allocation[cell];
    if (take <= 0) continue;
    const auto start = offsets[cell].start;
    const auto end = offsets[cell].end;
    std::vector<int> members;
    members.reserve(static_cast<std::size_t>(end - start));
    for (std::int64_t row = start; row < end; ++row) {
      members.push_back(static_cast<int>(row));
    }
    std::shuffle(members.begin(), members.end(), rng);
    if (take < static_cast<int>(members.size())) {
      members.resize(static_cast<std::size_t>(take));
    }
    for (const int row : members) {
      out.global_rows.push_back(row);
      out.cell_indices.push_back(static_cast<int>(cell));
    }
  }

  std::vector<std::size_t> order(out.global_rows.size());
  std::iota(order.begin(), order.end(), 0U);
  std::sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
    return out.global_rows[lhs] < out.global_rows[rhs];
  });
  auto rows = out.global_rows;
  auto cells = out.cell_indices;
  for (std::size_t i = 0; i < order.size(); ++i) {
    out.global_rows[i] = rows[order[i]];
    out.cell_indices[i] = cells[order[i]];
  }
  return out;
}

TrainingSelection sample_training_queries_from_row_lists(
    const std::vector<std::vector<int>>& rows_by_cell,
    const std::vector<std::string>& cell_strata,
    int max_rows,
    unsigned int seed,
    int min_molecules) {
  std::mt19937 rng(seed);
  const bool has_explicit_subset =
      !cell_strata.empty() &&
      std::any_of(cell_strata.begin(), cell_strata.end(), [](const std::string& value) {
        return value.empty();
      });

  std::vector<int> eligible_cells;
  for (std::size_t cell = 0; cell < rows_by_cell.size(); ++cell) {
    const int size = static_cast<int>(rows_by_cell[cell].size());
    if (size <= 0) continue;
    if (has_explicit_subset && cell_strata[cell].empty()) continue;
    if (size >= std::max(min_molecules, 1)) eligible_cells.push_back(static_cast<int>(cell));
  }
  if (eligible_cells.empty()) {
    if (has_explicit_subset) {
      throw std::runtime_error("explicit training strata did not label any cells above the training size threshold");
    }
    for (std::size_t cell = 0; cell < rows_by_cell.size(); ++cell) {
      if (!rows_by_cell[cell].empty()) eligible_cells.push_back(static_cast<int>(cell));
    }
  }

  int eligible_rows = 0;
  for (const int cell : eligible_cells) {
    eligible_rows += static_cast<int>(rows_by_cell[static_cast<std::size_t>(cell)].size());
  }
  if (eligible_rows == 0) return {};

  std::vector<int> allocation(rows_by_cell.size(), 0);
  if (max_rows <= 0 || max_rows >= eligible_rows) {
    for (const int cell : eligible_cells) {
      allocation[static_cast<std::size_t>(cell)] =
          static_cast<int>(rows_by_cell[static_cast<std::size_t>(cell)].size());
    }
  } else if (max_rows < static_cast<int>(eligible_cells.size())) {
    auto candidate_cells = eligible_cells;
    std::shuffle(candidate_cells.begin(), candidate_cells.end(), rng);
    candidate_cells.resize(static_cast<std::size_t>(max_rows));
    for (const int cell : candidate_cells) allocation[static_cast<std::size_t>(cell)] = 1;
  } else {
    std::vector<std::string> strata_order;
    std::unordered_map<std::string, int> stratum_lookup;
    std::vector<std::vector<int>> cells_by_stratum;
    for (const int cell : eligible_cells) {
      const std::string key = cell_strata[static_cast<std::size_t>(cell)];
      const auto it = stratum_lookup.find(key);
      if (it == stratum_lookup.end()) {
        const int stratum = static_cast<int>(strata_order.size());
        stratum_lookup.emplace(key, stratum);
        strata_order.push_back(key);
        cells_by_stratum.push_back({cell});
      } else {
        cells_by_stratum[static_cast<std::size_t>(it->second)].push_back(cell);
      }
    }

    int remaining = max_rows;
    for (const int cell : eligible_cells) {
      allocation[static_cast<std::size_t>(cell)] = 1;
      remaining -= 1;
      if (remaining == 0) break;
    }
    if (remaining > 0) {
      std::vector<double> stratum_weights(cells_by_stratum.size(), 0.0);
      std::vector<int> stratum_capacities(cells_by_stratum.size(), 0);
      std::vector<std::vector<int>> residual_by_cell(cells_by_stratum.size());
      for (std::size_t stratum = 0; stratum < cells_by_stratum.size(); ++stratum) {
        int total_residual = 0;
        for (const int cell : cells_by_stratum[stratum]) {
          const int size = static_cast<int>(rows_by_cell[static_cast<std::size_t>(cell)].size());
          const int residual = std::max(0, size - allocation[static_cast<std::size_t>(cell)]);
          residual_by_cell[stratum].push_back(residual);
          total_residual += residual;
        }
        stratum_capacities[stratum] = total_residual;
        stratum_weights[stratum] = total_residual > 0 ? std::sqrt(static_cast<double>(total_residual)) : 0.0;
      }
      const auto extra_by_stratum =
          weighted_capacity_allocation(stratum_weights, stratum_capacities, remaining, rng);
      for (std::size_t stratum = 0; stratum < cells_by_stratum.size(); ++stratum) {
        if (extra_by_stratum[stratum] <= 0) continue;
        std::vector<double> cell_weights(cells_by_stratum[stratum].size(), 0.0);
        std::vector<int> cell_capacities(cells_by_stratum[stratum].size(), 0);
        for (std::size_t local = 0; local < cells_by_stratum[stratum].size(); ++local) {
          const int residual = residual_by_cell[stratum][local];
          cell_capacities[local] = residual;
          cell_weights[local] = residual > 0 ? std::sqrt(static_cast<double>(residual)) : 0.0;
        }
        const auto extra_by_cell =
            weighted_capacity_allocation(cell_weights, cell_capacities, extra_by_stratum[stratum], rng);
        for (std::size_t local = 0; local < cells_by_stratum[stratum].size(); ++local) {
          allocation[static_cast<std::size_t>(cells_by_stratum[stratum][local])] += extra_by_cell[local];
        }
      }
    }
  }

  TrainingSelection out;
  for (std::size_t cell = 0; cell < rows_by_cell.size(); ++cell) {
    const int take = allocation[cell];
    if (take <= 0) continue;
    auto members = rows_by_cell[cell];
    std::shuffle(members.begin(), members.end(), rng);
    if (take < static_cast<int>(members.size())) members.resize(static_cast<std::size_t>(take));
    for (const int row : members) {
      out.global_rows.push_back(row);
      out.cell_indices.push_back(static_cast<int>(cell));
    }
  }

  std::vector<std::size_t> order(out.global_rows.size());
  std::iota(order.begin(), order.end(), 0U);
  std::sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
    return out.global_rows[lhs] < out.global_rows[rhs];
  });
  auto rows = out.global_rows;
  auto cells = out.cell_indices;
  for (std::size_t i = 0; i < order.size(); ++i) {
    out.global_rows[i] = rows[order[i]];
    out.cell_indices[i] = cells[order[i]];
  }
  return out;
}

std::vector<int> build_training_init_groups(
    const TrainingSelection& selection,
    const std::vector<std::string>& cell_strata,
    const std::string& nmf_init) {
  if (selection.global_rows.empty()) return {};
  if (nmf_init != "auto" && nmf_init != "random" && nmf_init != "cluster") {
    throw std::runtime_error("nmf_init must be one of auto, random, or cluster");
  }
  if (nmf_init == "random") return {};

  std::unordered_map<std::string, int> group_by_label;
  std::vector<int> out(selection.global_rows.size(), -1);
  int next_group = 0;
  for (std::size_t i = 0; i < selection.global_rows.size(); ++i) {
    const std::string& label = cell_strata[static_cast<std::size_t>(selection.cell_indices[i])];
    if (label.empty() || label == "__all_cells__") continue;
    const auto [it, inserted] = group_by_label.emplace(label, next_group);
    if (inserted) ++next_group;
    out[i] = it->second;
  }
  return next_group >= 2 ? out : std::vector<int>{};
}

struct LocalTrainingTable {
  TranscriptTable table;
  std::vector<int> query_indices;
};

LocalTrainingTable load_training_table(
    const std::string& store_dir,
    const std::vector<std::string>& genes,
    const CellTable& cells,
    const std::vector<InputStoreCellOffset>& offsets,
    const TrainingSelection& selection) {
  std::vector<char> selected_cell(cells.size(), 0);
  std::vector<int> global_to_local_cell(cells.size(), -1);
  for (const int cell : selection.cell_indices) {
    selected_cell[static_cast<std::size_t>(cell)] = 1;
  }
  LocalTrainingTable out;
  out.table.genes = genes;
  for (std::size_t cell = 0; cell < selected_cell.size(); ++cell) {
    if (!selected_cell[cell]) continue;
    global_to_local_cell[cell] = static_cast<int>(out.table.cells.size());
    out.table.cells.push_back(cells.cell_ids[cell]);
  }

  std::unordered_map<int, int> global_row_to_local;
  global_row_to_local.reserve(selection.global_rows.size());
  constexpr std::int64_t max_span = 512 * 1024;
  constexpr std::int64_t max_gap = 32 * 1024;
  std::size_t pos = 0;
  while (pos < selection.cell_indices.size()) {
    const int first_cell = selection.cell_indices[pos];
    std::int64_t range_start = offsets[static_cast<std::size_t>(first_cell)].start;
    std::int64_t range_end = offsets[static_cast<std::size_t>(first_cell)].end;
    std::size_t next = pos + 1;
    while (next < selection.cell_indices.size()) {
      const int cell = selection.cell_indices[next];
      const std::int64_t candidate_start = offsets[static_cast<std::size_t>(cell)].start;
      const std::int64_t candidate_end = offsets[static_cast<std::size_t>(cell)].end;
      if (candidate_start - range_end > max_gap) break;
      if (candidate_end - range_start > max_span) break;
      range_end = std::max(range_end, candidate_end);
      ++next;
    }

    const auto block = load_input_store_molecule_range(store_dir, range_start, range_end);
    for (std::size_t row = 0; row < block.size(); ++row) {
      const int global_cell = block.cell_idx[row];
      if (global_cell < 0 ||
          global_cell >= static_cast<int>(selected_cell.size()) ||
          !selected_cell[static_cast<std::size_t>(global_cell)]) {
        continue;
      }
      const int local_cell = global_to_local_cell[static_cast<std::size_t>(global_cell)];
      const int local_row = static_cast<int>(out.table.x.size());
      out.table.x.push_back(block.x[row]);
      out.table.y.push_back(block.y[row]);
      out.table.z.push_back(block.z[row]);
      out.table.gene_index.push_back(block.gene_idx[row]);
      out.table.cell_index.push_back(local_cell);
      global_row_to_local.emplace(static_cast<int>(block.row_index[row]), local_row);
    }
    pos = next;
  }

  out.query_indices.reserve(selection.global_rows.size());
  for (const int global : selection.global_rows) {
    const auto it = global_row_to_local.find(global);
    if (it == global_row_to_local.end()) {
      throw std::runtime_error("sampled training row was not found in loaded store molecules");
    }
    out.query_indices.push_back(it->second);
  }
  return out;
}

LocalTrainingTable load_training_table_from_row_lists(
    const std::string& store_dir,
    const std::vector<std::string>& genes,
    const CellTable& cells,
    const std::vector<int>& local_cell_by_global,
    const std::vector<std::vector<int>>& rows_by_cell,
    const TrainingSelection& selection) {
  std::vector<char> selected_cell(cells.size(), 0);
  std::vector<int> scope_to_training_cell(cells.size(), -1);
  for (const int cell : selection.cell_indices) {
    selected_cell[static_cast<std::size_t>(cell)] = 1;
  }

  LocalTrainingTable out;
  out.table.genes = genes;
  for (std::size_t cell = 0; cell < selected_cell.size(); ++cell) {
    if (!selected_cell[cell]) continue;
    scope_to_training_cell[cell] = static_cast<int>(out.table.cells.size());
    out.table.cells.push_back(cells.cell_ids[cell]);
  }

  std::vector<int> wanted_rows;
  for (std::size_t cell = 0; cell < rows_by_cell.size(); ++cell) {
    if (!selected_cell[cell]) continue;
    wanted_rows.insert(wanted_rows.end(), rows_by_cell[cell].begin(), rows_by_cell[cell].end());
  }
  std::sort(wanted_rows.begin(), wanted_rows.end());
  wanted_rows.erase(std::unique(wanted_rows.begin(), wanted_rows.end()), wanted_rows.end());
  std::unordered_set<int> wanted_lookup(wanted_rows.begin(), wanted_rows.end());
  std::unordered_map<int, int> global_row_to_local;
  global_row_to_local.reserve(selection.global_rows.size());

  constexpr std::int64_t max_span = 512 * 1024;
  constexpr std::int64_t max_gap = 32 * 1024;
  std::size_t pos = 0;
  while (pos < wanted_rows.size()) {
    std::int64_t range_start = wanted_rows[pos];
    std::int64_t range_end = range_start + 1;
    std::size_t next = pos + 1;
    while (next < wanted_rows.size()) {
      const std::int64_t candidate = wanted_rows[next];
      if (candidate - range_end > max_gap) break;
      if (candidate + 1 - range_start > max_span) break;
      range_end = candidate + 1;
      ++next;
    }

    const auto block = load_input_store_molecule_range(store_dir, range_start, range_end);
    for (std::size_t row = 0; row < block.size(); ++row) {
      const int global_row = static_cast<int>(block.row_index[row]);
      if (wanted_lookup.find(global_row) == wanted_lookup.end()) continue;
      const int global_cell = block.cell_idx[row];
      if (global_cell < 0 || global_cell >= static_cast<int>(local_cell_by_global.size())) {
        throw std::runtime_error("input-store molecule cell_idx is out of range");
      }
      const int scope_cell = local_cell_by_global[static_cast<std::size_t>(global_cell)];
      if (scope_cell < 0 || !selected_cell[static_cast<std::size_t>(scope_cell)]) continue;
      const int training_cell = scope_to_training_cell[static_cast<std::size_t>(scope_cell)];
      const int local_row = static_cast<int>(out.table.x.size());
      out.table.x.push_back(block.x[row]);
      out.table.y.push_back(block.y[row]);
      out.table.z.push_back(block.z[row]);
      out.table.gene_index.push_back(block.gene_idx[row]);
      out.table.cell_index.push_back(training_cell);
      global_row_to_local.emplace(global_row, local_row);
    }
    pos = next;
  }

  out.query_indices.reserve(selection.global_rows.size());
  for (const int global : selection.global_rows) {
    const auto it = global_row_to_local.find(global);
    if (it == global_row_to_local.end()) {
      throw std::runtime_error("sampled training row was not found in loaded crop molecules");
    }
    out.query_indices.push_back(it->second);
  }
  return out;
}

InputStoreMoleculeBlock filter_block_for_scope(
    const InputStoreMoleculeBlock& block,
    const std::optional<std::string>& analysis_crop,
    const std::vector<int>& local_cell_by_global) {
  if (!analysis_crop.has_value()) {
    return block;
  }
  if (block.crop_id.size() != block.size()) {
    throw std::runtime_error("analysis_crop requires crop_id metadata in the input store");
  }
  InputStoreMoleculeBlock out;
  out.row_index.reserve(block.size());
  out.obs_id.reserve(block.size());
  out.cell_idx.reserve(block.size());
  out.gene_idx.reserve(block.size());
  out.x.reserve(block.size());
  out.y.reserve(block.size());
  out.z.reserve(block.size());
  if (!block.qv.empty()) out.qv.reserve(block.size());
  if (!block.overlaps_nucleus.empty()) out.overlaps_nucleus.reserve(block.size());
  if (!block.nucleus_distance.empty()) out.nucleus_distance.reserve(block.size());
  out.crop_id.reserve(block.size());
  for (std::size_t row = 0; row < block.size(); ++row) {
    if (block.crop_id[row] != *analysis_crop) continue;
    const int global_cell = block.cell_idx[row];
    if (global_cell < 0 || global_cell >= static_cast<int>(local_cell_by_global.size())) {
      throw std::runtime_error("input-store molecule cell_idx is out of range");
    }
    const int local_cell = local_cell_by_global[static_cast<std::size_t>(global_cell)];
    if (local_cell < 0) continue;
    out.row_index.push_back(block.row_index[row]);
    out.obs_id.push_back(block.obs_id[row]);
    out.cell_idx.push_back(local_cell);
    out.gene_idx.push_back(block.gene_idx[row]);
    out.x.push_back(block.x[row]);
    out.y.push_back(block.y[row]);
    out.z.push_back(block.z[row]);
    if (!block.qv.empty()) out.qv.push_back(block.qv[row]);
    if (!block.overlaps_nucleus.empty()) out.overlaps_nucleus.push_back(block.overlaps_nucleus[row]);
    if (!block.nucleus_distance.empty()) out.nucleus_distance.push_back(block.nucleus_distance[row]);
    out.crop_id.push_back(block.crop_id[row]);
  }
  return out;
}

TranscriptTable make_local_table_for_scoped_block(
    const InputStoreMoleculeBlock& block,
    const std::vector<std::string>& genes,
    const CellTable& cells) {
  TranscriptTable out;
  out.genes = genes;
  std::unordered_map<int, int> block_cell_by_scope_cell;
  out.x = block.x;
  out.y = block.y;
  out.z = block.z;
  out.qv = block.qv;
  out.overlaps_nucleus = block.overlaps_nucleus;
  out.nucleus_distance = block.nucleus_distance;
  out.gene_index = block.gene_idx;
  out.cell_index.reserve(block.cell_idx.size());
  for (const int scope_cell : block.cell_idx) {
    auto it = block_cell_by_scope_cell.find(scope_cell);
    if (it == block_cell_by_scope_cell.end()) {
      const int block_cell = static_cast<int>(out.cells.size());
      it = block_cell_by_scope_cell.emplace(scope_cell, block_cell).first;
      out.cells.push_back(cells.cell_ids[static_cast<std::size_t>(scope_cell)]);
    }
    out.cell_index.push_back(it->second);
  }
  return out;
}

std::vector<double> compute_factor_margin(const DenseMatrix& scores) {
  std::vector<double> out(static_cast<std::size_t>(scores.rows()), 0.0);
  for (int row = 0; row < scores.rows(); ++row) {
    double best = -std::numeric_limits<double>::infinity();
    double second = -std::numeric_limits<double>::infinity();
    for (int col = 0; col < scores.cols(); ++col) {
      const double value = scores(row, col);
      if (value >= best) {
        second = best;
        best = value;
      } else if (value > second) {
        second = value;
      }
    }
    if (!std::isfinite(second)) second = 0.0;
    out[static_cast<std::size_t>(row)] = best - second;
  }
  return out;
}

int tile_index(double value, double tile_size) {
  return static_cast<int>(std::floor(value / tile_size));
}

class MoleculeParquetWriter {
 public:
  MoleculeParquetWriter(
      const std::filesystem::path& path,
      bool has_qv,
      bool has_overlaps_nucleus,
      bool has_nucleus_distance,
      int row_group_size)
      : has_qv_(has_qv),
        has_overlaps_nucleus_(has_overlaps_nucleus),
        has_nucleus_distance_(has_nucleus_distance),
        row_group_size_(std::max(row_group_size, 1)) {
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
        arrow::field("factor_margin", arrow::float64())};
    if (has_qv_) fields.push_back(arrow::field("qv", arrow::float64()));
    if (has_overlaps_nucleus_) fields.push_back(arrow::field("overlaps_nucleus", arrow::int32()));
    if (has_nucleus_distance_) fields.push_back(arrow::field("nucleus_distance", arrow::float64()));
    schema_ = std::make_shared<arrow::Schema>(fields);
    sink_ = arrow_unwrap(arrow::io::FileOutputStream::Open(path.string()), "Open molecule parquet output");
    writer_ = arrow_unwrap(
        parquet::arrow::FileWriter::Open(*schema_, arrow::default_memory_pool(), sink_),
        "Open molecule parquet writer");
  }

  void write_block(
      const InputStoreMoleculeBlock& block,
      const std::vector<int>& labels,
      const std::vector<double>& margin,
      const std::unordered_map<std::string, int>& crop_lookup,
      double tile_size) {
    if (labels.size() != block.size() || margin.size() != block.size()) {
      throw std::runtime_error("block labels/margins must match molecule block size");
    }
    std::vector<int> crop_idx;
    std::vector<int> tile_x;
    std::vector<int> tile_y;
    std::vector<int> tile_z;
    crop_idx.reserve(block.size());
    tile_x.reserve(block.size());
    tile_y.reserve(block.size());
    tile_z.reserve(block.size());
    for (std::size_t i = 0; i < block.size(); ++i) {
      int crop = -1;
      if (!block.crop_id.empty()) {
        const auto it = crop_lookup.find(block.crop_id[i]);
        if (it != crop_lookup.end()) crop = it->second;
      }
      crop_idx.push_back(crop);
      tile_x.push_back(tile_index(block.x[i], tile_size));
      tile_y.push_back(tile_index(block.y[i], tile_size));
      tile_z.push_back(tile_index(block.z[i], tile_size));
    }
    std::vector<std::shared_ptr<arrow::Array>> arrays = {
        build_int64_array(block.row_index),
        build_int32_array(crop_idx),
        build_int32_array(tile_x),
        build_int32_array(tile_y),
        build_int32_array(tile_z),
        build_double_array(block.x),
        build_double_array(block.y),
        build_double_array(block.z),
        build_int32_array(block.gene_idx),
        build_int32_array(block.cell_idx),
        build_int32_array(labels),
        build_double_array(margin)};
    if (has_qv_) arrays.push_back(build_double_array(block.qv));
    if (has_overlaps_nucleus_) arrays.push_back(build_int32_array(block.overlaps_nucleus));
    if (has_nucleus_distance_) arrays.push_back(build_double_array(block.nucleus_distance));
    const auto table = arrow::Table::Make(schema_, arrays);
    arrow_check(writer_->WriteTable(*table, row_group_size_), "Write molecule parquet block");
  }

  void close() {
    if (writer_) {
      arrow_check(writer_->Close(), "Close molecule parquet writer");
      writer_.reset();
    }
    if (sink_) {
      arrow_check(sink_->Close(), "Close molecule parquet sink");
      sink_.reset();
    }
  }

  ~MoleculeParquetWriter() {
    try {
      close();
    } catch (...) {
    }
  }

 private:
  bool has_qv_ = false;
  bool has_overlaps_nucleus_ = false;
  bool has_nucleus_distance_ = false;
  int row_group_size_ = 65536;
  std::shared_ptr<arrow::Schema> schema_;
  std::shared_ptr<arrow::io::FileOutputStream> sink_;
  std::unique_ptr<parquet::arrow::FileWriter> writer_;
};

void write_factors_parquet(
    const std::filesystem::path& path,
    const std::vector<std::string>& genes,
    const DenseMatrix& h,
    int row_group_size) {
  std::vector<int> factor_id;
  std::vector<int> gene_idx;
  std::vector<std::string> gene;
  std::vector<double> loading;
  const std::size_t total = static_cast<std::size_t>(h.rows()) * static_cast<std::size_t>(h.cols());
  factor_id.reserve(total);
  gene_idx.reserve(total);
  gene.reserve(total);
  loading.reserve(total);
  for (int factor = 0; factor < h.rows(); ++factor) {
    for (int gene_i = 0; gene_i < h.cols(); ++gene_i) {
      factor_id.push_back(factor + 1);
      gene_idx.push_back(gene_i);
      gene.push_back(genes[static_cast<std::size_t>(gene_i)]);
      loading.push_back(h(factor, gene_i));
    }
  }
  write_parquet_table(
      arrow::Table::Make(
          arrow::schema({
              arrow::field("factor_id", arrow::int32()),
              arrow::field("gene_idx", arrow::int32()),
              arrow::field("gene", arrow::utf8()),
              arrow::field("loading", arrow::float64())}),
          {
              build_int32_array(factor_id),
              build_int32_array(gene_idx),
              build_string_array(gene),
              build_double_array(loading)}),
      path,
      row_group_size);
}

void write_training_rows_parquet(
    const std::filesystem::path& path,
    const std::vector<int>& training_query_indices,
    int row_group_size) {
  std::vector<int> training_rank;
  std::vector<std::int64_t> obs_id;
  training_rank.reserve(training_query_indices.size());
  obs_id.reserve(training_query_indices.size());
  for (std::size_t i = 0; i < training_query_indices.size(); ++i) {
    training_rank.push_back(static_cast<int>(i));
    obs_id.push_back(static_cast<std::int64_t>(training_query_indices[i]));
  }
  write_parquet_table(
      arrow::Table::Make(
          arrow::schema({
              arrow::field("training_rank", arrow::int32()),
              arrow::field("obs_id", arrow::int64())}),
          {build_int32_array(training_rank), build_int64_array(obs_id)}),
      path,
      row_group_size);
}

void write_cells_parquet(
    const std::filesystem::path& path,
    const CellTable& cells,
    const std::vector<int>& transcript_counts,
    const std::vector<int>& factor_counts,
    int rank,
    int row_group_size) {
  std::vector<int> cell_idx;
  std::vector<int> dominant_factor;
  std::vector<double> dominant_fraction;
  std::vector<std::vector<double>> factor_fraction(static_cast<std::size_t>(rank));
  cell_idx.reserve(cells.size());
  dominant_factor.reserve(cells.size());
  dominant_fraction.reserve(cells.size());
  for (auto& values : factor_fraction) values.reserve(cells.size());

  for (std::size_t cell = 0; cell < cells.size(); ++cell) {
    cell_idx.push_back(static_cast<int>(cell));
    const int total = transcript_counts[cell];
    int best_factor = -1;
    int best_count = -1;
    for (int factor = 0; factor < rank; ++factor) {
      const int count = factor_counts[cell * static_cast<std::size_t>(rank) + static_cast<std::size_t>(factor)];
      if (count > best_count) {
        best_count = count;
        best_factor = factor;
      }
      factor_fraction[static_cast<std::size_t>(factor)].push_back(
          total > 0 ? static_cast<double>(count) / static_cast<double>(total) : 0.0);
    }
    dominant_factor.push_back(best_factor >= 0 ? best_factor + 1 : -1);
    dominant_fraction.push_back(
        best_factor >= 0 && total > 0
        ? static_cast<double>(best_count) / static_cast<double>(total)
        : 0.0);
  }

  std::vector<std::shared_ptr<arrow::Field>> fields = {
      arrow::field("cell_idx", arrow::int32()),
      arrow::field("cell_id", arrow::utf8()),
      arrow::field("x", arrow::float64()),
      arrow::field("y", arrow::float64()),
      arrow::field("z", arrow::float64()),
      arrow::field("transcript_count", arrow::int32()),
      arrow::field("dominant_factor", arrow::int32()),
      arrow::field("dominant_fraction", arrow::float64())};
  std::vector<std::shared_ptr<arrow::Array>> arrays = {
      build_int32_array(cell_idx),
      build_string_array(cells.cell_ids),
      build_double_array(cells.centroid_x),
      build_double_array(cells.centroid_y),
      build_double_array(cells.centroid_z),
      build_int32_array(transcript_counts),
      build_int32_array(dominant_factor),
      build_double_array(dominant_fraction)};
  if (!cells.cell_types.empty()) {
    fields.push_back(arrow::field("cell_type", arrow::utf8()));
    arrays.push_back(build_string_array(cells.cell_types));
  }
  if (!cells.sample_ids.empty()) {
    fields.push_back(arrow::field("sample_id", arrow::utf8()));
    arrays.push_back(build_string_array(cells.sample_ids));
  }
  if (!cells.fov_ids.empty()) {
    fields.push_back(arrow::field("fov_id", arrow::utf8()));
    arrays.push_back(build_string_array(cells.fov_ids));
  }
  for (int factor = 0; factor < rank; ++factor) {
    fields.push_back(arrow::field("factor_" + std::to_string(factor + 1) + "_fraction", arrow::float64()));
    arrays.push_back(build_double_array(factor_fraction[static_cast<std::size_t>(factor)]));
  }
  write_parquet_table(
      arrow::Table::Make(std::make_shared<arrow::Schema>(fields), arrays),
      path,
      row_group_size);
}

}  // namespace

StorePipelineResult run_basic_pipeline_store(
    const std::string& store_dir,
    const RunSourceInfo& source,
    const BasicPipelineOptions& options,
    const RunStorageOptions& storage_options,
    const std::string& out_dir,
    const std::optional<std::string>& analysis_crop,
    bool report_ncv_umap,
    const CellCountMatrix* preloaded_counts,
    bool verbose) {
  const auto pipeline_start = std::chrono::steady_clock::now();
  auto stage_start = pipeline_start;

  const auto store_manifest = read_input_store_manifest(store_dir);
  if (!store_manifest.has_molecule_rows || !store_manifest.has_cell_offsets) {
    throw std::runtime_error("store-backed fit requires molecule rows and cell offsets");
  }
  CellCountMatrix owned_counts;
  if (preloaded_counts == nullptr) {
    owned_counts = load_input_store_counts(store_dir);
    preloaded_counts = &owned_counts;
  }
  const CellCountMatrix& full_counts = *preloaded_counts;
  const auto offsets = load_input_store_cell_offsets(store_dir);
  if (offsets.size() != full_counts.cells.size()) {
    throw std::runtime_error("input-store cell offsets do not match cells table");
  }
  stage_start = std::chrono::steady_clock::now();
  StoreFitScope scope = analysis_crop.has_value()
      ? build_crop_fit_scope(store_dir, full_counts, offsets, *analysis_crop)
      : build_full_fit_scope(full_counts);
  if (scope.n_molecules == 0 || scope.counts.cells.size() == 0) {
    throw std::runtime_error("No molecules found in the requested input-store analysis scope");
  }
  const CellCountMatrix& counts = scope.counts;
  auto cell_strata =
      resolve_scope_cell_strata(scope, options.training_cell_strata, full_counts.cells.size());
  const std::size_t training_scope_cells = filter_training_strata_to_cell_types(
      counts.cells,
      cell_strata,
      options.training_scope_cell_types);
  if (verbose) {
    std::ostringstream message;
    message << "Loaded store metadata: " << scope.n_molecules
            << " molecules, " << counts.cells.size()
            << " cells, " << counts.genes.size() << " genes";
    if (analysis_crop.has_value()) message << " [analysis_crop=" << *analysis_crop << "]";
    if (!options.training_scope_cell_types.empty()) {
      message << " [training_scope_cell_types="
              << join_strings(unique_nonempty(options.training_scope_cell_types), ",")
              << "; selected_cells=" << training_scope_cells << "]";
    }
    emit_info(
        pipeline_start,
        message.str(),
        std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());
  }

  StorePipelineResult result;
  result.manifest.paths = make_run_paths(out_dir);
  ensure_directory(result.manifest.paths.root_dir);
  ensure_directory(result.manifest.paths.scores_dir);
  ensure_directory(result.manifest.paths.corrected_dir);

  NcvOptions ncv_options;
  ncv_options.k = options.ncv_k + 1;
  ncv_options.include_self = true;
  ncv_options.within_cell = true;

  stage_start = std::chrono::steady_clock::now();
  const int training_min_molecules = std::max(options.nmf_min_molecules, options.ncv_k + 1);
  auto selection = analysis_crop.has_value()
      ? sample_training_queries_from_row_lists(
            scope.crop_rows_by_cell,
            cell_strata,
            options.nmf_train_max_rows,
            options.seed,
            training_min_molecules)
      : sample_training_queries_from_offsets(
            offsets,
            cell_strata,
            options.nmf_train_max_rows,
            options.seed,
            training_min_molecules);
  if (selection.global_rows.empty()) {
    throw std::runtime_error("No training queries selected from input store; check crop size and nmf_min_molecules/ncv_k");
  }
  result.timing.training_query_sampling_sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count();
  if (verbose) {
    emit_info(
        pipeline_start,
        "Selected store-backed training queries: " + std::to_string(selection.global_rows.size()) + " rows",
        result.timing.training_query_sampling_sec);
  }

  stage_start = std::chrono::steady_clock::now();
  auto training_table = analysis_crop.has_value()
      ? load_training_table_from_row_lists(
            store_dir,
            counts.genes,
            counts.cells,
            scope.local_cell_by_global,
            scope.crop_rows_by_cell,
            selection)
      : load_training_table(store_dir, counts.genes, counts.cells, offsets, selection);
  std::vector<int> training_local_query_indices = training_table.query_indices;
  NcvOptions train_ncv_options = ncv_options;
  train_ncv_options.query_indices = training_local_query_indices;
  const auto raw_training_ncv = build_sparse_ncv_matrix(training_table.table, train_ncv_options);
  const auto compact_training = compact_training_ncv(raw_training_ncv);
  if (compact_training.matrix.rows() == 0 || compact_training.matrix.cols() == 0) {
    throw std::runtime_error("No usable sampled NCV rows/columns remained after training-matrix filtering");
  }
  if (!compact_training.kept_rows.empty() &&
      compact_training.kept_rows.size() != selection.global_rows.size()) {
    selection.global_rows = subset_int_vector(selection.global_rows, compact_training.kept_rows);
    selection.cell_indices = subset_int_vector(selection.cell_indices, compact_training.kept_rows);
    training_local_query_indices = subset_int_vector(training_local_query_indices, compact_training.kept_rows);
  }
  result.timing.training_ncv_sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count();
  if (verbose) {
    emit_info(pipeline_start, "Built store-backed sampled NCV matrix", result.timing.training_ncv_sec);
  }
  const auto nmf_init_groups = build_training_init_groups(selection, cell_strata, options.nmf_init);
  stage_start = std::chrono::steady_clock::now();
  NcvFeatureTransform compact_transform;
  if (is_weighted_ls_variant(options.nmf_variant)) {
    const auto column_weights = default_column_weights(compact_training.matrix);
    WeightedNmfOptions nmf_options;
    nmf_options.rank = options.rank;
    nmf_options.max_iterations = options.nmf_iterations;
    nmf_options.seed = options.seed;
    nmf_options.num_threads = options.num_threads;
    nmf_options.n_runs = options.nmf_n_runs;
    nmf_options.init_groups = nmf_init_groups;
    result.nmf = weighted_to_sparse_nmf_result(
        weighted_nmf(compact_training.matrix, column_weights, nmf_options));
    compact_transform.mode = options.nmf_variant;
    compact_transform.gene_weights = column_weights;
  } else {
    compact_transform = make_ncv_feature_transform(compact_training.matrix, options.nmf_variant);
    const SparseRowMatrix transformed_training_ncv =
        transform_sparse_ncv_matrix(compact_training.matrix, compact_transform);
    SparseNmfOptions nmf_options;
    nmf_options.rank = options.rank;
    nmf_options.max_iterations = options.nmf_iterations;
    nmf_options.seed = options.seed;
    nmf_options.num_threads = options.num_threads;
    nmf_options.n_runs = options.nmf_n_runs;
    nmf_options.loss_mode = "kl";
    nmf_options.init_groups = nmf_init_groups;
    result.nmf = sparse_nmf(transformed_training_ncv, {}, nmf_options);
    result.nmf.w = project_rows_to_factors_kl(transformed_training_ncv, result.nmf.h, options.num_threads);
  }
  order_nmf_factors_by_training_importance(result.nmf);
  NcvFeatureTransform full_transform = compact_transform;
  if (!compact_training.kept_cols.empty() &&
      static_cast<int>(compact_training.kept_cols.size()) != static_cast<int>(counts.genes.size())) {
    result.nmf.h = expand_h_to_full_genes(
        result.nmf.h,
        static_cast<int>(counts.genes.size()),
        compact_training.kept_cols);
    full_transform = expand_ncv_feature_transform(
        compact_transform,
        static_cast<int>(counts.genes.size()),
        compact_training.kept_cols);
  }
  result.timing.nmf_fit_sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count();
  if (verbose) {
    std::ostringstream message;
    message << (is_weighted_ls_variant(options.nmf_variant)
        ? "Fit store-backed weighted LS-NMF"
        : "Fit store-backed sparse KL NMF");
    if (options.nmf_variant != "kl" && !is_weighted_ls_variant(options.nmf_variant)) {
      message << " [" << options.nmf_variant << "]";
    }
    if (!nmf_init_groups.empty()) message << " [cluster init]";
    message << " [ordered by training W mass]";
    if (options.nmf_n_runs > 1) {
      message << " (selected run " << (result.nmf.selected_run + 1)
              << "/" << options.nmf_n_runs
              << ", seed=" << result.nmf.selected_seed
              << ", objective_mean=" << std::setprecision(6) << result.nmf.candidate_final_objective_mean
              << ", objective_sd=" << result.nmf.candidate_final_objective_sd
              << ", mean_best_component_cor=" << std::setprecision(3)
              << result.nmf.candidate_best_match_correlation_mean << ")";
    }
    emit_info(pipeline_start, message.str(), result.timing.nmf_fit_sec);
  }

  write_factors_parquet(
      result.manifest.paths.factors_parquet,
      counts.genes,
      result.nmf.h,
      storage_options.parquet_row_group_size);
  write_training_rows_parquet(
      result.manifest.paths.training_rows_parquet,
      selection.global_rows,
      storage_options.parquet_row_group_size);

  std::vector<std::string> crop_ids = analysis_crop.has_value()
      ? std::vector<std::string>{*analysis_crop}
      : unique_nonempty(counts.crop_ids);
  std::unordered_map<std::string, int> crop_lookup;
  for (std::size_t i = 0; i < crop_ids.size(); ++i) {
    crop_lookup.emplace(crop_ids[i], static_cast<int>(i));
  }
  const std::string molecule_scoring = resolve_molecule_scoring(options);

  MoleculeParquetWriter molecule_writer(
      result.manifest.paths.molecules_parquet,
      store_manifest.has_qv,
      store_manifest.has_overlaps_nucleus,
      store_manifest.has_nucleus_distance,
      storage_options.parquet_row_group_size);
  std::vector<int> transcript_counts(counts.cells.size(), 0);
  std::vector<int> factor_counts(counts.cells.size() * static_cast<std::size_t>(options.rank), 0);
  std::unordered_map<std::int64_t, int> training_rank_by_obs;
  std::vector<int> report_labels;
  DenseMatrix report_factor_scores;
  if (report_ncv_umap) {
    training_rank_by_obs.reserve(selection.global_rows.size());
    for (std::size_t rank = 0; rank < selection.global_rows.size(); ++rank) {
      training_rank_by_obs.emplace(static_cast<std::int64_t>(selection.global_rows[rank]), static_cast<int>(rank));
    }
    report_labels.assign(training_table.table.size(), 0);
    report_factor_scores = DenseMatrix(
        static_cast<int>(selection.global_rows.size()),
        result.nmf.h.rows(),
        0.0);
  }

  stage_start = std::chrono::steady_clock::now();
  double block_load_sec = 0.0;
  double block_projection_sec = 0.0;
  double block_smoothing_sec = 0.0;
  double block_write_sec = 0.0;
  constexpr std::int64_t target_block_rows = 512 * 1024;
  int cell = 0;
  while (cell < static_cast<int>(offsets.size())) {
    while (cell < static_cast<int>(offsets.size()) && offsets[static_cast<std::size_t>(cell)].start < 0) {
      ++cell;
    }
    if (cell >= static_cast<int>(offsets.size())) break;
    const int first_cell = cell;
    std::int64_t range_start = offsets[static_cast<std::size_t>(cell)].start;
    std::int64_t range_end = offsets[static_cast<std::size_t>(cell)].end;
    ++cell;
    while (cell < static_cast<int>(offsets.size())) {
      if (offsets[static_cast<std::size_t>(cell)].start < 0) {
        ++cell;
        continue;
      }
      const std::int64_t candidate_end = offsets[static_cast<std::size_t>(cell)].end;
      if (candidate_end - range_start > target_block_rows && range_end > range_start) {
        break;
      }
      range_end = candidate_end;
      ++cell;
    }
    const int last_cell = cell;
    auto block_stage_start = std::chrono::steady_clock::now();
    const auto raw_block = load_input_store_molecule_range(store_dir, range_start, range_end);
    block_load_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - block_stage_start).count();
    const auto block = filter_block_for_scope(raw_block, analysis_crop, scope.local_cell_by_global);
    if (block.size() == 0) continue;
    const auto local_table = make_local_table_for_scoped_block(block, counts.genes, counts.cells);
    block_stage_start = std::chrono::steady_clock::now();
    DenseMatrix scores = molecule_scoring == "gene_loadings"
        ? project_gene_loadings_to_factors(
            local_table,
            result.nmf.h,
            ncv_options,
            options.num_threads)
        : is_weighted_ls_variant(options.nmf_variant)
            ? project_ncv_to_factors(
                local_table,
                result.nmf.h,
                ncv_options,
                options.num_threads)
            : project_ncv_to_factors_kl(
            local_table,
            result.nmf.h,
            ncv_options,
            options.num_threads,
            10,
            1e-4,
            1e-10,
            full_transform);
    block_projection_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - block_stage_start).count();
    block_stage_start = std::chrono::steady_clock::now();
    const auto labels = assign_factors_per_cell(
        local_table,
        scores,
        options.graph_k,
        options.same_label_ratio,
        20,
        options.num_threads);
    block_smoothing_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - block_stage_start).count();
    const auto margins = compute_factor_margin(scores);
    for (std::size_t row = 0; row < block.size(); ++row) {
      const int global_cell = block.cell_idx[row];
      const int label = labels[row];
      transcript_counts[static_cast<std::size_t>(global_cell)] += 1;
      if (label >= 0 && label < options.rank) {
        factor_counts[static_cast<std::size_t>(global_cell) * static_cast<std::size_t>(options.rank) +
            static_cast<std::size_t>(label)] += 1;
      }
      if (report_ncv_umap) {
        const auto rank_it = training_rank_by_obs.find(block.row_index[row]);
        if (rank_it != training_rank_by_obs.end()) {
          const int training_rank = rank_it->second;
          const int local_query = training_local_query_indices[static_cast<std::size_t>(training_rank)];
          report_labels[static_cast<std::size_t>(local_query)] = label;
          for (int factor = 0; factor < scores.cols(); ++factor) {
            report_factor_scores(training_rank, factor) = scores(static_cast<int>(row), factor);
          }
        }
      }
    }
    block_stage_start = std::chrono::steady_clock::now();
    molecule_writer.write_block(block, labels, margins, crop_lookup, storage_options.tile_size);
    block_write_sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - block_stage_start).count();
  }
  molecule_writer.close();
  const double block_total_sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count();
  result.timing.projection_sec = block_projection_sec;
  result.timing.label_smoothing_sec = block_smoothing_sec;
  if (verbose) {
    std::ostringstream message;
    message << "Processed molecule blocks"
            << " [molecule_scoring=" << molecule_scoring << "]"
            << " [load=" << std::fixed << std::setprecision(3) << block_load_sec
            << "s, project=" << block_projection_sec
            << "s, smooth=" << block_smoothing_sec
            << "s, write=" << block_write_sec << "s]";
    emit_info(pipeline_start, message.str(), block_total_sec);
  }

  write_cells_parquet(
      result.manifest.paths.cells_parquet,
      counts.cells,
      transcript_counts,
      factor_counts,
      options.rank,
      storage_options.parquet_row_group_size);

  if (report_ncv_umap) {
    stage_start = std::chrono::steady_clock::now();
    BasicPipelineResult report_fit;
    report_fit.nmf = result.nmf;
    report_fit.training_query_indices = training_local_query_indices;
    report_fit.factor_scores = std::move(report_factor_scores);
    report_fit.labels = std::move(report_labels);
    std::vector<std::int64_t> training_obs_ids;
    training_obs_ids.reserve(selection.global_rows.size());
    for (const int row : selection.global_rows) {
      training_obs_ids.push_back(static_cast<std::int64_t>(row));
    }
    write_training_molecule_umap_report(
        result.manifest.paths.root_dir,
        training_table.table,
        report_fit,
        options.ncv_k,
        /*umap_neighbors=*/15,
        /*umap_epochs=*/200,
        options.seed,
        /*normalization_scale=*/5000.0,
        /*pca_dims=*/30,
        options.num_threads,
        &training_obs_ids);
    if (verbose) {
      emit_info(
          pipeline_start,
          "Built sampled NCV UMAP report sidecar from store-backed training subset",
          std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count());
    }
  }

  result.manifest.source = source;
  result.manifest.pipeline_options = options;
  result.manifest.storage_options = storage_options;
  result.manifest.analysis_crop = analysis_crop;
  result.manifest.genes = counts.genes;
  result.manifest.nmf_gene_weights = full_transform.gene_weights;
  result.manifest.nmf_diagnostics.final_objective = result.nmf.final_objective;
  result.manifest.nmf_diagnostics.selected_seed = result.nmf.selected_seed;
  result.manifest.nmf_diagnostics.selected_run = result.nmf.selected_run;
  result.manifest.nmf_diagnostics.candidate_final_objectives = result.nmf.candidate_final_objectives;
  result.manifest.nmf_diagnostics.candidate_best_match_correlations = result.nmf.candidate_best_match_correlations;
  result.manifest.nmf_diagnostics.selected_factor_stability = result.nmf.selected_factor_stability;
  result.manifest.nmf_diagnostics.candidate_final_objective_mean = result.nmf.candidate_final_objective_mean;
  result.manifest.nmf_diagnostics.candidate_final_objective_sd = result.nmf.candidate_final_objective_sd;
  result.manifest.nmf_diagnostics.candidate_best_match_correlation_mean =
      result.nmf.candidate_best_match_correlation_mean;
  result.manifest.nmf_transform_target_row_sum = full_transform.target_row_sum;
  result.manifest.crop_ids = crop_ids;
  result.manifest.n_transcripts = scope.n_molecules;
  result.manifest.n_cells = counts.cells.size();
  result.manifest.n_factors = static_cast<std::size_t>(result.nmf.h.rows());
  result.manifest.n_training_rows = selection.global_rows.size();
  result.manifest.has_z = store_manifest.has_z;
  result.manifest.has_qv = store_manifest.has_qv;
  result.manifest.has_transcript_id = false;
  result.manifest.has_cell_type = !counts.cells.cell_types.empty();
  result.manifest.has_sample_id = !counts.cells.sample_ids.empty();
  result.manifest.has_fov_id = !counts.cells.fov_ids.empty();
  result.manifest.has_nucleus_id = false;
  result.manifest.has_overlaps_nucleus = store_manifest.has_overlaps_nucleus;
  result.manifest.has_nucleus_distance = store_manifest.has_nucleus_distance;
  write_manifest_json(result.manifest);
  result.timing.total_sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - pipeline_start).count();
  if (verbose) {
    emit_info(pipeline_start, "Finished cell-blocked store-backed fit", result.timing.total_sec);
  }
  return result;
}

}  // namespace celladmix
