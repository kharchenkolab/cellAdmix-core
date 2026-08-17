#include "celladmix/pipeline.hpp"
#include "celladmix/nmf_euclidean.hpp"
#include "celladmix/pipeline_common.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace celladmix {

// Core end-to-end pipeline orchestration from NCV sampling through labeling.

namespace {

// Hold the compacted sampled NCV matrix and maps back to original rows/genes.
struct CompressedTrainingNcv {
  SparseRowMatrix matrix;
  std::vector<int> kept_rows;
  std::vector<int> kept_cols;
};

// Drop empty sampled NCV rows/columns before NMF and record the retained indices.
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
  kept_rows.reserve(static_cast<std::size_t>(x.rows()));
  std::vector<int> indptr;
  std::vector<int> indices;
  std::vector<double> values;
  indptr.reserve(static_cast<std::size_t>(x.rows() + 1));
  indptr.push_back(0);

  for (int row = 0; row < x.rows(); ++row) {
    const int start = x.indptr()[static_cast<std::size_t>(row)];
    const int end = x.indptr()[static_cast<std::size_t>(row + 1)];
    const int row_begin = static_cast<int>(indices.size());
    for (int p = start; p < end; ++p) {
      const int original_col = x.indices()[static_cast<std::size_t>(p)];
      const int mapped_col = col_map[static_cast<std::size_t>(original_col)];
      if (mapped_col < 0) {
        continue;
      }
      indices.push_back(mapped_col);
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

// Project an integer vector through a retained-row index set.
std::vector<int> subset_int_vector(const std::vector<int>& values, const std::vector<int>& keep) {
  std::vector<int> out;
  out.reserve(keep.size());
  for (const int index : keep) {
    out.push_back(values[static_cast<std::size_t>(index)]);
  }
  return out;
}

// Emit one timestamped stage log line.
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

// Build the default query index vector covering all transcript rows.
std::vector<int> all_query_indices(std::size_t n) {
  std::vector<int> out(n);
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = static_cast<int>(i);
  }
  return out;
}

// Infer per-cell training strata from cell types when explicit strata were not provided.
std::vector<std::string> infer_training_cell_strata(
    const TranscriptTable& table,
    const std::vector<std::string>& provided) {
  if (!provided.empty()) {
    if (provided.size() != table.num_cells()) {
      throw std::runtime_error("training_cell_strata must match the number of cells");
    }
    return provided;
  }

  std::vector<std::string> strata(table.num_cells(), "__all_cells__");
  if (!table.has_cell_types()) {
    return strata;
  }

  strata.assign(table.num_cells(), "");
  for (std::size_t i = 0; i < table.size(); ++i) {
    const std::size_t cell = static_cast<std::size_t>(table.cell_index[i]);
    if (!strata[cell].empty()) {
      continue;
    }
    strata[cell] = table.cell_types[i].empty() ? "__unlabeled__" : table.cell_types[i];
  }
  for (auto& value : strata) {
    if (value.empty()) {
      value = "__unlabeled__";
    }
  }
  return strata;
}

// Map sampled training rows back to cluster or cell-type groups for KL initialization.
std::vector<int> build_training_init_groups(
    const TranscriptTable& table,
    const std::vector<int>& training_query_indices,
    const std::vector<std::string>& provided_cell_strata,
    const std::string& nmf_init) {
  if (training_query_indices.empty()) {
    return {};
  }
  if (nmf_init != "auto" && nmf_init != "random" && nmf_init != "cluster") {
    throw std::runtime_error("nmf_init must be one of auto, random, or cluster");
  }
  if (nmf_init == "random") {
    return {};
  }

  const auto cell_strata = infer_training_cell_strata(table, provided_cell_strata);
  std::unordered_map<std::string, int> group_by_label;
  std::vector<int> out(training_query_indices.size(), -1);
  int next_group = 0;
  for (std::size_t row = 0; row < training_query_indices.size(); ++row) {
    const int query = training_query_indices[row];
    const int cell = table.cell_index[static_cast<std::size_t>(query)];
    const std::string& label = cell_strata[static_cast<std::size_t>(cell)];
    if (label.empty() || label == "__all_cells__") {
      continue;
    }
    const auto [it, inserted] = group_by_label.emplace(label, next_group);
    if (inserted) {
      ++next_group;
    }
    out[row] = it->second;
  }
  return next_group >= 2 ? out : std::vector<int>{};
}

// Allocate an integer sampling budget proportionally while respecting capacities.
std::vector<int> weighted_capacity_allocation(
    const std::vector<double>& weights,
    const std::vector<int>& capacities,
    int budget,
    std::mt19937& rng) {
  if (weights.size() != capacities.size()) {
    throw std::runtime_error("weights and capacities must match");
  }

  std::vector<int> allocation(weights.size(), 0);
  if (budget <= 0 || weights.empty()) {
    return allocation;
  }

  std::vector<std::size_t> active;
  active.reserve(weights.size());
  for (std::size_t i = 0; i < capacities.size(); ++i) {
    if (capacities[i] > 0) {
      active.push_back(i);
    }
  }

  int remaining = budget;
  while (remaining > 0 && !active.empty()) {
    double total_weight = 0.0;
    for (const auto idx : active) {
      total_weight += std::max(weights[idx], 0.0);
    }
    if (total_weight <= 0.0) {
      total_weight = static_cast<double>(active.size());
    }

    std::vector<double> fractional(weights.size(), 0.0);
    int assigned = 0;
    for (const auto idx : active) {
      const int remaining_capacity = capacities[idx] - allocation[idx];
      if (remaining_capacity <= 0) {
        continue;
      }
      const double effective_weight =
          total_weight > 0.0 ? std::max(weights[idx], 0.0) : 1.0;
      const double expected =
          static_cast<double>(remaining) * effective_weight / total_weight;
      const int add = std::min<int>(
          remaining_capacity,
          static_cast<int>(std::floor(expected)));
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
      if (leftover == 0) {
        break;
      }
      if (capacities[idx] - allocation[idx] <= 0) {
        continue;
      }
      allocation[idx] += 1;
      assigned += 1;
      leftover -= 1;
    }

    if (assigned == 0) {
      break;
    }
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

// Sample training query molecules across labeled cells and strata for NMF fitting.
std::vector<int> sample_training_queries(
    const TranscriptTable& table,
    int max_rows,
    unsigned int seed,
    int min_molecules,
    const std::vector<std::string>& provided_cell_strata) {
  std::mt19937 rng(seed);
  const auto by_cell = table.transcripts_by_cell();
  const auto cell_strata = infer_training_cell_strata(table, provided_cell_strata);
  const bool has_explicit_subset =
      !provided_cell_strata.empty() &&
      std::any_of(
          cell_strata.begin(),
          cell_strata.end(),
          [](const std::string& value) { return value.empty(); });

  std::vector<int> nonempty_cells;
  nonempty_cells.reserve(by_cell.size());
  for (std::size_t cell = 0; cell < by_cell.size(); ++cell) {
    if (!by_cell[cell].empty()) {
      nonempty_cells.push_back(static_cast<int>(cell));
    }
  }
  if (nonempty_cells.empty()) {
    return {};
  }

  std::vector<int> eligible_cells;
  eligible_cells.reserve(nonempty_cells.size());
  const int effective_min_molecules = std::max(min_molecules, 1);
  for (const int cell : nonempty_cells) {
    if (has_explicit_subset && cell_strata[static_cast<std::size_t>(cell)].empty()) {
      continue;
    }
    if (static_cast<int>(by_cell[static_cast<std::size_t>(cell)].size()) >= effective_min_molecules) {
      eligible_cells.push_back(cell);
    }
  }
  if (eligible_cells.empty()) {
    if (has_explicit_subset) {
      throw std::runtime_error(
          "explicit training_cell_strata did not label any cells above the training size threshold");
    }
    eligible_cells = nonempty_cells;
  }

  int eligible_rows = 0;
  for (const int cell : eligible_cells) {
    eligible_rows += static_cast<int>(by_cell[static_cast<std::size_t>(cell)].size());
  }
  if (max_rows <= 0 || max_rows >= eligible_rows) {
    std::vector<int> sampled;
    sampled.reserve(static_cast<std::size_t>(eligible_rows));
    for (const int cell : eligible_cells) {
      const auto& members = by_cell[static_cast<std::size_t>(cell)];
      sampled.insert(sampled.end(), members.begin(), members.end());
    }
    std::sort(sampled.begin(), sampled.end());
    return sampled;
  }
  const int n_eligible_cells = static_cast<int>(eligible_cells.size());

  std::vector<std::string> strata_order;
  std::unordered_map<std::string, int> stratum_lookup;
  std::vector<std::vector<int>> cells_by_stratum;
  for (const int cell : eligible_cells) {
    const std::string key = cell_strata[static_cast<std::size_t>(cell)];
    const auto it = stratum_lookup.find(key);
    if (it == stratum_lookup.end()) {
      const int stratum_idx = static_cast<int>(strata_order.size());
      stratum_lookup.emplace(key, stratum_idx);
      strata_order.push_back(key);
      cells_by_stratum.push_back({cell});
    } else {
      cells_by_stratum[static_cast<std::size_t>(it->second)].push_back(cell);
    }
  }

  std::vector<int> allocation(by_cell.size(), 0);

  if (max_rows < n_eligible_cells) {
    std::vector<int> stratum_cell_quota(cells_by_stratum.size(), 0);
    std::vector<std::size_t> stratum_order(cells_by_stratum.size());
    std::iota(stratum_order.begin(), stratum_order.end(), 0U);
    std::shuffle(stratum_order.begin(), stratum_order.end(), rng);
    std::stable_sort(stratum_order.begin(), stratum_order.end(), [&](std::size_t lhs, std::size_t rhs) {
      return cells_by_stratum[lhs].size() < cells_by_stratum[rhs].size();
    });

    const int guaranteed = std::min<int>(max_rows, static_cast<int>(stratum_order.size()));
    for (int i = 0; i < guaranteed; ++i) {
      stratum_cell_quota[stratum_order[static_cast<std::size_t>(i)]] = 1;
    }

    int remaining_cells = max_rows - guaranteed;
    if (remaining_cells > 0) {
      std::vector<double> stratum_weights(cells_by_stratum.size(), 0.0);
      std::vector<int> stratum_capacities(cells_by_stratum.size(), 0);
      for (std::size_t stratum = 0; stratum < cells_by_stratum.size(); ++stratum) {
        const int remaining_capacity =
            static_cast<int>(cells_by_stratum[stratum].size()) - stratum_cell_quota[stratum];
        stratum_capacities[stratum] = std::max(remaining_capacity, 0);
        stratum_weights[stratum] =
            remaining_capacity > 0 ? std::sqrt(static_cast<double>(remaining_capacity)) : 0.0;
      }
      const auto extra_by_stratum = weighted_capacity_allocation(
          stratum_weights,
          stratum_capacities,
          remaining_cells,
          rng);
      for (std::size_t stratum = 0; stratum < cells_by_stratum.size(); ++stratum) {
        stratum_cell_quota[stratum] += extra_by_stratum[stratum];
      }
    }

    for (std::size_t stratum = 0; stratum < cells_by_stratum.size(); ++stratum) {
      auto candidates = cells_by_stratum[stratum];
      std::shuffle(candidates.begin(), candidates.end(), rng);
      candidates.resize(static_cast<std::size_t>(stratum_cell_quota[stratum]));
      for (const int cell : candidates) {
        allocation[static_cast<std::size_t>(cell)] = 1;
      }
    }
  } else {
    int remaining = max_rows;
    for (const int cell : eligible_cells) {
      allocation[static_cast<std::size_t>(cell)] = 1;
      remaining -= 1;
    }

    if (remaining > 0) {
      std::vector<double> stratum_weights(cells_by_stratum.size(), 0.0);
      std::vector<int> stratum_capacities(cells_by_stratum.size(), 0);
      std::vector<std::vector<int>> residual_by_cell(cells_by_stratum.size());
      for (std::size_t stratum = 0; stratum < cells_by_stratum.size(); ++stratum) {
        int total_residual = 0;
        residual_by_cell[stratum].reserve(cells_by_stratum[stratum].size());
        for (const int cell : cells_by_stratum[stratum]) {
          const int residual =
              std::max<int>(0, static_cast<int>(by_cell[static_cast<std::size_t>(cell)].size()) - 1);
          residual_by_cell[stratum].push_back(residual);
          total_residual += residual;
        }
        stratum_capacities[stratum] = total_residual;
        stratum_weights[stratum] =
            total_residual > 0 ? std::sqrt(static_cast<double>(total_residual)) : 0.0;
      }

      const auto extra_by_stratum = weighted_capacity_allocation(
          stratum_weights,
          stratum_capacities,
          remaining,
          rng);

      for (std::size_t stratum = 0; stratum < cells_by_stratum.size(); ++stratum) {
        if (extra_by_stratum[stratum] <= 0) {
          continue;
        }
        std::vector<double> cell_weights(cells_by_stratum[stratum].size(), 0.0);
        std::vector<int> cell_capacities(cells_by_stratum[stratum].size(), 0);
        for (std::size_t local = 0; local < cells_by_stratum[stratum].size(); ++local) {
          const int residual = residual_by_cell[stratum][local];
          cell_capacities[local] = residual;
          cell_weights[local] =
              residual > 0 ? std::sqrt(static_cast<double>(residual)) : 0.0;
        }
        const auto extra_by_cell = weighted_capacity_allocation(
            cell_weights,
            cell_capacities,
            extra_by_stratum[stratum],
            rng);
        for (std::size_t local = 0; local < cells_by_stratum[stratum].size(); ++local) {
          const int cell = cells_by_stratum[stratum][local];
          allocation[static_cast<std::size_t>(cell)] += extra_by_cell[local];
        }
      }
    }
  }

  std::vector<int> sampled;
  sampled.reserve(static_cast<std::size_t>(max_rows));
  for (std::size_t cell = 0; cell < by_cell.size(); ++cell) {
    auto members = by_cell[cell];
    if (members.empty() || allocation[cell] == 0) {
      continue;
    }
    std::shuffle(members.begin(), members.end(), rng);
    if (allocation[cell] < static_cast<int>(members.size())) {
      members.resize(static_cast<std::size_t>(allocation[cell]));
    }
    sampled.insert(sampled.end(), members.begin(), members.end());
  }

  if (static_cast<int>(sampled.size()) > max_rows) {
    std::shuffle(sampled.begin(), sampled.end(), rng);
    sampled.resize(static_cast<std::size_t>(max_rows));
  }
  std::sort(sampled.begin(), sampled.end());
  return sampled;
}

}  // namespace

// Run the full transcript-level factorization, projection, and smoothing pipeline.
BasicPipelineResult run_basic_pipeline(
    const TranscriptTable& table,
    const BasicPipelineOptions& input_options,
    bool verbose) {
  const auto pipeline_start = std::chrono::steady_clock::now();
  BasicPipelineOptions options = input_options;
  if (options.ncv_k <= 0) {
    std::vector<char> gene_present(table.num_genes(), 0);
    for (const int gene : table.gene_index) {
      if (gene >= 0) {
        gene_present[static_cast<std::size_t>(gene)] = 1;
      }
    }
    const int genes_present = static_cast<int>(
        std::count(gene_present.begin(), gene_present.end(), 1));
    std::vector<int> cell_totals(static_cast<std::size_t>(table.num_cells()), 0);
    for (const int cell : table.cell_index) {
      if (cell >= 0) {
        ++cell_totals[static_cast<std::size_t>(cell)];
      }
    }
    cell_totals.erase(
        std::remove(cell_totals.begin(), cell_totals.end(), 0),
        cell_totals.end());
    double median_cell = 0.0;
    if (!cell_totals.empty()) {
      const auto mid = cell_totals.begin() +
          static_cast<std::ptrdiff_t>(cell_totals.size() / 2);
      std::nth_element(cell_totals.begin(), mid, cell_totals.end());
      median_cell = static_cast<double>(*mid);
    }
    options.ncv_k = resolve_auto_ncv_k(genes_present, median_cell);
    if (verbose) {
      emit_info(
          pipeline_start,
          "Resolved ncv_k=" + std::to_string(options.ncv_k) + " automatically (" +
              std::to_string(genes_present) + " genes present, median " +
              std::to_string(static_cast<long long>(median_cell)) +
              " molecules/cell)",
          0.0);
    }
  }
  NcvOptions ncv_options;
  ncv_options.k = options.ncv_k + 1;
  ncv_options.include_self = true;
  ncv_options.within_cell = true;

  BasicPipelineResult result;
  result.resolved_options = options;
  auto stage_start = std::chrono::steady_clock::now();
  const int training_min_molecules = std::max(options.nmf_min_molecules, options.ncv_k + 1);
  result.training_query_indices =
      sample_training_queries(
          table,
          options.nmf_train_max_rows,
          options.seed,
          training_min_molecules,
          options.training_cell_strata);
  result.timing.training_query_sampling_sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count();
  if (verbose) {
    emit_info(
        pipeline_start,
        "Selected training queries: " +
            std::to_string(result.training_query_indices.size()) + " rows",
        result.timing.training_query_sampling_sec);
  }

  NcvOptions train_ncv_options = ncv_options;
  train_ncv_options.query_indices = result.training_query_indices;
  stage_start = std::chrono::steady_clock::now();
  const SparseRowMatrix raw_training_ncv = build_sparse_ncv_matrix(table, train_ncv_options);
  const auto compact_training = compact_training_ncv(raw_training_ncv);
  const SparseRowMatrix& training_ncv = compact_training.matrix;
  if (training_ncv.rows() == 0 || training_ncv.cols() == 0) {
    throw std::runtime_error("No usable sampled NCV rows/columns remained after training-matrix filtering");
  }
  if (!compact_training.kept_rows.empty() &&
      compact_training.kept_rows.size() != result.training_query_indices.size()) {
    result.training_query_indices = subset_int_vector(
        result.training_query_indices,
        compact_training.kept_rows);
    train_ncv_options.query_indices = result.training_query_indices;
  }
  result.timing.training_ncv_sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count();
  if (verbose) {
    emit_info(pipeline_start, "Built sampled NCV matrix", result.timing.training_ncv_sec);
  }
  const auto nmf_init_groups = build_training_init_groups(
      table,
      result.training_query_indices,
      options.training_cell_strata,
      options.nmf_init);
  stage_start = std::chrono::steady_clock::now();
  NcvFeatureTransform compact_transform;
  if (is_weighted_ls_variant(options.nmf_variant)) {
    const auto column_weights = default_column_weights(training_ncv);
    WeightedNmfOptions nmf_options;
    nmf_options.rank = options.rank;
    nmf_options.max_iterations = options.nmf_iterations;
    nmf_options.seed = options.seed;
    nmf_options.num_threads = options.num_threads;
    nmf_options.n_runs = options.nmf_n_runs;
    nmf_options.init_groups = nmf_init_groups;
    result.nmf = weighted_to_sparse_nmf_result(
        weighted_nmf(training_ncv, column_weights, nmf_options));
    compact_transform.mode = options.nmf_variant;
    compact_transform.gene_weights = column_weights;
  } else {
    compact_transform = make_ncv_feature_transform(training_ncv, options.nmf_variant);
    const SparseRowMatrix transformed_training_ncv =
        transform_sparse_ncv_matrix(training_ncv, compact_transform);
    SparseNmfOptions nmf_options;
    nmf_options.rank = options.rank;
    nmf_options.max_iterations = options.nmf_iterations;
    nmf_options.seed = options.seed;
    nmf_options.num_threads = options.num_threads;
    nmf_options.n_runs = options.nmf_n_runs;
    nmf_options.loss_mode = "kl";
    nmf_options.init_groups = nmf_init_groups;
    result.nmf = sparse_nmf(transformed_training_ncv, {}, nmf_options);
    result.nmf.w = project_rows_to_factors_kl(
        transformed_training_ncv,
        result.nmf.h,
        options.num_threads);
  }
  order_nmf_factors_by_training_importance(result.nmf);
  result.nmf_transform = compact_transform;
  if (!compact_training.kept_cols.empty() &&
      static_cast<int>(compact_training.kept_cols.size()) != static_cast<int>(table.num_genes())) {
    result.nmf.h = expand_h_to_full_genes(
        result.nmf.h,
        static_cast<int>(table.num_genes()),
        compact_training.kept_cols);
    for (auto& candidate : result.nmf.candidate_h) {
      candidate = expand_h_to_full_genes(
          candidate,
          static_cast<int>(table.num_genes()),
          compact_training.kept_cols);
    }
    result.nmf_transform = expand_ncv_feature_transform(
        compact_transform,
        static_cast<int>(table.num_genes()),
        compact_training.kept_cols);
  }
  result.timing.nmf_fit_sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count();
  if (verbose) {
    std::ostringstream nmf_message;
    nmf_message << (is_weighted_ls_variant(options.nmf_variant)
        ? "Fit weighted LS-NMF"
        : "Fit sparse KL NMF");
    if (options.nmf_variant != "kl" && !is_weighted_ls_variant(options.nmf_variant)) {
      nmf_message << " [" << options.nmf_variant << "]";
    }
    if (!nmf_init_groups.empty()) {
      nmf_message << " [cluster init]";
    }
    nmf_message << " [ordered by training W mass]";
    if (options.nmf_n_runs > 1) {
      nmf_message << " (selected run " << (result.nmf.selected_run + 1)
                  << "/" << options.nmf_n_runs
                  << ", seed=" << result.nmf.selected_seed
                  << ", objective_mean=" << std::setprecision(6) << result.nmf.candidate_final_objective_mean
                  << ", objective_sd=" << result.nmf.candidate_final_objective_sd
                  << ", mean_matched_ownership_cor=" << std::setprecision(3)
                  << result.nmf.candidate_best_match_correlation_mean
                  << ", stable_factors=" << result.nmf.stable_factor_count
                  << "/" << result.nmf.h.rows() << ")";
    }
    emit_info(pipeline_start, nmf_message.str(), result.timing.nmf_fit_sec);
  }

  stage_start = std::chrono::steady_clock::now();
  if (options.return_ncv) {
    result.ncv = build_ncv_matrix(table, ncv_options);
  }
  const std::string molecule_scoring = resolve_molecule_scoring(options);
  if (molecule_scoring == "gene_loadings") {
    result.factor_scores = project_gene_loadings_to_factors(
        table,
        result.nmf.h,
        ncv_options,
        options.num_threads);
  } else if (is_weighted_ls_variant(options.nmf_variant)) {
    result.factor_scores = project_ncv_to_factors(
        table,
        result.nmf.h,
        ncv_options,
        options.num_threads);
  } else {
    result.factor_scores = project_ncv_to_factors_kl(
        table,
        result.nmf.h,
        ncv_options,
        options.num_threads,
        10,
        1e-4,
        1e-10,
        result.nmf_transform);
  }
  result.timing.projection_sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count();
  if (verbose) {
    emit_info(
        pipeline_start,
        "Projected all transcripts to factors [" + molecule_scoring + "]",
        result.timing.projection_sec);
  }

  stage_start = std::chrono::steady_clock::now();
  result.labels = assign_factors_per_cell(
      table,
      result.factor_scores,
      options.graph_k,
      options.same_label_ratio,
      20,
      options.num_threads);
  result.timing.label_smoothing_sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start).count();
  if (verbose) {
    emit_info(
        pipeline_start,
        "Smoothed and assigned transcript labels",
        result.timing.label_smoothing_sec);
  }
  result.timing.total_sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - pipeline_start).count();
  return result;
}

// Build a keep mask that removes one factor from one target cell type.
std::vector<bool> apply_removal_rule(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    int factor,
    const std::string& target_cell_type) {
  if (!table.has_cell_types()) {
    throw std::runtime_error("TranscriptTable must contain cell_types for apply_removal_rule");
  }
  if (labels.size() != table.size()) {
    throw std::runtime_error("labels length must match TranscriptTable size");
  }

  std::vector<bool> keep(table.size(), true);
  for (std::size_t i = 0; i < table.size(); ++i) {
    if (labels[i] == factor && table.cell_types[i] == target_cell_type) {
      keep[i] = false;
    }
  }
  return keep;
}

// Aggregate transcript counts into a dense cell-by-gene matrix.
DenseMatrix build_cell_gene_counts(
    const TranscriptTable& table,
    const std::vector<bool>& keep_mask) {
  if (!keep_mask.empty() && keep_mask.size() != table.size()) {
    throw std::runtime_error("keep_mask length must match TranscriptTable size");
  }

  DenseMatrix counts(static_cast<int>(table.num_cells()), static_cast<int>(table.num_genes()), 0.0);
  for (std::size_t i = 0; i < table.size(); ++i) {
    if (!keep_mask.empty() && !keep_mask[i]) {
      continue;
    }
    counts(table.cell_index[i], table.gene_index[i]) += 1.0;
  }
  return counts;
}

// Pick the factor with the largest total loading on a supplied marker set.
int select_factor_by_marker_sum(
    const DenseMatrix& h,
    const std::vector<std::string>& genes,
    const std::vector<std::string>& marker_genes) {
  std::unordered_set<std::string> markers(marker_genes.begin(), marker_genes.end());
  int best_factor = 0;
  double best_sum = -1.0;
  for (int factor = 0; factor < h.rows(); ++factor) {
    double score = 0.0;
    for (int gene = 0; gene < h.cols(); ++gene) {
      if (markers.count(genes[static_cast<std::size_t>(gene)]) != 0U) {
        score += h(factor, gene);
      }
    }
    if (score > best_sum) {
      best_sum = score;
      best_factor = factor;
    }
  }
  return best_factor;
}

}  // namespace celladmix
