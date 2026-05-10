#include "celladmix/ncv.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "celladmix/spatial.hpp"
#include "subpar/range.hpp"

namespace celladmix {

// NCV construction and fixed-H projection helpers for transcript-level features.

namespace {

// Reusable sparse accumulation buffers for one NCV row.
struct NcvScratch {
  std::vector<int> marks;
  std::vector<int> touched;
  std::vector<double> counts;
  int current_mark = 0;

  explicit NcvScratch(int n_genes)
      : marks(static_cast<std::size_t>(n_genes), 0),
        counts(static_cast<std::size_t>(n_genes), 0.0) {}

  void begin_row() {
    touched.clear();
    if (current_mark == std::numeric_limits<int>::max()) {
      std::fill(marks.begin(), marks.end(), 0);
      current_mark = 1;
    } else {
      current_mark += 1;
    }
  }

  void add(int gene, double value) {
    if (marks[static_cast<std::size_t>(gene)] != current_mark) {
      marks[static_cast<std::size_t>(gene)] = current_mark;
      counts[static_cast<std::size_t>(gene)] = value;
      touched.push_back(gene);
      return;
    }
    counts[static_cast<std::size_t>(gene)] += value;
  }
};

// Default to all transcript rows when no explicit query subset is supplied.
std::vector<int> normalize_queries(const TranscriptTable& table, const NcvOptions& options) {
  std::vector<int> queries = options.query_indices;
  if (queries.empty()) {
    queries.resize(table.size());
    for (std::size_t i = 0; i < table.size(); ++i) {
      queries[i] = static_cast<int>(i);
    }
  }
  return queries;
}

// Clamp the requested worker count to the amount of row-wise work.
int effective_threads(int requested, int num_tasks) {
  if (num_tasks <= 0) {
    return 1;
  }
  return std::max(1, subpar::sanitize_num_workers(requested, num_tasks));
}

bool is_identity_transform(const NcvFeatureTransform& transform) {
  return transform.mode.empty() || transform.mode == "kl";
}

double transform_raw_value(double value, int gene, const NcvFeatureTransform& transform) {
  if (transform.mode == "sqrt_kl") {
    return std::sqrt(std::max(value, 0.0));
  }
  if (transform.mode == "invsqrt_kl") {
    if (gene < 0 || gene >= static_cast<int>(transform.gene_weights.size())) {
      return 0.0;
    }
    return value * transform.gene_weights[static_cast<std::size_t>(gene)];
  }
  return value;
}

double transformed_scratch_row_sum(
    const NcvScratch& scratch,
    const NcvFeatureTransform& transform) {
  double total = 0.0;
  for (const int gene : scratch.touched) {
    total += transform_raw_value(scratch.counts[static_cast<std::size_t>(gene)], gene, transform);
  }
  return total;
}

double transformed_sparse_row_sum(
    const SparseRowMatrix& x,
    int row,
    const NcvFeatureTransform& transform) {
  double total = 0.0;
  for (int p = x.indptr()[static_cast<std::size_t>(row)];
       p < x.indptr()[static_cast<std::size_t>(row + 1)];
       ++p) {
    total += transform_raw_value(
        x.values()[static_cast<std::size_t>(p)],
        x.indices()[static_cast<std::size_t>(p)],
        transform);
  }
  return total;
}

double transform_row_scale(double transformed_sum, const NcvFeatureTransform& transform) {
  if (is_identity_transform(transform) || transform.target_row_sum <= 0.0 || transformed_sum <= 0.0) {
    return 1.0;
  }
  return transform.target_row_sum / transformed_sum;
}

// Score one implicit NCV row by direct multiplication against H and row normalization.
void fill_normalized_score_row(
    const NcvScratch& scratch,
    const DenseMatrix& h,
    std::vector<double>& score_data,
    std::size_t row_offset) {
  const auto& h_data = h.data();
  const int n_factors = h.rows();
  const int n_genes = h.cols();
  double row_sum = 0.0;
  for (int factor = 0; factor < n_factors; ++factor) {
    const std::size_t h_offset = static_cast<std::size_t>(factor * n_genes);
    double total = 0.0;
    for (const int gene : scratch.touched) {
      total += scratch.counts[static_cast<std::size_t>(gene)] *
          h_data[h_offset + static_cast<std::size_t>(gene)];
    }
    score_data[row_offset + static_cast<std::size_t>(factor)] = total;
    row_sum += total;
  }
  if (row_sum <= 1e-12) {
    return;
  }
  for (int factor = 0; factor < n_factors; ++factor) {
    score_data[row_offset + static_cast<std::size_t>(factor)] /= row_sum;
  }
}

// Reproject one implicit NCV row with a fixed-H KL row solver.
void fill_kl_projected_score_row(
    const NcvScratch& scratch,
    const DenseMatrix& h,
    const std::vector<double>& h_sums,
    std::vector<double>& score_data,
    std::size_t row_offset,
    int max_iterations,
    double tolerance,
    double update_epsilon,
    const NcvFeatureTransform& transform) {
  const int n_factors = h.rows();
  std::vector<double> weights(static_cast<std::size_t>(n_factors), 0.0);
  std::vector<double> numerators(static_cast<std::size_t>(n_factors), 0.0);

  const double transformed_sum = transformed_scratch_row_sum(scratch, transform);
  const double row_scale = transform_row_scale(transformed_sum, transform);
  double row_sum = 0.0;
  double init_sum = 0.0;
  for (const int gene : scratch.touched) {
    row_sum += row_scale * transform_raw_value(
        scratch.counts[static_cast<std::size_t>(gene)],
        gene,
        transform);
  }
  for (int factor = 0; factor < n_factors; ++factor) {
    double total = 0.0;
    for (const int gene : scratch.touched) {
      const double value = row_scale * transform_raw_value(
          scratch.counts[static_cast<std::size_t>(gene)],
          gene,
          transform);
      total += value * h(factor, gene);
    }
    weights[static_cast<std::size_t>(factor)] = std::max(total, update_epsilon);
    init_sum += weights[static_cast<std::size_t>(factor)];
  }
  if (init_sum <= update_epsilon) {
    const double uniform = std::max(row_sum / std::max(n_factors, 1), update_epsilon);
    std::fill(weights.begin(), weights.end(), uniform);
  } else {
    const double scale = std::max(row_sum, 1.0) / init_sum;
    for (double& value : weights) {
      value = std::max(value * scale, update_epsilon);
    }
  }

  for (int iter = 0; iter < max_iterations; ++iter) {
    std::fill(numerators.begin(), numerators.end(), 0.0);
    for (const int gene : scratch.touched) {
      const double value = row_scale * transform_raw_value(
          scratch.counts[static_cast<std::size_t>(gene)],
          gene,
          transform);
      double fitted = update_epsilon;
      for (int factor = 0; factor < n_factors; ++factor) {
        fitted += weights[static_cast<std::size_t>(factor)] * h(factor, gene);
      }
      const double ratio = value / std::max(fitted, update_epsilon);
      for (int factor = 0; factor < n_factors; ++factor) {
        numerators[static_cast<std::size_t>(factor)] += h(factor, gene) * ratio;
      }
    }

    double max_relative_change = 0.0;
    for (int factor = 0; factor < n_factors; ++factor) {
      const double old_value = weights[static_cast<std::size_t>(factor)];
      const double updated = std::max(
          old_value * numerators[static_cast<std::size_t>(factor)] /
              std::max(h_sums[static_cast<std::size_t>(factor)], update_epsilon),
          update_epsilon);
      weights[static_cast<std::size_t>(factor)] = updated;
      max_relative_change = std::max(
          max_relative_change,
          std::abs(updated - old_value) / std::max(old_value, update_epsilon));
    }
    if (max_relative_change < tolerance) {
      break;
    }
  }

  double total = 0.0;
  for (const double value : weights) {
    total += value;
  }
  if (total <= update_epsilon) {
    return;
  }
  for (int factor = 0; factor < n_factors; ++factor) {
    score_data[row_offset + static_cast<std::size_t>(factor)] =
        weights[static_cast<std::size_t>(factor)] / total;
  }
}

// Reproject one already materialized sparse NCV row with a fixed-H KL row solver.
void fill_kl_projected_sparse_row(
    const SparseRowMatrix& x,
    int row,
    const DenseMatrix& h,
    const std::vector<double>& h_sums,
    std::vector<double>& score_data,
    std::size_t row_offset,
    int max_iterations,
    double tolerance,
    double update_epsilon,
    const NcvFeatureTransform& transform) {
  const int n_factors = h.rows();
  std::vector<double> weights(static_cast<std::size_t>(n_factors), 0.0);
  std::vector<double> numerators(static_cast<std::size_t>(n_factors), 0.0);

  const int row_start = x.indptr()[static_cast<std::size_t>(row)];
  const int row_end = x.indptr()[static_cast<std::size_t>(row + 1)];
  const double transformed_sum = transformed_sparse_row_sum(x, row, transform);
  const double row_scale = transform_row_scale(transformed_sum, transform);
  double row_sum = 0.0;
  for (int p = row_start; p < row_end; ++p) {
    row_sum += row_scale * transform_raw_value(
        x.values()[static_cast<std::size_t>(p)],
        x.indices()[static_cast<std::size_t>(p)],
        transform);
  }
  double init_sum = 0.0;
  for (int factor = 0; factor < n_factors; ++factor) {
    double total = 0.0;
    for (int p = row_start; p < row_end; ++p) {
      const int gene = x.indices()[static_cast<std::size_t>(p)];
      const double value = row_scale * transform_raw_value(
          x.values()[static_cast<std::size_t>(p)],
          gene,
          transform);
      total += value * h(factor, gene);
    }
    weights[static_cast<std::size_t>(factor)] = std::max(total, update_epsilon);
    init_sum += weights[static_cast<std::size_t>(factor)];
  }
  if (init_sum <= update_epsilon) {
    const double uniform = std::max(row_sum / std::max(n_factors, 1), update_epsilon);
    std::fill(weights.begin(), weights.end(), uniform);
  } else {
    const double scale = std::max(row_sum, 1.0) / init_sum;
    for (double& value : weights) {
      value = std::max(value * scale, update_epsilon);
    }
  }

  for (int iter = 0; iter < max_iterations; ++iter) {
    std::fill(numerators.begin(), numerators.end(), 0.0);
    for (int p = row_start; p < row_end; ++p) {
      const int gene = x.indices()[static_cast<std::size_t>(p)];
      const double value = row_scale * transform_raw_value(
          x.values()[static_cast<std::size_t>(p)],
          gene,
          transform);
      double fitted = update_epsilon;
      for (int factor = 0; factor < n_factors; ++factor) {
        fitted += weights[static_cast<std::size_t>(factor)] * h(factor, gene);
      }
      const double ratio = value / std::max(fitted, update_epsilon);
      for (int factor = 0; factor < n_factors; ++factor) {
        numerators[static_cast<std::size_t>(factor)] += h(factor, gene) * ratio;
      }
    }

    double max_relative_change = 0.0;
    for (int factor = 0; factor < n_factors; ++factor) {
      const double old_value = weights[static_cast<std::size_t>(factor)];
      const double updated = std::max(
          old_value * numerators[static_cast<std::size_t>(factor)] /
              std::max(h_sums[static_cast<std::size_t>(factor)], update_epsilon),
          update_epsilon);
      weights[static_cast<std::size_t>(factor)] = updated;
      max_relative_change = std::max(
          max_relative_change,
          std::abs(updated - old_value) / std::max(old_value, update_epsilon));
    }
    if (max_relative_change < tolerance) {
      break;
    }
  }

  double total = 0.0;
  for (const double value : weights) {
    total += value;
  }
  if (total <= update_epsilon) {
    return;
  }
  for (int factor = 0; factor < n_factors; ++factor) {
    score_data[row_offset + static_cast<std::size_t>(factor)] =
        weights[static_cast<std::size_t>(factor)] / total;
  }
}

// Visit NCV rows in query order and expose each row through a sparse scratch buffer.
template <class Callback>
void for_each_ncv_row(
    const TranscriptTable& table,
    const NcvOptions& options,
    Callback&& callback) {
  const auto queries = normalize_queries(table, options);
  NcvScratch scratch(static_cast<int>(table.num_genes()));
  const auto by_cell = table.transcripts_by_cell();

  if (!options.within_cell) {
    const SpatialKnnIndex global_index(table);
    for (std::size_t qi = 0; qi < queries.size(); ++qi) {
      scratch.begin_row();
      const auto hits =
          global_index.query(queries[static_cast<std::size_t>(qi)], options.k, options.include_self);
      for (const auto& hit : hits) {
        scratch.add(table.gene_index[static_cast<std::size_t>(hit.index)], 1.0);
      }
      callback(static_cast<int>(qi), scratch);
    }
    return;
  }

  std::vector<std::vector<int>> query_positions_by_cell(table.num_cells());
  for (std::size_t qi = 0; qi < queries.size(); ++qi) {
    const int query = queries[qi];
    const int cell = table.cell_index[static_cast<std::size_t>(query)];
    query_positions_by_cell[static_cast<std::size_t>(cell)].push_back(static_cast<int>(qi));
  }

  for (std::size_t cell = 0; cell < by_cell.size(); ++cell) {
    const auto& members = by_cell[cell];
    const auto& query_positions = query_positions_by_cell[cell];
    if (members.empty() || query_positions.empty()) {
      continue;
    }

    const SpatialKnnIndex cell_index(table, members);
    for (const int qi : query_positions) {
      scratch.begin_row();
      const int query = queries[static_cast<std::size_t>(qi)];
      const auto hits = cell_index.query(query, options.k, options.include_self);
      for (const auto& hit : hits) {
        scratch.add(table.gene_index[static_cast<std::size_t>(hit.index)], 1.0);
      }
      callback(qi, scratch);
    }
  }
}

}  // namespace

// Materialize a dense NCV matrix for the requested transcript queries.
DenseMatrix build_ncv_matrix(const TranscriptTable& table, const NcvOptions& options) {
  if (table.gene_index.empty() || table.cell_index.empty()) {
    throw std::runtime_error("TranscriptTable must be finalized before NCV construction");
  }

  const auto queries = normalize_queries(table, options);
  DenseMatrix out(static_cast<int>(queries.size()), static_cast<int>(table.num_genes()), 0.0);
  for_each_ncv_row(table, options, [&](int qi, const NcvScratch& scratch) {
    for (const int gene : scratch.touched) {
      out(qi, gene) = scratch.counts[static_cast<std::size_t>(gene)];
    }
  });

  return out;
}

// Materialize a sparse NCV matrix for the requested transcript queries.
SparseRowMatrix build_sparse_ncv_matrix(const TranscriptTable& table, const NcvOptions& options) {
  if (table.gene_index.empty() || table.cell_index.empty()) {
    throw std::runtime_error("TranscriptTable must be finalized before NCV construction");
  }

  const auto queries = normalize_queries(table, options);
  std::vector<std::vector<int>> row_indices(static_cast<std::size_t>(queries.size()));
  std::vector<std::vector<double>> row_values(static_cast<std::size_t>(queries.size()));
  std::size_t total_nnz = 0;
  for_each_ncv_row(table, options, [&](int qi, const NcvScratch& scratch) {
    auto touched = scratch.touched;
    std::sort(touched.begin(), touched.end());
    auto& indices = row_indices[static_cast<std::size_t>(qi)];
    auto& values = row_values[static_cast<std::size_t>(qi)];
    indices.reserve(touched.size());
    values.reserve(touched.size());
    for (const int gene : touched) {
      indices.push_back(gene);
      values.push_back(scratch.counts[static_cast<std::size_t>(gene)]);
    }
    total_nnz += touched.size();
  });

  std::vector<int> indptr;
  std::vector<int> indices;
  std::vector<double> values;
  indptr.reserve(static_cast<std::size_t>(queries.size() + 1));
  indices.reserve(total_nnz);
  values.reserve(total_nnz);
  indptr.push_back(0);
  for (std::size_t row = 0; row < row_indices.size(); ++row) {
    indices.insert(
        indices.end(),
        row_indices[row].begin(),
        row_indices[row].end());
    values.insert(
        values.end(),
        row_values[row].begin(),
        row_values[row].end());
    indptr.push_back(static_cast<int>(indices.size()));
  }

  return SparseRowMatrix(
      static_cast<int>(queries.size()),
      static_cast<int>(table.num_genes()),
      std::move(indptr),
      std::move(indices),
      std::move(values));
}

NcvFeatureTransform make_ncv_feature_transform(
    const SparseRowMatrix& x,
    const std::string& mode) {
  if (mode != "kl" && mode != "invsqrt_kl" && mode != "sqrt_kl") {
    throw std::runtime_error("nmf_variant must be one of kl, invsqrt_kl, or sqrt_kl");
  }
  NcvFeatureTransform transform;
  transform.mode = mode;
  if (mode == "kl") {
    return transform;
  }

  std::vector<double> row_sums;
  row_sums.reserve(static_cast<std::size_t>(x.rows()));
  for (int row = 0; row < x.rows(); ++row) {
    double total = 0.0;
    for (int p = x.indptr()[static_cast<std::size_t>(row)];
         p < x.indptr()[static_cast<std::size_t>(row + 1)];
         ++p) {
      total += x.values()[static_cast<std::size_t>(p)];
    }
    if (total > 0.0) {
      row_sums.push_back(total);
    }
  }
  if (!row_sums.empty()) {
    const auto mid = row_sums.begin() + static_cast<std::ptrdiff_t>(row_sums.size() / 2);
    std::nth_element(row_sums.begin(), mid, row_sums.end());
    transform.target_row_sum = *mid;
  }

  if (mode == "invsqrt_kl") {
    transform.gene_weights.assign(static_cast<std::size_t>(x.cols()), 0.0);
    const auto col_sums = x.col_sums();
    double weight_sum = 0.0;
    int weight_count = 0;
    for (int col = 0; col < x.cols(); ++col) {
      if (col_sums[static_cast<std::size_t>(col)] <= 0.0) {
        continue;
      }
      const double weight = 1.0 / std::sqrt(col_sums[static_cast<std::size_t>(col)] + 1.0);
      transform.gene_weights[static_cast<std::size_t>(col)] = weight;
      weight_sum += weight;
      ++weight_count;
    }
    const double mean_weight = weight_count > 0 ? weight_sum / static_cast<double>(weight_count) : 1.0;
    if (mean_weight > 0.0) {
      for (double& weight : transform.gene_weights) {
        if (weight > 0.0) {
          weight /= mean_weight;
        }
      }
    }
  }
  return transform;
}

SparseRowMatrix transform_sparse_ncv_matrix(
    const SparseRowMatrix& x,
    const NcvFeatureTransform& transform) {
  if (is_identity_transform(transform)) {
    return x;
  }
  std::vector<int> indptr;
  std::vector<int> indices;
  std::vector<double> values;
  indptr.reserve(static_cast<std::size_t>(x.rows() + 1));
  indices.reserve(x.indices().size());
  values.reserve(x.values().size());
  indptr.push_back(0);
  for (int row = 0; row < x.rows(); ++row) {
    const double transformed_sum = transformed_sparse_row_sum(x, row, transform);
    const double row_scale = transform_row_scale(transformed_sum, transform);
    for (int p = x.indptr()[static_cast<std::size_t>(row)];
         p < x.indptr()[static_cast<std::size_t>(row + 1)];
         ++p) {
      const int gene = x.indices()[static_cast<std::size_t>(p)];
      const double value = row_scale * transform_raw_value(
          x.values()[static_cast<std::size_t>(p)],
          gene,
          transform);
      if (value <= 0.0) {
        continue;
      }
      indices.push_back(gene);
      values.push_back(value);
    }
    indptr.push_back(static_cast<int>(indices.size()));
  }
  return SparseRowMatrix(x.rows(), x.cols(), std::move(indptr), std::move(indices), std::move(values));
}

NcvFeatureTransform expand_ncv_feature_transform(
    const NcvFeatureTransform& compact,
    int full_cols,
    const std::vector<int>& kept_cols) {
  if (compact.gene_weights.empty() || static_cast<int>(kept_cols.size()) == full_cols) {
    return compact;
  }
  NcvFeatureTransform out = compact;
  out.gene_weights.assign(static_cast<std::size_t>(full_cols), 0.0);
  for (std::size_t local_col = 0; local_col < kept_cols.size(); ++local_col) {
    const int full_col = kept_cols[local_col];
    if (full_col >= 0 && full_col < full_cols && local_col < compact.gene_weights.size()) {
      out.gene_weights[static_cast<std::size_t>(full_col)] = compact.gene_weights[local_col];
    }
  }
  return out;
}

// Project transcript NCVs with the legacy normalized dot-product rule.
DenseMatrix project_ncv_to_factors(
    const TranscriptTable& table,
    const DenseMatrix& h,
    const NcvOptions& options,
    int num_threads) {
  if (table.gene_index.empty() || table.cell_index.empty()) {
    throw std::runtime_error("TranscriptTable must be finalized before NCV projection");
  }
  if (h.cols() != static_cast<int>(table.num_genes())) {
    throw std::runtime_error("h columns must match TranscriptTable genes");
  }

  const auto queries = normalize_queries(table, options);
  DenseMatrix scores(static_cast<int>(queries.size()), h.rows(), 0.0);
  if (num_threads <= 1 || !options.within_cell) {
    for_each_ncv_row(table, options, [&](int qi, const NcvScratch& scratch) {
      for (const int gene : scratch.touched) {
        const double count = scratch.counts[static_cast<std::size_t>(gene)];
        for (int factor = 0; factor < h.rows(); ++factor) {
          scores(qi, factor) += count * h(factor, gene);
        }
      }
    });
    scores.normalize_rows();
    return scores;
  }

  const auto by_cell = table.transcripts_by_cell();
  std::vector<std::vector<int>> query_positions_by_cell(table.num_cells());
  for (std::size_t qi = 0; qi < queries.size(); ++qi) {
    const int query = queries[qi];
    const int cell = table.cell_index[static_cast<std::size_t>(query)];
    query_positions_by_cell[static_cast<std::size_t>(cell)].push_back(static_cast<int>(qi));
  }

  auto& score_data = scores.data();
  const int workers = effective_threads(num_threads, static_cast<int>(by_cell.size()));
  subpar::parallelize_range<true>(workers, static_cast<int>(by_cell.size()), [&](int, int start, int length) {
    NcvScratch scratch(static_cast<int>(table.num_genes()));
    for (int cell = start; cell < start + length; ++cell) {
      const auto& members = by_cell[static_cast<std::size_t>(cell)];
      const auto& query_positions = query_positions_by_cell[static_cast<std::size_t>(cell)];
      if (members.empty() || query_positions.empty()) {
        continue;
      }
      const SpatialKnnIndex cell_index(table, members);
      for (const int qi : query_positions) {
        scratch.begin_row();
        const int query = queries[static_cast<std::size_t>(qi)];
        const auto hits = cell_index.query(query, options.k, options.include_self);
        for (const auto& hit : hits) {
          scratch.add(table.gene_index[static_cast<std::size_t>(hit.index)], 1.0);
        }
        const std::size_t row_offset = static_cast<std::size_t>(qi * h.rows());
        fill_normalized_score_row(scratch, h, score_data, row_offset);
      }
    }
  });
  return scores;
}

// Assign transcript factor potentials from the observed gene's normalized H row.
DenseMatrix project_gene_loadings_to_factors(
    const TranscriptTable& table,
    const DenseMatrix& h,
    const NcvOptions& options,
    int num_threads) {
  if (table.gene_index.empty()) {
    throw std::runtime_error("TranscriptTable must be finalized before gene-loading projection");
  }
  if (h.cols() != static_cast<int>(table.num_genes())) {
    throw std::runtime_error("h columns must match TranscriptTable genes");
  }

  const auto queries = normalize_queries(table, options);
  const int n_factors = h.rows();
  DenseMatrix scores(static_cast<int>(queries.size()), n_factors, 0.0);
  auto& score_data = scores.data();
  const auto& h_data = h.data();
  const int n_genes = h.cols();
  const double uniform = n_factors > 0 ? 1.0 / static_cast<double>(n_factors) : 0.0;

  const int workers = effective_threads(num_threads, static_cast<int>(queries.size()));
  subpar::parallelize_range<true>(workers, static_cast<int>(queries.size()), [&](int, int start, int length) {
    for (int row = start; row < start + length; ++row) {
      const int query = queries[static_cast<std::size_t>(row)];
      const int gene = table.gene_index[static_cast<std::size_t>(query)];
      const std::size_t row_offset = static_cast<std::size_t>(row * n_factors);
      double total = 0.0;
      if (gene >= 0 && gene < n_genes) {
        for (int factor = 0; factor < n_factors; ++factor) {
          const double value = std::max(
              0.0,
              h_data[static_cast<std::size_t>(factor * n_genes + gene)]);
          score_data[row_offset + static_cast<std::size_t>(factor)] = value;
          total += value;
        }
      }
      if (total <= 1e-12) {
        for (int factor = 0; factor < n_factors; ++factor) {
          score_data[row_offset + static_cast<std::size_t>(factor)] = uniform;
        }
        continue;
      }
      for (int factor = 0; factor < n_factors; ++factor) {
        score_data[row_offset + static_cast<std::size_t>(factor)] /= total;
      }
    }
  });
  return scores;
}

// Project transcript NCVs with the fixed-H KL row solver.
DenseMatrix project_ncv_to_factors_kl(
    const TranscriptTable& table,
    const DenseMatrix& h,
    const NcvOptions& options,
    int num_threads,
    int max_iterations,
    double tolerance,
    double update_epsilon,
    NcvFeatureTransform transform) {
  if (table.gene_index.empty() || table.cell_index.empty()) {
    throw std::runtime_error("TranscriptTable must be finalized before NCV projection");
  }
  if (h.cols() != static_cast<int>(table.num_genes())) {
    throw std::runtime_error("h columns must match TranscriptTable genes");
  }

  const auto queries = normalize_queries(table, options);
  DenseMatrix scores(static_cast<int>(queries.size()), h.rows(), 0.0);
  const auto h_sums = h.row_sums();

  if (num_threads <= 1 || !options.within_cell) {
    auto& score_data = scores.data();
    for_each_ncv_row(table, options, [&](int qi, const NcvScratch& scratch) {
      const std::size_t row_offset = static_cast<std::size_t>(qi * h.rows());
      fill_kl_projected_score_row(
          scratch,
          h,
          h_sums,
          score_data,
          row_offset,
          max_iterations,
          tolerance,
          update_epsilon,
          transform);
    });
    return scores;
  }

  const auto by_cell = table.transcripts_by_cell();
  std::vector<std::vector<int>> query_positions_by_cell(table.num_cells());
  for (std::size_t qi = 0; qi < queries.size(); ++qi) {
    const int query = queries[qi];
    const int cell = table.cell_index[static_cast<std::size_t>(query)];
    query_positions_by_cell[static_cast<std::size_t>(cell)].push_back(static_cast<int>(qi));
  }

  auto& score_data = scores.data();
  const int workers = effective_threads(num_threads, static_cast<int>(by_cell.size()));
  subpar::parallelize_range<true>(workers, static_cast<int>(by_cell.size()), [&](int, int start, int length) {
    NcvScratch scratch(static_cast<int>(table.num_genes()));
    for (int cell = start; cell < start + length; ++cell) {
      const auto& members = by_cell[static_cast<std::size_t>(cell)];
      const auto& query_positions = query_positions_by_cell[static_cast<std::size_t>(cell)];
      if (members.empty() || query_positions.empty()) {
        continue;
      }
      const SpatialKnnIndex cell_index(table, members);
      for (const int qi : query_positions) {
        scratch.begin_row();
        const int query = queries[static_cast<std::size_t>(qi)];
        const auto hits = cell_index.query(query, options.k, options.include_self);
        for (const auto& hit : hits) {
          scratch.add(table.gene_index[static_cast<std::size_t>(hit.index)], 1.0);
        }
        const std::size_t row_offset = static_cast<std::size_t>(qi * h.rows());
        fill_kl_projected_score_row(
            scratch,
            h,
            h_sums,
            score_data,
            row_offset,
            max_iterations,
            tolerance,
            update_epsilon,
            transform);
      }
    }
  });
  return scores;
}

// Project prebuilt sparse NCV rows with the same fixed-H KL row solver.
DenseMatrix project_rows_to_factors_kl(
    const SparseRowMatrix& x,
    const DenseMatrix& h,
    int num_threads,
    int max_iterations,
    double tolerance,
    double update_epsilon,
    NcvFeatureTransform transform) {
  if (x.cols() != h.cols()) {
    throw std::runtime_error("x columns must match h columns for KL row projection");
  }
  DenseMatrix scores(x.rows(), h.rows(), 0.0);
  const auto h_sums = h.row_sums();
  auto& score_data = scores.data();
  const int workers = effective_threads(num_threads, x.rows());
  if (workers <= 1) {
    for (int row = 0; row < x.rows(); ++row) {
      const std::size_t row_offset = static_cast<std::size_t>(row * h.rows());
      fill_kl_projected_sparse_row(
          x,
          row,
          h,
          h_sums,
          score_data,
          row_offset,
          max_iterations,
          tolerance,
          update_epsilon,
          transform);
    }
    return scores;
  }

  subpar::parallelize_range<true>(workers, x.rows(), [&](int, int start, int length) {
    for (int row = start; row < start + length; ++row) {
      const std::size_t row_offset = static_cast<std::size_t>(row * h.rows());
      fill_kl_projected_sparse_row(
          x,
          row,
          h,
          h_sums,
          score_data,
          row_offset,
          max_iterations,
          tolerance,
          update_epsilon,
          transform);
    }
  });
  return scores;
}

}  // namespace celladmix
