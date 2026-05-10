#include "celladmix/report_embedding.hpp"

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <cctype>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "celladmix/clustering.hpp"
#include "celladmix/ncv.hpp"
#include "irlba/irlba.hpp"

namespace celladmix {
namespace {

// Report-side embedding helpers for cell and sampled-NCV visualization artifacts.

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

// Create an output directory if it does not already exist.
void ensure_directory(const std::filesystem::path& path) {
  std::error_code ec;
  std::filesystem::create_directories(path, ec);
  if (ec) {
    throw std::runtime_error(
        "Could not create directory '" + path.string() + "': " + ec.message());
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

// Locate the report sidecar directory for a run.
std::filesystem::path report_dir(const std::string& run_path_or_dir) {
  return resolve_run_dir(run_path_or_dir) / "report";
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

// Write one Arrow table to parquet for later notebook/report use.
void write_parquet_table(
    const std::shared_ptr<arrow::Table>& table,
    const std::filesystem::path& path) {
  auto sink = arrow_unwrap(
      arrow::io::FileOutputStream::Open(path.string()),
      "Open parquet output");
  auto writer = arrow_unwrap(
      parquet::arrow::FileWriter::Open(*table->schema(), arrow::default_memory_pool(), sink),
      "Open parquet writer");
  arrow_check(writer->WriteTable(*table, std::max<int64_t>(1, std::min<int64_t>(table->num_rows(), 65536))),
      "Write parquet table");
  arrow_check(writer->Close(), "Close parquet writer");
  arrow_check(sink->Close(), "Close parquet sink");
}

// Read one Arrow table from parquet.
std::shared_ptr<arrow::Table> read_parquet_table(const std::filesystem::path& path) {
  auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path.string()), "Open parquet input");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open parquet reader");
  auto reader = arrow_unwrap(builder.Build(), "Build parquet reader");
  std::shared_ptr<arrow::Table> table;
  arrow_check(reader->ReadTable(&table), "Read parquet table");
  return table;
}

// Compute the top-two margin per row in a factor score matrix.
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
    if (!std::isfinite(second)) {
      second = 0.0;
    }
    out[static_cast<std::size_t>(row)] = best - second;
  }
  return out;
}

// Row-normalize component weights before visualization or comparison.
DenseMatrix normalized_component_weights(const DenseMatrix& weights) {
  DenseMatrix out = weights;
  out.normalize_rows();
  return out;
}

// Extract a row subset and row-normalize it for report overlays.
DenseMatrix subset_and_normalize_rows(const DenseMatrix& matrix, const std::vector<int>& row_indices) {
  if (matrix.rows() == static_cast<int>(row_indices.size())) {
    DenseMatrix out = matrix;
    out.normalize_rows();
    return out;
  }
  DenseMatrix out(static_cast<int>(row_indices.size()), matrix.cols(), 0.0);
  for (std::size_t i = 0; i < row_indices.size(); ++i) {
    const int source_row = row_indices[i];
    for (int col = 0; col < matrix.cols(); ++col) {
      out(static_cast<int>(i), col) = matrix(source_row, col);
    }
  }
  out.normalize_rows();
  return out;
}

// Drop all-zero NCV columns before PCA/UMAP feature construction.
SparseRowMatrix compact_sparse_columns(const SparseRowMatrix& x) {
  std::vector<double> col_sums = x.col_sums();
  std::vector<int> keep_cols;
  keep_cols.reserve(static_cast<std::size_t>(x.cols()));
  for (int j = 0; j < x.cols(); ++j) {
    if (col_sums[static_cast<std::size_t>(j)] > 0.0) {
      keep_cols.push_back(j);
    }
  }

  if (static_cast<int>(keep_cols.size()) == x.cols()) {
    return x;
  }

  std::vector<int> new_col(x.cols(), -1);
  for (std::size_t j = 0; j < keep_cols.size(); ++j) {
    new_col[static_cast<std::size_t>(keep_cols[j])] = static_cast<int>(j);
  }

  std::vector<int> indptr;
  std::vector<int> indices;
  std::vector<double> values;
  indptr.reserve(static_cast<std::size_t>(x.rows() + 1));
  indptr.push_back(0);
  for (int i = 0; i < x.rows(); ++i) {
    for (int p = x.indptr()[static_cast<std::size_t>(i)];
         p < x.indptr()[static_cast<std::size_t>(i + 1)];
         ++p) {
      const int mapped = new_col[static_cast<std::size_t>(x.indices()[static_cast<std::size_t>(p)])];
      if (mapped >= 0) {
        indices.push_back(mapped);
        values.push_back(x.values()[static_cast<std::size_t>(p)]);
      }
    }
    indptr.push_back(static_cast<int>(indices.size()));
  }

  return SparseRowMatrix(
      x.rows(),
      static_cast<int>(keep_cols.size()),
      std::move(indptr),
      std::move(indices),
      std::move(values));
}

// Convert sparse NCV counts into log-scaled features for embedding.
Eigen::MatrixXd sparse_ncv_features(
    const SparseRowMatrix& x,
    double normalization_scale) {
  Eigen::MatrixXd out(x.cols(), x.rows());
  out.setZero();
  for (int i = 0; i < x.rows(); ++i) {
    double denom = 0.0;
    for (int p = x.indptr()[static_cast<std::size_t>(i)];
         p < x.indptr()[static_cast<std::size_t>(i + 1)];
         ++p) {
      denom += x.values()[static_cast<std::size_t>(p)];
    }
    denom = std::max(denom, 1.0);
    for (int p = x.indptr()[static_cast<std::size_t>(i)];
         p < x.indptr()[static_cast<std::size_t>(i + 1)];
         ++p) {
      const int j = x.indices()[static_cast<std::size_t>(p)];
      out(j, i) = std::log1p(x.values()[static_cast<std::size_t>(p)] * normalization_scale / denom);
    }
  }
  return out;
}

// Run a lightweight IRLBA PCA step before NCV UMAP construction.
Eigen::MatrixXd run_pca_for_embedding(Eigen::MatrixXd features, int pca_dims, unsigned int seed) {
  if (features.cols() == 0 || features.rows() == 0) {
    return features;
  }
  const Eigen::VectorXd mean = features.rowwise().mean();
  features.colwise() -= mean;
  const int k = std::max<int>(
      1,
      std::min<int>(
          pca_dims,
          std::min<int>(features.rows(), features.cols())));
  if (k >= std::min<int>(features.rows(), features.cols())) {
    return features;
  }
  irlba::Options options;
  options.seed = static_cast<std::uint64_t>(seed);
  options.extra_work = std::max(10, std::min(30, k));
  Eigen::MatrixXd left;
  Eigen::MatrixXd right;
  Eigen::VectorXd singular;
  const auto status = irlba::compute(
      features,
      static_cast<Eigen::Index>(k),
      left,
      right,
      singular,
      options);
  if (!status.first) {
    throw std::runtime_error("IRLBA PCA for NCV UMAP did not converge within the iteration limit");
  }
  return singular.asDiagonal() * right.transpose();
}

// Convert DenseMatrix rows into Eigen columns for UMAP input.
Eigen::MatrixXd dense_rows_to_eigen_cols(const DenseMatrix& matrix) {
  Eigen::MatrixXd out(matrix.cols(), matrix.rows());
  for (int row = 0; row < matrix.rows(); ++row) {
    for (int col = 0; col < matrix.cols(); ++col) {
      out(col, row) = matrix(row, col);
    }
  }
  return out;
}

// Assemble the sampled-training NCV UMAP sidecar from one in-memory fit result.
TrainingMoleculeUmapData build_training_molecule_umap_data(
    const TranscriptTable& table,
    const BasicPipelineResult& fit,
    int ncv_k,
    int umap_neighbors,
    int umap_epochs,
    unsigned int seed,
    double normalization_scale,
    int pca_dims,
    int num_threads,
    const std::vector<std::int64_t>* obs_id_override) {
  if (fit.training_query_indices.size() != static_cast<std::size_t>(fit.nmf.w.rows())) {
    throw std::runtime_error("training_query_indices must align with NMF W rows");
  }
  if (obs_id_override != nullptr && obs_id_override->size() != fit.training_query_indices.size()) {
    throw std::runtime_error("obs_id_override must align with training_query_indices");
  }

  TrainingMoleculeUmapData out;
  const int n = static_cast<int>(fit.training_query_indices.size());
  const int rank = fit.nmf.w.cols();
  if (n == 0 || rank == 0) {
    return out;
  }

  const DenseMatrix components = normalized_component_weights(fit.nmf.w);
  const DenseMatrix projected = subset_and_normalize_rows(fit.factor_scores, fit.training_query_indices);
  const std::vector<double> component_margin = compute_factor_margin(components);
  const std::vector<double> projected_margin = compute_factor_margin(projected);
  NcvOptions ncv_options;
  ncv_options.k = ncv_k + 1;
  ncv_options.include_self = true;
  ncv_options.within_cell = true;
  ncv_options.query_indices = fit.training_query_indices;
  const SparseRowMatrix training_ncv_raw = build_sparse_ncv_matrix(table, ncv_options);
  const SparseRowMatrix training_ncv = compact_sparse_columns(training_ncv_raw);
  Eigen::MatrixXd umap_input = sparse_ncv_features(training_ncv, normalization_scale);
  if (umap_input.rows() > std::max(1, pca_dims)) {
    umap_input = run_pca_for_embedding(std::move(umap_input), pca_dims, seed);
  }
  Eigen::MatrixXd umap = Eigen::MatrixXd::Zero(2, n);
  if (n > 1 && rank > 0) {
    umap = umap_embed(
        umap_input,
        2,
        std::min(umap_neighbors, n - 1),
        umap_epochs,
        static_cast<int>(seed),
        num_threads,
        true);
  }

  out.training_rank.reserve(static_cast<std::size_t>(n));
  out.obs_id.reserve(static_cast<std::size_t>(n));
  out.cell_id.reserve(static_cast<std::size_t>(n));
  out.gene.reserve(static_cast<std::size_t>(n));
  out.x.reserve(static_cast<std::size_t>(n));
  out.y.reserve(static_cast<std::size_t>(n));
  out.z.reserve(static_cast<std::size_t>(n));
  out.factor_label.reserve(static_cast<std::size_t>(n));
  out.factor_margin.reserve(static_cast<std::size_t>(n));
  out.dominant_component.reserve(static_cast<std::size_t>(n));
  out.dominant_component_margin.reserve(static_cast<std::size_t>(n));
  out.umap_x.reserve(static_cast<std::size_t>(n));
  out.umap_y.reserve(static_cast<std::size_t>(n));
  out.component_weights = components;
  out.projected_component_weights = projected;

  const std::vector<double> final_margin = compute_factor_margin(fit.factor_scores);

  for (int training_rank = 0; training_rank < n; ++training_rank) {
    const int obs_id = fit.training_query_indices[static_cast<std::size_t>(training_rank)];
    const int margin_row = final_margin.size() == static_cast<std::size_t>(n) ? training_rank : obs_id;
    out.training_rank.push_back(training_rank);
    out.obs_id.push_back(
        obs_id_override == nullptr
        ? static_cast<std::int64_t>(obs_id)
        : (*obs_id_override)[static_cast<std::size_t>(training_rank)]);
    out.cell_id.push_back(table.cells[static_cast<std::size_t>(table.cell_index[static_cast<std::size_t>(obs_id)])]);
    out.gene.push_back(table.genes[static_cast<std::size_t>(table.gene_index[static_cast<std::size_t>(obs_id)])]);
    out.x.push_back(table.x[static_cast<std::size_t>(obs_id)]);
    out.y.push_back(table.y[static_cast<std::size_t>(obs_id)]);
    out.z.push_back(table.z[static_cast<std::size_t>(obs_id)]);
    out.factor_label.push_back(fit.labels[static_cast<std::size_t>(obs_id)] + 1);
    out.factor_margin.push_back(final_margin[static_cast<std::size_t>(margin_row)]);
    out.dominant_component.push_back(components.row_argmax(training_rank) + 1);
    out.dominant_component_margin.push_back(component_margin[static_cast<std::size_t>(training_rank)]);
    out.projected_component.push_back(projected.row_argmax(training_rank) + 1);
    out.projected_component_margin.push_back(projected_margin[static_cast<std::size_t>(training_rank)]);
    out.umap_x.push_back(umap(0, training_rank));
    out.umap_y.push_back(umap(1, training_rank));
  }

  return out;
}

// Rehydrate the report sidecar from a parquet table.
TrainingMoleculeUmapData table_to_training_molecule_umap(
    const std::shared_ptr<arrow::Table>& table) {
  TrainingMoleculeUmapData out;
  if (table->num_rows() == 0) {
    return out;
  }

  auto read_int32 = [&](const std::string& name) {
    std::vector<int> values(static_cast<std::size_t>(table->num_rows()));
    const int idx = table->schema()->GetFieldIndex(name);
    const auto chunked = table->column(idx);
    int64_t offset = 0;
    for (int c = 0; c < chunked->num_chunks(); ++c) {
      const auto chunk = std::static_pointer_cast<arrow::Int32Array>(chunked->chunk(c));
      for (int64_t i = 0; i < chunk->length(); ++i) {
        values[static_cast<std::size_t>(offset + i)] = chunk->Value(i);
      }
      offset += chunk->length();
    }
    return values;
  };
  auto read_int64 = [&](const std::string& name) {
    std::vector<std::int64_t> values(static_cast<std::size_t>(table->num_rows()));
    const int idx = table->schema()->GetFieldIndex(name);
    const auto chunked = table->column(idx);
    int64_t offset = 0;
    for (int c = 0; c < chunked->num_chunks(); ++c) {
      const auto chunk = std::static_pointer_cast<arrow::Int64Array>(chunked->chunk(c));
      for (int64_t i = 0; i < chunk->length(); ++i) {
        values[static_cast<std::size_t>(offset + i)] = chunk->Value(i);
      }
      offset += chunk->length();
    }
    return values;
  };
  auto read_double = [&](const std::string& name) {
    std::vector<double> values(static_cast<std::size_t>(table->num_rows()));
    const int idx = table->schema()->GetFieldIndex(name);
    const auto chunked = table->column(idx);
    int64_t offset = 0;
    for (int c = 0; c < chunked->num_chunks(); ++c) {
      const auto chunk = std::static_pointer_cast<arrow::DoubleArray>(chunked->chunk(c));
      for (int64_t i = 0; i < chunk->length(); ++i) {
        values[static_cast<std::size_t>(offset + i)] = chunk->Value(i);
      }
      offset += chunk->length();
    }
    return values;
  };
  auto read_string = [&](const std::string& name) {
    std::vector<std::string> values(static_cast<std::size_t>(table->num_rows()));
    const int idx = table->schema()->GetFieldIndex(name);
    const auto chunked = table->column(idx);
    int64_t offset = 0;
    for (int c = 0; c < chunked->num_chunks(); ++c) {
      const auto chunk = std::static_pointer_cast<arrow::StringArray>(chunked->chunk(c));
      for (int64_t i = 0; i < chunk->length(); ++i) {
        values[static_cast<std::size_t>(offset + i)] = chunk->GetString(i);
      }
      offset += chunk->length();
    }
    return values;
  };

  out.training_rank = read_int32("training_rank");
  out.obs_id = read_int64("obs_id");
  out.cell_id = read_string("cell_id");
  out.gene = read_string("gene");
  out.x = read_double("x");
  out.y = read_double("y");
  out.z = read_double("z");
  out.factor_label = read_int32("factor_label");
  out.factor_margin = read_double("factor_margin");
  out.dominant_component = read_int32("dominant_component");
  out.dominant_component_margin = read_double("dominant_component_margin");
  out.projected_component = read_int32("projected_component");
  out.projected_component_margin = read_double("projected_component_margin");
  out.umap_x = read_double("umap_x");
  out.umap_y = read_double("umap_y");

  auto has_numeric_suffix = [](const std::string& name, const std::string& prefix) {
    if (name.rfind(prefix, 0) != 0 || name.size() <= prefix.size()) {
      return false;
    }
    return std::all_of(
        name.begin() + static_cast<std::ptrdiff_t>(prefix.size()),
        name.end(),
        [](unsigned char ch) { return std::isdigit(ch) != 0; });
  };

  std::vector<std::string> component_names;
  std::vector<std::string> projected_component_names;
  for (int i = 0; i < table->num_columns(); ++i) {
    const auto name = table->schema()->field(i)->name();
    if (has_numeric_suffix(name, "component_")) {
      component_names.push_back(name);
    } else if (has_numeric_suffix(name, "projected_component_")) {
      projected_component_names.push_back(name);
    }
  }
  std::sort(
      component_names.begin(),
      component_names.end(),
      [](const std::string& lhs, const std::string& rhs) {
        return std::stoi(lhs.substr(10)) < std::stoi(rhs.substr(10));
      });
  std::sort(
      projected_component_names.begin(),
      projected_component_names.end(),
      [](const std::string& lhs, const std::string& rhs) {
        return std::stoi(lhs.substr(20)) < std::stoi(rhs.substr(20));
      });
  out.component_weights = DenseMatrix(
      static_cast<int>(table->num_rows()),
      static_cast<int>(component_names.size()),
      0.0);
  for (std::size_t j = 0; j < component_names.size(); ++j) {
    const auto values = read_double(component_names[j]);
    for (std::size_t i = 0; i < values.size(); ++i) {
      out.component_weights(static_cast<int>(i), static_cast<int>(j)) = values[i];
    }
  }
  out.projected_component_weights = DenseMatrix(
      static_cast<int>(table->num_rows()),
      static_cast<int>(projected_component_names.size()),
      0.0);
  for (std::size_t j = 0; j < projected_component_names.size(); ++j) {
    const auto values = read_double(projected_component_names[j]);
    for (std::size_t i = 0; i < values.size(); ++i) {
      out.projected_component_weights(static_cast<int>(i), static_cast<int>(j)) = values[i];
    }
  }
  return out;
}

}  // namespace

// Return the on-disk path of the sampled-training NCV UMAP sidecar.
std::string training_molecule_umap_parquet_path(const std::string& run_path_or_dir) {
  return (report_dir(run_path_or_dir) / "training_molecule_umap.parquet").string();
}

// Build and persist the sampled-training NCV UMAP report sidecar.
void write_training_molecule_umap_report(
    const std::string& run_root_dir,
    const TranscriptTable& table,
    const BasicPipelineResult& fit,
    int ncv_k,
    int umap_neighbors,
    int umap_epochs,
    unsigned int seed,
    double normalization_scale,
    int pca_dims,
    int num_threads,
    const std::vector<std::int64_t>* obs_id_override) {
  const auto data = build_training_molecule_umap_data(
      table,
      fit,
      ncv_k,
      umap_neighbors,
      umap_epochs,
      seed,
      normalization_scale,
      pca_dims,
      num_threads,
      obs_id_override);
  if (data.obs_id.empty()) {
    return;
  }

  const auto output_path = std::filesystem::path(training_molecule_umap_parquet_path(run_root_dir));
  ensure_directory(output_path.parent_path());

  std::vector<std::shared_ptr<arrow::Field>> fields = {
      arrow::field("training_rank", arrow::int32()),
      arrow::field("obs_id", arrow::int64()),
      arrow::field("cell_id", arrow::utf8()),
      arrow::field("gene", arrow::utf8()),
      arrow::field("x", arrow::float64()),
      arrow::field("y", arrow::float64()),
      arrow::field("z", arrow::float64()),
      arrow::field("factor_label", arrow::int32()),
      arrow::field("factor_margin", arrow::float64()),
      arrow::field("dominant_component", arrow::int32()),
      arrow::field("dominant_component_margin", arrow::float64()),
      arrow::field("projected_component", arrow::int32()),
      arrow::field("projected_component_margin", arrow::float64()),
      arrow::field("umap_x", arrow::float64()),
      arrow::field("umap_y", arrow::float64()),
  };
  std::vector<std::shared_ptr<arrow::Array>> arrays = {
      build_int32_array(data.training_rank),
      build_int64_array(data.obs_id),
      build_string_array(data.cell_id),
      build_string_array(data.gene),
      build_double_array(data.x),
      build_double_array(data.y),
      build_double_array(data.z),
      build_int32_array(data.factor_label),
      build_double_array(data.factor_margin),
      build_int32_array(data.dominant_component),
      build_double_array(data.dominant_component_margin),
      build_int32_array(data.projected_component),
      build_double_array(data.projected_component_margin),
      build_double_array(data.umap_x),
      build_double_array(data.umap_y),
  };

  for (int component = 0; component < data.component_weights.cols(); ++component) {
    std::vector<double> values(static_cast<std::size_t>(data.component_weights.rows()));
    for (int row = 0; row < data.component_weights.rows(); ++row) {
      values[static_cast<std::size_t>(row)] = data.component_weights(row, component);
    }
    fields.push_back(arrow::field("component_" + std::to_string(component + 1), arrow::float64()));
    arrays.push_back(build_double_array(values));
  }
  for (int component = 0; component < data.projected_component_weights.cols(); ++component) {
    std::vector<double> values(static_cast<std::size_t>(data.projected_component_weights.rows()));
    for (int row = 0; row < data.projected_component_weights.rows(); ++row) {
      values[static_cast<std::size_t>(row)] = data.projected_component_weights(row, component);
    }
    fields.push_back(arrow::field("projected_component_" + std::to_string(component + 1), arrow::float64()));
    arrays.push_back(build_double_array(values));
  }

  write_parquet_table(
      arrow::Table::Make(std::make_shared<arrow::Schema>(fields), arrays),
      output_path);
}

// Load the sampled-training NCV UMAP report sidecar back into memory.
TrainingMoleculeUmapData load_training_molecule_umap_report(const std::string& run_path_or_dir) {
  const auto path = std::filesystem::path(training_molecule_umap_parquet_path(run_path_or_dir));
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error("Training molecule UMAP report not found: " + path.string());
  }
  return table_to_training_molecule_umap(read_parquet_table(path));
}

}  // namespace celladmix
