#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "celladmix/bridge.hpp"
#include "celladmix/input_store.hpp"
#include "celladmix/membrane.hpp"
#include "celladmix/pipeline_store.hpp"
#include "celladmix/run_store.hpp"
#include "celladmix/workflow.hpp"
#include "celladmix/xenium.hpp"

namespace py = pybind11;

namespace {

// The Python API receives plain dict/list objects rather than exposing C++
// structs. These helpers keep serialization rules centralized and stable.
py::dict paths_to_dict(const celladmix::RunPaths& paths) {
  py::dict out;
  out["root_dir"] = paths.root_dir;
  out["run_json"] = paths.run_json;
  out["factors_parquet"] = paths.factors_parquet;
  out["cells_parquet"] = paths.cells_parquet;
  out["molecules_parquet"] = paths.molecules_parquet;
  out["training_rows_parquet"] = paths.training_rows_parquet;
  out["scores_dir"] = paths.scores_dir;
  out["corrected_dir"] = paths.corrected_dir;
  return out;
}

py::dict source_to_dict(const celladmix::RunSourceInfo& source) {
  py::dict out;
  out["type"] = source.type;
  out["path"] = source.path;
  out["used_parquet"] = source.used_parquet;
  py::list files;
  for (const auto& file : source.files) {
    py::dict f;
    f["role"] = file.role;
    f["path"] = file.path;
    f["exists"] = file.exists;
    files.append(f);
  }
  out["files"] = files;
  return out;
}

py::dict pipeline_options_to_dict(const celladmix::BasicPipelineOptions& options) {
  py::dict out;
  out["ncv_k"] = options.ncv_k;
  out["rank"] = options.rank;
  out["graph_k"] = options.graph_k;
  out["same_label_ratio"] = options.same_label_ratio;
  out["nmf_iterations"] = options.nmf_iterations;
  out["nmf_n_runs"] = options.nmf_n_runs;
  out["nmf_init"] = options.nmf_init;
  out["nmf_variant"] = options.nmf_variant;
  out["molecule_scoring"] = options.molecule_scoring;
  out["nmf_train_max_rows"] = options.nmf_train_max_rows;
  out["nmf_min_molecules"] = options.nmf_min_molecules;
  out["num_threads"] = options.num_threads;
  out["seed"] = options.seed;
  out["training_scope_cell_types"] = options.training_scope_cell_types;
  return out;
}

py::dict nmf_diagnostics_to_dict(const celladmix::NmfRunDiagnostics& diagnostics) {
  py::dict out;
  out["final_objective"] = diagnostics.final_objective;
  out["selected_seed"] = diagnostics.selected_seed;
  out["selected_run"] = diagnostics.selected_run + 1;
  out["candidate_final_objectives"] = diagnostics.candidate_final_objectives;
  out["candidate_best_match_correlations"] = diagnostics.candidate_best_match_correlations;
  out["selected_factor_stability"] = diagnostics.selected_factor_stability;
  out["candidate_final_objective_mean"] = diagnostics.candidate_final_objective_mean;
  out["candidate_final_objective_sd"] = diagnostics.candidate_final_objective_sd;
  out["candidate_best_match_correlation_mean"] =
      diagnostics.candidate_best_match_correlation_mean;
  out["stability_metric"] = diagnostics.stability_metric;
  out["stability_comparison_runs"] = diagnostics.stability_comparison_runs;
  out["stable_factor_count"] = diagnostics.stable_factor_count;
  out["stability_threshold"] = diagnostics.stability_threshold;
  return out;
}

py::dict manifest_to_dict(const celladmix::RunManifest& manifest) {
  py::dict out;
  out["format_version"] = manifest.format_version;
  out["run_type"] = manifest.run_type;
  out["source"] = source_to_dict(manifest.source);
  out["paths"] = paths_to_dict(manifest.paths);
  out["pipeline_options"] = pipeline_options_to_dict(manifest.pipeline_options);
  out["genes"] = manifest.genes;
  out["nmf_gene_weights"] = manifest.nmf_gene_weights;
  out["nmf_diagnostics"] = nmf_diagnostics_to_dict(manifest.nmf_diagnostics);
  out["crop_ids"] = manifest.crop_ids;
  out["n_transcripts"] = manifest.n_transcripts;
  out["n_cells"] = manifest.n_cells;
  out["n_factors"] = manifest.n_factors;
  out["n_training_rows"] = manifest.n_training_rows;
  out["has_z"] = manifest.has_z;
  out["has_qv"] = manifest.has_qv;
  out["has_transcript_id"] = manifest.has_transcript_id;
  out["has_cell_type"] = manifest.has_cell_type;
  out["has_sample_id"] = manifest.has_sample_id;
  out["has_fov_id"] = manifest.has_fov_id;
  out["has_nucleus_id"] = manifest.has_nucleus_id;
  out["has_overlaps_nucleus"] = manifest.has_overlaps_nucleus;
  out["has_nucleus_distance"] = manifest.has_nucleus_distance;
  if (manifest.analysis_crop.has_value()) {
    out["analysis_crop"] = *manifest.analysis_crop;
  } else {
    out["analysis_crop"] = py::none();
  }
  if (manifest.parent_run.has_value()) {
    out["parent_run"] = *manifest.parent_run;
  } else {
    out["parent_run"] = py::none();
  }
  return out;
}

py::dict input_manifest_to_dict(const celladmix::InputStoreManifest& manifest) {
  py::dict out;
  out["format_version"] = manifest.format_version;
  out["root_dir"] = manifest.root_dir;
  out["source_type"] = manifest.source_type;
  out["store_mode"] = manifest.store_mode;
  out["source_path"] = manifest.source_path;
  out["source_fingerprint"] = manifest.source_fingerprint;
  out["filter_signature"] = manifest.filter_signature;
  out["has_molecules"] = manifest.has_molecules;
  out["has_molecule_rows"] = manifest.has_molecule_rows;
  out["has_cell_gene_counts"] = manifest.has_cell_gene_counts;
  out["has_cell_offsets"] = manifest.has_cell_offsets;
  out["n_molecules"] = manifest.n_molecules;
  out["n_cells"] = manifest.n_cells;
  out["n_genes"] = manifest.n_genes;
  out["has_qv"] = manifest.has_qv;
  out["has_transcript_id"] = manifest.has_transcript_id;
  out["has_cell_type"] = manifest.has_cell_type;
  out["has_sample_id"] = manifest.has_sample_id;
  out["has_fov_id"] = manifest.has_fov_id;
  out["has_overlaps_nucleus"] = manifest.has_overlaps_nucleus;
  out["has_nucleus_distance"] = manifest.has_nucleus_distance;
  return out;
}

py::dict dense_matrix_to_dict(const celladmix::DenseMatrix& matrix) {
  py::dict out;
  out["rows"] = matrix.rows();
  out["cols"] = matrix.cols();
  out["data"] = matrix.data();
  return out;
}

celladmix::RunSourceInfo source_from_store_manifest(
    const std::string& store_dir,
    const celladmix::InputStoreManifest& manifest) {
  celladmix::RunSourceInfo source;
  source.type = manifest.source_type.empty() ? "input_store" : manifest.source_type;
  source.path = manifest.source_path.empty() ? store_dir : manifest.source_path;
  source.used_parquet = true;

  const auto paths = celladmix::make_input_store_paths(store_dir);
  auto add_file = [&](const std::string& role, const std::string& path) {
    celladmix::RunSourceFile file;
    file.role = role;
    file.path = path;
    file.exists = std::filesystem::exists(path);
    source.files.push_back(file);
  };
  add_file("input_store_manifest", paths.manifest_json);
  add_file("input_store_genes", paths.genes_parquet);
  add_file("input_store_cells", paths.cells_parquet);
  add_file("input_store_counts", paths.counts_parquet);
  add_file("input_store_molecules", paths.molecules_parquet);
  add_file("input_store_cell_offsets", paths.cell_offsets_parquet);
  return source;
}

std::function<void(const std::string&, double)> make_progress_logger(bool verbose) {
  const auto start = std::make_shared<std::chrono::steady_clock::time_point>(
      std::chrono::steady_clock::now());
  return [verbose, start](const std::string& message, double stage_sec) {
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
         << std::chrono::duration<double>(std::chrono::steady_clock::now() - *start).count()
         << "s] " << message;
    if (stage_sec >= 0.0) {
      line << " (" << std::fixed << std::setprecision(3) << stage_sec << "s)";
    }
    std::cout << line.str() << std::endl;
  };
}

std::unordered_map<std::string, std::string> make_annotation_map(
    const std::vector<std::string>& cell_ids,
    const std::vector<std::string>& labels) {
  if (cell_ids.size() != labels.size()) {
    throw std::runtime_error("annotation cell_ids and labels must have the same length");
  }
  std::unordered_map<std::string, std::string> out;
  out.reserve(cell_ids.size());
  for (std::size_t i = 0; i < cell_ids.size(); ++i) {
    if (!cell_ids[i].empty() && !labels[i].empty()) {
      out.emplace(cell_ids[i], labels[i]);
    }
  }
  return out;
}

void apply_annotation_override(
    celladmix::RunData& run_data,
    const std::vector<std::string>& cell_ids,
    const std::vector<std::string>& labels) {
  const auto mapping = make_annotation_map(cell_ids, labels);
  if (mapping.empty()) return;

  // Python callers pass the active annotation separately from the stored run.
  // Apply it to both cell-level and transcript-level tables before scoring.
  run_data.cells.cell_types.assign(run_data.cells.cell_ids.size(), "");
  for (std::size_t i = 0; i < run_data.cells.cell_ids.size(); ++i) {
    const auto found = mapping.find(run_data.cells.cell_ids[i]);
    if (found != mapping.end()) {
      run_data.cells.cell_types[i] = found->second;
    }
  }

  run_data.transcripts.cell_types.assign(run_data.transcripts.size(), "");
  for (std::size_t i = 0; i < run_data.transcripts.size(); ++i) {
    std::string cell_id = run_data.transcripts.orig_cell_id[i];
    if (!run_data.transcripts.cell_index.empty() &&
        run_data.transcripts.cell_index[i] >= 0 &&
        static_cast<std::size_t>(run_data.transcripts.cell_index[i]) <
            run_data.transcripts.cells.size()) {
      cell_id = run_data.transcripts.cells[static_cast<std::size_t>(run_data.transcripts.cell_index[i])];
    }
    const auto found = mapping.find(cell_id);
    if (found != mapping.end()) {
      run_data.transcripts.cell_types[i] = found->second;
    }
  }
}

std::vector<std::string> labels_for_cells(
    const celladmix::CellTable& cells,
    const std::vector<std::string>& cell_ids,
    const std::vector<std::string>& labels) {
  const auto mapping = make_annotation_map(cell_ids, labels);
  std::vector<std::string> out;
  out.reserve(cells.cell_ids.size());
  for (const auto& cell_id : cells.cell_ids) {
    const auto found = mapping.find(cell_id);
    out.push_back(found == mapping.end() ? "" : found->second);
  }
  return out;
}

py::list membrane_pair_scores_to_list(
    const celladmix::TranscriptTable& table,
    const celladmix::MembraneTestResult& result) {
  py::list out;
  for (const auto& row : result.pair_scores) {
    py::dict d;
    d["target_cell"] = row.target_cell >= 0 &&
            static_cast<std::size_t>(row.target_cell) < table.cells.size()
        ? table.cells[static_cast<std::size_t>(row.target_cell)]
        : "";
    d["source_cell"] = row.source_cell >= 0 &&
            static_cast<std::size_t>(row.source_cell) < table.cells.size()
        ? table.cells[static_cast<std::size_t>(row.source_cell)]
        : "";
    d["target_cell_type"] = row.target_type >= 0 &&
            static_cast<std::size_t>(row.target_type) < result.cell_types.size()
        ? result.cell_types[static_cast<std::size_t>(row.target_type)]
        : "";
    d["source_cell_type"] = row.source_type >= 0 &&
            static_cast<std::size_t>(row.source_type) < result.cell_types.size()
        ? result.cell_types[static_cast<std::size_t>(row.source_type)]
        : "";
    d["factor"] = row.factor + 1;
    d["factor_count"] = row.factor_count;
    d["scored_molecules"] = row.scored_molecules;
    d["mean_score"] = row.mean_score;
    d["fraction_positive"] = row.fraction_positive;
    d["mean_directional_weight"] = row.mean_directional_weight;
    d["used_in_summary"] = row.used_in_summary;
    out.append(d);
  }
  return out;
}

py::list membrane_summary_to_list(const celladmix::MembraneTestResult& result) {
  py::list out;
  for (const auto& row : result.summaries) {
    py::dict d;
    d["target_cell_type"] = row.target_type >= 0 &&
            static_cast<std::size_t>(row.target_type) < result.cell_types.size()
        ? result.cell_types[static_cast<std::size_t>(row.target_type)]
        : "";
    d["source_cell_type"] = row.source_type >= 0 &&
            static_cast<std::size_t>(row.source_type) < result.cell_types.size()
        ? result.cell_types[static_cast<std::size_t>(row.source_type)]
        : "";
    d["factor"] = row.factor + 1;
    d["n_pairs"] = row.n_pairs;
    d["mean_score"] = row.mean_score;
    d["q75_score"] = row.q75_score;
    d["fraction_positive_pairs"] = row.fraction_positive_pairs;
    d["p_value"] = row.p_value;
    d["neg_log10_p"] = row.neg_log10_p;
    out.append(d);
  }
  return out;
}

py::list bridge_pair_scores_to_list(
    const celladmix::TranscriptTable& table,
    const celladmix::BridgeTestResult& result) {
  py::list out;
  for (const auto& row : result.pair_scores) {
    py::dict d;
    d["target_cell"] = row.target_cell >= 0 &&
            static_cast<std::size_t>(row.target_cell) < table.cells.size()
        ? table.cells[static_cast<std::size_t>(row.target_cell)]
        : "";
    d["source_cell"] = row.source_cell >= 0 &&
            static_cast<std::size_t>(row.source_cell) < table.cells.size()
        ? table.cells[static_cast<std::size_t>(row.source_cell)]
        : "";
    d["target_cell_type"] = row.target_type >= 0 &&
            static_cast<std::size_t>(row.target_type) < result.cell_types.size()
        ? result.cell_types[static_cast<std::size_t>(row.target_type)]
        : "";
    d["source_cell_type"] = row.source_type >= 0 &&
            static_cast<std::size_t>(row.source_type) < result.cell_types.size()
        ? result.cell_types[static_cast<std::size_t>(row.source_type)]
        : "";
    d["factor"] = row.factor + 1;
    d["factor_count"] = row.factor_count;
    d["total_crossing_count"] = row.total_crossing_count;
    d["target_fraction"] = row.target_fraction;
    d["source_fraction"] = row.source_fraction;
    d["crossing_fraction"] = row.crossing_fraction;
    d["score"] = row.score;
    d["used_in_summary"] = row.used_in_summary;
    out.append(d);
  }
  return out;
}

py::list bridge_summary_to_list(const celladmix::BridgeTestResult& result) {
  py::list out;
  for (const auto& row : result.summaries) {
    py::dict d;
    d["target_cell_type"] = row.target_type >= 0 &&
            static_cast<std::size_t>(row.target_type) < result.cell_types.size()
        ? result.cell_types[static_cast<std::size_t>(row.target_type)]
        : "";
    d["source_cell_type"] = row.source_type >= 0 &&
            static_cast<std::size_t>(row.source_type) < result.cell_types.size()
        ? result.cell_types[static_cast<std::size_t>(row.source_type)]
        : "";
    d["factor"] = row.factor + 1;
    d["n_pairs"] = row.n_pairs;
    d["mean_score"] = row.mean_score;
    d["mean_null_score"] = row.mean_null_score;
    d["q75_score"] = row.q75_score;
    d["p_value"] = row.p_value;
    d["neg_log10_p"] = row.neg_log10_p;
    out.append(d);
  }
  return out;
}

py::dict counts_to_dict(const celladmix::CellCountMatrix& counts) {
  py::dict out;
  out["indptr"] = counts.indptr;
  out["indices"] = counts.indices;
  out["data"] = counts.values;
  out["genes"] = counts.genes;
  out["cells"] = counts.cells.cell_ids;
  out["cell_type"] = counts.cells.cell_types;
  out["transcript_counts"] = counts.transcript_counts;
  out["detected_genes"] = counts.detected_genes;
  return out;
}

py::dict clustering_to_dict(const celladmix::CellClusteringResult& result) {
  py::dict out;
  out["cell_id"] = result.cells.cell_ids;
  out["cell_type"] = result.cells.cell_types.empty()
      ? std::vector<std::string>(result.cells.size(), "")
      : result.cells.cell_types;
  out["cluster"] = result.clusters;
  out["transcript_count"] = result.transcript_counts;
  out["detected_genes"] = result.detected_genes;
  out["x"] = result.cells.centroid_x;
  out["y"] = result.cells.centroid_y;
  out["z"] = result.cells.centroid_z;
  out["analysis_crop"] = result.crop_ids;
  out["variable_genes"] = result.variable_genes;
  out["pca_variance_explained"] = result.pca_variance_explained;
  out["has_umap"] = result.has_umap;

  py::list pcs;
  for (int row = 0; row < result.pcs.rows(); ++row) {
    py::list values;
    for (int col = 0; col < result.pcs.cols(); ++col) {
      values.append(result.pcs(row, col));
    }
    pcs.append(values);
  }
  out["pcs"] = pcs;

  std::vector<double> umap_1(result.cells.size(), std::numeric_limits<double>::quiet_NaN());
  std::vector<double> umap_2(result.cells.size(), std::numeric_limits<double>::quiet_NaN());
  if (result.has_umap) {
    for (std::size_t i = 0; i < result.cells.size(); ++i) {
      if (result.umap.cols() > 0) {
        umap_1[i] = result.umap(static_cast<int>(i), 0);
      }
      if (result.umap.cols() > 1) {
        umap_2[i] = result.umap(static_cast<int>(i), 1);
      }
    }
  }
  out["umap_1"] = umap_1;
  out["umap_2"] = umap_2;
  return out;
}

celladmix::CellClusteringOptions make_clustering_options(
    int min_molecules,
    int min_genes,
    int cells_max,
    int n_variable_genes,
    int pca_dims,
    int graph_k,
    double cluster_resolution,
    bool compute_umap,
    int umap_neighbors,
    int umap_epochs,
    int num_threads,
    bool umap_parallel_optimization,
    double normalization_scale,
    unsigned int seed) {
  celladmix::CellClusteringOptions options;
  options.min_molecules = min_molecules;
  options.min_genes = min_genes;
  options.cells_max = cells_max;
  options.n_variable_genes = n_variable_genes;
  options.pca_dims = pca_dims;
  options.graph_k = graph_k;
  options.cluster_resolution = cluster_resolution;
  options.compute_umap = compute_umap;
  options.umap_neighbors = umap_neighbors;
  options.umap_epochs = umap_epochs;
  options.num_threads = num_threads;
  options.umap_parallel_optimization = umap_parallel_optimization;
  options.normalization_scale = normalization_scale;
  options.seed = seed;
  return options;
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "Python bindings for the C++ cellAdmix core";

  // Long-running C++ kernels release the GIL so progress logging and other
  // Python threads are not blocked while native work is running.
  m.def(
      "build_xenium_store",
      [](const std::string& bundle_dir,
         const std::string& store_dir,
         const std::vector<std::string>& cell_filter,
         const std::vector<std::string>& gene_filter,
         double min_qv,
         bool keep_unassigned,
         bool keep_non_gene,
         bool read_cells,
         bool prefer_parquet,
         bool materialize_molecules,
         bool force,
         int num_threads,
         int parquet_row_group_size) {
        celladmix::XeniumLoadOptions load_options;
        load_options.cell_filter = cell_filter;
        load_options.gene_filter = gene_filter;
        load_options.min_qv = min_qv;
        load_options.keep_unassigned = keep_unassigned;
        load_options.keep_non_gene = keep_non_gene;
        load_options.read_cells = read_cells;
        load_options.prefer_parquet = prefer_parquet;

        celladmix::InputStoreBuildOptions store_options;
        store_options.store_dir = store_dir;
        store_options.materialize_molecules = materialize_molecules;
        store_options.force = force;
        store_options.num_threads = num_threads;
        store_options.parquet_row_group_size = parquet_row_group_size;

        celladmix::InputStoreManifest manifest;
        {
          py::gil_scoped_release release;
          manifest =
              celladmix::build_xenium_input_store(bundle_dir, load_options, store_options);
        }
        return input_manifest_to_dict(manifest);
      },
      py::arg("bundle_dir"),
      py::arg("store_dir"),
      py::arg("cell_filter") = std::vector<std::string>{},
      py::arg("gene_filter") = std::vector<std::string>{},
      py::arg("min_qv") = -1.0,
      py::arg("keep_unassigned") = false,
      py::arg("keep_non_gene") = false,
      py::arg("read_cells") = true,
      py::arg("prefer_parquet") = true,
      py::arg("materialize_molecules") = true,
      py::arg("force") = false,
      py::arg("num_threads") = 1,
      py::arg("parquet_row_group_size") = 65536);

  m.def(
      "read_input_store_manifest",
      [](const std::string& store_dir) {
        return input_manifest_to_dict(celladmix::read_input_store_manifest(store_dir));
      },
      py::arg("store_dir"));

  m.def(
      "fit_store",
      [](const std::string& store_dir,
         const std::string& out_dir,
         int rank,
         int ncv_k,
         int graph_k,
         double same_label_ratio,
         int nmf_iterations,
         const std::string& nmf_init,
         const std::string& nmf_variant,
         const std::string& molecule_scoring,
         int nmf_n_runs,
         int nmf_train_max_rows,
         int nmf_min_molecules,
         int num_threads,
         unsigned int seed,
         const std::vector<std::string>& training_cell_strata,
         const std::vector<std::string>& training_scope_cell_types,
         bool report_ncv_umap,
         bool verbose) {
        celladmix::BasicPipelineOptions pipeline_options;
        pipeline_options.ncv_k = ncv_k;
        pipeline_options.rank = rank;
        pipeline_options.graph_k = graph_k;
        pipeline_options.same_label_ratio = same_label_ratio;
        pipeline_options.nmf_iterations = nmf_iterations;
        pipeline_options.nmf_init = nmf_init;
        pipeline_options.nmf_variant = nmf_variant;
        pipeline_options.molecule_scoring = molecule_scoring;
        pipeline_options.nmf_n_runs = nmf_n_runs;
        pipeline_options.nmf_train_max_rows = nmf_train_max_rows;
        pipeline_options.nmf_min_molecules = nmf_min_molecules;
        pipeline_options.num_threads = num_threads;
        pipeline_options.seed = seed;
        pipeline_options.return_ncv = false;
        pipeline_options.training_cell_strata = training_cell_strata;
        pipeline_options.training_scope_cell_types = training_scope_cell_types;

        celladmix::RunStorageOptions storage_options;
        const auto store_manifest = celladmix::read_input_store_manifest(store_dir);
        const auto source = source_from_store_manifest(store_dir, store_manifest);
        const auto counts = celladmix::load_input_store_counts(store_dir);

        celladmix::StorePipelineResult fit;
        {
          // NMF fitting can run for seconds to minutes; keep Python responsive.
          py::gil_scoped_release release;
          fit = celladmix::run_basic_pipeline_store(
              store_dir,
              source,
              pipeline_options,
              storage_options,
              out_dir,
              std::nullopt,
              report_ncv_umap,
              &counts,
              verbose);
        }
        auto out = manifest_to_dict(fit.manifest);
        out["h"] = dense_matrix_to_dict(fit.nmf.h);
        out["nmf_final_objective"] = fit.nmf.final_objective;
        return out;
      },
      py::arg("store_dir"),
      py::arg("out_dir"),
      py::arg("rank"),
      py::arg("ncv_k") = 0,
      py::arg("graph_k") = 10,
      py::arg("same_label_ratio") = 5.0,
      py::arg("nmf_iterations") = 150,
      py::arg("nmf_init") = "auto",
      py::arg("nmf_variant") = "invsqrt_kl",
      py::arg("molecule_scoring") = "gene_loadings",
      py::arg("nmf_n_runs") = 1,
      py::arg("nmf_train_max_rows") = 10000,
      py::arg("nmf_min_molecules") = 10,
      py::arg("num_threads") = 1,
      py::arg("seed") = 1U,
      py::arg("training_cell_strata") = std::vector<std::string>{},
      py::arg("training_scope_cell_types") = std::vector<std::string>{},
      py::arg("report_ncv_umap") = false,
      py::arg("verbose") = false);

  m.def(
      "read_run_manifest",
      [](const std::string& path) {
        return manifest_to_dict(celladmix::read_run_manifest(path));
      },
      py::arg("path"));

  m.def(
      "load_run_cells_path",
      [](const std::string& path) {
        return celladmix::read_run_manifest(path).paths.cells_parquet;
      },
      py::arg("path"));

  m.def(
      "score_membrane",
      [](const std::string& run_path,
         const std::string& image_path,
         double pixel_size,
         double x_offset,
         double y_offset,
         double epsilon,
         int cell_candidate_k,
         int candidate_pairs_per_type_pair,
         double cell_candidate_halo,
         int min_factor_molecules,
         int min_pairs,
         int max_cells_per_type_pair,
         double control_distance_fraction,
         int line_samples,
         int num_threads,
         unsigned int seed,
         const std::vector<std::string>& annotation_cell_ids,
         const std::vector<std::string>& annotation_labels,
         bool verbose) {
        celladmix::MembraneImageOptions image_options;
        image_options.image_path = image_path;
        image_options.pixel_size = pixel_size;
        image_options.x_offset = x_offset;
        image_options.y_offset = y_offset;
        image_options.epsilon = epsilon;

        celladmix::MembraneTestOptions options;
        options.cell_candidate_k = cell_candidate_k;
        options.candidate_pairs_per_type_pair = candidate_pairs_per_type_pair;
        options.cell_candidate_halo = cell_candidate_halo;
        options.min_factor_molecules = min_factor_molecules;
        options.min_pairs = min_pairs;
        options.max_cells_per_type_pair = max_cells_per_type_pair;
        options.control_distance_fraction = control_distance_fraction;
        options.line_samples = line_samples;
        options.num_threads = num_threads;
        options.seed = seed;
        options.progress = make_progress_logger(verbose);

        auto run_data = celladmix::load_run_data(run_path);
        apply_annotation_override(run_data, annotation_cell_ids, annotation_labels);

        celladmix::MembraneTestResult result;
        {
          // Scoring reads image tiles and loops over many candidate pairs.
          py::gil_scoped_release release;
          result = celladmix::run_membrane_test(
              run_data.transcripts, run_data.cells, run_data.labels, image_options, options);
        }
        py::dict out;
        out["cell_types"] = result.cell_types;
        out["scores"] = membrane_pair_scores_to_list(run_data.transcripts, result);
        out["summary"] = membrane_summary_to_list(result);
        out["image_path"] = image_path;
        out["pixel_size"] = pixel_size;
        return out;
      },
      py::arg("run_path"),
      py::arg("image_path"),
      py::arg("pixel_size"),
      py::arg("x_offset") = 0.0,
      py::arg("y_offset") = 0.0,
      py::arg("epsilon") = 1.0,
      py::arg("cell_candidate_k") = 50,
      py::arg("candidate_pairs_per_type_pair") = 400,
      py::arg("cell_candidate_halo") = -1.0,
      py::arg("min_factor_molecules") = 5,
      py::arg("min_pairs") = 5,
      py::arg("max_cells_per_type_pair") = 400,
      py::arg("control_distance_fraction") = 0.025,
      py::arg("line_samples") = 16,
      py::arg("num_threads") = 1,
      py::arg("seed") = 1U,
      py::arg("annotation_cell_ids") = std::vector<std::string>{},
      py::arg("annotation_labels") = std::vector<std::string>{},
      py::arg("verbose") = false);

  m.def(
      "score_bridge",
      [](const std::string& run_path,
         const std::string& candidate_mode,
         int candidate_k,
         int cell_candidate_k,
         int candidate_pairs_per_type_pair,
         int crossing_k,
         int min_type_pair_contacts,
         int min_factor_molecules,
         int min_pairs,
         int max_cells_per_type_pair,
         int null_iterations,
         int null_max_iterations,
         double cell_candidate_halo,
         double null_step,
         double null_pool_fraction,
         int num_threads,
         unsigned int seed,
         bool compute_null,
         bool fast_null_crossing,
         const std::vector<std::string>& annotation_cell_ids,
         const std::vector<std::string>& annotation_labels,
         bool verbose) {
        celladmix::BridgeTestOptions options;
        options.candidate_mode = candidate_mode;
        options.candidate_k = candidate_k;
        options.cell_candidate_k = cell_candidate_k;
        options.candidate_pairs_per_type_pair = candidate_pairs_per_type_pair;
        options.crossing_k = crossing_k;
        options.min_type_pair_contacts = min_type_pair_contacts;
        options.min_factor_molecules = min_factor_molecules;
        options.min_pairs = min_pairs;
        options.max_cells_per_type_pair = max_cells_per_type_pair;
        options.null_iterations = null_iterations;
        options.null_max_iterations = null_max_iterations;
        options.cell_candidate_halo = cell_candidate_halo;
        options.null_step = null_step;
        options.null_pool_fraction = null_pool_fraction;
        options.num_threads = num_threads;
        options.seed = seed;
        options.compute_null = compute_null;
        options.fast_null_crossing = fast_null_crossing;
        options.progress = make_progress_logger(verbose);

        auto run_data = celladmix::load_run_data(run_path);
        apply_annotation_override(run_data, annotation_cell_ids, annotation_labels);

        celladmix::BridgeTestResult result;
        {
          // Bridge scoring is CPU-heavy and does not need Python objects.
          py::gil_scoped_release release;
          result =
              celladmix::run_bridge_test(run_data.transcripts, run_data.cells, run_data.labels, options);
        }
        py::dict out;
        out["cell_types"] = result.cell_types;
        out["scores"] = bridge_pair_scores_to_list(run_data.transcripts, result);
        out["summary"] = bridge_summary_to_list(result);
        return out;
      },
      py::arg("run_path"),
      py::arg("candidate_mode") = "cell_center",
      py::arg("candidate_k") = 10,
      py::arg("cell_candidate_k") = 50,
      py::arg("candidate_pairs_per_type_pair") = 1000,
      py::arg("crossing_k") = 20,
      py::arg("min_type_pair_contacts") = 5,
      py::arg("min_factor_molecules") = 5,
      py::arg("min_pairs") = 5,
      py::arg("max_cells_per_type_pair") = 400,
      py::arg("null_iterations") = 3,
      py::arg("null_max_iterations") = 50,
      py::arg("cell_candidate_halo") = -1.0,
      py::arg("null_step") = 0.02,
      py::arg("null_pool_fraction") = 0.1,
      py::arg("num_threads") = 1,
      py::arg("seed") = 1U,
      py::arg("compute_null") = true,
      py::arg("fast_null_crossing") = true,
      py::arg("annotation_cell_ids") = std::vector<std::string>{},
      py::arg("annotation_labels") = std::vector<std::string>{},
      py::arg("verbose") = false);

  m.def(
      "correct_run",
      [](const std::string& run_path,
         const std::string& out_dir,
         const std::vector<int>& factors,
         const std::vector<std::string>& target_cell_types,
         const std::vector<std::string>& annotation_cell_ids,
         const std::vector<std::string>& annotation_labels) {
        if (factors.size() != target_cell_types.size()) {
          throw std::runtime_error("factors and target_cell_types must have the same length");
        }
        auto run_data = celladmix::load_run_data(run_path);
        apply_annotation_override(run_data, annotation_cell_ids, annotation_labels);

        // Correction rules are Python-facing one-based factors; the core label
        // vector is zero-based.
        std::vector<bool> keep_mask(run_data.transcripts.size(), true);
        for (std::size_t rule = 0; rule < factors.size(); ++rule) {
          const int zero_based_factor = factors[rule] - 1;
          std::vector<bool> rule_keep(run_data.transcripts.size(), true);
          if (run_data.transcripts.cell_types.size() == run_data.transcripts.size()) {
            for (std::size_t i = 0; i < run_data.transcripts.size(); ++i) {
              if (run_data.labels[i] == zero_based_factor &&
                  run_data.transcripts.cell_types[i] == target_cell_types[rule]) {
                rule_keep[i] = false;
              }
            }
          } else {
            rule_keep = celladmix::apply_removal_rule(
                run_data.transcripts,
                run_data.labels,
                zero_based_factor,
                target_cell_types[rule]);
          }
          for (std::size_t i = 0; i < keep_mask.size(); ++i) {
            keep_mask[i] = keep_mask[i] && rule_keep[i];
          }
        }

        celladmix::RunManifest manifest;
        {
          // Writing corrected parquet files can be slow on large runs.
          py::gil_scoped_release release;
          manifest = celladmix::write_corrected_run(out_dir, run_path, run_data, keep_mask);
        }
        int kept = 0;
        for (const bool keep : keep_mask) kept += keep ? 1 : 0;
        auto out = manifest_to_dict(manifest);
        out["n_removed"] = static_cast<int>(keep_mask.size()) - kept;
        out["correction_summary_parquet"] =
            celladmix::correction_summary_parquet_path(manifest.paths.root_dir);
        return out;
      },
      py::arg("run_path"),
      py::arg("out_dir"),
      py::arg("factors"),
      py::arg("target_cell_types"),
      py::arg("annotation_cell_ids") = std::vector<std::string>{},
      py::arg("annotation_labels") = std::vector<std::string>{});

  m.def(
      "collect_counts",
      [](const std::string& run_path) {
        celladmix::CellCountMatrix counts;
        {
          // Count collection can scan corrected run output.
          py::gil_scoped_release release;
          counts = celladmix::collect_run_counts(run_path);
        }
        return counts_to_dict(counts);
      },
      py::arg("run_path"));

  m.def(
      "cluster_store_counts",
      [](const std::string& store_dir,
         int min_molecules,
         int min_genes,
         int cells_max,
         int n_variable_genes,
         int pca_dims,
         int graph_k,
         double cluster_resolution,
         bool compute_umap,
         int umap_neighbors,
         int umap_epochs,
         int num_threads,
         bool umap_parallel_optimization,
         double normalization_scale,
         unsigned int seed) {
        auto options = make_clustering_options(
            min_molecules,
            min_genes,
            cells_max,
            n_variable_genes,
            pca_dims,
            graph_k,
            cluster_resolution,
            compute_umap,
            umap_neighbors,
            umap_epochs,
            num_threads,
            umap_parallel_optimization,
            normalization_scale,
            seed);
        celladmix::CellClusteringResult result;
        {
          // Load counts while the GIL is released; conversion back to dict
          // happens after the native result is complete.
          py::gil_scoped_release release;
          const auto counts = celladmix::load_input_store_counts(store_dir);
          result = celladmix::cluster_cell_counts(counts, options);
        }
        return clustering_to_dict(result);
      },
      py::arg("store_dir"),
      py::arg("min_molecules") = 10,
      py::arg("min_genes") = 5,
      py::arg("cells_max") = -1,
      py::arg("n_variable_genes") = 1000,
      py::arg("pca_dims") = 30,
      py::arg("graph_k") = 15,
      py::arg("cluster_resolution") = 1.0,
      py::arg("compute_umap") = true,
      py::arg("umap_neighbors") = 15,
      py::arg("umap_epochs") = 200,
      py::arg("num_threads") = 1,
      py::arg("umap_parallel_optimization") = true,
      py::arg("normalization_scale") = 5000.0,
      py::arg("seed") = 1U);

  m.def(
      "cluster_run_counts",
      [](const std::string& run_path,
         int min_molecules,
         int min_genes,
         int cells_max,
         int n_variable_genes,
         int pca_dims,
         int graph_k,
         double cluster_resolution,
         bool compute_umap,
         int umap_neighbors,
         int umap_epochs,
         int num_threads,
         bool umap_parallel_optimization,
         double normalization_scale,
         unsigned int seed) {
        auto options = make_clustering_options(
            min_molecules,
            min_genes,
            cells_max,
            n_variable_genes,
            pca_dims,
            graph_k,
            cluster_resolution,
            compute_umap,
            umap_neighbors,
            umap_epochs,
            num_threads,
            umap_parallel_optimization,
            normalization_scale,
            seed);
        celladmix::CellClusteringResult result;
        {
          // The native clustering path owns all temporary numeric buffers.
          py::gil_scoped_release release;
          const auto counts = celladmix::collect_run_counts(run_path);
          result = celladmix::cluster_cell_counts(counts, options);
        }
        return clustering_to_dict(result);
      },
      py::arg("run_path"),
      py::arg("min_molecules") = 10,
      py::arg("min_genes") = 5,
      py::arg("cells_max") = -1,
      py::arg("n_variable_genes") = 1000,
      py::arg("pca_dims") = 30,
      py::arg("graph_k") = 15,
      py::arg("cluster_resolution") = 1.0,
      py::arg("compute_umap") = true,
      py::arg("umap_neighbors") = 15,
      py::arg("umap_epochs") = 200,
      py::arg("num_threads") = 1,
      py::arg("umap_parallel_optimization") = true,
      py::arg("normalization_scale") = 5000.0,
      py::arg("seed") = 1U);

  m.def(
      "training_strata_for_cells",
      [](const std::string& store_dir,
         const std::vector<std::string>& annotation_cell_ids,
         const std::vector<std::string>& annotation_labels) {
        const auto counts = celladmix::load_input_store_counts(store_dir);
        return labels_for_cells(counts.cells, annotation_cell_ids, annotation_labels);
      },
      py::arg("store_dir"),
      py::arg("annotation_cell_ids"),
      py::arg("annotation_labels"));
}
