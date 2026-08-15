#include "celladmix/run_store.hpp"
#include "celladmix/simulation.hpp"
#include "test_framework.hpp"

#include <filesystem>
#include <numeric>

using namespace celladmix;

namespace {

std::filesystem::path unique_temp_dir(const std::string& stem) {
  auto path = std::filesystem::temp_directory_path() /
      (stem + "_" + std::to_string(std::rand()));
  std::filesystem::remove_all(path);
  std::filesystem::create_directories(path);
  return path;
}

}  // namespace

TEST_CASE("Run store writes and reloads a basic fit run") {
  SimulationParams params;
  params.transcripts_per_cell = 16;
  params.admixture_per_target_cell = 4;
  params.seed = 13U;
  const auto sim = simulate_nsclc_admixture(params);

  BasicPipelineOptions options;
  options.ncv_k = 6;
  options.rank = 2;
  options.graph_k = 4;
  options.nmf_iterations = 40;
  options.nmf_n_runs = 3;
  options.nmf_train_max_rows = 40;
  options.return_ncv = false;
  options.seed = 13U;
  const auto fit = run_basic_pipeline(sim.transcripts, options);

  const auto dir = unique_temp_dir("celladmix_run_store_fit");
  RunSourceInfo source;
  source.type = "synthetic";
  source.path = "simulation";

  RunStorageOptions storage;
  storage.tile_size = 10.0;
  storage.parquet_row_group_size = 16;

  const auto manifest = write_basic_run(
      dir.string(),
      source,
      sim.transcripts,
      fit,
      options,
      storage,
      std::nullopt,
      nullptr,
      &sim.cells);

  REQUIRE_EQ(manifest.n_transcripts, sim.transcripts.size());
  REQUIRE_EQ(manifest.n_factors, static_cast<std::size_t>(fit.nmf.h.rows()));
  REQUIRE_EQ(manifest.n_training_rows, fit.training_query_indices.size());

  const auto loaded_manifest = read_run_manifest(dir.string());
  REQUIRE_EQ(loaded_manifest.n_transcripts, manifest.n_transcripts);
  REQUIRE_EQ(loaded_manifest.n_factors, manifest.n_factors);
  REQUIRE_EQ(loaded_manifest.n_training_rows, manifest.n_training_rows);

  // Multirun stability diagnostics round-trip through the manifest JSON.
  REQUIRE_EQ(manifest.nmf_diagnostics.stability_comparison_runs, 2);
  REQUIRE_EQ(loaded_manifest.nmf_diagnostics.stability_metric, std::string{"ownership_matched"});
  REQUIRE_EQ(
      loaded_manifest.nmf_diagnostics.stability_comparison_runs,
      manifest.nmf_diagnostics.stability_comparison_runs);
  REQUIRE_EQ(
      loaded_manifest.nmf_diagnostics.stable_factor_count,
      manifest.nmf_diagnostics.stable_factor_count);
  REQUIRE_NEAR(
      loaded_manifest.nmf_diagnostics.stability_threshold,
      manifest.nmf_diagnostics.stability_threshold,
      1e-12);
  REQUIRE_EQ(
      loaded_manifest.nmf_diagnostics.selected_factor_stability.size(),
      manifest.nmf_diagnostics.selected_factor_stability.size());
  REQUIRE_EQ(
      loaded_manifest.nmf_diagnostics.candidate_best_match_correlations.size(),
      static_cast<std::size_t>(3));

  const auto training_obs_ids = load_run_training_obs_ids(dir.string());
  REQUIRE_EQ(training_obs_ids.size(), fit.training_query_indices.size());
  for (std::size_t i = 0; i < training_obs_ids.size(); ++i) {
    REQUIRE_EQ(training_obs_ids[i], static_cast<std::int64_t>(fit.training_query_indices[i]));
  }

  const auto run_data = load_run_data(dir.string());
  REQUIRE_EQ(run_data.transcripts.size(), sim.transcripts.size());
  REQUIRE_EQ(run_data.labels.size(), sim.transcripts.size());
  REQUIRE_EQ(run_data.cells.size(), sim.cells.size());

  const auto sparse_counts = collect_run_counts(dir.string());
  REQUIRE_EQ(sparse_counts.cells.size(), sim.cells.size());
  REQUIRE_EQ(sparse_counts.genes.size(), sim.transcripts.genes.size());
  REQUIRE_EQ(sparse_counts.indptr.size(), sim.cells.size() + 1U);
  REQUIRE_EQ(
      std::accumulate(sparse_counts.transcript_counts.begin(), sparse_counts.transcript_counts.end(), 0),
      static_cast<int>(sim.transcripts.size()));
  const auto dense_counts = build_cell_gene_counts(run_data.transcripts);
  for (std::size_t cell = 0; cell < sparse_counts.cells.size(); ++cell) {
    for (int pos = sparse_counts.indptr[cell]; pos < sparse_counts.indptr[cell + 1U]; ++pos) {
      const int gene = sparse_counts.indices[static_cast<std::size_t>(pos)];
      REQUIRE_EQ(
          sparse_counts.values[static_cast<std::size_t>(pos)],
          dense_counts(static_cast<int>(cell), gene));
    }
  }
}

TEST_CASE("Run store can write a corrected child run") {
  SimulationParams params;
  params.transcripts_per_cell = 16;
  params.admixture_per_target_cell = 4;
  params.seed = 17U;
  const auto sim = simulate_nsclc_admixture(params);

  BasicPipelineOptions options;
  options.ncv_k = 6;
  options.rank = 2;
  options.graph_k = 4;
  options.nmf_iterations = 40;
  options.nmf_train_max_rows = 40;
  options.return_ncv = false;
  options.seed = 17U;
  const auto fit = run_basic_pipeline(sim.transcripts, options);

  const auto parent_dir = unique_temp_dir("celladmix_run_store_parent");
  RunSourceInfo source;
  source.type = "synthetic";
  source.path = "simulation";

  write_basic_run(
      parent_dir.string(),
      source,
      sim.transcripts,
      fit,
      options,
      RunStorageOptions{},
      std::nullopt,
      nullptr,
      &sim.cells);

  const auto parent_run = load_run_data(parent_dir.string());
  std::vector<bool> keep_mask(parent_run.transcripts.size(), true);
  for (std::size_t i = 0; i < keep_mask.size(); ++i) {
    keep_mask[i] = parent_run.labels[i] != 0;
  }

  const auto corrected_dir = unique_temp_dir("celladmix_run_store_corrected");
  const auto corrected_manifest = write_corrected_run(
      corrected_dir.string(),
      parent_dir.string(),
      parent_run,
      keep_mask);

  REQUIRE_EQ(corrected_manifest.run_type, std::string("corrected"));
  REQUIRE(corrected_manifest.parent_run.has_value());
  REQUIRE_EQ(corrected_manifest.n_training_rows, fit.training_query_indices.size());

  const auto corrected_run = load_run_data(corrected_dir.string());
  REQUIRE(corrected_run.transcripts.size() <= parent_run.transcripts.size());
  REQUIRE_EQ(load_run_cells(corrected_dir.string()).size(), parent_run.cells.size());
  REQUIRE_EQ(load_run_training_obs_ids(corrected_dir.string()).size(), fit.training_query_indices.size());
  REQUIRE(std::filesystem::exists(correction_summary_parquet_path(corrected_dir.string())));

  const auto corrected_counts = collect_run_counts(corrected_dir.string());
  REQUIRE_EQ(corrected_counts.cells.size(), parent_run.cells.size());
  REQUIRE_EQ(
      std::accumulate(corrected_counts.transcript_counts.begin(), corrected_counts.transcript_counts.end(), 0),
      static_cast<int>(corrected_run.transcripts.size()));
}
