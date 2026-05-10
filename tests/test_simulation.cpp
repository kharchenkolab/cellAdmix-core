#include <algorithm>
#include <unordered_set>
#include <vector>

#include "celladmix/bridge.hpp"
#include "celladmix/pipeline.hpp"
#include "celladmix/simulation.hpp"
#include "test_framework.hpp"

using namespace celladmix;

TEST_CASE("Integrated simulation reduces malignant markers in fibroblasts after correction") {
  SimulationParams params;
  params.transcripts_per_cell = 80;
  params.admixture_per_target_cell = 24;
  params.seed = 17;

  const auto sim = simulate_nsclc_admixture(params);

  BasicPipelineOptions options;
  options.ncv_k = 12;
  options.rank = 2;
  options.graph_k = 6;
  options.same_label_ratio = 4.0;
  options.nmf_iterations = 250;
  options.seed = 17;

  const auto fit = run_basic_pipeline(sim.transcripts, options);
  const int malignant_factor = select_factor_by_marker_sum(
      fit.nmf.h, sim.transcripts.genes, sim.malignant_markers);
  const auto keep_mask = apply_removal_rule(
      sim.transcripts, fit.labels, malignant_factor, "fibroblast");

  const auto raw_counts = build_cell_gene_counts(sim.transcripts);
  const auto corrected_counts = build_cell_gene_counts(sim.transcripts, keep_mask);

  std::vector<int> malignant_gene_indices;
  for (const auto& marker : sim.malignant_markers) {
    const auto it = std::find(sim.transcripts.genes.begin(), sim.transcripts.genes.end(), marker);
    REQUIRE(it != sim.transcripts.genes.end());
    malignant_gene_indices.push_back(static_cast<int>(std::distance(sim.transcripts.genes.begin(), it)));
  }

  double fibro_raw = 0.0;
  double fibro_corrected = 0.0;
  for (int cell = 0; cell < static_cast<int>(sim.transcripts.num_cells()); ++cell) {
    if (sim.cells.cell_types[static_cast<std::size_t>(cell)] != "fibroblast") {
      continue;
    }
    for (const int gene : malignant_gene_indices) {
      fibro_raw += raw_counts(cell, gene);
      fibro_corrected += corrected_counts(cell, gene);
    }
  }

  REQUIRE_GT(fibro_raw, fibro_corrected);
  REQUIRE_LT(fibro_corrected / fibro_raw, 0.75);

  const int fibro_factor = malignant_factor == 0 ? 1 : 0;
  const auto bridge = compute_bridge_evidence(
      sim.transcripts,
      fit.labels,
      sim.fibro_cell,
      sim.malignant_cell,
      malignant_factor,
      5);
  const auto fibro_bridge = compute_bridge_evidence(
      sim.transcripts,
      fit.labels,
      sim.fibro_cell,
      sim.malignant_cell,
      fibro_factor,
      5);
  REQUIRE_GE(bridge.score, fibro_bridge.score);

  std::vector<std::pair<double, double>> admixture_points;
  std::vector<std::pair<double, double>> control_points;
  for (std::size_t i = 0; i < sim.transcripts.size(); ++i) {
    if (sim.transcripts.orig_cell_id[i] != "fibroblast_1") {
      continue;
    }
    if (sim.true_admixture[i] && admixture_points.size() < 8U) {
      admixture_points.emplace_back(sim.transcripts.x[i], sim.transcripts.y[i]);
    } else if (!sim.true_admixture[i] && sim.transcripts.gene_key[i] == "COL1A1" &&
               sim.transcripts.x[i] > 60.0 &&
               control_points.size() < 8U) {
      control_points.emplace_back(sim.transcripts.x[i], sim.transcripts.y[i]);
    }
  }
  REQUIRE_EQ(admixture_points.size(), control_points.size());
  const double membrane = membrane_score(
      sim.membrane_image,
      admixture_points,
      control_points,
      {sim.cells.centroid_x[static_cast<std::size_t>(sim.fibro_cell)],
       sim.cells.centroid_y[static_cast<std::size_t>(sim.fibro_cell)]},
      0.5);
  REQUIRE_GT(membrane, 0.0);
}

TEST_CASE("Sampled pipeline training still produces full transcript labels without materializing NCV") {
  SimulationParams params;
  params.transcripts_per_cell = 60;
  params.admixture_per_target_cell = 18;
  params.seed = 23;

  const auto sim = simulate_nsclc_admixture(params);

  BasicPipelineOptions options;
  options.ncv_k = 10;
  options.rank = 2;
  options.graph_k = 5;
  options.same_label_ratio = 4.0;
  options.nmf_iterations = 200;
  options.nmf_train_max_rows = 80;
  options.return_ncv = false;
  options.seed = 23;

  const auto fit = run_basic_pipeline(sim.transcripts, options);
  REQUIRE_EQ(fit.training_query_indices.size(), 80U);
  REQUIRE_EQ(fit.ncv.rows(), 0);
  REQUIRE_EQ(fit.ncv.cols(), 0);
  REQUIRE_EQ(fit.factor_scores.rows(), static_cast<int>(sim.transcripts.size()));
  REQUIRE_EQ(fit.factor_scores.cols(), options.rank);
  REQUIRE_EQ(fit.labels.size(), sim.transcripts.size());

  const int malignant_factor = select_factor_by_marker_sum(
      fit.nmf.h, sim.transcripts.genes, sim.malignant_markers);
  const auto keep_mask = apply_removal_rule(
      sim.transcripts, fit.labels, malignant_factor, "fibroblast");
  const auto raw_counts = build_cell_gene_counts(sim.transcripts);
  const auto corrected_counts = build_cell_gene_counts(sim.transcripts, keep_mask);

  double fibro_raw = 0.0;
  double fibro_corrected = 0.0;
  for (int cell = 0; cell < static_cast<int>(sim.transcripts.num_cells()); ++cell) {
    if (sim.cells.cell_types[static_cast<std::size_t>(cell)] != "fibroblast") {
      continue;
    }
    for (int gene = 0; gene < static_cast<int>(sim.transcripts.num_genes()); ++gene) {
      const auto& name = sim.transcripts.genes[static_cast<std::size_t>(gene)];
      if (std::find(sim.malignant_markers.begin(), sim.malignant_markers.end(), name) ==
          sim.malignant_markers.end()) {
        continue;
      }
      fibro_raw += raw_counts(cell, gene);
      fibro_corrected += corrected_counts(cell, gene);
    }
  }

  REQUIRE_GT(fibro_raw, fibro_corrected);
}

TEST_CASE("Cell-centric sampled training preserves cell type coverage under a tight row budget") {
  SimulationParams params;
  params.transcripts_per_cell = 40;
  params.admixture_per_target_cell = 8;
  params.seed = 31;

  const auto sim = simulate_nsclc_admixture(params);

  BasicPipelineOptions options;
  options.ncv_k = 8;
  options.rank = 2;
  options.graph_k = 4;
  options.nmf_iterations = 80;
  options.nmf_train_max_rows = 2;
  options.return_ncv = false;
  options.seed = 31;

  const auto fit = run_basic_pipeline(sim.transcripts, options);
  REQUIRE_EQ(fit.training_query_indices.size(), 2U);

  std::unordered_set<std::string> sampled_cell_types;
  for (const int idx : fit.training_query_indices) {
    sampled_cell_types.insert(sim.transcripts.cell_types[static_cast<std::size_t>(idx)]);
  }
  REQUIRE(sampled_cell_types.count("malignant") == 1U);
  REQUIRE(sampled_cell_types.count("fibroblast") == 1U);
}

TEST_CASE("Minimum training molecules excludes undersized cells from NMF training") {
  SimulationParams params;
  params.transcripts_per_cell = 40;
  params.admixture_per_target_cell = 8;
  params.seed = 37;

  const auto sim = simulate_nsclc_admixture(params);

  BasicPipelineOptions options;
  options.ncv_k = 8;
  options.rank = 2;
  options.graph_k = 4;
  options.nmf_iterations = 80;
  options.nmf_train_max_rows = 16;
  options.nmf_min_molecules = 45;
  options.return_ncv = false;
  options.seed = 37;

  const auto fit = run_basic_pipeline(sim.transcripts, options);
  REQUIRE_EQ(fit.training_query_indices.size(), 16U);
  for (const int idx : fit.training_query_indices) {
    REQUIRE_EQ(sim.transcripts.cell_types[static_cast<std::size_t>(idx)], std::string("fibroblast"));
  }
}
