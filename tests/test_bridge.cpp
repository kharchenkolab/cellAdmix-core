#include "celladmix/bridge.hpp"
#include "test_framework.hpp"

using namespace celladmix;

TEST_CASE("Bridge evidence is stronger for admixture factor than native factor") {
  TranscriptTable table;
  table.gene_key = {
      "COL1A1", "COL1A1", "COL3A1", "KRT19", "KRT19", "KRT8",
      "KRT19", "KRT8", "KRT17", "KRT19", "KRT8", "KRT17"};
  table.orig_cell_id = {
      "fib", "fib", "fib", "fib", "fib", "fib",
      "mal", "mal", "mal", "mal", "mal", "mal"};
  table.cell_types = {
      "fibroblast", "fibroblast", "fibroblast", "fibroblast", "fibroblast", "fibroblast",
      "malignant", "malignant", "malignant", "malignant", "malignant", "malignant"};
  table.x = {60, 65, 70, 46, 47, 48, 40, 38, 36, 34, 32, 30};
  table.y = {20, 25, 30, 22, 28, 35, 20, 24, 28, 32, 36, 40};
  table.z = std::vector<double>(12, 0.0);
  table.finalize();

  const std::vector<int> labels = {
      0, 0, 0, 1, 1, 1,
      1, 1, 1, 1, 1, 1};

  const auto malignant_bridge = compute_bridge_evidence(table, labels, 0, 1, 1, 3);
  const auto fibro_bridge = compute_bridge_evidence(table, labels, 0, 1, 0, 3);

  REQUIRE_GT(malignant_bridge.crossing_fraction, 0.6);
  REQUIRE_GT(malignant_bridge.score, fibro_bridge.score);
  REQUIRE_GT(malignant_bridge.score, 0.2);
}

TEST_CASE("Adjacent cell discovery uses transcript-level cross-cell neighbors") {
  TranscriptTable table;
  table.gene_key = {"A", "A", "A", "A", "A", "A", "A", "A"};
  table.orig_cell_id = {"c1", "c1", "c2", "c2", "c3", "c3", "c3", "c3"};
  table.x = {0.0, 0.5, 2.0, 2.5, 20.0, 20.3, 20.6, 20.9};
  table.y = {0.0, 0.2, 0.0, 0.2, 0.0, 0.2, 0.0, 0.2};
  table.z = std::vector<double>(8, 0.0);
  table.finalize();

  const auto pairs = discover_adjacent_cell_pairs(table, 3);
  REQUIRE_EQ(pairs.size(), 1U);
  REQUIRE_EQ(pairs[0].cell_a, 0);
  REQUIRE_EQ(pairs[0].cell_b, 1);
}

TEST_CASE("Original-style bridge test returns ordered cell-type summaries") {
  TranscriptTable table;
  table.gene_key = {
      "COL1A1", "COL1A1", "KRT19", "KRT8",
      "KRT19", "KRT8", "KRT17", "KRT19"};
  table.orig_cell_id = {
      "fib1", "fib1", "fib1", "fib1",
      "mal1", "mal1", "mal1", "mal1"};
  table.cell_types = {
      "fibroblast", "fibroblast", "fibroblast", "fibroblast",
      "malignant", "malignant", "malignant", "malignant"};
  table.x = {60, 65, 46, 47, 40, 38, 36, 34};
  table.y = {20, 25, 22, 28, 20, 24, 28, 32};
  table.z = std::vector<double>(8, 0.0);
  table.finalize();

  CellTable cells;
  cells.cell_ids = table.cells;
  cells.cell_types = {"fibroblast", "malignant"};

  const std::vector<int> labels = {0, 0, 1, 1, 1, 1, 1, 1};
  BridgeTestOptions options;
  options.candidate_k = 3;
  options.crossing_k = 3;
  options.min_type_pair_contacts = 1;
  options.min_factor_molecules = 1;
  options.min_pairs = 1;
  options.null_iterations = 0;
  options.compute_null = false;

  const auto result = run_bridge_test(table, cells, labels, options);
  REQUIRE(!result.candidates.empty());
  REQUIRE(!result.pair_scores.empty());
  REQUIRE(!result.summaries.empty());

  const auto it = std::find_if(result.summaries.begin(), result.summaries.end(), [&](const auto& row) {
    return result.cell_types[static_cast<std::size_t>(row.target_type)] == "fibroblast" &&
        result.cell_types[static_cast<std::size_t>(row.source_type)] == "malignant" &&
        row.factor == 1;
  });
  REQUIRE(it != result.summaries.end());
  REQUIRE_GE(it->q75_score, 0.0);

  options.candidate_mode = "cell_center";
  options.cell_candidate_k = 1;
  options.candidate_pairs_per_type_pair = 10;
  const auto cell_center_result = run_bridge_test(table, cells, labels, options);
  REQUIRE(!cell_center_result.candidates.empty());
  REQUIRE(!cell_center_result.pair_scores.empty());
  REQUIRE(!cell_center_result.summaries.empty());
}

TEST_CASE("Fast bridge null crossing matches combined-kNN null crossing") {
  TranscriptTable table;
  table.gene_key = {
      "A", "A", "B", "B", "C", "C",
      "A", "B", "B", "C", "C", "C",
      "A", "A", "B", "B", "C", "C",
      "A", "B", "B", "C", "C", "C"};
  table.orig_cell_id = {
      "target1", "target1", "target1", "target1", "target1", "target1",
      "source1", "source1", "source1", "source1", "source1", "source1",
      "source2", "source2", "source2", "source2", "source2", "source2",
      "source3", "source3", "source3", "source3", "source3", "source3"};
  table.cell_types = {
      "target", "target", "target", "target", "target", "target",
      "source", "source", "source", "source", "source", "source",
      "source", "source", "source", "source", "source", "source",
      "source", "source", "source", "source", "source", "source"};
  table.x = {
      0.0, 0.2, 0.4, 0.1, 0.3, 0.5,
      1.0, 1.2, 1.4, 1.1, 1.3, 1.5,
      8.0, 8.2, 8.4, 8.1, 8.3, 8.5,
      -8.0, -8.2, -8.4, -8.1, -8.3, -8.5};
  table.y = {
      0.0, 0.1, 0.0, 0.4, 0.5, 0.4,
      0.0, 0.1, 0.0, 0.4, 0.5, 0.4,
      0.0, 0.1, 0.0, 0.4, 0.5, 0.4,
      0.0, 0.1, 0.0, 0.4, 0.5, 0.4};
  table.z = std::vector<double>(24, 0.0);
  table.finalize();

  CellTable cells;
  cells.cell_ids = table.cells;
  cells.cell_types = {"target", "source", "source", "source"};

  const std::vector<int> labels = {
      0, 0, 1, 1, 1, 0,
      0, 1, 1, 1, 0, 0,
      0, 0, 1, 1, 0, 1,
      0, 1, 1, 0, 0, 1};

  BridgeTestOptions options;
  options.candidate_k = 4;
  options.crossing_k = 4;
  options.min_type_pair_contacts = 1;
  options.min_factor_molecules = 1;
  options.min_pairs = 1;
  options.max_cells_per_type_pair = 20;
  options.null_iterations = 2;
  options.null_max_iterations = 12;
  options.num_threads = 1;
  options.seed = 5;
  options.compute_null = true;

  auto old_options = options;
  old_options.fast_null_crossing = false;
  auto fast_options = options;
  fast_options.fast_null_crossing = true;

  const auto old_result = run_bridge_test(table, cells, labels, old_options);
  const auto fast_result = run_bridge_test(table, cells, labels, fast_options);

  REQUIRE_EQ(old_result.summaries.size(), fast_result.summaries.size());
  for (std::size_t i = 0; i < old_result.summaries.size(); ++i) {
    const auto& old_row = old_result.summaries[i];
    const auto& fast_row = fast_result.summaries[i];
    REQUIRE_EQ(old_row.target_type, fast_row.target_type);
    REQUIRE_EQ(old_row.source_type, fast_row.source_type);
    REQUIRE_EQ(old_row.factor, fast_row.factor);
    REQUIRE_EQ(old_row.n_pairs, fast_row.n_pairs);
    REQUIRE_NEAR(old_row.mean_score, fast_row.mean_score, 1e-12);
    REQUIRE_NEAR(old_row.mean_null_score, fast_row.mean_null_score, 1e-12);
    REQUIRE_NEAR(old_row.q75_score, fast_row.q75_score, 1e-12);
    REQUIRE_NEAR(old_row.p_value, fast_row.p_value, 1e-12);
  }
}
