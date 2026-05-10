#include "celladmix/graph.hpp"
#include "test_framework.hpp"

using namespace celladmix;

TEST_CASE("ICM smoothing removes a local factor flip in a chain graph") {
  DenseMatrix scores(5, 2, 0.0);
  scores(0, 0) = 0.9; scores(0, 1) = 0.1;
  scores(1, 0) = 0.8; scores(1, 1) = 0.2;
  scores(2, 0) = 0.3; scores(2, 1) = 0.7;
  scores(3, 0) = 0.85; scores(3, 1) = 0.15;
  scores(4, 0) = 0.8; scores(4, 1) = 0.2;

  KnnGraph graph;
  graph.n_nodes = 5;
  graph.indptr = {0, 1, 3, 5, 7, 8};
  graph.indices = {
      1,
      0, 2,
      1, 3,
      2, 4,
      3};

  const auto labels = smooth_labels_icm(scores, graph, 5.0, 20);
  REQUIRE_EQ(labels.size(), 5U);
  for (const int label : labels) {
    REQUIRE_EQ(label, 0);
  }
}

TEST_CASE("Per-cell factor assignment matches between serial and threaded paths") {
  TranscriptTable table;
  table.gene_key = {"A", "A", "B", "B", "A", "B", "A", "B"};
  table.orig_cell_id = {"cell1", "cell1", "cell1", "cell1", "cell2", "cell2", "cell2", "cell2"};
  table.x = {0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0};
  table.y = {0.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0, 5.0};
  table.z = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  table.finalize();

  DenseMatrix scores(static_cast<int>(table.size()), 2, 0.0);
  for (int i = 0; i < scores.rows(); ++i) {
    scores(i, 0) = 0.85;
    scores(i, 1) = 0.15;
  }
  scores(2, 0) = 0.35; scores(2, 1) = 0.65;
  scores(6, 0) = 0.40; scores(6, 1) = 0.60;

  const auto serial = assign_factors_per_cell(table, scores, 2, 5.0, 20, 1);
  const auto threaded = assign_factors_per_cell(table, scores, 2, 5.0, 20, 4);
  REQUIRE_EQ(serial.size(), threaded.size());
  for (std::size_t i = 0; i < serial.size(); ++i) {
    REQUIRE_EQ(serial[i], threaded[i]);
  }
}
