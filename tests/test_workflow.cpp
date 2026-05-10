#include "celladmix/workflow.hpp"
#include "test_framework.hpp"

using namespace celladmix;

TEST_CASE("Workflow transcript subsetting honors crop ids and keep masks") {
  TranscriptTable table;
  table.gene_key = {"A", "B", "C", "D"};
  table.orig_cell_id = {"c1", "c1", "c2", "c2"};
  table.cell_types = {"fib", "fib", "mal", "mal"};
  table.transcript_ids = {"t1", "t2", "t3", "t4"};
  table.x = {0.0, 1.0, 10.0, 11.0};
  table.y = {0.0, 1.0, 10.0, 11.0};
  table.z = {0.0, 0.0, 0.0, 0.0};
  table.finalize();

  const std::vector<std::string> crop_ids = {"left", "left", "right", "right"};
  const std::vector<bool> keep_mask = {true, false, true, true};

  const auto subset = subset_transcripts(table, std::string("right"), &crop_ids, &keep_mask);
  REQUIRE_EQ(subset.size(), 2U);
  REQUIRE_EQ(subset.orig_cell_id[0], std::string("c2"));
  REQUIRE_EQ(subset.orig_cell_id[1], std::string("c2"));
  REQUIRE_EQ(subset.transcript_ids[0], std::string("t3"));
  REQUIRE_EQ(subset.transcript_ids[1], std::string("t4"));
}

TEST_CASE("Workflow cell summaries aggregate dominant factor by cell") {
  TranscriptTable table;
  table.gene_key = {"A", "A", "B", "B", "C", "C"};
  table.orig_cell_id = {"c1", "c1", "c1", "c2", "c2", "c2"};
  table.cell_types = {"fib", "fib", "fib", "mal", "mal", "mal"};
  table.x = {0.0, 1.0, 2.0, 10.0, 11.0, 12.0};
  table.y = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  table.z = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  table.finalize();

  const std::vector<int> labels = {0, 0, 1, 1, 1, 0};
  const auto summaries = summarize_cells(table, &labels);
  REQUIRE_EQ(summaries.size(), 2U);
  REQUIRE_EQ(summaries[0].cell_id, std::string("c1"));
  REQUIRE_EQ(summaries[0].cell_type, std::string("fib"));
  REQUIRE_EQ(summaries[0].dominant_factor, 0);
  REQUIRE_NEAR(summaries[0].dominant_fraction, 2.0 / 3.0, 1e-9);
  REQUIRE_EQ(summaries[1].cell_id, std::string("c2"));
  REQUIRE_EQ(summaries[1].dominant_factor, 1);
  REQUIRE_NEAR(summaries[1].dominant_fraction, 2.0 / 3.0, 1e-9);
}
