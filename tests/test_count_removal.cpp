#include "celladmix/count_removal.hpp"

#include <string>
#include <vector>

#include "test_framework.hpp"

using namespace celladmix;

namespace {

// Two cells, two genes: cell c1 has five molecules of g1 and two of g2,
// cell c2 has three molecules of g1.
TranscriptTable make_table() {
  TranscriptTable t;
  const std::vector<std::pair<std::string, std::string>> rows = {
      {"c1", "g1"}, {"c1", "g1"}, {"c1", "g1"}, {"c1", "g1"}, {"c1", "g1"},
      {"c1", "g2"}, {"c1", "g2"}, {"c2", "g1"}, {"c2", "g1"}, {"c2", "g1"}};
  for (const auto& [cell, gene] : rows) {
    t.orig_cell_id.push_back(cell);
    t.gene_key.push_back(gene);
    t.x.push_back(0.0);
    t.y.push_back(0.0);
  }
  return t;
}

int removed_in(const CountRemovalResult& r, std::size_t from, std::size_t to) {
  int n = 0;
  for (std::size_t i = from; i < to; ++i) n += r.keep_mask[i] ? 0 : 1;
  return n;
}

}  // namespace

TEST_CASE("A drawn removal takes the integer amounts and stays within each cell and gene") {
  const auto t = make_table();
  const CountRemoval amounts{{"c1", "c2", "c1"}, {"g1", "g1", "g2"}, {3.0, 10.0, 0.0}};
  const auto r = draw_removal(t, amounts, 7);
  REQUIRE_EQ(removed_in(r, 0, 5), 3);   // c1 g1: three of five
  REQUIRE_EQ(removed_in(r, 5, 7), 0);   // c1 g2: amount zero
  REQUIRE_EQ(removed_in(r, 7, 10), 3);  // c2 g1: all three (amount above the count)
  REQUIRE_EQ(r.removed, 6);
}

TEST_CASE("A fractional amount is rounded to its expectation over draws") {
  const auto t = make_table();
  const CountRemoval amounts{{"c1"}, {"g1"}, {2.25}};
  int total = 0;
  for (unsigned int seed = 0; seed < 400; ++seed) {
    const auto r = draw_removal(t, amounts, seed);
    const int k = removed_in(r, 0, 5);
    REQUIRE(k == 2 || k == 3);
    total += k;
  }
  const double mean = total / 400.0;
  REQUIRE(mean > 2.15 && mean < 2.35);
}

TEST_CASE("A limit restores the removed molecules with the fewest votes beyond the amount") {
  const auto t = make_table();
  // the vote removed four molecules of c1 g1 and both of c2's first two
  std::vector<bool> keep(10, true);
  std::vector<int> votes(10, 0);
  for (std::size_t i = 0; i < 4; ++i) {
    keep[i] = false;
    votes[i] = static_cast<int>(i) + 1;  // molecule 3 has the most votes
  }
  keep[7] = false;
  keep[8] = false;
  votes[7] = votes[8] = 2;
  const CountRemoval amounts{{"c1"}, {"g1"}, {2.0}};
  const auto r = limit_removal(t, votes, keep, amounts, 1.0, 3);
  // c1 g1 keeps its two most-voted removals (molecules 2 and 3)
  REQUIRE(r.keep_mask[0] && r.keep_mask[1]);
  REQUIRE(!r.keep_mask[2] && !r.keep_mask[3]);
  // c2 g1 has no amount: its removals are restored
  REQUIRE(r.keep_mask[7] && r.keep_mask[8]);
  REQUIRE_EQ(r.restored, 4);
  REQUIRE_EQ(r.removed, 2);
  // a scale of two lets all four of c1 g1 stand
  const auto r2 = limit_removal(t, votes, keep, amounts, 2.0, 3);
  REQUIRE_EQ(removed_in(r2, 0, 5), 4);
}
