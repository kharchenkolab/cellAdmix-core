// Molecule-level removal from count-level amounts.
//
// A count-level correction (the generative model) says how many molecules
// of each gene to remove from each cell, as a fractional amount. Two
// operations turn such amounts into a keep mask over a run's molecules:
// limiting an existing molecule-level removal (the factor-based cleanup's
// vote) so that no cell and gene loses more than its amount, and drawing a
// removal of exactly the amount at random from each cell's molecules of the
// gene. Both round a fractional amount stochastically (the integer part
// plus one with probability equal to the fraction), so the expected number
// removed equals the amount.

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "celladmix/types.hpp"

namespace celladmix {

// Removal amounts per (cell identifier, gene name); entries absent from
// the lists have an amount of zero.
struct CountRemoval {
  std::vector<std::string> cell_ids;
  std::vector<std::string> genes;
  std::vector<double> amounts;
};

struct CountRemovalResult {
  std::vector<bool> keep_mask;
  // The summed fractional amounts over the cells and genes the operation
  // touched, the rounded limits, the molecules restored by a limit (zero
  // for a draw), and the molecules removed in the end.
  double amount_molecules = 0.0;
  std::int64_t limit_molecules = 0;
  std::int64_t restored = 0;
  std::int64_t removed = 0;
};

// Limit a removal to the amounts times scale: within each cell and gene the
// molecules the removal took (keep_mask false) beyond the rounded limit are
// restored, those with the fewest votes first and at random among equal
// votes. votes may be empty (all equal).
CountRemovalResult limit_removal(const TranscriptTable& transcripts,
                                 const std::vector<int>& votes,
                                 const std::vector<bool>& keep_mask,
                                 const CountRemoval& amounts,
                                 double scale,
                                 unsigned int seed);

// Remove the rounded amount of every cell and gene, drawn at random from
// the cell's molecules of that gene (all of them when the amount exceeds
// their number).
CountRemovalResult draw_removal(const TranscriptTable& transcripts,
                                const CountRemoval& amounts,
                                unsigned int seed);

}  // namespace celladmix
