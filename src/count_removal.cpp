#include "celladmix/count_removal.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_map>

namespace celladmix {

namespace {

std::string molecule_cell(const TranscriptTable& t, std::size_t i) {
  if (!t.cell_index.empty() && t.cell_index[i] >= 0 &&
      static_cast<std::size_t>(t.cell_index[i]) < t.cells.size()) {
    return t.cells[static_cast<std::size_t>(t.cell_index[i])];
  }
  return t.orig_cell_id[i];
}

std::string molecule_gene(const TranscriptTable& t, std::size_t i) {
  if (!t.gene_key.empty()) return t.gene_key[i];
  return t.genes[static_cast<std::size_t>(t.gene_index[i])];
}

// The amount of every molecule's cell and gene (negative when absent).
std::vector<double> molecule_amounts(const TranscriptTable& t, const CountRemoval& a,
                                     std::vector<std::int64_t>& key_of) {
  if (a.cell_ids.size() != a.genes.size() || a.genes.size() != a.amounts.size()) {
    throw std::runtime_error("count removal: cell_ids, genes and amounts must have equal length");
  }
  std::unordered_map<std::string, int> cell_code;
  std::unordered_map<std::string, int> gene_code;
  for (const auto& c : a.cell_ids) cell_code.emplace(c, static_cast<int>(cell_code.size()));
  for (const auto& g : a.genes) gene_code.emplace(g, static_cast<int>(gene_code.size()));
  const std::int64_t n_genes = static_cast<std::int64_t>(gene_code.size());
  std::unordered_map<std::int64_t, double> amount_of;
  amount_of.reserve(a.amounts.size() * 2);
  for (std::size_t r = 0; r < a.amounts.size(); ++r) {
    const std::int64_t key = static_cast<std::int64_t>(cell_code[a.cell_ids[r]]) * n_genes +
        gene_code[a.genes[r]];
    amount_of[key] += std::max(a.amounts[r], 0.0);
  }
  const std::size_t n = t.size();
  std::vector<double> out(n, -1.0);
  key_of.assign(n, -1);
  for (std::size_t i = 0; i < n; ++i) {
    const auto c = cell_code.find(molecule_cell(t, i));
    if (c == cell_code.end()) continue;
    const auto g = gene_code.find(molecule_gene(t, i));
    if (g == gene_code.end()) continue;
    const std::int64_t key = static_cast<std::int64_t>(c->second) * n_genes + g->second;
    const auto found = amount_of.find(key);
    if (found == amount_of.end()) continue;
    out[i] = found->second;
    key_of[i] = key;
  }
  return out;
}

// Molecule indices grouped by key, in a random order within each group.
std::unordered_map<std::int64_t, std::vector<std::size_t>> groups_of(
    const std::vector<std::int64_t>& key_of, const std::vector<char>& member, std::mt19937_64& rng) {
  std::vector<std::size_t> order(key_of.size());
  std::iota(order.begin(), order.end(), 0);
  std::shuffle(order.begin(), order.end(), rng);
  std::unordered_map<std::int64_t, std::vector<std::size_t>> groups;
  for (const std::size_t i : order) {
    if (member[i] && key_of[i] >= 0) groups[key_of[i]].push_back(i);
  }
  return groups;
}

std::int64_t rounded(double amount, std::mt19937_64& rng) {
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  const double base = std::floor(amount);
  return static_cast<std::int64_t>(base) + (unit(rng) < amount - base ? 1 : 0);
}

}  // namespace

CountRemovalResult limit_removal(const TranscriptTable& transcripts,
                                 const std::vector<int>& votes,
                                 const std::vector<bool>& keep_mask,
                                 const CountRemoval& amounts,
                                 double scale,
                                 unsigned int seed) {
  const std::size_t n = transcripts.size();
  if (keep_mask.size() != n || (!votes.empty() && votes.size() != n)) {
    throw std::runtime_error("limit_removal: keep_mask and votes must match the molecules");
  }
  if (!(scale >= 0.0)) throw std::runtime_error("limit_removal: scale must be non-negative");
  std::vector<std::int64_t> key_of;
  const auto amount = molecule_amounts(transcripts, amounts, key_of);
  CountRemovalResult result;
  result.keep_mask = keep_mask;
  std::vector<char> removed(n, 0);
  for (std::size_t i = 0; i < n; ++i) removed[i] = keep_mask[i] ? 0 : 1;
  // molecules removed in a cell and gene without an amount have a limit of
  // zero; give them a (negative) key of their own per cell and gene
  std::unordered_map<std::string, std::int64_t> extra_key;
  for (std::size_t i = 0; i < n; ++i) {
    if (removed[i] && key_of[i] < 0) {
      const std::string k = molecule_cell(transcripts, i) + '\t' + molecule_gene(transcripts, i);
      const auto it = extra_key.emplace(k, -2 - static_cast<std::int64_t>(extra_key.size())).first;
      key_of[i] = it->second;
    }
  }
  std::mt19937_64 rng(seed);
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::shuffle(order.begin(), order.end(), rng);
  std::unordered_map<std::int64_t, std::vector<std::size_t>> groups;
  for (const std::size_t i : order) {
    if (removed[i]) groups[key_of[i]].push_back(i);
  }
  for (auto& [key, idx] : groups) {
    double a = 0.0;
    if (key >= 0) a = amount[idx.front()] * scale;
    const std::int64_t limit = key >= 0 ? rounded(a, rng) : 0;
    result.amount_molecules += a;
    result.limit_molecules += limit;
    if (static_cast<std::int64_t>(idx.size()) <= limit) continue;
    // most votes first; the shuffle above breaks ties at random
    if (!votes.empty()) {
      std::stable_sort(idx.begin(), idx.end(), [&](std::size_t x, std::size_t y) {
        return votes[x] > votes[y];
      });
    }
    for (std::size_t r = static_cast<std::size_t>(std::max<std::int64_t>(limit, 0)); r < idx.size(); ++r) {
      result.keep_mask[idx[r]] = true;
      ++result.restored;
    }
  }
  for (std::size_t i = 0; i < n; ++i) result.removed += result.keep_mask[i] ? 0 : 1;
  return result;
}

CountRemovalResult draw_removal(const TranscriptTable& transcripts,
                                const CountRemoval& amounts,
                                unsigned int seed) {
  const std::size_t n = transcripts.size();
  std::vector<std::int64_t> key_of;
  const auto amount = molecule_amounts(transcripts, amounts, key_of);
  std::mt19937_64 rng(seed);
  std::vector<char> member(n, 1);
  auto groups = groups_of(key_of, member, rng);
  CountRemovalResult result;
  result.keep_mask.assign(n, true);
  for (auto& [key, idx] : groups) {
    const double a = amount[idx.front()];
    if (a <= 0.0) continue;
    const std::int64_t count = std::min<std::int64_t>(rounded(a, rng), static_cast<std::int64_t>(idx.size()));
    result.amount_molecules += a;
    result.limit_molecules += count;
    for (std::int64_t r = 0; r < count; ++r) result.keep_mask[idx[static_cast<std::size_t>(r)]] = false;
    result.removed += count;
  }
  return result;
}

}  // namespace celladmix
