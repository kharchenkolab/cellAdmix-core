#include "celladmix/generative.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <random>
#include <string>
#include <vector>

#include "test_framework.hpp"

using namespace celladmix;

namespace {

// Two cell types on a split grid: type 0 (source) on the left, type 1
// (target) on the right. Target cells express their own genes; source-owned
// gene content is planted into target cells in proportion to how close they
// sit to the boundary, sampling the source profile (transfer), except one
// low-share gene that receives a large disproportionate extra (induction).
struct Fixture {
  std::vector<int> indptr{0};
  std::vector<int> indices;
  std::vector<double> values;
  std::vector<std::string> cell_ids;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<int> type_codes;
  int n_genes = 8;
  double planted_transfer = 0.0;
  double planted_induced = 0.0;
  std::vector<double> target_exposure;  // per target cell, in column order
  GenerativePairSpec pair;
};

Fixture make_fixture(bool with_induction) {
  Fixture fx;
  // genes 0..3 source-owned (decreasing profile), genes 4..6 target-owned,
  // gene 3 is the low-share source gene that carries the induction.
  const double src_profile[4] = {400.0, 200.0, 100.0, 20.0};
  std::mt19937 rng(7);
  const int side = 32;
  for (int iy = 0; iy < side; ++iy) {
    for (int ix = 0; ix < side; ++ix) {
      const int c = static_cast<int>(fx.cell_ids.size());
      fx.cell_ids.push_back("c" + std::to_string(c));
      fx.x.push_back(ix * 30.0);
      fx.y.push_back(iy * 30.0);
      const bool source = ix < side / 2;
      fx.type_codes.push_back(source ? 0 : 1);
      std::vector<double> gene_counts(static_cast<std::size_t>(fx.n_genes), 0.0);
      if (source) {
        for (int g = 0; g < 4; ++g) gene_counts[static_cast<std::size_t>(g)] = src_profile[g] / 4.0;
      } else {
        for (int g = 4; g < 7; ++g) gene_counts[static_cast<std::size_t>(g)] = 60.0;
        // transfer: decays with distance from the boundary column
        const int dist = ix - side / 2;
        const double dose = std::max(0.0, 0.12 - 0.05 * dist);
        for (int g = 0; g < 4; ++g) {
          const double add = dose * src_profile[g];
          gene_counts[static_cast<std::size_t>(g)] += add;
          fx.planted_transfer += add;
        }
        if (with_induction && dist <= 1) {
          gene_counts[3] += 40.0;  // far beyond gene 3's transfer share
          fx.planted_induced += 40.0;
        }
      }
      // shared background gene 7 everywhere
      gene_counts[7] = 30.0;
      for (int g = 0; g < fx.n_genes; ++g) {
        if (gene_counts[static_cast<std::size_t>(g)] <= 0.0) continue;
        fx.indices.push_back(g);
        fx.values.push_back(gene_counts[static_cast<std::size_t>(g)]);
      }
      fx.indptr.push_back(static_cast<int>(fx.values.size()));
    }
  }
  fx.pair.source_type = 0;
  fx.pair.target_type = 1;
  fx.pair.pool = {0, 1, 2, 3};
  fx.pair.strict = {0, 1, 2};
  return fx;
}

double removed_from_targets(const Fixture& fx, const GenerativeResult& res,
                            int gene_lo, int gene_hi) {
  double out = 0.0;
  const int n_cells = static_cast<int>(fx.cell_ids.size());
  for (int c = 0; c < n_cells; ++c) {
    if (fx.type_codes[static_cast<std::size_t>(c)] != 1) continue;
    for (int p = fx.indptr[c]; p < fx.indptr[c + 1]; ++p) {
      const int g = fx.indices[static_cast<std::size_t>(p)];
      if (g >= gene_lo && g <= gene_hi) out += res.removed[static_cast<std::size_t>(p)];
    }
  }
  return out;
}

GenerativeOptions test_options() {
  GenerativeOptions opt;
  opt.neighbor_k = 14;
  opt.near_min_molecules = 100.0;
  opt.far_min_molecules = 1e18;    // force the exposure-0 ambient fallback
  opt.induced_min_excess = 20.0;
  opt.num_threads = 2;
  return opt;
}

}  // namespace

TEST_CASE("PAVA produces a monotone weighted fit") {
  const std::vector<double> y{1.0, 0.5, 2.0, 1.5, 3.0};
  const std::vector<double> w{1.0, 1.0, 1.0, 1.0, 1.0};
  const auto fit = pava_increasing(y, w);
  for (std::size_t i = 1; i < fit.size(); ++i) REQUIRE_GE(fit[i], fit[i - 1]);
  REQUIRE_NEAR(fit[0], 0.75, 1e-9);
  REQUIRE_NEAR(fit[1], 0.75, 1e-9);
  REQUIRE_NEAR(fit[4], 3.0, 1e-9);
}

TEST_CASE("Generative fit removes planted transfer, never target-owned genes") {
  const auto fx = make_fixture(false);
  const auto res = fit_generative(
      fx.indptr, fx.indices, fx.values, fx.n_genes, fx.cell_ids, fx.x, fx.y,
      fx.type_codes, 2, "", "", {fx.pair}, {}, {}, {}, test_options());
  const double removed = removed_from_targets(fx, res, 0, 3);
  std::cerr << "[debug] transfer removed " << removed << " of "
            << fx.planted_transfer << "\n";
  REQUIRE_GT(removed, 0.6 * fx.planted_transfer);
  // target-owned genes are structurally untouchable
  REQUIRE_NEAR(removed_from_targets(fx, res, 4, 6), 0.0, 1e-9);
  REQUIRE_EQ(static_cast<int>(res.pairs.size()), 1);
  REQUIRE_GT(res.pairs[0].posterior_molecules, 0.0);
}

TEST_CASE("Generative fit flags and retains disproportionate induction") {
  const auto fx = make_fixture(true);
  const auto res = fit_generative(
      fx.indptr, fx.indices, fx.values, fx.n_genes, fx.cell_ids, fx.x, fx.y,
      fx.type_codes, 2, "", "", {fx.pair}, {}, {}, {}, test_options());
  bool flagged = false;
  for (const auto& row : res.induced) {
    std::cerr << "[debug] induced gene " << row.gene << " excess " << row.excess
              << " expected " << row.expected << " z " << row.z << "\n";
    flagged = flagged || row.gene == 3;
  }
  REQUIRE(flagged);
  // the induced gene keeps most of its content; without the induced term it
  // would be removed like the other source genes
  const double removed_ind = removed_from_targets(fx, res, 3, 3);
  auto opt_noind = test_options();
  opt_noind.use_induced = false;
  const auto res_noind = fit_generative(
      fx.indptr, fx.indices, fx.values, fx.n_genes, fx.cell_ids, fx.x, fx.y,
      fx.type_codes, 2, "", "", {fx.pair}, {}, {}, {}, opt_noind);
  const double removed_noind = removed_from_targets(fx, res_noind, 3, 3);
  REQUIRE_LT(removed_ind, 0.5 * removed_noind);
}

TEST_CASE("Cluster and explicit program initializations behave") {
  const auto fx = make_fixture(false);
  auto opt = test_options();
  opt.init_mode = "clusters";
  opt.n_programs = 3;
  const auto res = fit_generative(
      fx.indptr, fx.indices, fx.values, fx.n_genes, fx.cell_ids, fx.x, fx.y,
      fx.type_codes, 2, "", "", {fx.pair}, {}, {}, {}, opt);
  REQUIRE_GT(removed_from_targets(fx, res, 0, 3), 0.6 * fx.planted_transfer);
  REQUIRE_NEAR(removed_from_targets(fx, res, 4, 6), 0.0, 1e-9);
  REQUIRE_GE(res.n_programs_used[1], 1);
  // explicit programs: the true target profile, one program
  std::vector<double> prog(static_cast<std::size_t>(2 * fx.n_genes), 0.0);
  for (int g = 4; g < 7; ++g) prog[static_cast<std::size_t>(fx.n_genes + g)] = 60.0;
  prog[static_cast<std::size_t>(fx.n_genes + 7)] = 30.0;
  for (int g = 0; g < 4; ++g) prog[static_cast<std::size_t>(g)] = 1.0;
  prog[7] = 30.0;
  const auto res2 = fit_generative(
      fx.indptr, fx.indices, fx.values, fx.n_genes, fx.cell_ids, fx.x, fx.y,
      fx.type_codes, 2, "", "", {fx.pair}, {}, prog, {0, 1}, test_options());
  REQUIRE_GT(removed_from_targets(fx, res2, 0, 3), 0.6 * fx.planted_transfer);
}

TEST_CASE("The no-retention split removes the induced share") {
  const auto fx = make_fixture(true);
  const auto res = fit_generative(
      fx.indptr, fx.indices, fx.values, fx.n_genes, fx.cell_ids, fx.x, fx.y,
      fx.type_codes, 2, "", "", {fx.pair}, {}, {}, {}, test_options());
  double kept3 = 0.0;
  double kept3_strict = 0.0;
  const int n_cells = static_cast<int>(fx.cell_ids.size());
  for (int c = 0; c < n_cells; ++c) {
    if (fx.type_codes[static_cast<std::size_t>(c)] != 1) continue;
    for (int p = fx.indptr[c]; p < fx.indptr[c + 1]; ++p) {
      if (fx.indices[static_cast<std::size_t>(p)] != 3) continue;
      kept3 += fx.values[static_cast<std::size_t>(p)] - res.removed[static_cast<std::size_t>(p)];
      kept3_strict += fx.values[static_cast<std::size_t>(p)] -
          res.removed_without_retention[static_cast<std::size_t>(p)];
    }
  }
  REQUIRE_LT(kept3_strict, 0.6 * kept3);
}

TEST_CASE("Permuted exposure collapses generative removal") {
  auto fx = make_fixture(false);
  // exposure override: the true exposures, shuffled across target cells
  std::vector<double> expo;
  {
    // reproduce column order of target cells and their true exposures via a
    // second fit? Not needed: pass a constant permutation seed and shuffle
    // the boundary-distance-derived exposure levels.
    std::vector<int> target_cols;
    for (std::size_t c = 0; c < fx.cell_ids.size(); ++c) {
      if (fx.type_codes[c] == 1) target_cols.push_back(static_cast<int>(c));
    }
    expo.resize(target_cols.size());
    for (std::size_t i = 0; i < target_cols.size(); ++i) {
      const double xi = fx.x[static_cast<std::size_t>(target_cols[i])];
      expo[i] = xi < 17 * 30.0 ? 3.0 : (xi < 18 * 30.0 ? 1.0 : 0.0);
    }
    std::mt19937 rng(3);
    std::shuffle(expo.begin(), expo.end(), rng);
  }
  fx.pair.exposure = expo;
  auto opt = test_options();
  opt.use_ambient = false;  // control convention: no exposure-independent removal
  const auto res = fit_generative(
      fx.indptr, fx.indices, fx.values, fx.n_genes, fx.cell_ids, fx.x, fx.y,
      fx.type_codes, 2, "", "", {fx.pair}, {}, {}, {}, opt);
  const double removed = removed_from_targets(fx, res, 0, 3);
  std::cerr << "[debug] permuted removed " << removed << " of "
            << fx.planted_transfer << "\n";
  REQUIRE_LT(removed, 0.05 * fx.planted_transfer);
}
