#include "celladmix/generative.hpp"
#include "celladmix/clustering.hpp"

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
  opt.dose_ladder.clear();  // the override permutes exposures only; see the
                            // covariate control below for the ladder
  const auto res = fit_generative(
      fx.indptr, fx.indices, fx.values, fx.n_genes, fx.cell_ids, fx.x, fx.y,
      fx.type_codes, 2, "", "", {fx.pair}, {}, {}, {}, opt);
  const double removed = removed_from_targets(fx, res, 0, 3);
  std::cerr << "[debug] permuted removed " << removed << " of "
            << fx.planted_transfer << "\n";
  REQUIRE_LT(removed, 0.05 * fx.planted_transfer);
}

// ---- dose covariates: ladder strata and joint estimation -----------------------

TEST_CASE("Nearest-type rank reports the neighbor position of the first cell of a type") {
  // Five cells on a line: type 0 at 0 and 40, type 1 at 10, 20, 30.
  const std::vector<double> x{0.0, 10.0, 20.0, 30.0, 40.0};
  const std::vector<double> y{0.0, 0.0, 0.0, 0.0, 0.0};
  const std::vector<int> types{0, 1, 1, 1, 0};
  const auto rank = cell_nearest_type_rank(x, y, types, 2, 3);
  // Cell at 30: neighbors in order 20 (type 1), 40 (type 0), 10 (type 1).
  REQUIRE_NEAR(rank(3, 1), 1.0, 1e-12);
  REQUIRE_NEAR(rank(3, 0), 2.0, 1e-12);
  // Cell at 20: neighbors 10, 30 (type 1), then 0 or 40 (type 0) at rank 3.
  REQUIRE_NEAR(rank(2, 0), 3.0, 1e-12);
  // Cell at 0 with max_k = 1: its single neighbor is type 1, so type 0 is
  // beyond the window and reports max_k + 1.
  const auto rank1 = cell_nearest_type_rank(x, y, types, 2, 1);
  REQUIRE_NEAR(rank1(0, 0), 2.0, 1e-12);
  REQUIRE_NEAR(rank1(0, 1), 1.0, 1e-12);
}

TEST_CASE("Additive isotonic fit separates two correlated monotone components") {
  // Rows carry two ordered covariates whose levels are correlated (level b
  // tends to follow level a); the rate is an additive function of both.
  std::mt19937 rng(11);
  std::uniform_int_distribution<int> da(0, 3);
  std::uniform_real_distribution<double> noise(-0.005, 0.005);
  const double fa[4] = {0.0, 0.02, 0.05, 0.09};
  const double fb[3] = {0.0, 0.03, 0.03};
  std::vector<double> yv;
  std::vector<double> wv;
  std::vector<int> sa;
  std::vector<int> sb;
  for (int i = 0; i < 4000; ++i) {
    const int a = da(rng);
    // b correlates with a: mostly the same band, sometimes shifted
    int b = std::min(a, 2);
    if (rng() % 3 == 0) b = std::max(0, b - 1);
    yv.push_back(0.10 + fa[a] + fb[b] + noise(rng));
    wv.push_back(1.0);
    sa.push_back(a);
    sb.push_back(b);
  }
  const auto fit = additive_isotonic_fit(yv, wv, {sa, sb}, {4, 3}, {0, 0});
  REQUIRE_NEAR(fit.intercept, 0.10, 0.01);
  for (int a = 0; a < 4; ++a) REQUIRE_NEAR(fit.components[0][a], fa[a], 0.01);
  for (int b = 0; b < 3; ++b) REQUIRE_NEAR(fit.components[1][b], fb[b], 0.01);
  // The marginal fit of the same rate on covariate b alone credits b with
  // part of a's effect, which the joint fit does not.
  std::vector<double> mb(3, 0.0);
  std::vector<double> wb(3, 0.0);
  for (std::size_t i = 0; i < yv.size(); ++i) {
    mb[static_cast<std::size_t>(sb[i])] += yv[i];
    wb[static_cast<std::size_t>(sb[i])] += 1.0;
  }
  for (int b = 0; b < 3; ++b) mb[static_cast<std::size_t>(b)] /= wb[static_cast<std::size_t>(b)];
  const auto marginal_b = pava_increasing(mb, wb);
  REQUIRE_GT(marginal_b[2] - marginal_b[0], fb[2] + 0.02);
}

namespace {

// Three types on a grid: a pure source-A region on the left, a mixed zone in
// which each cell is A or B at random, and the target region on the right.
// Transfer into a target cell comes from B only, in proportion to the number
// of B cells among its nearest neighbors, and B's profile carries a gene that
// A owns. Across columns both A and B exposures rise toward the zone, so A
// exposure co-varies with B's transfer; within a column the two counts vary
// independently, which is what lets a joint fit tell them apart.
struct ConfoundedFixture {
  std::vector<int> indptr{0};
  std::vector<int> indices;
  std::vector<double> values;
  std::vector<std::string> cell_ids;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<int> type_codes;
  int n_genes = 9;
  double planted_b = 0.0;
  GenerativePairSpec pair_a;
  GenerativePairSpec pair_b;
};

ConfoundedFixture make_confounded_fixture(int neighbor_k) {
  ConfoundedFixture fx;
  const int side = 32;
  std::mt19937 rng(23);
  std::vector<int> col_of;
  for (int iy = 0; iy < side; ++iy) {
    for (int ix = 0; ix < side; ++ix) {
      const int c = static_cast<int>(fx.cell_ids.size());
      fx.cell_ids.push_back("c" + std::to_string(c));
      fx.x.push_back(ix * 30.0);
      fx.y.push_back(iy * 30.0);
      col_of.push_back(ix);
      int type = 2;
      if (ix < 10) {
        type = 0;
      } else if (ix < 14) {
        type = (rng() % 2 == 0) ? 0 : 1;
      }
      fx.type_codes.push_back(type);
    }
  }
  const int n = static_cast<int>(fx.cell_ids.size());
  // B exposure per cell: B cells among the neighbor_k nearest (brute force).
  std::vector<int> b_count(static_cast<std::size_t>(n), 0);
  for (int c = 0; c < n; ++c) {
    std::vector<std::pair<double, int>> d;
    d.reserve(static_cast<std::size_t>(n));
    for (int o = 0; o < n; ++o) {
      if (o == c) continue;
      const double dx = fx.x[static_cast<std::size_t>(o)] - fx.x[static_cast<std::size_t>(c)];
      const double dy = fx.y[static_cast<std::size_t>(o)] - fx.y[static_cast<std::size_t>(c)];
      d.emplace_back(dx * dx + dy * dy, o);
    }
    std::partial_sort(d.begin(), d.begin() + neighbor_k, d.end());
    for (int h = 0; h < neighbor_k; ++h) {
      if (fx.type_codes[static_cast<std::size_t>(d[static_cast<std::size_t>(h)].second)] == 1) {
        ++b_count[static_cast<std::size_t>(c)];
      }
    }
  }
  // genes 0,1: A-owned (gene 1 also expressed by B); genes 2,3: B-owned;
  // genes 4..6: target-owned; gene 7: shared background.
  for (int c = 0; c < n; ++c) {
    std::vector<double> gene_counts(static_cast<std::size_t>(fx.n_genes), 0.0);
    const int type = fx.type_codes[static_cast<std::size_t>(c)];
    if (type == 0) {
      gene_counts[0] = 120.0;
      gene_counts[1] = 120.0;
    } else if (type == 1) {
      gene_counts[1] = 40.0;
      gene_counts[2] = 100.0;
      gene_counts[3] = 100.0;
    } else {
      for (int g = 4; g < 7; ++g) gene_counts[static_cast<std::size_t>(g)] = 60.0;
      const double dose = 0.025 * b_count[static_cast<std::size_t>(c)];
      const double prof[3] = {40.0, 100.0, 100.0};
      const int gidx[3] = {1, 2, 3};
      for (int k = 0; k < 3; ++k) {
        const double add = dose * prof[k];
        gene_counts[static_cast<std::size_t>(gidx[k])] += add;
        fx.planted_b += add;
      }
    }
    gene_counts[7] = 30.0;
    for (int g = 0; g < fx.n_genes; ++g) {
      if (gene_counts[static_cast<std::size_t>(g)] <= 0.0) continue;
      fx.indices.push_back(g);
      fx.values.push_back(gene_counts[static_cast<std::size_t>(g)]);
    }
    fx.indptr.push_back(static_cast<int>(fx.values.size()));
  }
  fx.pair_a.source_type = 0;
  fx.pair_a.target_type = 2;
  fx.pair_a.pool = {0, 1};
  fx.pair_a.strict = {0, 1};
  fx.pair_b.source_type = 1;
  fx.pair_b.target_type = 2;
  fx.pair_b.pool = {2, 3};
  fx.pair_b.strict = {2, 3};
  return fx;
}

}  // namespace

TEST_CASE("Joint dose does not credit a co-located source with another's transfer") {
  auto opt = test_options();
  const auto fx = make_confounded_fixture(opt.neighbor_k);
  opt.use_induced = false;
  opt.dose_mode = "marginal";
  opt.dose_ladder.clear();
  const auto marginal = fit_generative(
      fx.indptr, fx.indices, fx.values, fx.n_genes, fx.cell_ids, fx.x, fx.y,
      fx.type_codes, 3, "", "", {fx.pair_a, fx.pair_b}, {}, {}, {}, opt);
  opt.dose_mode = "joint";
  const auto joint = fit_generative(
      fx.indptr, fx.indices, fx.values, fx.n_genes, fx.cell_ids, fx.x, fx.y,
      fx.type_codes, 3, "", "", {fx.pair_a, fx.pair_b}, {}, {}, {}, opt);
  REQUIRE_EQ(static_cast<int>(marginal.pairs.size()), 2);
  REQUIRE_EQ(static_cast<int>(joint.pairs.size()), 2);
  std::cerr << "[debug] A prior marginal " << marginal.pairs[0].prior_molecules
            << " joint " << joint.pairs[0].prior_molecules
            << "; B prior marginal " << marginal.pairs[1].prior_molecules
            << " joint " << joint.pairs[1].prior_molecules
            << "; A cross-attributed " << joint.pairs[0].cross_attributed_molecules << "\n";
  // The marginal fit credits A with transfer it never made; the joint fit
  // attributes that content to B's exposure instead.
  REQUIRE_GT(marginal.pairs[0].prior_molecules, 0.1 * marginal.pairs[1].prior_molecules);
  REQUIRE_LT(joint.pairs[0].prior_molecules, 0.25 * marginal.pairs[0].prior_molecules);
  REQUIRE_GT(joint.pairs[0].cross_attributed_molecules, 0.0);
  // B's own dose is essentially unchanged.
  REQUIRE_GT(joint.pairs[1].prior_molecules, 0.8 * marginal.pairs[1].prior_molecules);
  REQUIRE_NEAR(marginal.pairs[0].cross_attributed_molecules, 0.0, 1e-12);
}

TEST_CASE("A dose ladder reaches transfer in cells with no source among the nearest neighbors") {
  // The fixture's transfer decays over three columns from the boundary, so
  // the outermost planted cells can sit with no source cell among their
  // nearest neighbors; with neighbor_k small, the ladder must find them.
  auto fx = make_fixture(false);
  auto opt = test_options();
  opt.neighbor_k = 6;
  opt.use_induced = false;
  opt.dose_reference_min_molecules = 2000.0;
  opt.dose_ladder.clear();
  const auto plain = fit_generative(
      fx.indptr, fx.indices, fx.values, fx.n_genes, fx.cell_ids, fx.x, fx.y,
      fx.type_codes, 2, "", "", {fx.pair}, {}, {}, {}, opt);
  opt.dose_ladder = {12, 24, 48};
  const auto ladder = fit_generative(
      fx.indptr, fx.indices, fx.values, fx.n_genes, fx.cell_ids, fx.x, fx.y,
      fx.type_codes, 2, "", "", {fx.pair}, {}, {}, {}, opt);
  const double removed_plain = removed_from_targets(fx, plain, 0, 3);
  const double removed_ladder = removed_from_targets(fx, ladder, 0, 3);
  std::cerr << "[debug] planted " << fx.planted_transfer << " removed plain "
            << removed_plain << " ladder " << removed_ladder
            << " unexposed prior " << ladder.pairs[0].prior_unexposed_molecules << "\n";
  REQUIRE_NEAR(plain.pairs[0].prior_unexposed_molecules, 0.0, 1e-9);
  REQUIRE_GT(ladder.pairs[0].prior_unexposed_molecules, 0.0);
  REQUIRE_GT(removed_ladder, removed_plain);
  REQUIRE_LE(removed_ladder, 1.05 * fx.planted_transfer);
  REQUIRE_NEAR(removed_from_targets(fx, ladder, 4, 6), 0.0, 1e-9);
}

TEST_CASE("Permuting the covariates inside the fit collapses removal under joint dose and a ladder") {
  auto fx = make_fixture(false);
  auto opt = test_options();
  opt.use_ambient = false;
  opt.use_induced = false;
  opt.dose_mode = "joint";
  opt.dose_ladder = {20, 40};
  opt.dose_reference_min_molecules = 2000.0;
  const auto res = fit_generative(
      fx.indptr, fx.indices, fx.values, fx.n_genes, fx.cell_ids, fx.x, fx.y,
      fx.type_codes, 2, "", "", {fx.pair}, {}, {}, {}, opt);
  opt.permute_covariates_seed = 5;
  const auto control = fit_generative(
      fx.indptr, fx.indices, fx.values, fx.n_genes, fx.cell_ids, fx.x, fx.y,
      fx.type_codes, 2, "", "", {fx.pair}, {}, {}, {}, opt);
  const double removed = removed_from_targets(fx, res, 0, 3);
  const double removed_control = removed_from_targets(fx, control, 0, 3);
  std::cerr << "[debug] joint+ladder removed " << removed << " control " << removed_control
            << "; control prior " << control.pairs[0].prior_molecules
            << " (unexposed " << control.pairs[0].prior_unexposed_molecules << ")"
            << " posterior " << control.pairs[0].posterior_molecules
            << "; fit prior " << res.pairs[0].prior_molecules << "\n";
  REQUIRE_GT(removed, 0.6 * fx.planted_transfer);
  REQUIRE_LT(removed_control, 0.05 * fx.planted_transfer);
}
