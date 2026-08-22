#include "celladmix/generative.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <random>
#include <thread>
#include <unordered_map>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>

#include "celladmix/clustering.hpp"

namespace celladmix {

namespace {

// ---- small utilities --------------------------------------------------------

template <typename T>
T arrow_unwrap(arrow::Result<T>&& result, const char* context) {
  if (!result.ok()) {
    throw std::runtime_error(std::string(context) + ": " + result.status().ToString());
  }
  return std::move(result).ValueUnsafe();
}

void arrow_check(const arrow::Status& status, const char* context) {
  if (!status.ok()) {
    throw std::runtime_error(std::string(context) + ": " + status.ToString());
  }
}

std::shared_ptr<arrow::Table> read_parquet(const std::string& path) {
  auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path), "Open parquet input");
  parquet::arrow::FileReaderBuilder builder;
  arrow_check(builder.Open(input), "Open parquet reader");
  auto reader = arrow_unwrap(builder.Build(), "Build parquet reader");
  std::shared_ptr<arrow::Table> table;
  arrow_check(reader->ReadTable(&table), "Read parquet table");
  return table;
}

bool table_has_column(const arrow::Table& table, const std::string& name) {
  return table.schema()->GetFieldIndex(name) >= 0;
}

std::vector<int> extract_int_column(const arrow::Table& table, const std::string& name) {
  const auto column = table.GetColumnByName(name);
  if (!column) {
    throw std::runtime_error("generative: missing parquet column " + name);
  }
  std::vector<int> out;
  out.reserve(static_cast<std::size_t>(table.num_rows()));
  for (const auto& chunk : column->chunks()) {
    if (chunk->type_id() == arrow::Type::INT32) {
      const auto& arr = static_cast<const arrow::Int32Array&>(*chunk);
      for (std::int64_t i = 0; i < arr.length(); ++i) out.push_back(arr.Value(i));
    } else if (chunk->type_id() == arrow::Type::INT64) {
      const auto& arr = static_cast<const arrow::Int64Array&>(*chunk);
      for (std::int64_t i = 0; i < arr.length(); ++i) out.push_back(static_cast<int>(arr.Value(i)));
    } else {
      throw std::runtime_error("generative: unexpected type for column " + name);
    }
  }
  return out;
}

std::vector<std::string> extract_string_column(const arrow::Table& table, const std::string& name) {
  const auto column = table.GetColumnByName(name);
  if (!column) {
    throw std::runtime_error("generative: missing parquet column " + name);
  }
  std::vector<std::string> out;
  out.reserve(static_cast<std::size_t>(table.num_rows()));
  for (const auto& chunk : column->chunks()) {
    if (chunk->type_id() == arrow::Type::STRING) {
      const auto& arr = static_cast<const arrow::StringArray&>(*chunk);
      for (std::int64_t i = 0; i < arr.length(); ++i) out.push_back(arr.GetString(i));
    } else if (chunk->type_id() == arrow::Type::LARGE_STRING) {
      const auto& arr = static_cast<const arrow::LargeStringArray&>(*chunk);
      for (std::int64_t i = 0; i < arr.length(); ++i) out.push_back(arr.GetString(i));
    } else {
      throw std::runtime_error("generative: unexpected type for column " + name);
    }
  }
  return out;
}

// Split [0, n) into contiguous chunks and run fn(begin, end, thread_index)
// on each from its own thread.
void parallel_rows(int n, int num_threads, const std::function<void(int, int, int)>& fn) {
  const int workers = std::max(1, std::min(num_threads, n));
  if (workers == 1) {
    fn(0, n, 0);
    return;
  }
  std::vector<std::thread> threads;
  threads.reserve(static_cast<std::size_t>(workers));
  const int chunk = (n + workers - 1) / workers;
  for (int w = 0; w < workers; ++w) {
    const int begin = w * chunk;
    const int end = std::min(n, begin + chunk);
    if (begin >= end) break;
    threads.emplace_back([&fn, begin, end, w]() { fn(begin, end, w); });
  }
  for (auto& t : threads) t.join();
}

// ---- per-target-type working state ------------------------------------------

// Sparse block of the count matrix for one target type: rows are the type's
// cells (in column order), columns are genes; gpos maps each stored value
// back to its position in the input matrix's value array.
struct TypeBlock {
  std::vector<int> cols;        // global column index per row
  std::vector<int> rowptr;
  std::vector<int> gene;        // gene index per stored value
  std::vector<double> value;
  std::vector<std::int64_t> gpos;
  std::vector<double> totals;   // per-row molecule totals
};

struct PairState {
  int pair_index = 0;
  const GenerativePairSpec* spec = nullptr;
  std::vector<int> rows;        // row index within the block per target cell
  std::vector<double> exposure;
  std::vector<double> psi;      // G, decontaminated, zeroed on target genes
  std::vector<int> guide;
  double share_guide = 0.0;
  std::vector<double> lam_cell; // dose prior per row of the block
};

}  // namespace

// ---- exported helpers --------------------------------------------------------

std::vector<double> pava_increasing(
    const std::vector<double>& y, const std::vector<double>& w) {
  const std::size_t n = y.size();
  std::vector<double> val;
  std::vector<double> wt;
  std::vector<std::vector<std::size_t>> idx;
  for (std::size_t i = 0; i < n; ++i) {
    val.push_back(y[i]);
    wt.push_back(w[i]);
    idx.push_back({i});
    while (val.size() > 1 && val[val.size() - 2] > val.back() + 1e-15) {
      const std::size_t a = val.size() - 2;
      const std::size_t b = val.size() - 1;
      const double merged_w = wt[a] + wt[b];
      const double merged_v = (val[a] * wt[a] + val[b] * wt[b]) / merged_w;
      val[a] = merged_v;
      wt[a] = merged_w;
      idx[a].insert(idx[a].end(), idx[b].begin(), idx[b].end());
      val.pop_back();
      wt.pop_back();
      idx.pop_back();
    }
  }
  std::vector<double> out(n, 0.0);
  for (std::size_t j = 0; j < val.size(); ++j) {
    for (const std::size_t i : idx[j]) out[i] = val[j];
  }
  return out;
}

// ---- the fit -----------------------------------------------------------------

GenerativeResult fit_generative(
    const std::vector<int>& counts_indptr,
    const std::vector<int>& counts_indices,
    const std::vector<double>& counts_values,
    int n_genes,
    const std::vector<std::string>& cell_ids,
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<int>& type_codes,
    int n_types,
    const std::string& molecules_parquet,
    const std::string& cells_parquet,
    const std::vector<GenerativePairSpec>& pairs,
    const std::vector<int>& factor_to_type,
    const std::vector<double>& programs_flat,
    const std::vector<int>& program_type,
    const GenerativeOptions& opt) {
  const int n_cells = static_cast<int>(cell_ids.size());
  const int G = n_genes;
  if (static_cast<int>(counts_indptr.size()) != n_cells + 1) {
    throw std::runtime_error("generative: counts indptr does not match cells");
  }
  if (x.size() != cell_ids.size() || y.size() != cell_ids.size() ||
      type_codes.size() != cell_ids.size()) {
    throw std::runtime_error("generative: cell metadata does not match cells");
  }
  GenerativeResult result;
  result.removed.assign(counts_values.size(), 0.0);
  result.removed_without_retention.assign(counts_values.size(), 0.0);
  result.n_programs_used.assign(static_cast<std::size_t>(n_types), 0);
  result.ambient_scale.assign(static_cast<std::size_t>(n_cells), 0.0);
  result.pair_cells.resize(pairs.size());
  result.pair_dose.resize(pairs.size());
  result.pair_alpha.resize(pairs.size());
  result.pair_rho.resize(pairs.size());

  // Per-cell totals and per-type pseudobulk profiles (counts per million).
  std::vector<double> totals(static_cast<std::size_t>(n_cells), 0.0);
  for (int c = 0; c < n_cells; ++c) {
    for (int p = counts_indptr[c]; p < counts_indptr[c + 1]; ++p) {
      totals[static_cast<std::size_t>(c)] += counts_values[static_cast<std::size_t>(p)];
    }
  }
  DenseMatrix profiles_cpm(G, n_types, 0.0);
  {
    std::vector<double> type_totals(static_cast<std::size_t>(n_types), 0.0);
    for (int c = 0; c < n_cells; ++c) {
      const int t = type_codes[static_cast<std::size_t>(c)];
      if (t < 0) continue;
      type_totals[static_cast<std::size_t>(t)] += totals[static_cast<std::size_t>(c)];
      for (int p = counts_indptr[c]; p < counts_indptr[c + 1]; ++p) {
        profiles_cpm(counts_indices[static_cast<std::size_t>(p)], t) +=
            counts_values[static_cast<std::size_t>(p)];
      }
    }
    for (int t = 0; t < n_types; ++t) {
      const double denom = std::max(type_totals[static_cast<std::size_t>(t)], 1.0);
      for (int g = 0; g < G; ++g) profiles_cpm(g, t) *= 1e6 / denom;
    }
  }
  std::vector<int> top_type(static_cast<std::size_t>(G), -1);
  for (int g = 0; g < G; ++g) {
    double best = -1.0;
    for (int t = 0; t < n_types; ++t) {
      if (profiles_cpm(g, t) > best) {
        best = profiles_cpm(g, t);
        top_type[static_cast<std::size_t>(g)] = t;
      }
    }
  }

  // Spatial covariates: exposure counts and nearest-type distances.
  const DenseMatrix exposure_counts =
      cell_neighbor_type_counts(x, y, type_codes, n_types, opt.neighbor_k);
  const DenseMatrix nearest_dist = cell_nearest_type_distance(x, y, type_codes, n_types);

  // Cells of each type in column order.
  std::vector<std::vector<int>> type_cols(static_cast<std::size_t>(n_types));
  for (int c = 0; c < n_cells; ++c) {
    const int t = type_codes[static_cast<std::size_t>(c)];
    if (t >= 0) type_cols[static_cast<std::size_t>(t)].push_back(c);
  }

  // ---- molecule-derived preparation ----------------------------------------
  // Cytoplasmic gene x cell counts (for interface-local transfer profiles) and
  // per-type expression-program initializations from factor-labeled molecules.
  // Without a molecule table, whole-cell counts and pseudobulk programs serve.
  // Cytoplasmic counts are only ever summed over column subsets, so they are
  // stored per column as sparse (gene, count) pairs.
  std::vector<std::vector<std::pair<int, float>>> cyto(static_cast<std::size_t>(n_cells));
  const bool have_molecules = !molecules_parquet.empty();
  const bool explicit_programs = !programs_flat.empty();
  if (explicit_programs) {
    if (programs_flat.size() != program_type.size() * static_cast<std::size_t>(G)) {
      throw std::runtime_error("generative: programs_flat does not match program_type x genes");
    }
  }
  // Program source: explicit programs override the mode; "factors" needs a
  // molecule table and otherwise falls back to pseudobulk.
  const bool factor_programs = !explicit_programs && have_molecules &&
      opt.init_mode == "factors";
  const bool cluster_programs = !explicit_programs && opt.init_mode == "clusters";
  std::vector<DenseMatrix> program_init(static_cast<std::size_t>(n_types));
  result.whole_cell_profiles = true;
  std::vector<int> alignment = factor_to_type;
  if (have_molecules) {
    const auto cells_table = read_parquet(cells_parquet);
    const auto run_cell_ids = extract_string_column(*cells_table, "cell_id");
    std::vector<int> run_cell_idx;
    if (table_has_column(*cells_table, "cell_idx")) {
      run_cell_idx = extract_int_column(*cells_table, "cell_idx");
    } else {
      run_cell_idx.resize(run_cell_ids.size());
      for (std::size_t i = 0; i < run_cell_ids.size(); ++i) {
        run_cell_idx[i] = static_cast<int>(i);
      }
    }
    std::unordered_map<std::string, int> col_of_id;
    col_of_id.reserve(cell_ids.size() * 2);
    for (int c = 0; c < n_cells; ++c) col_of_id[cell_ids[static_cast<std::size_t>(c)]] = c;
    int max_run_idx = 0;
    for (const int i : run_cell_idx) max_run_idx = std::max(max_run_idx, i);
    std::vector<int> col_of_run(static_cast<std::size_t>(max_run_idx + 1), -1);
    for (std::size_t i = 0; i < run_cell_ids.size(); ++i) {
      const auto it = col_of_id.find(run_cell_ids[i]);
      if (it != col_of_id.end()) col_of_run[static_cast<std::size_t>(run_cell_idx[i])] = it->second;
    }

    const auto mol = read_parquet(molecules_parquet);
    const auto gene_idx = extract_int_column(*mol, "gene_idx");
    const auto cell_idx = extract_int_column(*mol, "cell_idx");
    const auto labels = extract_int_column(*mol, "factor_label");
    const bool has_nucleus = table_has_column(*mol, "overlaps_nucleus");
    std::vector<int> nucleus;
    if (has_nucleus) nucleus = extract_int_column(*mol, "overlaps_nucleus");
    result.whole_cell_profiles = !has_nucleus;

    int n_factors = 0;
    for (const int k : labels) n_factors = std::max(n_factors, k + 1);

    // One pass accumulates everything needed: cytoplasmic per-column counts,
    // per-factor global gene profiles (for alignment), per-(type, factor)
    // gene profiles (for program initializations), and per-type factor shares.
    std::vector<std::vector<double>> factor_gene(
        static_cast<std::size_t>(n_factors), std::vector<double>(static_cast<std::size_t>(G), 0.0));
    std::vector<DenseMatrix> type_factor_gene(static_cast<std::size_t>(n_types));
    for (int t = 0; t < n_types; ++t) type_factor_gene[static_cast<std::size_t>(t)] = DenseMatrix(n_factors, G, 0.0);
    DenseMatrix type_factor_share(n_types, n_factors, 0.0);
    {
      // temporary dense per-column gene tallies would be too large; build the
      // cytoplasmic counts by sorting (column, gene) keys.
      std::vector<std::int64_t> keys;
      keys.reserve(gene_idx.size());
      for (std::size_t i = 0; i < gene_idx.size(); ++i) {
        const int ci = cell_idx[i];
        if (ci < 0 || ci > max_run_idx) continue;
        const int col = col_of_run[static_cast<std::size_t>(ci)];
        if (col < 0) continue;
        const int g = gene_idx[i];
        const int k = labels[i];
        const int t = type_codes[static_cast<std::size_t>(col)];
        if (k >= 0) {
          factor_gene[static_cast<std::size_t>(k)][static_cast<std::size_t>(g)] += 1.0;
          if (t >= 0) {
            type_factor_gene[static_cast<std::size_t>(t)](k, g) += 1.0;
            type_factor_share(t, k) += 1.0;
          }
        }
        if (!has_nucleus || nucleus[i] == 0) {
          keys.push_back(static_cast<std::int64_t>(col) * G + g);
        }
      }
      std::sort(keys.begin(), keys.end());
      for (std::size_t i = 0; i < keys.size();) {
        std::size_t j = i;
        while (j < keys.size() && keys[j] == keys[i]) ++j;
        const int col = static_cast<int>(keys[i] / G);
        const int g = static_cast<int>(keys[i] % G);
        cyto[static_cast<std::size_t>(col)].push_back({g, static_cast<float>(j - i)});
        i = j;
      }
    }

    // Factor-to-type alignment: supplied, or derived by cosine similarity of
    // square-root profiles with a required margin over the second-best type.
    if (alignment.empty() && factor_programs) {
      alignment.assign(static_cast<std::size_t>(n_factors), -1);
      for (int k = 0; k < n_factors; ++k) {
        double total = 0.0;
        for (int g = 0; g < G; ++g) total += factor_gene[static_cast<std::size_t>(k)][static_cast<std::size_t>(g)];
        if (total < opt.program_min_molecules) continue;
        std::vector<double> fk(static_cast<std::size_t>(G));
        double fk_norm = 0.0;
        for (int g = 0; g < G; ++g) {
          fk[static_cast<std::size_t>(g)] = std::sqrt(
              factor_gene[static_cast<std::size_t>(k)][static_cast<std::size_t>(g)] / total);
          fk_norm += fk[static_cast<std::size_t>(g)] * fk[static_cast<std::size_t>(g)];
        }
        fk_norm = std::sqrt(std::max(fk_norm, 1e-12));
        double best = -1.0;
        double second = -1.0;
        int best_t = -1;
        for (int t = 0; t < n_types; ++t) {
          double dot = 0.0;
          double tn = 0.0;
          double tsum = 0.0;
          for (int g = 0; g < G; ++g) tsum += profiles_cpm(g, t);
          for (int g = 0; g < G; ++g) {
            const double tv = std::sqrt(profiles_cpm(g, t) / std::max(tsum, 1.0));
            dot += fk[static_cast<std::size_t>(g)] * tv;
            tn += tv * tv;
          }
          const double cosv = dot / (fk_norm * std::sqrt(std::max(tn, 1e-12)));
          if (cosv > best) {
            second = best;
            best = cosv;
            best_t = t;
          } else if (cosv > second) {
            second = cosv;
          }
        }
        if (best >= opt.align_min_cosine && best >= opt.align_margin * std::max(second, 1e-9)) {
          alignment[static_cast<std::size_t>(k)] = best_t;
        }
      }
    }
    result.factor_alignment = alignment;

    // Program initializations: the type's aligned factors plus unaligned
    // factors carrying at least program_share_min of the type's molecules;
    // pseudobulk fallback when none qualify.
    for (int t = 0; factor_programs && t < n_types; ++t) {
      double type_total = 0.0;
      for (int k = 0; k < n_factors; ++k) type_total += type_factor_share(t, k);
      std::vector<int> own;
      for (int k = 0; k < n_factors; ++k) {
        const bool aligned_here =
            k < static_cast<int>(alignment.size()) && alignment[static_cast<std::size_t>(k)] == t;
        const bool unaligned = k >= static_cast<int>(alignment.size()) ||
            alignment[static_cast<std::size_t>(k)] == -1;
        const double share = type_factor_share(t, k) / std::max(type_total, 1.0);
        if (aligned_here || (unaligned && share >= opt.program_share_min)) own.push_back(k);
      }
      std::vector<std::vector<double>> progs;
      for (const int k : own) {
        double s = 0.0;
        for (int g = 0; g < G; ++g) s += type_factor_gene[static_cast<std::size_t>(t)](k, g);
        if (s < opt.program_min_molecules) continue;
        std::vector<double> prog(static_cast<std::size_t>(G));
        for (int g = 0; g < G; ++g) prog[static_cast<std::size_t>(g)] = type_factor_gene[static_cast<std::size_t>(t)](k, g) / s;
        progs.push_back(std::move(prog));
      }
      if (progs.empty()) {
        std::vector<double> prog(static_cast<std::size_t>(G));
        double s = 0.0;
        for (int g = 0; g < G; ++g) s += profiles_cpm(g, t);
        for (int g = 0; g < G; ++g) prog[static_cast<std::size_t>(g)] = profiles_cpm(g, t) / std::max(s, 1.0);
        progs.push_back(std::move(prog));
      }
      DenseMatrix F(static_cast<int>(progs.size()), G, 0.0);
      for (std::size_t k = 0; k < progs.size(); ++k) {
        for (int g = 0; g < G; ++g) F(static_cast<int>(k), g) = progs[k][static_cast<std::size_t>(g)];
      }
      program_init[static_cast<std::size_t>(t)] = std::move(F);
    }
  } else {
    // No molecule table: whole-cell counts serve as the transfer profiles.
    for (int c = 0; c < n_cells; ++c) {
      for (int p = counts_indptr[c]; p < counts_indptr[c + 1]; ++p) {
        cyto[static_cast<std::size_t>(c)].push_back(
            {counts_indices[static_cast<std::size_t>(p)],
             static_cast<float>(counts_values[static_cast<std::size_t>(p)])});
      }
    }
  }
  if (explicit_programs) {
    for (int t = 0; t < n_types; ++t) {
      std::vector<int> rows;
      for (std::size_t r = 0; r < program_type.size(); ++r) {
        if (program_type[r] == t) rows.push_back(static_cast<int>(r));
      }
      if (rows.empty()) continue;
      DenseMatrix F(static_cast<int>(rows.size()), G, 0.0);
      for (std::size_t k = 0; k < rows.size(); ++k) {
        double sum = 0.0;
        for (int g = 0; g < G; ++g) {
          const double v = std::max(
              programs_flat[static_cast<std::size_t>(rows[k]) * G + g], 0.0);
          F(static_cast<int>(k), g) = v;
          sum += v;
        }
        for (int g = 0; g < G; ++g) F(static_cast<int>(k), g) /= std::max(sum, 1e-12);
      }
      program_init[static_cast<std::size_t>(t)] = std::move(F);
    }
  }
  // Any type still without programs (pseudobulk mode, "factors" without a
  // molecule table, explicit programs omitting a type, or a type whose
  // factor programs were all too small): one pooled profile. The "clusters"
  // mode replaces this per type inside the fit, where the dose weights are
  // available.
  for (int t = 0; t < n_types; ++t) {
    if (program_init[static_cast<std::size_t>(t)].rows() > 0) continue;
    DenseMatrix F(1, G, 0.0);
    double sum = 0.0;
    for (int g = 0; g < G; ++g) sum += profiles_cpm(g, t);
    for (int g = 0; g < G; ++g) F(0, g) = profiles_cpm(g, t) / std::max(sum, 1.0);
    program_init[static_cast<std::size_t>(t)] = std::move(F);
  }

  // Interface-local raw transfer profile of source S toward target type T:
  // cytoplasmic profile of S cells within near_um of the nearest T cell,
  // widening to 100 um and then all S cells when too few molecules.
  const auto near_profile = [&](int S, int T) {
    std::vector<double> prof(static_cast<std::size_t>(G), 0.0);
    const double radii[3] = {opt.near_um, 100.0, std::numeric_limits<double>::infinity()};
    for (const double r : radii) {
      std::fill(prof.begin(), prof.end(), 0.0);
      double total = 0.0;
      for (const int c : type_cols[static_cast<std::size_t>(S)]) {
        const double d = nearest_dist(c, T);
        if (!(d <= r)) continue;
        for (const auto& [g, v] : cyto[static_cast<std::size_t>(c)]) {
          prof[static_cast<std::size_t>(g)] += v;
          total += v;
        }
      }
      if (total >= opt.near_min_molecules || std::isinf(r)) {
        for (auto& v : prof) v /= std::max(total, 1.0);
        return prof;
      }
    }
    return prof;
  };
  // Cache raw profiles per pair (decontamination is applied per round).
  std::vector<std::vector<double>> psi_raw(pairs.size());
  for (std::size_t j = 0; j < pairs.size(); ++j) {
    psi_raw[j] = near_profile(pairs[j].source_type, pairs[j].target_type);
  }

  // Per-gene own-expression fraction of each type, from the previous outer
  // round, used to decontaminate the source profiles.
  std::vector<std::vector<double>> own_frac(static_cast<std::size_t>(n_types));

  for (int round = 0; round < opt.outer_rounds; ++round) {
    const bool last_round = round == opt.outer_rounds - 1;
    std::vector<std::vector<double>> own_frac_new(static_cast<std::size_t>(n_types));

    for (int T = 0; T < n_types; ++T) {
      // Pairs into this target type.
      std::vector<PairState> plist;
      for (std::size_t j = 0; j < pairs.size(); ++j) {
        if (pairs[j].target_type == T) {
          PairState ps;
          ps.pair_index = static_cast<int>(j);
          ps.spec = &pairs[j];
          plist.push_back(std::move(ps));
        }
      }
      if (plist.empty() || type_cols[static_cast<std::size_t>(T)].empty()) continue;

      // Sparse block for this type's cells.
      TypeBlock blk;
      blk.cols = type_cols[static_cast<std::size_t>(T)];
      const int n = static_cast<int>(blk.cols.size());
      blk.rowptr.assign(static_cast<std::size_t>(n + 1), 0);
      std::size_t nnz = 0;
      for (int i = 0; i < n; ++i) {
        const int c = blk.cols[static_cast<std::size_t>(i)];
        nnz += static_cast<std::size_t>(counts_indptr[c + 1] - counts_indptr[c]);
        blk.rowptr[static_cast<std::size_t>(i + 1)] = static_cast<int>(nnz);
      }
      blk.gene.resize(nnz);
      blk.value.resize(nnz);
      blk.gpos.resize(nnz);
      blk.totals.assign(static_cast<std::size_t>(n), 0.0);
      for (int i = 0; i < n; ++i) {
        const int c = blk.cols[static_cast<std::size_t>(i)];
        int q = blk.rowptr[static_cast<std::size_t>(i)];
        for (int p = counts_indptr[c]; p < counts_indptr[c + 1]; ++p, ++q) {
          blk.gene[static_cast<std::size_t>(q)] = counts_indices[static_cast<std::size_t>(p)];
          blk.value[static_cast<std::size_t>(q)] = counts_values[static_cast<std::size_t>(p)];
          blk.gpos[static_cast<std::size_t>(q)] = p;
          blk.totals[static_cast<std::size_t>(i)] += counts_values[static_cast<std::size_t>(p)];
        }
      }
      std::vector<int> row_of_col(static_cast<std::size_t>(n_cells), -1);
      for (int i = 0; i < n; ++i) row_of_col[static_cast<std::size_t>(blk.cols[static_cast<std::size_t>(i)])] = i;

      const int S_n = static_cast<int>(plist.size());
      if (S_n > 16 || program_init[static_cast<std::size_t>(T)].rows() > 16) {
        throw std::runtime_error("generative: more than 16 sources or programs per type");
      }

      // Decontaminated transfer profiles, zeroed on target-owned genes.
      for (auto& ps : plist) {
        ps.psi = psi_raw[static_cast<std::size_t>(ps.pair_index)];
        const auto& of = own_frac[static_cast<std::size_t>(ps.spec->source_type)];
        if (!of.empty()) {
          double s = 0.0;
          for (int g = 0; g < G; ++g) {
            ps.psi[static_cast<std::size_t>(g)] *= of[static_cast<std::size_t>(g)];
            s += ps.psi[static_cast<std::size_t>(g)];
          }
          for (auto& v : ps.psi) v /= std::max(s, 1e-12);
        }
        for (int g = 0; g < G; ++g) {
          if (top_type[static_cast<std::size_t>(g)] == T) ps.psi[static_cast<std::size_t>(g)] = 0.0;
        }
      }

      // Exposure and target-cell rows per pair (all target cells of the type).
      for (auto& ps : plist) {
        ps.rows.resize(static_cast<std::size_t>(n));
        ps.exposure.resize(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
          ps.rows[static_cast<std::size_t>(i)] = i;
          ps.exposure[static_cast<std::size_t>(i)] = ps.spec->exposure.empty()
              ? exposure_counts(blk.cols[static_cast<std::size_t>(i)], ps.spec->source_type)
              : ps.spec->exposure[static_cast<std::size_t>(i)];
        }
        ps.guide = ps.spec->guide.empty() ? ps.spec->pool : ps.spec->guide;
        ps.share_guide = 0.0;
        for (const int g : ps.guide) ps.share_guide += ps.psi[static_cast<std::size_t>(g)];
      }

      // Isotonic dose-response of guide-gene content on exposure, converted
      // to a whole-profile contamination fraction per exposure stratum.
      // values: current per-entry kept counts (raw counts on the first call).
      const auto dose_response = [&](const PairState& ps, const std::vector<double>& kept,
                                     const std::vector<int>& guide) {
        std::vector<char> in_guide(static_cast<std::size_t>(G), 0);
        for (const int g : guide) in_guide[static_cast<std::size_t>(g)] = 1;
        std::unordered_map<double, std::pair<double, double>> strata;  // e -> (m, t)
        std::vector<double> row_guide(static_cast<std::size_t>(n), 0.0);
        parallel_rows(n, opt.num_threads, [&](int begin, int end, int) {
          for (int i = begin; i < end; ++i) {
            double s = 0.0;
            for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
              if (in_guide[static_cast<std::size_t>(blk.gene[static_cast<std::size_t>(p)])]) {
                s += kept[static_cast<std::size_t>(p)];
              }
            }
            row_guide[static_cast<std::size_t>(i)] = s;
          }
        });
        for (int i = 0; i < n; ++i) {
          auto& slot = strata[ps.exposure[static_cast<std::size_t>(i)]];
          slot.first += row_guide[static_cast<std::size_t>(i)];
          slot.second += blk.totals[static_cast<std::size_t>(i)];
        }
        std::vector<double> levels;
        for (const auto& [e, mt] : strata) levels.push_back(e);
        std::sort(levels.begin(), levels.end());
        if (levels.empty() || levels.front() != 0.0) {
          return std::unordered_map<double, double>{};
        }
        std::vector<double> rate;
        std::vector<double> weight;
        for (const double e : levels) {
          const auto& mt = strata[e];
          rate.push_back(mt.first / std::max(mt.second, 1.0));
          weight.push_back(mt.second);
        }
        const auto fitted = pava_increasing(rate, weight);
        std::unordered_map<double, double> lam;
        const double share = std::max(ps.share_guide, 1e-6);
        for (std::size_t s = 0; s < levels.size(); ++s) {
          lam[levels[s]] = std::min(std::max(fitted[s] - fitted[0], 0.0) / share, opt.lambda_max);
        }
        return lam;
      };

      for (auto& ps : plist) {
        const auto lam = dose_response(ps, blk.value, ps.guide);
        ps.lam_cell.assign(static_cast<std::size_t>(n), 0.0);
        for (int i = 0; i < n; ++i) {
          const auto it = lam.find(ps.exposure[static_cast<std::size_t>(i)]);
          if (it != lam.end()) ps.lam_cell[static_cast<std::size_t>(i)] = it->second;
        }
      }

      // Ambient tier: union of strict genes; level and profile from target
      // cells far from every source of this type's pairs, falling back to
      // half the pooled zero-exposure rate.
      std::vector<char> strictU(static_cast<std::size_t>(G), 0);
      for (const auto& ps : plist) {
        for (const int g : ps.spec->strict) strictU[static_cast<std::size_t>(g)] = 1;
      }
      bool any_strict = std::any_of(strictU.begin(), strictU.end(), [](char v) { return v != 0; });
      std::vector<double> a_prof(static_cast<std::size_t>(G), 0.0);
      double lam_a = 0.0;
      const bool use_ambient = opt.use_ambient && any_strict;
      if (use_ambient) {
        std::vector<char> far(static_cast<std::size_t>(n), 1);
        for (int i = 0; i < n; ++i) {
          const int c = blk.cols[static_cast<std::size_t>(i)];
          for (const auto& ps : plist) {
            if (nearest_dist(c, ps.spec->source_type) <= opt.far_um) {
              far[static_cast<std::size_t>(i)] = 0;
              break;
            }
          }
        }
        double far_mol = 0.0;
        for (int i = 0; i < n; ++i) {
          if (far[static_cast<std::size_t>(i)]) far_mol += blk.totals[static_cast<std::size_t>(i)];
        }
        std::vector<double> rates(static_cast<std::size_t>(G), 0.0);
        if (far_mol >= opt.far_min_molecules) {
          for (int i = 0; i < n; ++i) {
            if (!far[static_cast<std::size_t>(i)]) continue;
            for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
              rates[static_cast<std::size_t>(blk.gene[static_cast<std::size_t>(p)])] += blk.value[static_cast<std::size_t>(p)];
            }
          }
          for (auto& v : rates) v /= far_mol;
        } else {
          std::vector<char> e0(static_cast<std::size_t>(n), 0);
          for (const auto& ps : plist) {
            for (int i = 0; i < n; ++i) {
              if (ps.exposure[static_cast<std::size_t>(i)] == 0.0) e0[static_cast<std::size_t>(i)] = 1;
            }
          }
          double t0 = 0.0;
          for (int i = 0; i < n; ++i) {
            if (!e0[static_cast<std::size_t>(i)]) continue;
            t0 += blk.totals[static_cast<std::size_t>(i)];
            for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
              rates[static_cast<std::size_t>(blk.gene[static_cast<std::size_t>(p)])] += blk.value[static_cast<std::size_t>(p)];
            }
          }
          for (auto& v : rates) v = 0.5 * v / std::max(t0, 1.0);
        }
        for (int g = 0; g < G; ++g) {
          if (strictU[static_cast<std::size_t>(g)]) lam_a += rates[static_cast<std::size_t>(g)];
        }
        if (lam_a > 0.0) {
          for (int g = 0; g < G; ++g) {
            a_prof[static_cast<std::size_t>(g)] =
                strictU[static_cast<std::size_t>(g)] ? rates[static_cast<std::size_t>(g)] / lam_a : 0.0;
          }
        }
      }

      // Own programs; the production configuration zeroes them on strict
      // genes (a gene the target does not express cannot be a program's).
      // The "clusters" mode derives the programs here, where the dose
      // priors are available: the type's cells are clustered on their gene
      // fractions by weighted k-means, with weights favoring lightly dosed
      // cells so that a contamination pattern cannot seed a program.
      DenseMatrix F = program_init[static_cast<std::size_t>(T)];
      if (cluster_programs && n >= 20) {
        std::vector<double> lam_tot(static_cast<std::size_t>(n), 0.0);
        for (const auto& ps : plist) {
          for (int i = 0; i < n; ++i) {
            lam_tot[static_cast<std::size_t>(i)] += ps.lam_cell[static_cast<std::size_t>(i)];
          }
        }
        std::vector<double> wgt(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
          wgt[static_cast<std::size_t>(i)] = blk.totals[static_cast<std::size_t>(i)] /
              (1.0 + opt.dose_weight * lam_tot[static_cast<std::size_t>(i)]);
        }
        const int Kc = std::min(opt.n_programs, n);
        // squared norms of the per-cell gene-fraction vectors
        std::vector<double> xnorm(static_cast<std::size_t>(n), 0.0);
        for (int i = 0; i < n; ++i) {
          const double ti = std::max(blk.totals[static_cast<std::size_t>(i)], 1.0);
          for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
            const double f = blk.value[static_cast<std::size_t>(p)] / ti;
            xnorm[static_cast<std::size_t>(i)] += f * f;
          }
        }
        // weighted k-means++ seeding, deterministic per type
        std::mt19937 rng(opt.seed + static_cast<unsigned int>(T) * 7919u);
        std::vector<std::vector<double>> cent;
        {
          std::discrete_distribution<int> first(wgt.begin(), wgt.end());
          std::vector<double> d2(static_cast<std::size_t>(n),
                                 std::numeric_limits<double>::max());
          int pick = first(rng);
          while (static_cast<int>(cent.size()) < Kc) {
            std::vector<double> c(static_cast<std::size_t>(G), 0.0);
            const double ti = std::max(blk.totals[static_cast<std::size_t>(pick)], 1.0);
            for (int p = blk.rowptr[static_cast<std::size_t>(pick)]; p < blk.rowptr[static_cast<std::size_t>(pick + 1)]; ++p) {
              c[static_cast<std::size_t>(blk.gene[static_cast<std::size_t>(p)])] =
                  blk.value[static_cast<std::size_t>(p)] / ti;
            }
            cent.push_back(std::move(c));
            if (static_cast<int>(cent.size()) == Kc) break;
            const auto& cc = cent.back();
            double cn = 0.0;
            for (const double v : cc) cn += v * v;
            std::vector<double> probs(static_cast<std::size_t>(n), 0.0);
            for (int i = 0; i < n; ++i) {
              double dot = 0.0;
              const double ti2 = std::max(blk.totals[static_cast<std::size_t>(i)], 1.0);
              for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
                dot += blk.value[static_cast<std::size_t>(p)] / ti2 *
                    cc[static_cast<std::size_t>(blk.gene[static_cast<std::size_t>(p)])];
              }
              const double d = std::max(xnorm[static_cast<std::size_t>(i)] - 2.0 * dot + cn, 0.0);
              d2[static_cast<std::size_t>(i)] = std::min(d2[static_cast<std::size_t>(i)], d);
              probs[static_cast<std::size_t>(i)] = d2[static_cast<std::size_t>(i)] * wgt[static_cast<std::size_t>(i)];
            }
            std::discrete_distribution<int> next(probs.begin(), probs.end());
            pick = next(rng);
          }
        }
        // Lloyd iterations with weights
        std::vector<int> assign(static_cast<std::size_t>(n), 0);
        for (int it = 0; it < 10; ++it) {
          std::vector<double> cn(cent.size(), 0.0);
          for (std::size_t k = 0; k < cent.size(); ++k) {
            for (const double v : cent[k]) cn[k] += v * v;
          }
          parallel_rows(n, opt.num_threads, [&](int begin, int end, int) {
            for (int i = begin; i < end; ++i) {
              const double ti = std::max(blk.totals[static_cast<std::size_t>(i)], 1.0);
              double best = std::numeric_limits<double>::max();
              int bk = 0;
              for (std::size_t k = 0; k < cent.size(); ++k) {
                double dot = 0.0;
                for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
                  dot += blk.value[static_cast<std::size_t>(p)] / ti *
                      cent[k][static_cast<std::size_t>(blk.gene[static_cast<std::size_t>(p)])];
                }
                const double d = xnorm[static_cast<std::size_t>(i)] - 2.0 * dot + cn[k];
                if (d < best) { best = d; bk = static_cast<int>(k); }
              }
              assign[static_cast<std::size_t>(i)] = bk;
            }
          });
          std::vector<std::vector<double>> nc(cent.size(),
              std::vector<double>(static_cast<std::size_t>(G), 0.0));
          std::vector<double> nw(cent.size(), 0.0);
          for (int i = 0; i < n; ++i) {
            const int k = assign[static_cast<std::size_t>(i)];
            const double w = wgt[static_cast<std::size_t>(i)];
            const double ti = std::max(blk.totals[static_cast<std::size_t>(i)], 1.0);
            nw[static_cast<std::size_t>(k)] += w;
            for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
              nc[static_cast<std::size_t>(k)][static_cast<std::size_t>(blk.gene[static_cast<std::size_t>(p)])] +=
                  w * blk.value[static_cast<std::size_t>(p)] / ti;
            }
          }
          for (std::size_t k = 0; k < cent.size(); ++k) {
            if (nw[k] <= 0.0) continue;
            for (int g = 0; g < G; ++g) nc[k][static_cast<std::size_t>(g)] /= nw[k];
            cent[k] = std::move(nc[k]);
          }
        }
        // programs: dose-weighted count sums per cluster, keeping clusters
        // with enough molecules and share
        std::vector<std::vector<double>> progs;
        double type_mol = 0.0;
        for (int i = 0; i < n; ++i) type_mol += blk.totals[static_cast<std::size_t>(i)];
        for (std::size_t k = 0; k < cent.size(); ++k) {
          std::vector<double> prog(static_cast<std::size_t>(G), 0.0);
          double mol = 0.0;
          for (int i = 0; i < n; ++i) {
            if (assign[static_cast<std::size_t>(i)] != static_cast<int>(k)) continue;
            mol += blk.totals[static_cast<std::size_t>(i)];
            const double w = 1.0 / (1.0 + opt.dose_weight * lam_tot[static_cast<std::size_t>(i)]);
            for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
              prog[static_cast<std::size_t>(blk.gene[static_cast<std::size_t>(p)])] +=
                  w * blk.value[static_cast<std::size_t>(p)];
            }
          }
          if (mol < opt.program_min_molecules || mol < opt.program_share_min * type_mol) {
            continue;
          }
          double sum = 0.0;
          for (const double v : prog) sum += v;
          for (auto& v : prog) v /= std::max(sum, 1e-12);
          progs.push_back(std::move(prog));
        }
        if (!progs.empty()) {
          F = DenseMatrix(static_cast<int>(progs.size()), G, 0.0);
          for (std::size_t k = 0; k < progs.size(); ++k) {
            for (int g = 0; g < G; ++g) F(static_cast<int>(k), g) = progs[k][static_cast<std::size_t>(g)];
          }
        }
      }
      result.n_programs_used[static_cast<std::size_t>(T)] = F.rows();
      const int K = F.rows();
      if (opt.use_ambient && any_strict) {
        for (int k = 0; k < K; ++k) {
          double s = 0.0;
          for (int g = 0; g < G; ++g) {
            if (strictU[static_cast<std::size_t>(g)]) F(k, g) = 0.0;
            s += F(k, g);
          }
          for (int g = 0; g < G; ++g) F(k, g) /= std::max(s, 1e-12);
        }
      }

      // Induced-expression support: per pair, genes owned by the source whose
      // exposure-linked excess is disproportionate to the transfer profile.
      DenseMatrix M(S_n, G, 0.0);
      std::vector<GenerativeInducedGene> induced_rows;
      if (opt.use_induced) {
        for (int j = 0; j < S_n; ++j) {
          const auto& ps = plist[static_cast<std::size_t>(j)];
          std::vector<double> c_exp(static_cast<std::size_t>(G), 0.0);
          std::vector<double> c_un(static_cast<std::size_t>(G), 0.0);
          double t_exp = 0.0;
          double t_un = 0.0;
          int n_exp = 0;
          int n_un = 0;
          for (int i = 0; i < n; ++i) {
            const bool exposed = ps.exposure[static_cast<std::size_t>(i)] > 0.0;
            auto& acc = exposed ? c_exp : c_un;
            (exposed ? t_exp : t_un) += blk.totals[static_cast<std::size_t>(i)];
            (exposed ? n_exp : n_un) += 1;
            for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
              acc[static_cast<std::size_t>(blk.gene[static_cast<std::size_t>(p)])] += blk.value[static_cast<std::size_t>(p)];
            }
          }
          if (n_exp < 50 || n_un < 50) continue;
          double dose_tot = 0.0;
          for (int i = 0; i < n; ++i) {
            dose_tot += blk.totals[static_cast<std::size_t>(i)] * ps.lam_cell[static_cast<std::size_t>(i)];
          }
          dose_tot = std::max(dose_tot, 1.0);
          // weighted proportional fit over source-owned genes
          double num = 0.0;
          double den = 0.0;
          std::vector<int> gset;
          std::vector<double> excess(static_cast<std::size_t>(G), 0.0);
          std::vector<double> var(static_cast<std::size_t>(G), 0.0);
          for (int g = 0; g < G; ++g) {
            if (top_type[static_cast<std::size_t>(g)] != ps.spec->source_type) continue;
            gset.push_back(g);
            const double ex = c_exp[static_cast<std::size_t>(g)] -
                c_un[static_cast<std::size_t>(g)] / std::max(t_un, 1.0) * t_exp;
            const double vr = c_exp[static_cast<std::size_t>(g)] +
                (t_exp / std::max(t_un, 1.0)) * (t_exp / std::max(t_un, 1.0)) *
                    c_un[static_cast<std::size_t>(g)] + 1.0;
            excess[static_cast<std::size_t>(g)] = ex;
            var[static_cast<std::size_t>(g)] = vr;
            const double p = ps.psi[static_cast<std::size_t>(g)];
            num += ex * p / vr;
            den += p * p / vr;
          }
          const double slope = std::max(num / std::max(den, 1e-12), 0.0);
          for (const int g : gset) {
            const double expected = slope * ps.psi[static_cast<std::size_t>(g)];
            const double z = (excess[static_cast<std::size_t>(g)] - expected) /
                std::sqrt(var[static_cast<std::size_t>(g)] +
                          (opt.profile_cv * expected) * (opt.profile_cv * expected));
            if (z > opt.induced_z && excess[static_cast<std::size_t>(g)] > opt.induced_min_excess) {
              M(j, g) = std::max(excess[static_cast<std::size_t>(g)] - expected, 0.0) / dose_tot;
              induced_rows.push_back({plist[static_cast<std::size_t>(j)].pair_index, g,
                                      excess[static_cast<std::size_t>(g)], expected, z});
            }
          }
        }
      }

      // ---- expectation-maximization on the nonzero pattern ------------------
      const std::size_t bn = blk.value.size();
      std::vector<double> Theta(static_cast<std::size_t>(n) * K, 0.0);
      std::vector<double> Alpha(static_cast<std::size_t>(n) * S_n, 0.0);
      std::vector<double> Lam(static_cast<std::size_t>(n) * S_n, 0.0);
      std::vector<double> U(static_cast<std::size_t>(n) * S_n, 0.0);
      std::vector<double> Rho(static_cast<std::size_t>(n) * S_n, 1.0);
      std::vector<double> beta(static_cast<std::size_t>(n), lam_a);
      for (int j = 0; j < S_n; ++j) {
        const auto& ps = plist[static_cast<std::size_t>(j)];
        for (int i = 0; i < n; ++i) {
          Lam[static_cast<std::size_t>(i) * S_n + j] = ps.lam_cell[static_cast<std::size_t>(i)];
        }
      }
      Alpha = Lam;
      U = Lam;
      const double eps_g = opt.floor_total / G;
      {
        // Theta init: own share of the cell spread over the programs by the
        // programs' share of the type's pooled expression.
        std::vector<double> colsum(static_cast<std::size_t>(G), 0.0);
        double total_sum = 0.0;
        for (std::size_t p = 0; p < bn; ++p) {
          colsum[static_cast<std::size_t>(blk.gene[p])] += blk.value[p];
          total_sum += blk.value[p];
        }
        std::vector<double> prog_share(static_cast<std::size_t>(K), 0.0);
        double ps_sum = 0.0;
        for (int k = 0; k < K; ++k) {
          double s = 0.0;
          for (int g = 0; g < G; ++g) s += F(k, g) * colsum[static_cast<std::size_t>(g)] / std::max(total_sum, 1.0);
          prog_share[static_cast<std::size_t>(k)] = std::max(s, 1e-3);
          ps_sum += prog_share[static_cast<std::size_t>(k)];
        }
        for (int k = 0; k < K; ++k) prog_share[static_cast<std::size_t>(k)] /= ps_sum;
        for (int i = 0; i < n; ++i) {
          double lam_tot = 0.0;
          for (int j = 0; j < S_n; ++j) lam_tot += Lam[static_cast<std::size_t>(i) * S_n + j];
          const double own0 = std::max(1.0 - lam_tot - lam_a, 0.2);
          for (int k = 0; k < K; ++k) {
            Theta[static_cast<std::size_t>(i) * K + k] = own0 * prog_share[static_cast<std::size_t>(k)];
          }
        }
      }
      std::vector<double> w_dose(static_cast<std::size_t>(n));
      for (int i = 0; i < n; ++i) {
        double lam_tot = 0.0;
        for (int j = 0; j < S_n; ++j) lam_tot += Lam[static_cast<std::size_t>(i) * S_n + j];
        w_dose[static_cast<std::size_t>(i)] = 1.0 / (1.0 + opt.dose_weight * lam_tot);
      }
      std::vector<double> share_psi(static_cast<std::size_t>(S_n), 0.0);
      for (int j = 0; j < S_n; ++j) {
        double s = 0.0;
        for (int g = 0; g < G; ++g) s += plist[static_cast<std::size_t>(j)].psi[static_cast<std::size_t>(g)];
        share_psi[static_cast<std::size_t>(j)] = std::max(s, 1e-6);
      }
      std::vector<double> M_rowsum(static_cast<std::size_t>(S_n), 0.0);
      const auto refresh_M_rowsum = [&]() {
        for (int j = 0; j < S_n; ++j) {
          double s = 0.0;
          for (int g = 0; g < G; ++g) s += M(j, g);
          M_rowsum[static_cast<std::size_t>(j)] = s;
        }
      };
      refresh_M_rowsum();
      const bool any_induced = opt.use_induced &&
          std::any_of(M_rowsum.begin(), M_rowsum.end(), [](double v) { return v > 0.0; });

      const int workers = std::max(1, std::min(opt.num_threads, n));
      std::vector<double> rate_nz(bn);
      std::vector<double> W(bn);
      const auto compute_rates = [&]() {
        parallel_rows(n, opt.num_threads, [&](int begin, int end, int) {
          for (int i = begin; i < end; ++i) {
            const double* th = &Theta[static_cast<std::size_t>(i) * K];
            const double* al = &Alpha[static_cast<std::size_t>(i) * S_n];
            const double* uu = &U[static_cast<std::size_t>(i) * S_n];
            const double* rr = &Rho[static_cast<std::size_t>(i) * S_n];
            const double bi = beta[static_cast<std::size_t>(i)];
            for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
              const int g = blk.gene[static_cast<std::size_t>(p)];
              double r = eps_g;
              for (int k = 0; k < K; ++k) r += th[k] * F(k, g);
              for (int j = 0; j < S_n; ++j) {
                r += al[j] * plist[static_cast<std::size_t>(j)].psi[static_cast<std::size_t>(g)];
              }
              if (lam_a > 0.0) r += bi * a_prof[static_cast<std::size_t>(g)];
              if (any_induced) {
                for (int j = 0; j < S_n; ++j) r += uu[j] * rr[j] * M(j, g);
              }
              rate_nz[static_cast<std::size_t>(p)] = r;
            }
          }
        });
      };
      const auto expected_total = [&]() {
        double tot = 0.0;
        std::vector<double> f_row(static_cast<std::size_t>(K), 0.0);
        for (int k = 0; k < K; ++k) {
          double s = 0.0;
          for (int g = 0; g < G; ++g) s += F(k, g);
          f_row[static_cast<std::size_t>(k)] = s;
        }
        double t_sum = 0.0;
        std::vector<double> th_sum(static_cast<std::size_t>(K), 0.0);
        std::vector<double> al_sum(static_cast<std::size_t>(S_n), 0.0);
        std::vector<double> ur_sum(static_cast<std::size_t>(S_n), 0.0);
        double tb = 0.0;
        for (int i = 0; i < n; ++i) {
          const double ti = blk.totals[static_cast<std::size_t>(i)];
          t_sum += ti;
          tb += ti * beta[static_cast<std::size_t>(i)];
          for (int k = 0; k < K; ++k) th_sum[static_cast<std::size_t>(k)] += ti * Theta[static_cast<std::size_t>(i) * K + k];
          for (int j = 0; j < S_n; ++j) {
            al_sum[static_cast<std::size_t>(j)] += ti * Alpha[static_cast<std::size_t>(i) * S_n + j];
            ur_sum[static_cast<std::size_t>(j)] += ti * U[static_cast<std::size_t>(i) * S_n + j] * Rho[static_cast<std::size_t>(i) * S_n + j];
          }
        }
        for (int k = 0; k < K; ++k) tot += th_sum[static_cast<std::size_t>(k)] * f_row[static_cast<std::size_t>(k)];
        for (int j = 0; j < S_n; ++j) tot += al_sum[static_cast<std::size_t>(j)] * share_psi[static_cast<std::size_t>(j)];
        if (lam_a > 0.0) tot += tb;  // a_prof sums to one
        if (any_induced) {
          for (int j = 0; j < S_n; ++j) tot += ur_sum[static_cast<std::size_t>(j)] * M_rowsum[static_cast<std::size_t>(j)];
        }
        return tot + t_sum * G * eps_g;
      };

      double ll_prev = -std::numeric_limits<double>::infinity();
      const auto em_pass = [&](int n_iter, bool update_F) {
        for (int it = 0; it < n_iter; ++it) {
          compute_rates();
          const bool check_ll = (it + 1) % 10 == 0 || it == n_iter - 1;
          const double exp_tot = check_ll ? expected_total() : 0.0;
          parallel_rows(n, opt.num_threads, [&](int begin, int end, int) {
            for (int i = begin; i < end; ++i) {
              const double ti = blk.totals[static_cast<std::size_t>(i)];
              for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
                W[static_cast<std::size_t>(p)] = blk.value[static_cast<std::size_t>(p)] /
                    std::max(ti * rate_nz[static_cast<std::size_t>(p)], 1e-300);
              }
            }
          });

          // Gene-side accumulations (per-thread buffers, reduced afterwards):
          // wsum_M[j][g] = sum_i t_i U_ij Rho_ij W_ig and
          // f_num[k][g] = sum_i Theta_ik t_i w_dose_i W_ig.
          const bool need_M = any_induced;
          std::vector<std::vector<double>> f_num_t;
          std::vector<std::vector<double>> wsum_t;
          if (update_F) f_num_t.assign(static_cast<std::size_t>(workers), std::vector<double>(static_cast<std::size_t>(K) * G, 0.0));
          if (need_M) wsum_t.assign(static_cast<std::size_t>(workers), std::vector<double>(static_cast<std::size_t>(S_n) * G, 0.0));

          parallel_rows(n, opt.num_threads, [&](int begin, int end, int tid) {
            std::vector<double>* f_acc = update_F ? &f_num_t[static_cast<std::size_t>(tid)] : nullptr;
            std::vector<double>* m_acc = need_M ? &wsum_t[static_cast<std::size_t>(tid)] : nullptr;
            for (int i = begin; i < end; ++i) {
              const double ti = blk.totals[static_cast<std::size_t>(i)];
              double* th = &Theta[static_cast<std::size_t>(i) * K];
              double* al = &Alpha[static_cast<std::size_t>(i) * S_n];
              double* uu = &U[static_cast<std::size_t>(i) * S_n];
              double* rr = &Rho[static_cast<std::size_t>(i) * S_n];
              // row-level dot products with the gene-side factors
              double wf[16];
              double wp[16];
              double wm[16];
              double wa = 0.0;
              for (int k = 0; k < K; ++k) wf[k] = 0.0;
              for (int j = 0; j < S_n; ++j) {
                wp[j] = 0.0;
                wm[j] = 0.0;
              }
              for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
                const int g = blk.gene[static_cast<std::size_t>(p)];
                const double wv = W[static_cast<std::size_t>(p)];
                for (int k = 0; k < K; ++k) wf[k] += wv * F(k, g);
                for (int j = 0; j < S_n; ++j) {
                  wp[j] += wv * plist[static_cast<std::size_t>(j)].psi[static_cast<std::size_t>(g)];
                }
                if (need_M) {
                  for (int j = 0; j < S_n; ++j) wm[j] += wv * M(j, g);
                }
                if (lam_a > 0.0) wa += wv * a_prof[static_cast<std::size_t>(g)];
              }
              // own-program weights (no prior); the F numerator below uses
              // the updated weights, matching the reference update order.
              const double tw = ti * w_dose[static_cast<std::size_t>(i)];
              for (int k = 0; k < K; ++k) th[k] *= wf[k];
              // contamination fractions: bounded evidence update
              for (int j = 0; j < S_n; ++j) {
                const double lam_ij = Lam[static_cast<std::size_t>(i) * S_n + j];
                const double ev = al[j] * wp[j] * ti / share_psi[static_cast<std::size_t>(j)];
                double a_new = (ev + opt.alpha_prior_strength * lam_ij) /
                    (ti + opt.alpha_prior_strength);
                al[j] = std::min(a_new, opt.alpha_cap * lam_ij);
              }
              // ambient scale
              if (lam_a > 0.0) {
                const double evb = beta[static_cast<std::size_t>(i)] * wa * ti;
                beta[static_cast<std::size_t>(i)] = std::min(
                    (evb + opt.ambient_prior_strength * lam_a) / (ti + opt.ambient_prior_strength),
                    opt.ambient_cap * lam_a);
              }
              // induced per-cell activity
              if (need_M) {
                for (int j = 0; j < S_n; ++j) {
                  const double zind = uu[j] * rr[j] * wm[j] * ti;
                  const double denom = ti * uu[j] * M_rowsum[static_cast<std::size_t>(j)];
                  rr[j] = std::min((opt.rho_shape + zind) / (opt.rho_shape + denom), opt.rho_cap);
                }
              }
              // gene-side accumulations use the UPDATED Theta/Rho, matching
              // the reference implementation's update order.
              if (update_F || need_M) {
                for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
                  const int g = blk.gene[static_cast<std::size_t>(p)];
                  const double wv = W[static_cast<std::size_t>(p)];
                  if (update_F) {
                    double* acc = &(*f_acc)[static_cast<std::size_t>(g)];
                    for (int k = 0; k < K; ++k) {
                      acc[static_cast<std::size_t>(k) * G] += th[k] * tw * wv;
                    }
                  }
                  if (need_M) {
                    double* acc = &(*m_acc)[static_cast<std::size_t>(g)];
                    for (int j = 0; j < S_n; ++j) {
                      acc[static_cast<std::size_t>(j) * G] += ti * uu[j] * rr[j] * wv;
                    }
                  }
                }
              }
            }
          });

          if (need_M) {
            std::vector<double> dsum(static_cast<std::size_t>(S_n), 0.0);
            for (int i = 0; i < n; ++i) {
              const double ti = blk.totals[static_cast<std::size_t>(i)];
              for (int j = 0; j < S_n; ++j) {
                dsum[static_cast<std::size_t>(j)] += ti * U[static_cast<std::size_t>(i) * S_n + j] * Rho[static_cast<std::size_t>(i) * S_n + j];
              }
            }
            for (int j = 0; j < S_n; ++j) {
              for (int g = 0; g < G; ++g) {
                double ws = 0.0;
                for (int w = 0; w < workers; ++w) {
                  ws += wsum_t[static_cast<std::size_t>(w)][static_cast<std::size_t>(j) * G + g];
                }
                M(j, g) *= ws / std::max(dsum[static_cast<std::size_t>(j)], 1e-9);
              }
            }
            refresh_M_rowsum();
          }
          if (update_F) {
            for (int k = 0; k < K; ++k) {
              double s = 0.0;
              for (int g = 0; g < G; ++g) {
                double acc = 0.0;
                for (int w = 0; w < workers; ++w) {
                  acc += f_num_t[static_cast<std::size_t>(w)][static_cast<std::size_t>(k) * G + g];
                }
                double v = F(k, g) * acc;
                if (opt.use_ambient && any_strict && strictU[static_cast<std::size_t>(g)]) v = 0.0;
                v += 1e-8;
                F(k, g) = v;
                s += v;
              }
              for (int g = 0; g < G; ++g) F(k, g) /= std::max(s, 1e-12);
            }
          }

          if (check_ll) {
            // log-likelihood at the pre-update rates, as in the reference.
            std::vector<double> partial(static_cast<std::size_t>(workers), 0.0);
            parallel_rows(n, opt.num_threads, [&](int begin, int end, int tid) {
              double local = 0.0;
              for (int i = begin; i < end; ++i) {
                const double ti = blk.totals[static_cast<std::size_t>(i)];
                for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
                  local += blk.value[static_cast<std::size_t>(p)] *
                      std::log(std::max(ti * rate_nz[static_cast<std::size_t>(p)], 1e-300));
                }
              }
              partial[static_cast<std::size_t>(tid)] = local;
            });
            double lsum = 0.0;
            for (const double v : partial) lsum += v;
            const double ll = lsum - exp_tot;
            const double gain = ll - ll_prev;
            ll_prev = ll;
            if (it > 10 && std::abs(gain) < 1e-7 * std::abs(ll)) break;
          }
        }
      };

      em_pass(opt.em_iterations, true);

      // Top-up passes: re-measure the residual exposure gradient of the guide
      // genes on the kept counts (excluding induced-support genes, whose
      // gradient is retained by design), add it to the dose prior, continue.
      std::vector<double> kept(bn);
      const auto compute_kept = [&]() {
        compute_rates();
        parallel_rows(n, opt.num_threads, [&](int begin, int end, int) {
          for (int i = begin; i < end; ++i) {
            const double* th = &Theta[static_cast<std::size_t>(i) * K];
            const double* uu = &U[static_cast<std::size_t>(i) * S_n];
            const double* rr = &Rho[static_cast<std::size_t>(i) * S_n];
            for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
              const int g = blk.gene[static_cast<std::size_t>(p)];
              double own = eps_g;
              for (int k = 0; k < K; ++k) own += th[k] * F(k, g);
              if (any_induced) {
                for (int j = 0; j < S_n; ++j) own += uu[j] * rr[j] * M(j, g);
              }
              kept[static_cast<std::size_t>(p)] = blk.value[static_cast<std::size_t>(p)] * own /
                  std::max(rate_nz[static_cast<std::size_t>(p)], 1e-300);
            }
          }
        });
      };

      double t_all = 0.0;
      for (int i = 0; i < n; ++i) t_all += blk.totals[static_cast<std::size_t>(i)];
      for (int topup = 0; topup < opt.topup_passes; ++topup) {
        compute_kept();
        double added = 0.0;
        for (int j = 0; j < S_n; ++j) {
          auto& ps = plist[static_cast<std::size_t>(j)];
          std::vector<int> guide2;
          for (const int g : ps.guide) {
            if (M(j, g) <= 0.0) guide2.push_back(g);
          }
          if (static_cast<int>(guide2.size()) < 3) continue;
          double share2 = 0.0;
          for (const int g : guide2) share2 += ps.psi[static_cast<std::size_t>(g)];
          if (share2 <= 1e-6) continue;
          PairState tmp = ps;
          tmp.share_guide = share2;
          const auto lam = dose_response(tmp, kept, guide2);
          for (int i = 0; i < n; ++i) {
            const auto it = lam.find(ps.exposure[static_cast<std::size_t>(i)]);
            const double extra = it == lam.end() ? 0.0 : it->second;
            double& l = Lam[static_cast<std::size_t>(i) * S_n + j];
            l = std::min(l + extra, opt.lambda_max);
            added += extra * blk.totals[static_cast<std::size_t>(i)];
          }
        }
        if (added < 1e-4 * t_all) break;
        for (std::size_t q = 0; q < Alpha.size(); ++q) Alpha[q] = std::max(Alpha[q], Lam[q]);
        em_pass(opt.topup_em_iterations, false);
      }

      // Final split, own fractions, removed mass, and per-pair outputs.
      compute_kept();
      std::vector<double> own_num(static_cast<std::size_t>(G), 0.0);
      std::vector<double> own_den(static_cast<std::size_t>(G), 0.0);
      for (int i = 0; i < n; ++i) {
        const double* th = &Theta[static_cast<std::size_t>(i) * K];
        for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
          const int g = blk.gene[static_cast<std::size_t>(p)];
          double own_rate = 0.0;
          for (int k = 0; k < K; ++k) own_rate += th[k] * F(k, g);
          own_num[static_cast<std::size_t>(g)] += blk.value[static_cast<std::size_t>(p)] * own_rate /
              std::max(rate_nz[static_cast<std::size_t>(p)], 1e-300);
          own_den[static_cast<std::size_t>(g)] += blk.value[static_cast<std::size_t>(p)];
        }
      }
      auto& of_new = own_frac_new[static_cast<std::size_t>(T)];
      of_new.assign(static_cast<std::size_t>(G), 0.0);
      for (int g = 0; g < G; ++g) {
        of_new[static_cast<std::size_t>(g)] = own_num[static_cast<std::size_t>(g)] /
            std::max(own_den[static_cast<std::size_t>(g)], 1.0);
      }

      if (last_round) {
        // Both removal policies from the same fit: with retention (the
        // induced share stays) and without (it is removed too).
        compute_rates();
        parallel_rows(n, opt.num_threads, [&](int begin, int end, int) {
          for (int i = begin; i < end; ++i) {
            const double* th = &Theta[static_cast<std::size_t>(i) * K];
            for (int p = blk.rowptr[static_cast<std::size_t>(i)]; p < blk.rowptr[static_cast<std::size_t>(i + 1)]; ++p) {
              const int g = blk.gene[static_cast<std::size_t>(p)];
              double own = eps_g;
              for (int k = 0; k < K; ++k) own += th[k] * F(k, g);
              const double keep2 = own / std::max(rate_nz[static_cast<std::size_t>(p)], 1e-300);
              result.removed_without_retention[static_cast<std::size_t>(blk.gpos[p])] =
                  blk.value[p] * (1.0 - keep2);
            }
          }
        });
        for (std::size_t p = 0; p < bn; ++p) {
          result.removed[static_cast<std::size_t>(blk.gpos[p])] = blk.value[p] - kept[p];
        }
        for (int i = 0; i < n; ++i) {
          result.ambient_scale[static_cast<std::size_t>(blk.cols[static_cast<std::size_t>(i)])] =
              lam_a > 0.0 ? beta[static_cast<std::size_t>(i)] : 0.0;
        }
        for (int j = 0; j < S_n; ++j) {
          const auto& ps = plist[static_cast<std::size_t>(j)];
          const int pj = ps.pair_index;
          auto& cells_out = result.pair_cells[static_cast<std::size_t>(pj)];
          auto& dose_out = result.pair_dose[static_cast<std::size_t>(pj)];
          auto& alpha_out = result.pair_alpha[static_cast<std::size_t>(pj)];
          auto& rho_out = result.pair_rho[static_cast<std::size_t>(pj)];
          cells_out.resize(static_cast<std::size_t>(n));
          dose_out.resize(static_cast<std::size_t>(n));
          alpha_out.resize(static_cast<std::size_t>(n));
          rho_out.resize(static_cast<std::size_t>(n));
          GenerativePairSummary summary;
          summary.pair = pj;
          double dose_sum = 0.0;
          int n_exposed = 0;
          for (int i = 0; i < n; ++i) {
            cells_out[static_cast<std::size_t>(i)] = blk.cols[static_cast<std::size_t>(i)];
            dose_out[static_cast<std::size_t>(i)] = Lam[static_cast<std::size_t>(i) * S_n + j];
            alpha_out[static_cast<std::size_t>(i)] = Alpha[static_cast<std::size_t>(i) * S_n + j];
            rho_out[static_cast<std::size_t>(i)] = Rho[static_cast<std::size_t>(i) * S_n + j];
            const double ti = blk.totals[static_cast<std::size_t>(i)];
            summary.prior_molecules += Lam[static_cast<std::size_t>(i) * S_n + j] * ti;
            summary.posterior_molecules += Alpha[static_cast<std::size_t>(i) * S_n + j] * ti;
            summary.induced_molecules += U[static_cast<std::size_t>(i) * S_n + j] *
                Rho[static_cast<std::size_t>(i) * S_n + j] * ti * M_rowsum[static_cast<std::size_t>(j)];
            if (ps.exposure[static_cast<std::size_t>(i)] > 0.0) {
              dose_sum += Lam[static_cast<std::size_t>(i) * S_n + j];
              ++n_exposed;
            }
          }
          summary.mean_dose_exposed = n_exposed > 0 ? dose_sum / n_exposed : 0.0;
          result.pairs.push_back(summary);
        }
        for (const auto& row : induced_rows) result.induced.push_back(row);
      }
    }
    own_frac = std::move(own_frac_new);
  }

  return result;
}

}  // namespace celladmix
