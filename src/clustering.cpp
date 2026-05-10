#include "celladmix/clustering.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "celladmix/hnsw.hpp"
#include "celladmix/workflow.hpp"
#include "celladmix/xenium.hpp"
#include "irlba/irlba.hpp"
#include "subpar/subpar.hpp"

#include "knncolle/knncolle.hpp"
#include "umappp/umappp.hpp"

namespace celladmix {

// Cell-level feature construction, Louvain clustering, and UMAP embedding.

namespace {

struct WeightedAdjList {
  std::vector<int> indptr;
  std::vector<int> indices;
  std::vector<double> weights;

  int n_nodes() const { return static_cast<int>(indptr.empty() ? 0 : indptr.size() - 1); }

  int neighbor_count(int i) const { return indptr[static_cast<std::size_t>(i + 1)] - indptr[static_cast<std::size_t>(i)]; }

  const int* neighbor_ids(int i) const { return indices.data() + indptr[static_cast<std::size_t>(i)]; }

  const double* neighbor_weights(int i) const { return weights.data() + indptr[static_cast<std::size_t>(i)]; }

  static WeightedAdjList from_edge_list(
      const std::vector<int>& src,
      const std::vector<int>& dst,
      const std::vector<double>& wt,
      int n_nodes) {
    if (src.size() != dst.size() || src.size() != wt.size()) {
      throw std::runtime_error("edge lists must have matching lengths");
    }
    WeightedAdjList out;
    out.indptr.assign(static_cast<std::size_t>(n_nodes + 1), 0);
    for (const int node : src) {
      out.indptr[static_cast<std::size_t>(node + 1)] += 1;
    }
    for (int i = 1; i <= n_nodes; ++i) {
      out.indptr[static_cast<std::size_t>(i)] += out.indptr[static_cast<std::size_t>(i - 1)];
    }
    out.indices.assign(src.size(), 0);
    out.weights.assign(src.size(), 0.0);
    std::vector<int> write_pos = out.indptr;
    for (std::size_t i = 0; i < src.size(); ++i) {
      const int pos = write_pos[static_cast<std::size_t>(src[i])]++;
      out.indices[static_cast<std::size_t>(pos)] = dst[i];
      out.weights[static_cast<std::size_t>(pos)] = wt[i];
    }
    return out;
  }
};

bool is_control_or_codeword_gene(const std::string& gene) {
  return xenium_feature_name_is_control(gene);
}

std::vector<unsigned char> biological_gene_mask(const std::vector<std::string>& genes) {
  std::vector<unsigned char> keep(genes.size(), 1);
  for (std::size_t i = 0; i < genes.size(); ++i) {
    keep[i] = is_control_or_codeword_gene(genes[i]) ? 0 : 1;
  }
  return keep;
}

int count_detected_biological_genes(
    const CellCountMatrix& counts,
    int row,
    const std::vector<unsigned char>& gene_keep) {
  int detected = 0;
  for (int p = counts.indptr[static_cast<std::size_t>(row)];
       p < counts.indptr[static_cast<std::size_t>(row + 1)];
       ++p) {
    const int gene = counts.indices[static_cast<std::size_t>(p)];
    if (gene >= 0 && gene < static_cast<int>(gene_keep.size()) &&
        gene_keep[static_cast<std::size_t>(gene)]) {
      ++detected;
    }
  }
  return detected;
}

struct LouvainGraph {
  WeightedAdjList adj;
  std::vector<double> degree;
  double total_weight_twice = 0.0;
};

struct CellNeighborGraph {
  WeightedAdjList louvain_graph;
  knncolle::NeighborList<int, double> umap_neighbors;
  int umap_neighbor_count = 0;
};

CellTable cell_table_from_transcripts(const TranscriptTable& table, const std::vector<int>& cell_indices);

// Allocate an integer sampling budget proportionally while respecting capacities.
std::vector<int> weighted_capacity_allocation(
    const std::vector<double>& weights,
    const std::vector<int>& capacities,
    int budget,
    std::mt19937& rng) {
  if (weights.size() != capacities.size()) {
    throw std::runtime_error("weights and capacities must have matching lengths");
  }

  std::vector<int> allocation(weights.size(), 0);
  if (budget <= 0 || weights.empty()) {
    return allocation;
  }

  std::vector<std::size_t> active;
  active.reserve(weights.size());
  for (std::size_t i = 0; i < capacities.size(); ++i) {
    if (capacities[i] > 0) {
      active.push_back(i);
    }
  }

  int remaining = budget;
  while (remaining > 0 && !active.empty()) {
    double total_weight = 0.0;
    for (const auto idx : active) {
      total_weight += std::max(weights[idx], 0.0);
    }
    if (total_weight <= 0.0) {
      total_weight = static_cast<double>(active.size());
    }

    std::vector<double> fractional(weights.size(), 0.0);
    int assigned = 0;
    for (const auto idx : active) {
      const int remaining_capacity = capacities[idx] - allocation[idx];
      if (remaining_capacity <= 0) {
        continue;
      }
      const double expected = static_cast<double>(remaining) *
          (total_weight > 0.0 ? std::max(weights[idx], 0.0) : 1.0) / total_weight;
      const int add = std::min<int>(
          remaining_capacity,
          static_cast<int>(std::floor(expected)));
      allocation[idx] += add;
      assigned += add;
      fractional[idx] = expected - std::floor(expected);
    }

    int leftover = remaining - assigned;
    std::shuffle(active.begin(), active.end(), rng);
    std::stable_sort(active.begin(), active.end(), [&](std::size_t lhs, std::size_t rhs) {
      if (fractional[lhs] == fractional[rhs]) {
        return (capacities[lhs] - allocation[lhs]) > (capacities[rhs] - allocation[rhs]);
      }
      return fractional[lhs] > fractional[rhs];
    });
    for (const auto idx : active) {
      if (leftover == 0) {
        break;
      }
      if (capacities[idx] - allocation[idx] <= 0) {
        continue;
      }
      allocation[idx] += 1;
      assigned += 1;
      leftover -= 1;
    }

    if (assigned == 0) {
      break;
    }
    remaining -= assigned;
    active.erase(
        std::remove_if(
            active.begin(),
            active.end(),
            [&](std::size_t idx) { return capacities[idx] - allocation[idx] <= 0; }),
        active.end());
  }

  return allocation;
}

// Recover a crop identifier for a cell from its member transcripts.
std::string crop_id_for_cell(
    const std::vector<std::string>* transcript_crop_ids,
    const TranscriptTable& table,
    int cell_index,
    const std::vector<std::vector<int>>& by_cell) {
  if (transcript_crop_ids == nullptr || transcript_crop_ids->size() != table.size()) {
    return {};
  }
  const auto& members = by_cell[static_cast<std::size_t>(cell_index)];
  if (members.empty()) {
    return {};
  }
  return (*transcript_crop_ids)[static_cast<std::size_t>(members.front())];
}

// Coarsen cell sizes into broad bins for balanced sampling.
int size_bin(int count) {
  if (count < 25) {
    return 0;
  }
  if (count < 100) {
    return 1;
  }
  if (count < 300) {
    return 2;
  }
  return 3;
}

// Downsample clustering cells across crop and size strata.
std::vector<int> sample_clustering_cells(
    const TranscriptTable& table,
    const std::vector<int>& eligible_cells,
    const std::vector<int>& cell_sizes,
    const std::vector<std::vector<int>>& by_cell,
    const std::vector<std::string>* transcript_crop_ids,
    int cells_max,
    unsigned int seed) {
  if (cells_max <= 0 || static_cast<int>(eligible_cells.size()) <= cells_max) {
    return eligible_cells;
  }

  std::mt19937 rng(seed);
  std::vector<std::string> stratum_order;
  std::unordered_map<std::string, int> stratum_lookup;
  std::vector<std::vector<int>> cells_by_stratum;
  for (const int cell : eligible_cells) {
    const std::string crop = crop_id_for_cell(transcript_crop_ids, table, cell, by_cell);
    const std::string key = (crop.empty() ? "__all__" : crop) + "|size" + std::to_string(size_bin(cell_sizes[static_cast<std::size_t>(cell)]));
    const auto it = stratum_lookup.find(key);
    if (it == stratum_lookup.end()) {
      const int idx = static_cast<int>(stratum_order.size());
      stratum_lookup.emplace(key, idx);
      stratum_order.push_back(key);
      cells_by_stratum.push_back({cell});
    } else {
      cells_by_stratum[static_cast<std::size_t>(it->second)].push_back(cell);
    }
  }

  std::vector<int> stratum_quota(cells_by_stratum.size(), 0);
  const int guaranteed = std::min<int>(cells_max, static_cast<int>(cells_by_stratum.size()));
  std::vector<std::size_t> order(cells_by_stratum.size());
  std::iota(order.begin(), order.end(), 0U);
  std::shuffle(order.begin(), order.end(), rng);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
    return cells_by_stratum[lhs].size() < cells_by_stratum[rhs].size();
  });
  for (int i = 0; i < guaranteed; ++i) {
    stratum_quota[order[static_cast<std::size_t>(i)]] = 1;
  }

  int remaining = cells_max - guaranteed;
  if (remaining > 0) {
    std::vector<double> weights(cells_by_stratum.size(), 0.0);
    std::vector<int> capacities(cells_by_stratum.size(), 0);
    for (std::size_t stratum = 0; stratum < cells_by_stratum.size(); ++stratum) {
      capacities[stratum] = std::max<int>(
          static_cast<int>(cells_by_stratum[stratum].size()) - stratum_quota[stratum],
          0);
      weights[stratum] = capacities[stratum] > 0
          ? std::sqrt(static_cast<double>(capacities[stratum]))
          : 0.0;
    }
    const auto extra = weighted_capacity_allocation(weights, capacities, remaining, rng);
    for (std::size_t stratum = 0; stratum < cells_by_stratum.size(); ++stratum) {
      stratum_quota[stratum] += extra[stratum];
    }
  }

  std::vector<int> out;
  out.reserve(static_cast<std::size_t>(cells_max));
  for (std::size_t stratum = 0; stratum < cells_by_stratum.size(); ++stratum) {
    auto pool = cells_by_stratum[stratum];
    std::shuffle(pool.begin(), pool.end(), rng);
    pool.resize(static_cast<std::size_t>(std::min<int>(
        stratum_quota[stratum],
        static_cast<int>(pool.size()))));
    out.insert(out.end(), pool.begin(), pool.end());
  }
  std::sort(out.begin(), out.end());
  return out;
}

// Build sparse cell-by-gene counts for the eligible clustering cells.
CellCountMatrix build_sparse_counts(
    const TranscriptTable& table,
    const std::vector<int>& eligible_cells,
    const std::vector<std::vector<int>>& by_cell,
    const std::vector<std::string>* transcript_crop_ids) {
  CellCountMatrix out;
  out.genes = table.genes;
  out.cells = cell_table_from_transcripts(table, eligible_cells);
  out.transcript_counts.reserve(eligible_cells.size());
  out.detected_genes.reserve(eligible_cells.size());
  out.crop_ids.reserve(eligible_cells.size());
  out.indptr.reserve(eligible_cells.size() + 1);
  out.indptr.push_back(0);
  const auto gene_keep = biological_gene_mask(table.genes);

  for (const int cell : eligible_cells) {
    const auto& members = by_cell[static_cast<std::size_t>(cell)];
    std::unordered_map<int, int> gene_counts;
    gene_counts.reserve(members.size());
    for (const int tx : members) {
      gene_counts[table.gene_index[static_cast<std::size_t>(tx)]] += 1;
    }

    std::vector<std::pair<int, int>> entries;
    entries.reserve(gene_counts.size());
    for (const auto& kv : gene_counts) {
      entries.emplace_back(kv.first, kv.second);
    }
    std::sort(entries.begin(), entries.end(), [](const auto& lhs, const auto& rhs) {
      return lhs.first < rhs.first;
    });

    for (const auto& entry : entries) {
      out.indices.push_back(entry.first);
      out.values.push_back(static_cast<double>(entry.second));
    }
    out.transcript_counts.push_back(static_cast<int>(members.size()));
    int detected = 0;
    for (const auto& entry : entries) {
      const int gene = entry.first;
      if (gene >= 0 && gene < static_cast<int>(gene_keep.size()) &&
          gene_keep[static_cast<std::size_t>(gene)]) {
        ++detected;
      }
    }
    out.detected_genes.push_back(detected);
    out.crop_ids.push_back(crop_id_for_cell(transcript_crop_ids, table, cell, by_cell));
    out.indptr.push_back(static_cast<int>(out.indices.size()));
  }

  return out;
}

// Rank genes by variance after per-cell normalization and log scaling.
std::vector<int> select_variable_genes(
    const CellCountMatrix& counts,
    int n_genes,
    int n_variable_genes,
    double normalization_scale) {
  const int n_cells = static_cast<int>(counts.transcript_counts.size());
  const auto gene_keep = biological_gene_mask(counts.genes);
  std::vector<double> sum(static_cast<std::size_t>(n_genes), 0.0);
  std::vector<double> sumsq(static_cast<std::size_t>(n_genes), 0.0);
  std::vector<int> detected(static_cast<std::size_t>(n_genes), 0);

  for (int row = 0; row < n_cells; ++row) {
    const double denom = std::max(counts.transcript_counts[static_cast<std::size_t>(row)], 1);
    for (int p = counts.indptr[static_cast<std::size_t>(row)];
         p < counts.indptr[static_cast<std::size_t>(row + 1)];
         ++p) {
      const int gene = counts.indices[static_cast<std::size_t>(p)];
      if (gene < 0 || gene >= static_cast<int>(gene_keep.size()) ||
          !gene_keep[static_cast<std::size_t>(gene)]) {
        continue;
      }
      const double normalized = std::log1p(
          counts.values[static_cast<std::size_t>(p)] * normalization_scale / denom);
      sum[static_cast<std::size_t>(gene)] += normalized;
      sumsq[static_cast<std::size_t>(gene)] += normalized * normalized;
      detected[static_cast<std::size_t>(gene)] += 1;
    }
  }

  std::vector<std::pair<double, int>> ranked;
  ranked.reserve(static_cast<std::size_t>(n_genes));
  for (int gene = 0; gene < n_genes; ++gene) {
    if (detected[static_cast<std::size_t>(gene)] < 2) {
      continue;
    }
    const double mean = sum[static_cast<std::size_t>(gene)] / std::max(n_cells, 1);
    double variance = sumsq[static_cast<std::size_t>(gene)] / std::max(n_cells, 1) - mean * mean;
    if (variance < 0.0) {
      variance = 0.0;
    }
    if (variance > 0.0) {
      ranked.emplace_back(variance, gene);
    }
  }

  if (ranked.empty()) {
    for (int gene = 0; gene < n_genes; ++gene) {
      if (gene < static_cast<int>(gene_keep.size()) &&
          gene_keep[static_cast<std::size_t>(gene)] &&
          detected[static_cast<std::size_t>(gene)] > 0) {
        ranked.emplace_back(1.0, gene);
      }
    }
  }
  if (ranked.empty()) {
    throw std::runtime_error("no informative genes available for clustering");
  }

  const int keep = std::min<int>(
      n_variable_genes > 0 ? n_variable_genes : static_cast<int>(ranked.size()),
      static_cast<int>(ranked.size()));
  std::partial_sort(
      ranked.begin(),
      ranked.begin() + keep,
      ranked.end(),
      [](const auto& lhs, const auto& rhs) {
        if (lhs.first == rhs.first) {
          return lhs.second < rhs.second;
        }
        return lhs.first > rhs.first;
      });

  std::vector<int> selected;
  selected.reserve(static_cast<std::size_t>(keep));
  for (int i = 0; i < keep; ++i) {
    selected.push_back(ranked[static_cast<std::size_t>(i)].second);
  }
  std::sort(selected.begin(), selected.end());
  return selected;
}

// Materialize the selected genes into a dense feature matrix for PCA.
Eigen::MatrixXd materialize_feature_matrix(
    const CellCountMatrix& counts,
    const std::vector<int>& selected_genes,
    double normalization_scale) {
  const int n_cells = static_cast<int>(counts.transcript_counts.size());
  const int n_features = static_cast<int>(selected_genes.size());
  std::unordered_map<int, int> gene_lookup;
  gene_lookup.reserve(selected_genes.size());
  for (int i = 0; i < n_features; ++i) {
    gene_lookup.emplace(selected_genes[static_cast<std::size_t>(i)], i);
  }

  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n_features, n_cells);
  for (int row = 0; row < n_cells; ++row) {
    const double denom = std::max(counts.transcript_counts[static_cast<std::size_t>(row)], 1);
    for (int p = counts.indptr[static_cast<std::size_t>(row)];
         p < counts.indptr[static_cast<std::size_t>(row + 1)];
         ++p) {
      const int gene = counts.indices[static_cast<std::size_t>(p)];
      const auto it = gene_lookup.find(gene);
      if (it == gene_lookup.end()) {
        continue;
      }
      out(it->second, row) = std::log1p(
          counts.values[static_cast<std::size_t>(p)] * normalization_scale / denom);
    }
  }
  return out;
}

// Hold PCA coordinates and explained-variance fractions.
struct PcaResult {
  Eigen::MatrixXd pcs;
  std::vector<double> variance_explained;
};

// Run truncated PCA on the cell feature matrix and keep the leading dimensions.
PcaResult run_pca(Eigen::MatrixXd features, int pca_dims, unsigned int seed) {
  if (features.cols() == 0 || features.rows() == 0) {
    throw std::runtime_error("empty feature matrix");
  }
  const Eigen::VectorXd mean = features.rowwise().mean();
  features.colwise() -= mean;
  const int k = std::max<int>(
      1,
      std::min<int>(
          pca_dims,
          std::min<int>(features.rows(), features.cols())));

  const double total_var = features.squaredNorm();
  if (total_var <= 0.0) {
    std::vector<double> variance(static_cast<std::size_t>(k), 0.0);
    return {
        Eigen::MatrixXd::Zero(k, features.cols()),
        std::move(variance)};
  }

  irlba::Options irlba_options;
  irlba_options.seed = static_cast<std::uint64_t>(seed);
  irlba_options.extra_work = std::max(10, std::min(30, k));
  Eigen::MatrixXd left;
  Eigen::MatrixXd right;
  Eigen::VectorXd singular;
  const auto irlba_status = irlba::compute(
      features,
      static_cast<Eigen::Index>(k),
      left,
      right,
      singular,
      irlba_options);
  if (!irlba_status.first) {
    throw std::runtime_error("IRLBA PCA did not converge within the iteration limit");
  }

  Eigen::MatrixXd pcs = singular.asDiagonal() * right.transpose();
  std::vector<double> variance;
  variance.reserve(static_cast<std::size_t>(k));
  for (int i = 0; i < k; ++i) {
    const double value = singular[static_cast<Eigen::Index>(i)];
    variance.push_back((value * value) / total_var);
  }
  return {std::move(pcs), std::move(variance)};
}

// L2-normalize column vectors before cosine-like neighbor search.
Eigen::MatrixXd normalize_columns_l2(const Eigen::MatrixXd& matrix) {
  Eigen::MatrixXd out = matrix;
  for (Eigen::Index col = 0; col < out.cols(); ++col) {
    const double norm = out.col(col).norm();
    if (norm > 1e-12) {
      out.col(col) /= norm;
    } else {
      out.col(col).setZero();
    }
  }
  return out;
}

// Build a Louvain-oriented HNSW graph and a matching UMAP neighbor list.
CellNeighborGraph build_cell_neighbor_graph(
    const Eigen::MatrixXd& vectors,
    int k,
    int umap_neighbors,
    int num_threads,
    unsigned int seed) {
  const int n = static_cast<int>(vectors.cols());
  CellNeighborGraph out;
  out.louvain_graph.indptr.assign(static_cast<std::size_t>(n + 1), 0);
  out.umap_neighbors.resize(static_cast<std::size_t>(n));
  if (n <= 1 || k <= 0) {
    return out;
  }

  const int graph_k = std::min(k, n - 1);
  const int umap_k = umap_neighbors > 0 ? std::min(umap_neighbors, n - 1) : 0;
  const int neighbor_k = std::max(graph_k, umap_k);
  out.umap_neighbor_count = umap_k;
  const Eigen::MatrixXd normalized = normalize_columns_l2(vectors);

  const int workers = std::max(1, num_threads);
  HnswNeighborOptions hnsw_options;
  hnsw_options.seed = seed;
  hnsw_options.exact_threshold = 128;
  const auto neighbors = hnsw_l2_neighbors_from_columns(
      normalized,
      neighbor_k,
      workers,
      hnsw_options);

  std::vector<std::vector<int>> per_cell_dst(static_cast<std::size_t>(n));
  std::vector<std::vector<double>> per_cell_wt(static_cast<std::size_t>(n));
  std::vector<std::vector<double>> per_cell_dist(static_cast<std::size_t>(n));
  subpar::parallelize_range(workers, n, [&](int, int start, int length) {
    for (int i = start, end = start + length; i < end; ++i) {
      auto& cell_dst = per_cell_dst[static_cast<std::size_t>(i)];
      auto& cell_wt = per_cell_wt[static_cast<std::size_t>(i)];
      auto& cell_dist = per_cell_dist[static_cast<std::size_t>(i)];
      cell_dst.reserve(static_cast<std::size_t>(neighbor_k));
      cell_wt.reserve(static_cast<std::size_t>(neighbor_k));
      cell_dist.reserve(static_cast<std::size_t>(neighbor_k));
      const auto& current = neighbors[static_cast<std::size_t>(i)];
      for (const auto& hit : current) {
        const int nb = hit.first;
        const double cosine = std::max(0.0, std::min(1.0, 1.0 - 0.5 * hit.second));
        if (cosine <= 0.0) {
          continue;
        }
        cell_dst.push_back(nb);
        cell_wt.push_back(cosine);
        cell_dist.push_back(std::max(0.0, 1.0 - cosine));
      }
    }
  });

  std::size_t total_edges = 0;
  for (const auto& edges : per_cell_dst) {
    total_edges += edges.size();
  }
  std::vector<int> src;
  std::vector<int> dst;
  std::vector<double> wt;
  src.reserve(total_edges);
  dst.reserve(total_edges);
  wt.reserve(total_edges);
  for (int i = 0; i < n; ++i) {
    const auto& cell_dst = per_cell_dst[static_cast<std::size_t>(i)];
    const auto& cell_wt = per_cell_wt[static_cast<std::size_t>(i)];
    const auto& cell_dist = per_cell_dist[static_cast<std::size_t>(i)];
    auto& umap_cell = out.umap_neighbors[static_cast<std::size_t>(i)];
    umap_cell.reserve(static_cast<std::size_t>(std::min<int>(umap_k, cell_dst.size())));
    for (std::size_t j = 0; j < cell_dst.size(); ++j) {
      if (static_cast<int>(j) < graph_k) {
        src.push_back(i);
        dst.push_back(cell_dst[j]);
        wt.push_back(cell_wt[j]);
      }
      if (static_cast<int>(j) < umap_k) {
        umap_cell.push_back({cell_dst[j], cell_dist[j]});
      }
    }
  }

  out.louvain_graph = WeightedAdjList::from_edge_list(src, dst, wt, n);
  return out;
}

// Relabel arbitrary community ids into a compact zero-based range.
std::vector<int> reindex_membership_zero_based(const std::vector<int>& membership, int* n_clusters_out = nullptr) {
  std::unordered_map<int, int> remap;
  remap.reserve(membership.size());
  std::vector<int> out(membership.size(), 0);
  int next_id = 0;
  for (std::size_t i = 0; i < membership.size(); ++i) {
    auto it = remap.find(membership[i]);
    if (it == remap.end()) {
      it = remap.emplace(membership[i], next_id++).first;
    }
    out[i] = it->second;
  }
  if (n_clusters_out != nullptr) {
    *n_clusters_out = next_id;
  }
  return out;
}

// Wrap the adjacency list with node degrees and total graph weight for Louvain.
LouvainGraph make_louvain_graph(const WeightedAdjList& graph) {
  LouvainGraph out;
  out.adj = graph;
  out.degree.resize(static_cast<std::size_t>(graph.n_nodes()), 0.0);
  for (int i = 0; i < graph.n_nodes(); ++i) {
    double degree = 0.0;
    for (int j = 0; j < graph.neighbor_count(i); ++j) {
      degree += graph.neighbor_weights(i)[j];
    }
    out.degree[static_cast<std::size_t>(i)] = degree;
    out.total_weight_twice += degree;
  }
  return out;
}

// Compute modularity gain for moving one node into one target community.
double community_edge_gain(
    int node,
    int target_comm,
    double resolution,
    const std::vector<double>& comm_tot,
    const std::vector<double>& neigh_w,
    const LouvainGraph& graph) {
  return neigh_w[static_cast<std::size_t>(target_comm)] -
      resolution * graph.degree[static_cast<std::size_t>(node)] *
          comm_tot[static_cast<std::size_t>(target_comm)] /
          std::max(graph.total_weight_twice, 1e-12);
}

// Run one Louvain pass over the current graph level.
bool louvain_one_level(
    const LouvainGraph& graph,
    double resolution,
    std::vector<int>& membership,
    int max_passes) {
  const int n = graph.adj.n_nodes();
  if (n == 0) {
    return false;
  }

  int n_comms = 0;
  membership = reindex_membership_zero_based(membership, &n_comms);
  std::vector<double> comm_tot(static_cast<std::size_t>(n_comms), 0.0);
  for (int i = 0; i < n; ++i) {
    comm_tot[static_cast<std::size_t>(membership[static_cast<std::size_t>(i)])] +=
        graph.degree[static_cast<std::size_t>(i)];
  }

  bool moved_any = false;
  std::vector<double> neigh_w(static_cast<std::size_t>(std::max(n_comms, 1)), 0.0);
  std::vector<int> touched;
  touched.reserve(64);

  for (int pass = 0; pass < max_passes; ++pass) {
    int moved = 0;
    for (int node = 0; node < n; ++node) {
      const int cur_comm = membership[static_cast<std::size_t>(node)];
      const int nc = graph.adj.neighbor_count(node);
      const int* nb_ids = graph.adj.neighbor_ids(node);
      const double* nb_wts = graph.adj.neighbor_weights(node);

      touched.clear();
      for (int ai = 0; ai < nc; ++ai) {
        const int nb = nb_ids[ai];
        if (nb == node) {
          continue;
        }
        const int comm = membership[static_cast<std::size_t>(nb)];
        if (neigh_w[static_cast<std::size_t>(comm)] == 0.0) {
          touched.push_back(comm);
        }
        neigh_w[static_cast<std::size_t>(comm)] += nb_wts[ai];
      }

      comm_tot[static_cast<std::size_t>(cur_comm)] -= graph.degree[static_cast<std::size_t>(node)];

      int best_comm = cur_comm;
      double best_gain = 0.0;
      for (const int comm : touched) {
        const double gain = community_edge_gain(node, comm, resolution, comm_tot, neigh_w, graph);
        if (gain > best_gain + 1e-12 ||
            (std::abs(gain - best_gain) <= 1e-12 && comm < best_comm)) {
          best_gain = gain;
          best_comm = comm;
        }
      }

      membership[static_cast<std::size_t>(node)] = best_comm;
      comm_tot[static_cast<std::size_t>(best_comm)] += graph.degree[static_cast<std::size_t>(node)];
      if (best_comm != cur_comm) {
        moved_any = true;
        ++moved;
      }

      for (const int comm : touched) {
        neigh_w[static_cast<std::size_t>(comm)] = 0.0;
      }
    }

    if (moved == 0) {
      break;
    }
  }

  return moved_any;
}

// Collapse a graph according to the current community assignment.
WeightedAdjList aggregate_graph(
    const LouvainGraph& graph,
    const std::vector<int>& membership,
    int n_comms) {
  auto key = [](int a, int b) -> std::uint64_t {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) |
        static_cast<std::uint32_t>(b);
  };

  std::unordered_map<std::uint64_t, double> edge_weight;
  edge_weight.reserve(graph.adj.indices.size());
  for (int i = 0; i < graph.adj.n_nodes(); ++i) {
    const int ci = membership[static_cast<std::size_t>(i)];
    for (int j = 0; j < graph.adj.neighbor_count(i); ++j) {
      const int nb = graph.adj.neighbor_ids(i)[j];
      const double w = graph.adj.neighbor_weights(i)[j];
      if (nb < i) {
        continue;
      }
      const int cj = membership[static_cast<std::size_t>(nb)];
      if (i == nb) {
        edge_weight[key(ci, ci)] += w;
      } else if (ci == cj) {
        edge_weight[key(ci, ci)] += 2.0 * w;
      } else {
        const int a = std::min(ci, cj);
        const int b = std::max(ci, cj);
        edge_weight[key(a, b)] += w;
      }
    }
  }

  std::vector<int> src;
  std::vector<int> dst;
  std::vector<double> wt;
  src.reserve(edge_weight.size() * 2);
  dst.reserve(edge_weight.size() * 2);
  wt.reserve(edge_weight.size() * 2);

  for (const auto& kv : edge_weight) {
    const int a = static_cast<int>(kv.first >> 32);
    const int b = static_cast<int>(kv.first & 0xffffffffu);
    src.push_back(a);
    dst.push_back(b);
    wt.push_back(kv.second);
    if (a != b) {
      src.push_back(b);
      dst.push_back(a);
      wt.push_back(kv.second);
    }
  }
  return WeightedAdjList::from_edge_list(src, dst, wt, n_comms);
}

// Run full Louvain clustering and return zero-based community ids for original nodes.
std::vector<int> run_louvain_zero_based(
    const WeightedAdjList& graph,
    double resolution,
    int max_passes) {
  LouvainGraph current = make_louvain_graph(graph);
  const int n0 = current.adj.n_nodes();
  std::vector<int> membership(static_cast<std::size_t>(n0));
  std::iota(membership.begin(), membership.end(), 0);
  std::vector<int> finest_to_current = membership;

  for (int level = 0; level < max_passes; ++level) {
    const bool moved = louvain_one_level(current, resolution, membership, max_passes);
    int n_comms = 0;
    membership = reindex_membership_zero_based(membership, &n_comms);
    for (int& id : finest_to_current) {
      id = membership[static_cast<std::size_t>(id)];
    }
    if (!moved || n_comms == current.adj.n_nodes()) {
      break;
    }
    current = make_louvain_graph(aggregate_graph(current, membership, n_comms));
    membership.resize(static_cast<std::size_t>(current.adj.n_nodes()));
    std::iota(membership.begin(), membership.end(), 0);
  }

  return reindex_membership_zero_based(finest_to_current);
}

Eigen::MatrixXd umap_embed_neighbors_impl(
    knncolle::NeighborList<int, double> neighbors,
    int nobs,
    int ndim_out,
    int n_epochs,
    int seed,
    int num_threads,
    bool parallel_optimization,
    const std::function<void(const std::string&, double)>* progress = nullptr) {
  auto stage_start = std::chrono::steady_clock::now();
  auto emit_progress = [&](const std::string& message) {
    if (progress == nullptr || !*progress) {
      return;
    }
    const auto now = std::chrono::steady_clock::now();
    (*progress)(message, std::chrono::duration<double>(now - stage_start).count());
    stage_start = now;
  };

  if (nobs == 0) {
    return Eigen::MatrixXd(ndim_out, 0);
  }

  std::vector<double> embedding(static_cast<std::size_t>(ndim_out * nobs));
  std::mt19937 rng(static_cast<std::mt19937::result_type>(seed));
  std::uniform_real_distribution<double> uniform(-10.0, 10.0);
  for (double& value : embedding) {
    value = uniform(rng);
  }

  umappp::Options options;
  options.num_epochs = n_epochs;
  options.seed = static_cast<std::uint64_t>(seed);
  options.spread = 1.0;
  options.min_dist = 0.1;
  options.initialize = umappp::InitializeMethod::NONE;
  options.num_threads = std::max(1, num_threads);
  options.parallel_optimization = parallel_optimization;

  auto status = umappp::initialize(std::move(neighbors), ndim_out, embedding.data(), options);
  emit_progress("Initialized UMAP layout");
  status.run();
  emit_progress(
      std::string("Optimized UMAP layout") +
      (parallel_optimization ? " with parallel optimization" : ""));
  return Eigen::Map<Eigen::MatrixXd>(embedding.data(), ndim_out, nobs);
}

// Run UMAP on an observation-by-feature matrix with HNSW neighbor search.
Eigen::MatrixXd umap_embed_impl(
    const Eigen::MatrixXd& data,
    int ndim_out,
    int n_neighbors,
    int n_epochs,
    int seed,
    int num_threads,
    bool parallel_optimization,
    const std::function<void(const std::string&, double)>* progress = nullptr) {
  auto stage_start = std::chrono::steady_clock::now();
  auto emit_progress = [&](const std::string& message) {
    if (progress == nullptr || !*progress) {
      return;
    }
    const auto now = std::chrono::steady_clock::now();
    (*progress)(message, std::chrono::duration<double>(now - stage_start).count());
    stage_start = now;
  };

  const int ndim_in = static_cast<int>(data.rows());
  const int nobs = static_cast<int>(data.cols());
  if (nobs == 0 || ndim_in == 0) {
    return Eigen::MatrixXd(ndim_out, 0);
  }
  n_neighbors = std::min(n_neighbors, nobs - 1);

  HnswNeighborOptions hnsw_options;
  hnsw_options.seed = static_cast<unsigned int>(std::max(1, seed));
  auto neighbors = hnsw_l2_neighbors_from_columns(
      data,
      n_neighbors,
      std::max(1, num_threads),
      hnsw_options);
  for (auto& current : neighbors) {
    for (auto& hit : current) {
      hit.second = std::sqrt(std::max(0.0, hit.second));
    }
  }
  emit_progress(
      "Built and queried HNSW UMAP neighbors using " + std::to_string(std::max(1, num_threads)) +
      " thread(s)");
  return umap_embed_neighbors_impl(
      std::move(neighbors),
      nobs,
      ndim_out,
      n_epochs,
      seed,
      num_threads,
      parallel_optimization,
      progress);
}

// Convert an Eigen row-major style embedding into the package DenseMatrix type.
DenseMatrix eigen_to_dense_rows(const Eigen::MatrixXd& matrix) {
  DenseMatrix out(static_cast<int>(matrix.rows()), static_cast<int>(matrix.cols()), 0.0);
  for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
    for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
      out(static_cast<int>(row), static_cast<int>(col)) = matrix(row, col);
    }
  }
  return out;
}

// Build a minimal cell table for the clustered cell subset.
CellTable cell_table_from_transcripts(const TranscriptTable& table, const std::vector<int>& cell_indices) {
  const std::vector<double> centroid_x = infer_cell_centroids_x(table);
  const std::vector<double> centroid_y = infer_cell_centroids_y(table);
  const std::vector<double> centroid_z = infer_cell_centroids_z(table);

  CellTable out;
  out.cell_ids.reserve(cell_indices.size());
  out.centroid_x.reserve(cell_indices.size());
  out.centroid_y.reserve(cell_indices.size());
  out.centroid_z.reserve(cell_indices.size());
  if (table.has_cell_types()) {
    out.cell_types.assign(cell_indices.size(), "");
  }

  for (const int cell_index : cell_indices) {
    out.cell_ids.push_back(table.cells[static_cast<std::size_t>(cell_index)]);
    out.centroid_x.push_back(centroid_x[static_cast<std::size_t>(cell_index)]);
    out.centroid_y.push_back(centroid_y[static_cast<std::size_t>(cell_index)]);
    out.centroid_z.push_back(centroid_z[static_cast<std::size_t>(cell_index)]);
  }

  if (table.has_cell_types()) {
    std::vector<std::string> cell_types(table.num_cells(), "");
    for (std::size_t i = 0; i < table.size(); ++i) {
      const std::size_t cell = static_cast<std::size_t>(table.cell_index[i]);
      if (cell_types[cell].empty()) {
        cell_types[cell] = table.cell_types[i];
      }
    }
    out.cell_types.clear();
    out.cell_types.reserve(cell_indices.size());
    for (const int cell_index : cell_indices) {
      out.cell_types.push_back(cell_types[static_cast<std::size_t>(cell_index)]);
    }
  }

  return out;
}

}  // namespace

// Public UMAP wrapper shared between clustering and report helpers.
Eigen::MatrixXd umap_embed(
    const Eigen::MatrixXd& data,
    int ndim_out,
    int n_neighbors,
    int n_epochs,
    int seed,
    int num_threads,
    bool parallel_optimization) {
  return umap_embed_impl(data, ndim_out, n_neighbors, n_epochs, seed, num_threads, parallel_optimization);
}

// Build cell features, cluster cells, and return the embedding and summaries.
CellClusteringResult cluster_cell_counts(
    const CellCountMatrix& counts_all,
    const CellClusteringOptions& options,
    const std::string& progress_prefix) {
  if (counts_all.transcript_counts.empty()) {
    throw std::runtime_error("cannot cluster an empty cell count matrix");
  }

  auto stage_start = std::chrono::steady_clock::now();
  auto emit_progress = [&](const std::string& message) {
    if (!options.progress) {
      return;
    }
    const auto now = std::chrono::steady_clock::now();
    options.progress(message, std::chrono::duration<double>(now - stage_start).count());
    stage_start = now;
  };

  std::vector<int> eligible_rows;
  eligible_rows.reserve(counts_all.transcript_counts.size());
  const auto gene_keep = biological_gene_mask(counts_all.genes);
  for (std::size_t row = 0; row < counts_all.transcript_counts.size(); ++row) {
    const int detected_genes = row < counts_all.detected_genes.size()
        ? counts_all.detected_genes[row]
        : count_detected_biological_genes(counts_all, static_cast<int>(row), gene_keep);
    if (counts_all.transcript_counts[row] >= std::max(options.min_molecules, 1) &&
        detected_genes >= std::max(options.min_genes, 0)) {
      eligible_rows.push_back(static_cast<int>(row));
    }
  }
  if (eligible_rows.empty()) {
    throw std::runtime_error("no cells meet the clustering min_molecules/min_genes thresholds");
  }
  emit_progress(
      progress_prefix + "Indexed cell counts: " + std::to_string(eligible_rows.size()) +
      " eligible cells (min_molecules=" + std::to_string(std::max(options.min_molecules, 1)) +
      ", min_genes=" + std::to_string(std::max(options.min_genes, 0)) + ")");

  if (options.cells_max > 0 && static_cast<int>(eligible_rows.size()) > options.cells_max) {
    std::mt19937 rng(options.seed);
    std::vector<std::string> stratum_order;
    std::unordered_map<std::string, int> stratum_lookup;
    std::vector<std::vector<int>> rows_by_stratum;
    for (const int row : eligible_rows) {
      const std::string crop = row < static_cast<int>(counts_all.crop_ids.size())
          ? counts_all.crop_ids[static_cast<std::size_t>(row)]
          : "";
      const std::string key = (crop.empty() ? "__all__" : crop) + "|size" +
          std::to_string(size_bin(counts_all.transcript_counts[static_cast<std::size_t>(row)]));
      const auto it = stratum_lookup.find(key);
      if (it == stratum_lookup.end()) {
        const int idx = static_cast<int>(stratum_order.size());
        stratum_lookup.emplace(key, idx);
        stratum_order.push_back(key);
        rows_by_stratum.push_back({row});
      } else {
        rows_by_stratum[static_cast<std::size_t>(it->second)].push_back(row);
      }
    }

    std::vector<int> quota(rows_by_stratum.size(), 0);
    const int guaranteed = std::min<int>(options.cells_max, static_cast<int>(rows_by_stratum.size()));
    std::vector<std::size_t> order(rows_by_stratum.size());
    std::iota(order.begin(), order.end(), 0U);
    std::shuffle(order.begin(), order.end(), rng);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
      return rows_by_stratum[lhs].size() < rows_by_stratum[rhs].size();
    });
    for (int i = 0; i < guaranteed; ++i) {
      quota[order[static_cast<std::size_t>(i)]] = 1;
    }
    const int remaining = options.cells_max - guaranteed;
    if (remaining > 0) {
      std::vector<double> weights(rows_by_stratum.size(), 0.0);
      std::vector<int> capacities(rows_by_stratum.size(), 0);
      for (std::size_t stratum = 0; stratum < rows_by_stratum.size(); ++stratum) {
        capacities[stratum] = std::max<int>(
            static_cast<int>(rows_by_stratum[stratum].size()) - quota[stratum],
            0);
        weights[stratum] = capacities[stratum] > 0 ? std::sqrt(static_cast<double>(capacities[stratum])) : 0.0;
      }
      const auto extra = weighted_capacity_allocation(weights, capacities, remaining, rng);
      for (std::size_t stratum = 0; stratum < rows_by_stratum.size(); ++stratum) {
        quota[stratum] += extra[stratum];
      }
    }
    eligible_rows.clear();
    for (std::size_t stratum = 0; stratum < rows_by_stratum.size(); ++stratum) {
      auto pool = rows_by_stratum[stratum];
      std::shuffle(pool.begin(), pool.end(), rng);
      pool.resize(static_cast<std::size_t>(std::min<int>(quota[stratum], static_cast<int>(pool.size()))));
      eligible_rows.insert(eligible_rows.end(), pool.begin(), pool.end());
    }
    std::sort(eligible_rows.begin(), eligible_rows.end());
  }
  emit_progress(
      progress_prefix + "Selected clustering cells: " + std::to_string(eligible_rows.size()) +
      " cells");

  CellCountMatrix counts_subset;
  counts_subset.genes = counts_all.genes;
  counts_subset.indptr.reserve(eligible_rows.size() + 1);
  counts_subset.indptr.push_back(0);
  counts_subset.cells.cell_ids.reserve(eligible_rows.size());
  counts_subset.cells.centroid_x.reserve(eligible_rows.size());
  counts_subset.cells.centroid_y.reserve(eligible_rows.size());
  counts_subset.cells.centroid_z.reserve(eligible_rows.size());
  if (!counts_all.cells.cell_types.empty()) counts_subset.cells.cell_types.reserve(eligible_rows.size());
  if (!counts_all.cells.sample_ids.empty()) counts_subset.cells.sample_ids.reserve(eligible_rows.size());
  if (!counts_all.cells.fov_ids.empty()) counts_subset.cells.fov_ids.reserve(eligible_rows.size());
  counts_subset.transcript_counts.reserve(eligible_rows.size());
  counts_subset.detected_genes.reserve(eligible_rows.size());
  counts_subset.crop_ids.reserve(eligible_rows.size());
  for (const int source_row : eligible_rows) {
    const auto row = static_cast<std::size_t>(source_row);
    counts_subset.cells.cell_ids.push_back(counts_all.cells.cell_ids[row]);
    counts_subset.cells.centroid_x.push_back(counts_all.cells.centroid_x[row]);
    counts_subset.cells.centroid_y.push_back(counts_all.cells.centroid_y[row]);
    counts_subset.cells.centroid_z.push_back(counts_all.cells.centroid_z[row]);
    if (!counts_all.cells.cell_types.empty()) counts_subset.cells.cell_types.push_back(counts_all.cells.cell_types[row]);
    if (!counts_all.cells.sample_ids.empty()) counts_subset.cells.sample_ids.push_back(counts_all.cells.sample_ids[row]);
    if (!counts_all.cells.fov_ids.empty()) counts_subset.cells.fov_ids.push_back(counts_all.cells.fov_ids[row]);
    counts_subset.transcript_counts.push_back(counts_all.transcript_counts[row]);
    counts_subset.detected_genes.push_back(row < counts_all.detected_genes.size()
        ? counts_all.detected_genes[row]
        : count_detected_biological_genes(counts_all, source_row, gene_keep));
    counts_subset.crop_ids.push_back(row < counts_all.crop_ids.size() ? counts_all.crop_ids[row] : "");
    for (int p = counts_all.indptr[row]; p < counts_all.indptr[row + 1]; ++p) {
      counts_subset.indices.push_back(counts_all.indices[static_cast<std::size_t>(p)]);
      counts_subset.values.push_back(counts_all.values[static_cast<std::size_t>(p)]);
    }
    counts_subset.indptr.push_back(static_cast<int>(counts_subset.indices.size()));
  }

  emit_progress(
      progress_prefix + "Loaded sparse cell-gene counts: " +
      std::to_string(counts_subset.values.size()) + " non-zero entries");
  const auto selected_genes = select_variable_genes(
      counts_subset,
      static_cast<int>(counts_subset.genes.size()),
      options.n_variable_genes,
      options.normalization_scale);
  emit_progress(
      progress_prefix + "Selected variable genes: " + std::to_string(selected_genes.size()) +
      " genes");
  Eigen::MatrixXd features = materialize_feature_matrix(
      counts_subset,
      selected_genes,
      options.normalization_scale);
  emit_progress(
      progress_prefix + "Materialized dense clustering matrix: " +
      std::to_string(features.rows()) + " x " + std::to_string(features.cols()));
  const auto pca = run_pca(std::move(features), options.pca_dims, options.seed);
  emit_progress(
      progress_prefix + "Computed cell PCA: " + std::to_string(pca.pcs.rows()) +
      " x " + std::to_string(pca.pcs.cols()));

  CellNeighborGraph neighbor_graph = build_cell_neighbor_graph(
      pca.pcs,
      options.graph_k,
      options.compute_umap ? options.umap_neighbors : 0,
      options.num_threads,
      options.seed);
  const WeightedAdjList& graph = neighbor_graph.louvain_graph;
  emit_progress(
      progress_prefix + "Built HNSW cell KNN graph: " + std::to_string(graph.indices.size()) +
      " directed edges using " + std::to_string(std::max(1, options.num_threads)) +
      " thread(s); reusing " + std::to_string(neighbor_graph.umap_neighbor_count) +
      " cosine-distance neighbors for UMAP");
  std::vector<int> clusters = run_louvain_zero_based(
      graph,
      options.cluster_resolution,
      100);
  emit_progress(progress_prefix + "Ran Louvain clustering");
  int n_clusters = 0;
  clusters = reindex_membership_zero_based(clusters, &n_clusters);
  for (int& value : clusters) {
    value += 1;
  }

  Eigen::MatrixXd umap;
  bool has_umap = false;
  if (options.compute_umap && pca.pcs.cols() > 1 && neighbor_graph.umap_neighbor_count > 0) {
    const std::function<void(const std::string&, double)> umap_progress =
        options.progress
            ? [&](const std::string& message, double seconds) {
                options.progress(progress_prefix + message, seconds);
              }
            : std::function<void(const std::string&, double)>();
    umap = umap_embed_neighbors_impl(
        std::move(neighbor_graph.umap_neighbors),
        static_cast<int>(pca.pcs.cols()),
        2,
        options.umap_epochs,
        static_cast<int>(options.seed),
        options.num_threads,
        options.umap_parallel_optimization,
        options.progress ? &umap_progress : nullptr);
    has_umap = true;
    emit_progress(progress_prefix + "Computed cell UMAP");
  } else {
    emit_progress(progress_prefix + "Skipped cell UMAP");
  }

  CellClusteringResult out;
  out.cells = counts_subset.cells;
  out.transcript_counts = counts_subset.transcript_counts;
  out.detected_genes = counts_subset.detected_genes;
  out.crop_ids = counts_subset.crop_ids;
  out.clusters = std::move(clusters);
  out.variable_genes.reserve(selected_genes.size());
  for (const int gene : selected_genes) {
    out.variable_genes.push_back(counts_subset.genes[static_cast<std::size_t>(gene)]);
  }
  out.pca_variance_explained = pca.variance_explained;

  Eigen::MatrixXd pcs_rows(pca.pcs.cols(), pca.pcs.rows());
  for (Eigen::Index col = 0; col < pca.pcs.cols(); ++col) {
    pcs_rows.row(col) = pca.pcs.col(col).transpose();
  }
  out.pcs = eigen_to_dense_rows(pcs_rows);
  out.has_umap = has_umap;
  if (has_umap) {
    Eigen::MatrixXd umap_rows(umap.cols(), umap.rows());
    for (Eigen::Index col = 0; col < umap.cols(); ++col) {
      umap_rows.row(col) = umap.col(col).transpose();
    }
    out.umap = eigen_to_dense_rows(umap_rows);
  }
  emit_progress(progress_prefix + "Assembled cell clustering result");
  return out;
}

CellClusteringResult cluster_cell_counts(
    const CellCountMatrix& counts,
    const CellClusteringOptions& options) {
  return cluster_cell_counts(counts, options, "");
}

// Build cell features, cluster cells, and return the embedding and summaries.
CellClusteringResult cluster_cells(
    const TranscriptTable& table,
    const CellClusteringOptions& options,
    const std::vector<std::string>* transcript_crop_ids) {
  if (table.size() == 0) {
    throw std::runtime_error("cannot cluster an empty transcript table");
  }

  auto stage_start = std::chrono::steady_clock::now();
  auto emit_progress = [&](const std::string& message) {
    if (!options.progress) {
      return;
    }
    const auto now = std::chrono::steady_clock::now();
    options.progress(message, std::chrono::duration<double>(now - stage_start).count());
    stage_start = now;
  };

  const auto by_cell = table.transcripts_by_cell();
  std::vector<int> cell_sizes(table.num_cells(), 0);
  std::vector<int> detected_genes(table.num_cells(), 0);
  std::vector<int> eligible_cells;
  eligible_cells.reserve(table.num_cells());
  const auto gene_keep = biological_gene_mask(table.genes);
  std::vector<int> seen_gene(table.num_genes(), -1);
  for (std::size_t cell = 0; cell < by_cell.size(); ++cell) {
    cell_sizes[cell] = static_cast<int>(by_cell[cell].size());
    for (const int tx : by_cell[cell]) {
      const int gene = table.gene_index[static_cast<std::size_t>(tx)];
      if (gene < 0 || gene >= static_cast<int>(gene_keep.size()) ||
          !gene_keep[static_cast<std::size_t>(gene)] ||
          seen_gene[static_cast<std::size_t>(gene)] == static_cast<int>(cell)) {
        continue;
      }
      seen_gene[static_cast<std::size_t>(gene)] = static_cast<int>(cell);
      detected_genes[cell] += 1;
    }
    if (cell_sizes[cell] >= std::max(options.min_molecules, 1) &&
        detected_genes[cell] >= std::max(options.min_genes, 0)) {
      eligible_cells.push_back(static_cast<int>(cell));
    }
  }
  if (eligible_cells.empty()) {
    throw std::runtime_error("no cells meet the clustering min_molecules/min_genes thresholds");
  }
  emit_progress(
      "Indexed transcripts by cell: " + std::to_string(eligible_cells.size()) +
      " eligible cells (min_molecules=" + std::to_string(std::max(options.min_molecules, 1)) +
      ", min_genes=" + std::to_string(std::max(options.min_genes, 0)) + ")");

  const auto clustering_cells = sample_clustering_cells(
      table,
      eligible_cells,
      cell_sizes,
      by_cell,
      transcript_crop_ids,
      options.cells_max,
      options.seed);
  emit_progress(
      "Selected clustering cells: " + std::to_string(clustering_cells.size()) +
      " cells");
  const auto counts_subset = build_sparse_counts(table, clustering_cells, by_cell, transcript_crop_ids);
  emit_progress(
      "Built sparse cell-gene counts: " + std::to_string(counts_subset.values.size()) +
      " non-zero entries");
  return cluster_cell_counts(counts_subset, options, "");
}

}  // namespace celladmix
