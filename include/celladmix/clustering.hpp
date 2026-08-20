// Cell-level clustering and UMAP embedding utilities used for training-strata
// discovery and notebook visualization.

#pragma once

#include <functional>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "celladmix/matrix.hpp"
#include "celladmix/types.hpp"

namespace celladmix {

struct CellClusteringOptions {
  int min_molecules = 10;
  int min_genes = 5;
  int cells_max = -1;
  int n_variable_genes = 1000;
  int pca_dims = 30;
  int graph_k = 15;
  double cluster_resolution = 1.0;
  bool compute_umap = true;
  int umap_neighbors = 15;
  int umap_epochs = 200;
  int num_threads = 1;
  bool umap_parallel_optimization = true;
  double normalization_scale = 5000.0;
  unsigned int seed = 1;
  std::function<void(const std::string&, double)> progress;
};

// Persistable cell-level clustering outputs and derived embeddings.
struct CellClusteringResult {
  CellTable cells;
  std::vector<int> transcript_counts;
  std::vector<int> detected_genes;
  std::vector<int> clusters;
  std::vector<std::string> crop_ids;
  std::vector<std::string> variable_genes;
  std::vector<double> pca_variance_explained;
  DenseMatrix pcs;
  bool has_umap = false;
  DenseMatrix umap;
};

// Sparse cell-by-gene counts used by store-backed clustering without requiring
// transcript-level materialization.
struct CellCountMatrix {
  CellTable cells;
  std::vector<std::string> genes;
  std::vector<int> transcript_counts;
  std::vector<int> detected_genes;
  std::vector<std::string> crop_ids;
  std::vector<int> indptr;
  std::vector<int> indices;
  std::vector<double> values;
};

// Cluster cells from transcript-level data and return PCA/UMAP summaries.
CellClusteringResult cluster_cells(
    const TranscriptTable& table,
    const CellClusteringOptions& options = {},
    const std::vector<std::string>* transcript_crop_ids = nullptr);

// Cluster cells from precomputed sparse cell-by-gene counts.
CellClusteringResult cluster_cell_counts(
    const CellCountMatrix& counts,
    const CellClusteringOptions& options = {});

// Run a 2D or higher-dimensional UMAP embedding on column-major input data.
Eigen::MatrixXd umap_embed(
    const Eigen::MatrixXd& data,
    int ndim_out,
    int n_neighbors,
    int n_epochs,
    int seed,
    int num_threads = 1,
    bool parallel_optimization = false);


// Count cell types among each cell's k nearest neighbor cells (2D centroids).
// type_codes uses -1 for unlabeled cells; unlabeled cells occupy space (they
// are valid neighbors) but contribute no counts. Returns a dense cells x
// n_types matrix.
DenseMatrix cell_neighbor_type_counts(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<int>& type_codes,
    int n_types,
    int k);

// Distance from every cell to the nearest cell of each type, excluding the
// cell itself. Returns a dense cells x n_types matrix; infinity for types
// with no cells.
DenseMatrix cell_nearest_type_distance(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<int>& type_codes,
    int n_types);

}  // namespace celladmix
