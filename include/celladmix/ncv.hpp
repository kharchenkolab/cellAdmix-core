// Neighborhood composition vector construction and fixed-H projection helpers.

#pragma once

#include <string>
#include <vector>

#include "celladmix/matrix.hpp"
#include "celladmix/types.hpp"

namespace celladmix {

struct NcvOptions {
  int k = 20;
  bool include_self = true;
  bool within_cell = true;
  std::vector<int> query_indices;
};

struct NcvFeatureTransform {
  std::string mode = "kl";
  std::vector<double> gene_weights;
  double target_row_sum = 0.0;
};

// Build dense NCV rows for the requested transcript queries.
DenseMatrix build_ncv_matrix(const TranscriptTable& table, const NcvOptions& options = {});
// Build sparse NCV rows for the requested transcript queries.
SparseRowMatrix build_sparse_ncv_matrix(const TranscriptTable& table, const NcvOptions& options = {});

// Build and apply the NCV feature transform used before KL-NMF.
NcvFeatureTransform make_ncv_feature_transform(
    const SparseRowMatrix& x,
    const std::string& mode);
SparseRowMatrix transform_sparse_ncv_matrix(
    const SparseRowMatrix& x,
    const NcvFeatureTransform& transform);
NcvFeatureTransform expand_ncv_feature_transform(
    const NcvFeatureTransform& compact,
    int full_cols,
    const std::vector<int>& kept_cols);

// Project NCVs to factor space using simple normalized dot products.
DenseMatrix project_ncv_to_factors(
    const TranscriptTable& table,
    const DenseMatrix& h,
    const NcvOptions& options = {},
    int num_threads = 1);

// Assign molecule-level factor potentials directly from gene loadings.
DenseMatrix project_gene_loadings_to_factors(
    const TranscriptTable& table,
    const DenseMatrix& h,
    const NcvOptions& options = {},
    int num_threads = 1);

// Project NCVs to factor space with fixed-H KL multiplicative updates.
DenseMatrix project_ncv_to_factors_kl(
    const TranscriptTable& table,
    const DenseMatrix& h,
    const NcvOptions& options = {},
    int num_threads = 1,
    int max_iterations = 10,
    double tolerance = 1e-4,
    double update_epsilon = 1e-10,
    NcvFeatureTransform transform = {});

// Reproject already-built sparse NCV rows with the fixed-H KL solver.
DenseMatrix project_rows_to_factors_kl(
    const SparseRowMatrix& x,
    const DenseMatrix& h,
    int num_threads = 1,
    int max_iterations = 10,
    double tolerance = 1e-4,
    double update_epsilon = 1e-10,
    NcvFeatureTransform transform = {});

}  // namespace celladmix
