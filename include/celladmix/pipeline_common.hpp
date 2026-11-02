// Shared helpers used by both in-memory and store-backed fit pipelines.

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "celladmix/matrix.hpp"
#include "celladmix/nmf_euclidean.hpp"
#include "celladmix/nmf_kl.hpp"

#include "celladmix/pipeline.hpp"

namespace celladmix {

// Expand compact factor loadings back onto the full gene universe.
DenseMatrix expand_h_to_full_genes(
    const DenseMatrix& compact_h,
    int full_cols,
    const std::vector<int>& kept_cols);

// Return true for the weighted least-squares NMF compatibility mode.
// Place profiles given over `profile_genes` (row-major, rank x genes) onto
// `target_genes`: columns are matched by gene name, genes absent from the
// target are dropped, and target genes absent from the profiles get zero
// loadings. `matched`, when supplied, receives the number of genes found.
std::vector<double> align_profiles_to_genes(
    const std::vector<double>& profiles,
    int rank,
    const std::vector<std::string>& profile_genes,
    const std::vector<std::string>& target_genes,
    int* matched = nullptr);

bool is_weighted_ls_variant(const std::string& mode);

// Resolve the automatic NCV neighborhood size for a panel and cell-size
// profile: grows with sqrt(genes_present) from the 400-gene/k=20 anchor,
// capped at half the median cell's molecule count, floored at 20.
int resolve_auto_ncv_k(int genes_present, double median_cell_molecules);

// Resolve the molecule node-potential mode requested by pipeline options.
std::string resolve_molecule_scoring(const BasicPipelineOptions& options);

// Mean and sample standard deviation for small diagnostic vectors.
std::pair<double, double> mean_and_sd(const std::vector<double>& values);

// Convert weighted LS-NMF output into the common NMF result shape.
SparseNmfResult weighted_to_sparse_nmf_result(const WeightedNmfResult& fit);

}  // namespace celladmix
