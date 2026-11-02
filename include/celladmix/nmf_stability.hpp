// Multirun factor-stability diagnostics shared by the KL and weighted-Euclidean
// NMF solvers.
//
// Stability is measured on gene-ownership profiles: each gene's loadings are
// rescaled to sum to one across factors, so a factor is described by which
// genes it owns rather than by raw loading magnitudes. This cancels gene
// abundance per gene, which makes the null similarity between unrelated
// factors close to zero and the values comparable across NMF variants and
// ranks. Factors of the selected run are paired one-to-one with factors of
// each comparison run by a maximum-similarity (Hungarian) assignment, and the
// per-factor stability is the mean matched Pearson correlation across
// comparison runs.

#pragma once

#include <vector>

#include "celladmix/matrix.hpp"

namespace celladmix {

// Threshold above which a factor counts as reproducibly re-found.
inline constexpr double kStableFactorThreshold = 0.3;

// Per-gene ownership shares: each column of H rescaled to sum to one across
// factors. Columns with no loading stay zero.
DenseMatrix ownership_profiles(const DenseMatrix& h);

// Pearson correlations between every row of lhs and every row of rhs.
// Rows must have equal width; zero-variance rows correlate at zero.
DenseMatrix row_correlation_matrix(const DenseMatrix& lhs, const DenseMatrix& rhs);

// Maximum-total-similarity one-to-one assignment on a square matrix.
// Returns match such that row i is paired with column match[i].
std::vector<int> hungarian_match(const DenseMatrix& similarity);

struct FactorMatchResult {
  // For each row of `profiles`: the matched row of `reference` (-1 when the
  // reference has fewer rows) and the ownership correlation of the pair.
  std::vector<int> match;
  std::vector<double> correlation;
};

// One-to-one matching of the factors of `profiles` to those of `reference`
// by gene-ownership correlation, the restart-stability metric. Rectangular
// inputs are handled by padding the assignment with zero-similarity slots.
FactorMatchResult match_factor_profiles(
    const DenseMatrix& profiles,
    const DenseMatrix& reference);

struct FactorClusterResult {
  // Cluster index per input row (0-based), numbered in order of first
  // appearance along the rows.
  std::vector<int> cluster;
  // Mean ownership correlation of each row with the other rows of its
  // cluster; NaN for a row alone in its cluster.
  std::vector<double> within_correlation;
  int n_clusters = 0;
};

// Average-linkage agglomerative clustering of pooled factor profiles (one row
// per factor of one factorization) on gene-ownership correlation distance
// (1 - correlation). `group` names, per row, the factorization the row comes
// from; ownership shares are normalized within each group, so a gene's
// loading is divided among the factors of one factorization. Clusters merge
// while their average pairwise distance stays at or below
// 1 - min_correlation, and never merge two clusters that both hold a factor
// of the same factorization: one factorization finds a program once, so a
// cluster is the one-to-one correspondence of a program across
// factorizations. Deterministic: ties resolve by row order.
FactorClusterResult cluster_factor_profiles(
    const DenseMatrix& profiles,
    const std::vector<int>& group,
    double min_correlation);

struct FactorStabilityResult {
  // Mean matched ownership correlation per selected-run factor, averaged over
  // the comparison runs. Filled with 1.0 when there are no comparison runs.
  std::vector<double> factor_stability;
  // Mean matched ownership correlation per candidate run against the selected
  // run (1.0 for the selected run itself).
  std::vector<double> run_matched_means;
  // Number of runs the per-factor stability was averaged over.
  int comparison_runs = 0;
  // Number of factors at or above `threshold`.
  int stable_factor_count = 0;
  double threshold = kStableFactorThreshold;
};

// Compute matched ownership stability of the selected run's factors.
// `comparison_runs` lists the run indices to average the per-factor stability
// over (typically all independent random restarts except the selected run);
// `run_matched_means` is still reported for every candidate run.
FactorStabilityResult matched_ownership_stability(
    const std::vector<const DenseMatrix*>& candidate_h,
    int selected_run,
    const std::vector<int>& comparison_runs,
    double stable_threshold = kStableFactorThreshold);

}  // namespace celladmix
