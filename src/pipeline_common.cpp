// Shared helpers used by both in-memory and store-backed fit pipelines.

#include "celladmix/pipeline_common.hpp"

#include <cmath>
#include <numeric>
#include <stdexcept>

namespace celladmix {

DenseMatrix expand_h_to_full_genes(
    const DenseMatrix& compact_h,
    int full_cols,
    const std::vector<int>& kept_cols) {
  DenseMatrix full_h(compact_h.rows(), full_cols, 0.0);
  for (int factor = 0; factor < compact_h.rows(); ++factor) {
    for (std::size_t local_col = 0; local_col < kept_cols.size(); ++local_col) {
      full_h(factor, kept_cols[local_col]) =
          compact_h(factor, static_cast<int>(local_col));
    }
  }
  return full_h;
}

bool is_weighted_ls_variant(const std::string& mode) {
  return mode == "ls_nmf";
}

std::string resolve_molecule_scoring(const BasicPipelineOptions& options) {
  if (options.molecule_scoring == "auto") {
    return is_weighted_ls_variant(options.nmf_variant)
        ? "gene_loadings"
        : "ncv_projection";
  }
  if (options.molecule_scoring == "gene_loadings" ||
      options.molecule_scoring == "ncv_projection") {
    return options.molecule_scoring;
  }
  throw std::runtime_error(
      "molecule_scoring must be one of auto, gene_loadings, or ncv_projection");
}

std::pair<double, double> mean_and_sd(const std::vector<double>& values) {
  if (values.empty()) {
    return {0.0, 0.0};
  }
  const double mean = std::accumulate(values.begin(), values.end(), 0.0) /
      static_cast<double>(values.size());
  double ss = 0.0;
  for (const double value : values) {
    const double delta = value - mean;
    ss += delta * delta;
  }
  const double sd = values.size() > 1
      ? std::sqrt(ss / static_cast<double>(values.size() - 1))
      : 0.0;
  return {mean, sd};
}

SparseNmfResult weighted_to_sparse_nmf_result(const WeightedNmfResult& fit) {
  SparseNmfResult out;
  out.w = fit.w;
  out.h = fit.h;
  out.losses = fit.losses;
  out.final_objective = fit.losses.empty() ? 0.0 : fit.losses.back();
  out.selected_seed = fit.selected_seed;
  out.selected_run = fit.selected_run;
  out.candidate_final_objectives = fit.candidate_final_losses;
  out.candidate_best_match_correlations = fit.candidate_matched_correlations;
  out.selected_factor_stability = fit.selected_factor_stability;
  const auto [mean, sd] = mean_and_sd(out.candidate_final_objectives);
  out.candidate_final_objective_mean = mean;
  out.candidate_final_objective_sd = sd;
  out.candidate_best_match_correlation_mean =
      mean_and_sd(out.candidate_best_match_correlations).first;
  out.stability_comparison_runs = fit.stability_comparison_runs;
  out.stable_factor_count = fit.stable_factor_count;
  out.stability_threshold = fit.stability_threshold;
  return out;
}

}  // namespace celladmix
