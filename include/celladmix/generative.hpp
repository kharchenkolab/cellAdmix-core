// Generative admixture model: per-target-type Poisson mixture over genes.
//
// For each annotated target cell type, the observed molecule counts of each
// cell are modeled as a mixture of the cell's own expression programs,
// contamination from each detected source type (proportional to the source's
// interface-local expression profile, with a per-cell fraction anchored to
// the measured exposure dose), an ambient background on the strict genes,
// a sparse induced-expression term on genes whose exposure-linked excess is
// disproportionate to the transfer expectation, and a small uniform floor.
// Each observed count is split among the components in proportion to their
// fitted rates; the removed mass is the contamination + ambient share.
// docs/generative.md describes the model and its validation.

#pragma once

#include <string>
#include <vector>

#include "celladmix/matrix.hpp"

namespace celladmix {

struct GenerativeOptions {
  // Contamination fraction: prior strength (pseudo-molecules), cap on the
  // per-cell update relative to the dose prior, and cap on the prior itself.
  double alpha_prior_strength = 300.0;
  double alpha_cap = 5.0;
  double lambda_max = 0.6;
  // Ambient tier (strict genes): prior strength and per-cell scale cap.
  double ambient_prior_strength = 300.0;
  double ambient_cap = 5.0;
  // Uniform floor as a fraction of cell content spread over all genes.
  double floor_total = 0.002;
  // Induced-expression screen: deviation threshold in combined-noise units,
  // minimum excess molecules, and the relative profile-uncertainty allowance.
  double induced_z = 8.0;
  double induced_min_excess = 50.0;
  double profile_cv = 0.15;
  // Per-cell induced activity: Gamma prior shape and cap.
  double rho_shape = 2.0;
  double rho_cap = 30.0;
  // Interface-local source profiles: radius (um) and the minimum molecule
  // count below which the radius widens (100 um, then all source cells).
  double near_um = 30.0;
  double near_min_molecules = 30000.0;
  // Own-program learning down-weights each cell by
  // 1 / (1 + dose_weight * prior contamination fraction).
  double dose_weight = 20.0;
  // Ambient reference: minimum lateral distance (um) from any source cell
  // and the minimum molecule count for the far-cell estimate; otherwise the
  // level falls back to half the pooled zero-exposure strict rate.
  double far_um = 250.0;
  double far_min_molecules = 20000.0;
  // Expectation-maximization schedule.
  int em_iterations = 40;
  int topup_em_iterations = 15;
  int topup_passes = 3;
  int outer_rounds = 2;
  // Exposure: number of nearest cells defining each pair's dose covariate.
  int neighbor_k = 15;
  // Expression-program initialization from factor-labeled molecules:
  // minimum molecules for a program and the minimum share of a type's
  // molecules for an unaligned factor to join its programs.
  double program_min_molecules = 5000.0;
  double program_share_min = 0.05;
  // Factor-to-type alignment derivation (used when no explicit alignment is
  // supplied): minimum cosine and required margin over the second-best type.
  double align_min_cosine = 0.3;
  double align_margin = 1.15;
  // Component switches: the production configuration has both on. Turning
  // the ambient tier off gives the exposure-permutation control its meaning
  // (an exposure-independent removal would register as power under any
  // permutation); turning the induced term off removes retention.
  bool use_ambient = true;
  bool use_induced = true;
  int num_threads = 1;
};

// One source -> target pair, from the admixture audit. Gene sets are indices
// into the count matrix rows. guide defaults to pool when empty; a held-out
// evaluation can restrict it (and set excluded_genes) without touching the
// model itself.
struct GenerativePairSpec {
  int source_type = -1;
  int target_type = -1;
  std::vector<int> pool;
  std::vector<int> strict;
  std::vector<int> guide;
  // Optional per-target-cell exposure override, aligned with the cells of
  // the target type in column order. Empty -> computed from positions.
  std::vector<double> exposure;
};

struct GenerativeInducedGene {
  int pair = 0;
  int gene = 0;
  double excess = 0.0;
  double expected = 0.0;
  double z = 0.0;
};

struct GenerativePairSummary {
  int pair = 0;
  double prior_molecules = 0.0;
  double posterior_molecules = 0.0;
  double induced_molecules = 0.0;
  double mean_dose_exposed = 0.0;
};

struct GenerativeResult {
  // Removed molecule mass on the pattern of the input counts: value p of
  // the result aligns with value p of the input matrix.
  std::vector<double> removed;
  // Per pair: the target cells (columns of the count matrix), the prior
  // contamination fraction from the dose, the fitted per-cell fraction, and
  // the per-cell induced activity.
  std::vector<std::vector<int>> pair_cells;
  std::vector<std::vector<double>> pair_dose;
  std::vector<std::vector<double>> pair_alpha;
  std::vector<std::vector<double>> pair_rho;
  // Per cell: fitted ambient scale (zero for cells of types without one).
  std::vector<double> ambient_scale;
  std::vector<GenerativeInducedGene> induced;
  std::vector<GenerativePairSummary> pairs;
  // Factor-to-type alignment used for the program initializations
  // (factor index -> type code, -1 unaligned).
  std::vector<int> factor_alignment;
  // True when the run carries no nucleus-overlap flag and whole-cell
  // profiles were used in place of cytoplasmic ones.
  bool whole_cell_profiles = false;
};

// Fit the model on a fitted run. counts is the gene x cell matrix in
// compressed-sparse-column form (indptr over cells); x, y, type_codes are
// aligned with its columns (type_codes use -1 for unannotated cells);
// cell_ids are the column names used to match the run's cell table.
// molecules_parquet may be empty, in which case expression programs are
// initialized from type pseudobulk profiles and transfer profiles from
// whole-cell counts (used by unit tests and runs without factor labels).
// factor_to_type supplies an explicit factor alignment (factor index ->
// type code, -1 unaligned); empty derives it from the data.
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
    const GenerativeOptions& options);

// Weighted pool-adjacent-violators: monotone non-decreasing fit of y with
// weights w. Exposed for tests.
std::vector<double> pava_increasing(
    const std::vector<double>& y, const std::vector<double>& w);

}  // namespace celladmix
