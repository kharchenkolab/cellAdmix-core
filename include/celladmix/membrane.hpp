// Membrane-image utilities and stain-based factor scoring.

#pragma once

#include <functional>
#include <utility>
#include <string>
#include <vector>

#include "celladmix/types.hpp"

namespace celladmix {

struct Image2D {
  int width = 0;
  int height = 0;
  double pixel_size = 1.0;
  double x_origin = 0.0;
  double y_origin = 0.0;
  std::vector<double> values;

  // Sample the image at physical coordinates using bilinear interpolation.
  double sample_bilinear(double x, double y) const;
};

// Robustly rescale image values to [0, 1] using low/high quantiles.
std::pair<double, double> normalize_image_quantile(
    Image2D& image,
    double low_quantile = 0.05,
    double high_quantile = 0.995);

// Sample a line segment and return the maximum membrane intensity encountered.
double line_sample_max(
    const Image2D& image,
    double x0,
    double y0,
    double x1,
    double y1,
    double step = 0.5);

// Compare membrane signal along real and control directions from one cell.
double membrane_score(
    const Image2D& image,
    const std::vector<std::pair<double, double>>& real_points,
    const std::vector<std::pair<double, double>>& control_points,
    std::pair<double, double> self_centroid,
    double step = 0.5);

// Configuration for reading a single-channel Xenium stain image. Coordinates
// are interpreted in the same physical units as molecule coordinates.
struct MembraneImageOptions {
  std::string image_path;
  double pixel_size = 1.0;
  double x_offset = 0.0;
  double y_offset = 0.0;
  double epsilon = 1.0;
};

// Controls stain-based scoring of factor-assigned molecules in target cells.
struct MembraneTestOptions {
  int cell_candidate_k = 50;
  int candidate_pairs_per_type_pair = 400;
  double cell_candidate_halo = -1.0;
  int min_factor_molecules = 5;
  int min_pairs = 5;
  int max_cells_per_type_pair = 200;
  double control_distance_fraction = 0.025;
  int line_samples = 16;
  int num_threads = 1;
  unsigned int seed = 1U;
  std::function<void(const std::string&, double)> progress;
};

// Per-factor membrane evidence for one ordered target/source cell pair.
struct MembranePairScore {
  int target_cell = -1;
  int source_cell = -1;
  int target_type = -1;
  int source_type = -1;
  int factor = -1;
  int factor_count = 0;
  int scored_molecules = 0;
  double mean_score = 0.0;
  double fraction_positive = 0.0;
  double mean_directional_weight = 0.0;
  bool used_in_summary = false;
};

// Summary of membrane evidence for one target/source cell-type pair and factor.
struct MembraneSummary {
  int target_type = -1;
  int source_type = -1;
  int factor = -1;
  int n_pairs = 0;
  double mean_score = 0.0;
  double q75_score = 0.0;
  double fraction_positive_pairs = 0.0;
  double p_value = 1.0;
  double neg_log10_p = 0.0;
};

struct MembraneTestResult {
  std::vector<std::string> cell_types;
  std::vector<MembranePairScore> pair_scores;
  std::vector<MembraneSummary> summaries;
};

// Read a physical-coordinate stain-image window, downsampling if needed.
Image2D read_stain_image_crop(
    const MembraneImageOptions& image_options,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    int max_pixels_per_axis = 512);

// Run membrane-stain scoring over fitted molecule factor labels.
MembraneTestResult run_membrane_test(
    const TranscriptTable& table,
    const CellTable& cells,
    const std::vector<int>& labels,
    const MembraneImageOptions& image_options,
    const MembraneTestOptions& options = {});

}  // namespace celladmix
