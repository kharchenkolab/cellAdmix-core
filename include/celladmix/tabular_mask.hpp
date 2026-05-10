// TIFF segmentation-mask helpers shared by tabular loading and input-store
// streaming.

#pragma once

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

namespace celladmix {

struct TabularTiffMask {
  std::uint32_t width = 0;
  std::uint32_t height = 0;
  std::vector<std::uint32_t> labels;
  bool is_binary = false;
  std::size_t n_nonzero_labels = 0;
};

TabularTiffMask read_labeled_tiff_mask(const std::string& path);

std::string assign_cell_id_from_mask(double x, double y, const TabularTiffMask& mask);

std::vector<std::string> assign_cell_ids_from_mask(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const TabularTiffMask& mask);

}  // namespace celladmix
