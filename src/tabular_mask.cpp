// TIFF segmentation-mask loading and molecule-to-mask assignment.

#include "celladmix/tabular_mask.hpp"

#include <tiffio.h>

#include <cmath>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <unordered_set>

namespace celladmix {

TabularTiffMask read_labeled_tiff_mask(const std::string& path) {
  TIFFSetWarningHandler(nullptr);
  TIFF* tif = TIFFOpen(path.c_str(), "r");
  if (tif == nullptr) {
    throw std::runtime_error("Could not open TIFF segmentation mask: " + path);
  }

  std::uint32_t width = 0;
  std::uint32_t height = 0;
  std::uint16_t bits_per_sample = 8;
  std::uint16_t samples_per_pixel = 1;
  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
  TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
  TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);
  if (samples_per_pixel != 1) {
    TIFFClose(tif);
    throw std::runtime_error("Only single-channel labeled TIFF masks are supported");
  }

  TabularTiffMask mask;
  mask.width = width;
  mask.height = height;
  mask.labels.assign(static_cast<std::size_t>(width) * static_cast<std::size_t>(height), 0);
  bool saw_nonzero = false;
  std::uint32_t first_nonzero = 0;
  bool has_multiple_nonzero_values = false;

  const tmsize_t scanline_size = TIFFScanlineSize(tif);
  std::vector<std::uint8_t> row_buf(static_cast<std::size_t>(scanline_size));

  for (std::uint32_t row = 0; row < height; ++row) {
    if (TIFFReadScanline(tif, row_buf.data(), row, 0) < 0) {
      TIFFClose(tif);
      throw std::runtime_error("Error reading TIFF scanline " + std::to_string(row));
    }
    auto* dst = mask.labels.data() + static_cast<std::size_t>(row) * static_cast<std::size_t>(width);
    if (bits_per_sample == 8) {
      for (std::uint32_t col = 0; col < width; ++col) {
        dst[col] = row_buf[static_cast<std::size_t>(col)];
        if (dst[col] > 0) {
          if (!saw_nonzero) {
            first_nonzero = dst[col];
            saw_nonzero = true;
          } else if (dst[col] != first_nonzero) {
            has_multiple_nonzero_values = true;
          }
        }
      }
    } else if (bits_per_sample == 16) {
      const auto* src = reinterpret_cast<const std::uint16_t*>(row_buf.data());
      for (std::uint32_t col = 0; col < width; ++col) {
        dst[col] = src[col];
        if (dst[col] > 0) {
          if (!saw_nonzero) {
            first_nonzero = dst[col];
            saw_nonzero = true;
          } else if (dst[col] != first_nonzero) {
            has_multiple_nonzero_values = true;
          }
        }
      }
    } else if (bits_per_sample == 32) {
      const auto* src = reinterpret_cast<const std::uint32_t*>(row_buf.data());
      for (std::uint32_t col = 0; col < width; ++col) {
        dst[col] = src[col];
        if (dst[col] > 0) {
          if (!saw_nonzero) {
            first_nonzero = dst[col];
            saw_nonzero = true;
          } else if (dst[col] != first_nonzero) {
            has_multiple_nonzero_values = true;
          }
        }
      }
    } else {
      TIFFClose(tif);
      throw std::runtime_error(
          "Unsupported TIFF bits/sample for labeled mask: " + std::to_string(bits_per_sample));
    }
  }

  TIFFClose(tif);
  mask.is_binary = saw_nonzero && !has_multiple_nonzero_values;

  if (mask.is_binary) {
    std::cerr << "[INFO] Binary TIFF mask detected (single nonzero value "
              << first_nonzero << "); running connected-component labeling\n";
    std::vector<std::uint8_t> visited(mask.labels.size(), 0);
    std::queue<std::size_t> queue;
    std::uint32_t next_label = 1;

    auto enqueue_if_foreground = [&](int row, int col) {
      if (row < 0 || col < 0 ||
          row >= static_cast<int>(mask.height) ||
          col >= static_cast<int>(mask.width)) {
        return;
      }
      const std::size_t idx = static_cast<std::size_t>(row) * static_cast<std::size_t>(mask.width) +
                              static_cast<std::size_t>(col);
      if (visited[idx] || mask.labels[idx] == 0) {
        return;
      }
      visited[idx] = 1;
      queue.push(idx);
    };

    for (std::size_t start = 0; start < mask.labels.size(); ++start) {
      if (mask.labels[start] == 0 || visited[start]) {
        continue;
      }
      visited[start] = 1;
      queue.push(start);
      while (!queue.empty()) {
        const std::size_t idx = queue.front();
        queue.pop();
        mask.labels[idx] = next_label;
        const int row = static_cast<int>(idx / static_cast<std::size_t>(mask.width));
        const int col = static_cast<int>(idx % static_cast<std::size_t>(mask.width));
        enqueue_if_foreground(row - 1, col);
        enqueue_if_foreground(row + 1, col);
        enqueue_if_foreground(row, col - 1);
        enqueue_if_foreground(row, col + 1);
      }
      ++next_label;
    }
    mask.n_nonzero_labels = static_cast<std::size_t>(next_label - 1);
    std::cerr << "[INFO] Connected-component labeling produced "
              << mask.n_nonzero_labels << " mask labels\n";
  } else {
    std::unordered_set<std::uint32_t> labels_seen;
    for (const std::uint32_t value : mask.labels) {
      if (value > 0) {
        labels_seen.insert(value);
      }
    }
    mask.n_nonzero_labels = labels_seen.size();
    std::cerr << "[INFO] Multi-label TIFF mask detected with "
              << mask.n_nonzero_labels << " distinct nonzero labels\n";
  }
  return mask;
}

std::string assign_cell_id_from_mask(double x, double y, const TabularTiffMask& mask) {
  const int col = static_cast<int>(std::llround(x)) - 1;
  const int row = static_cast<int>(std::llround(y)) - 1;
  if (row < 0 || col < 0 ||
      row >= static_cast<int>(mask.height) ||
      col >= static_cast<int>(mask.width)) {
    return {};
  }
  const std::uint32_t label = mask.labels[static_cast<std::size_t>(row) * static_cast<std::size_t>(mask.width) +
                                         static_cast<std::size_t>(col)];
  return label > 0 ? std::to_string(label) : std::string();
}

std::vector<std::string> assign_cell_ids_from_mask(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const TabularTiffMask& mask) {
  std::vector<std::string> out(x.size(), "");
  for (std::size_t i = 0; i < x.size(); ++i) {
    out[i] = assign_cell_id_from_mask(x[i], y[i], mask);
  }
  return out;
}

}  // namespace celladmix
