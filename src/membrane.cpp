#include "celladmix/membrane.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <atomic>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include <openjpeg.h>
#include <tiffio.h>

#include "celladmix/bridge.hpp"
#include "celladmix/spatial.hpp"
#include "subpar/range.hpp"

namespace celladmix {

// Image sampling and membrane-enrichment helpers for simulation and diagnostics.

// Sample the membrane image at continuous coordinates with bilinear interpolation.
double Image2D::sample_bilinear(double x, double y) const {
  if (width <= 0 || height <= 0 || values.size() != static_cast<std::size_t>(width * height)) {
    throw std::runtime_error("Image2D is not initialized");
  }
  const double px = (x - x_origin) / pixel_size;
  const double py = (y - y_origin) / pixel_size;
  const double clamped_x = std::clamp(px, 0.0, static_cast<double>(width - 1));
  const double clamped_y = std::clamp(py, 0.0, static_cast<double>(height - 1));

  const int x0 = static_cast<int>(std::floor(clamped_x));
  const int y0 = static_cast<int>(std::floor(clamped_y));
  const int x1 = std::min(x0 + 1, width - 1);
  const int y1 = std::min(y0 + 1, height - 1);
  const double tx = clamped_x - static_cast<double>(x0);
  const double ty = clamped_y - static_cast<double>(y0);

  auto at = [&](int xi, int yi) {
    return values[static_cast<std::size_t>(yi * width + xi)];
  };

  const double v00 = at(x0, y0);
  const double v10 = at(x1, y0);
  const double v01 = at(x0, y1);
  const double v11 = at(x1, y1);

  const double v0 = (1.0 - tx) * v00 + tx * v10;
  const double v1 = (1.0 - tx) * v01 + tx * v11;
  return (1.0 - ty) * v0 + ty * v1;
}

std::pair<double, double> normalize_image_quantile(
    Image2D& image,
    double low_quantile,
    double high_quantile) {
  if (image.values.empty()) {
    return {0.0, 1.0};
  }
  if (!(low_quantile >= 0.0 && low_quantile <= 1.0 &&
        high_quantile >= 0.0 && high_quantile <= 1.0 &&
        high_quantile > low_quantile)) {
    throw std::runtime_error("image normalization quantiles must satisfy 0 <= low < high <= 1");
  }

  std::vector<double> sorted = image.values;
  std::sort(sorted.begin(), sorted.end());
  auto quantile = [&](double q) {
    const double pos = q * static_cast<double>(sorted.size() - 1U);
    const auto lo = static_cast<std::size_t>(std::floor(pos));
    const auto hi = static_cast<std::size_t>(std::ceil(pos));
    if (lo == hi) {
      return sorted[lo];
    }
    const double frac = pos - static_cast<double>(lo);
    return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
  };
  const double low = quantile(low_quantile);
  const double high = quantile(high_quantile);
  const double denom = std::max(high - low, 1e-12);
  for (double& value : image.values) {
    value = std::clamp((value - low) / denom, 0.0, 1.0);
  }
  return {low, high};
}

// Evaluate the strongest membrane signal encountered along a line segment.
double line_sample_max(
    const Image2D& image,
    double x0,
    double y0,
    double x1,
    double y1,
    double step) {
  const double dx = x1 - x0;
  const double dy = y1 - y0;
  const double length = std::sqrt(dx * dx + dy * dy);
  if (length < 1e-12) {
    return image.sample_bilinear(x0, y0);
  }
  const int n_steps = std::max(1, static_cast<int>(std::ceil(length / std::max(step, 1e-6))));
  double best = 0.0;
  for (int i = 0; i <= n_steps; ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(n_steps);
    const double x = x0 + t * dx;
    const double y = y0 + t * dy;
    best = std::max(best, image.sample_bilinear(x, y));
  }
  return best;
}

// Compare real boundary rays against control rays to score membrane enrichment.
double membrane_score(
    const Image2D& image,
    const std::vector<std::pair<double, double>>& real_points,
    const std::vector<std::pair<double, double>>& control_points,
    std::pair<double, double> self_centroid,
    double step) {
  if (real_points.size() != control_points.size()) {
    throw std::runtime_error("real_points and control_points must have the same length");
  }
  if (real_points.empty()) {
    return 0.0;
  }

  double total = 0.0;
  for (std::size_t i = 0; i < real_points.size(); ++i) {
    const double real_value = line_sample_max(
        image,
        real_points[i].first,
        real_points[i].second,
        self_centroid.first,
        self_centroid.second,
        step);
    const double control_value = line_sample_max(
        image,
        control_points[i].first,
        control_points[i].second,
        self_centroid.first,
        self_centroid.second,
        step);
    total += std::log((real_value + 1e-6) / (control_value + 1e-6));
  }
  return total / static_cast<double>(real_points.size());
}

namespace {

int effective_threads(int requested, int num_tasks) {
  if (num_tasks <= 0) {
    return 1;
  }
  return std::max(1, subpar::sanitize_num_workers(requested, num_tasks));
}

double elapsed_seconds(const std::chrono::steady_clock::time_point& start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

void emit_membrane_progress(
    const MembraneTestOptions& options,
    const std::string& message,
    const std::chrono::steady_clock::time_point& start) {
  if (options.progress) {
    options.progress(message, elapsed_seconds(start));
  }
}

std::uint64_t pair_key(int first, int second) {
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(first)) << 32) |
      static_cast<std::uint32_t>(second);
}

double mean_value(const std::vector<double>& values) {
  if (values.empty()) {
    return 0.0;
  }
  return std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
}

double quantile_75(std::vector<double> values) {
  if (values.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  std::sort(values.begin(), values.end());
  const double pos = 0.75 * static_cast<double>(values.size() - 1U);
  const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
  const std::size_t hi = static_cast<std::size_t>(std::ceil(pos));
  if (lo == hi) {
    return values[lo];
  }
  const double frac = pos - static_cast<double>(lo);
  return values[lo] * (1.0 - frac) + values[hi] * frac;
}

double one_sample_wilcoxon_greater(const std::vector<double>& values) {
  std::vector<std::pair<double, double>> abs_and_sign;
  abs_and_sign.reserve(values.size());
  for (const double value : values) {
    if (std::abs(value) <= 1e-12 || !std::isfinite(value)) {
      continue;
    }
    abs_and_sign.emplace_back(std::abs(value), value > 0.0 ? 1.0 : -1.0);
  }
  if (abs_and_sign.empty()) {
    return 1.0;
  }
  std::sort(abs_and_sign.begin(), abs_and_sign.end(), [](const auto& lhs, const auto& rhs) {
    return lhs.first < rhs.first;
  });

  double w_plus = 0.0;
  std::size_t i = 0;
  while (i < abs_and_sign.size()) {
    std::size_t j = i + 1U;
    while (j < abs_and_sign.size() && abs_and_sign[j].first == abs_and_sign[i].first) {
      ++j;
    }
    const double rank = (static_cast<double>(i + 1U) + static_cast<double>(j)) / 2.0;
    for (std::size_t k = i; k < j; ++k) {
      if (abs_and_sign[k].second > 0.0) {
        w_plus += rank;
      }
    }
    i = j;
  }

  const double n = static_cast<double>(abs_and_sign.size());
  const double mean = n * (n + 1.0) / 4.0;
  const double variance = n * (n + 1.0) * (2.0 * n + 1.0) / 24.0;
  if (variance <= 0.0) {
    return 1.0;
  }
  const double z = (w_plus - mean - 0.5) / std::sqrt(variance);
  return std::min(1.0, std::max(0.0, 0.5 * std::erfc(z / std::sqrt(2.0))));
}

std::vector<std::string> resolve_cell_types(
    const TranscriptTable& table,
    const CellTable& cells) {
  std::vector<std::string> out(table.num_cells(), "");

  if (cells.cell_types.size() == cells.cell_ids.size() && !cells.cell_types.empty()) {
    std::unordered_map<std::string, std::string> by_id;
    by_id.reserve(cells.cell_ids.size());
    for (std::size_t i = 0; i < cells.cell_ids.size(); ++i) {
      by_id.emplace(cells.cell_ids[i], cells.cell_types[i]);
    }
    for (std::size_t cell = 0; cell < table.cells.size(); ++cell) {
      const auto it = by_id.find(table.cells[cell]);
      if (it != by_id.end()) {
        out[cell] = it->second;
      }
    }
  }

  if (std::any_of(out.begin(), out.end(), [](const std::string& value) { return value.empty(); }) &&
      table.cell_types.size() == table.size()) {
    for (std::size_t i = 0; i < table.size(); ++i) {
      auto& value = out[static_cast<std::size_t>(table.cell_index[i])];
      if (value.empty()) {
        value = table.cell_types[i];
      }
    }
  }

  return out;
}

double cell_radius_xy(
    const TranscriptTable& table,
    const std::vector<int>& members,
    double cx,
    double cy) {
  double out = 0.0;
  for (const int row : members) {
    const double dx = table.x[static_cast<std::size_t>(row)] - cx;
    const double dy = table.y[static_cast<std::size_t>(row)] - cy;
    out = std::max(out, std::sqrt(dx * dx + dy * dy));
  }
  return out;
}

struct MembranePreparedData {
  std::vector<std::string> type_names;
  std::vector<int> cell_type_index;
  std::vector<std::vector<int>> members_by_cell;
  std::vector<int> factor_counts;
  std::vector<double> factor_fractions;
  std::vector<double> centroid_x;
  std::vector<double> centroid_y;
  std::vector<double> centroid_z;
  int n_factors = 0;
};

MembranePreparedData prepare_membrane_data(
    const TranscriptTable& table,
    const CellTable& cells,
    const std::vector<int>& labels) {
  if (table.cell_index.empty()) {
    throw std::runtime_error("TranscriptTable must provide cell indices before membrane scoring");
  }
  if (labels.size() != table.size()) {
    throw std::runtime_error("labels must match TranscriptTable size");
  }

  MembranePreparedData data;
  const auto cell_types = resolve_cell_types(table, cells);
  data.cell_type_index.assign(table.num_cells(), -1);
  std::unordered_map<std::string, int> type_lookup;
  for (std::size_t cell = 0; cell < cell_types.size(); ++cell) {
    const auto& type = cell_types[cell];
    if (type.empty()) {
      continue;
    }
    auto it = type_lookup.find(type);
    if (it == type_lookup.end()) {
      const int idx = static_cast<int>(data.type_names.size());
      it = type_lookup.emplace(type, idx).first;
      data.type_names.push_back(type);
    }
    data.cell_type_index[cell] = it->second;
  }

  data.n_factors = labels.empty() ? 0 : (*std::max_element(labels.begin(), labels.end()) + 1);
  if (data.n_factors <= 0) {
    throw std::runtime_error("Membrane scoring requires non-negative factor labels");
  }

  data.members_by_cell = table.transcripts_by_cell();
  const int n_cells = static_cast<int>(table.num_cells());
  data.factor_counts.assign(static_cast<std::size_t>(n_cells * data.n_factors), 0);
  data.factor_fractions.assign(static_cast<std::size_t>(n_cells * data.n_factors), 0.0);
  for (std::size_t i = 0; i < table.size(); ++i) {
    const int factor = labels[i];
    if (factor < 0 || factor >= data.n_factors) {
      continue;
    }
    const int cell = table.cell_index[i];
    data.factor_counts[static_cast<std::size_t>(cell * data.n_factors + factor)] += 1;
  }
  for (int cell = 0; cell < n_cells; ++cell) {
    const double denom = static_cast<double>(std::max<std::size_t>(
        data.members_by_cell[static_cast<std::size_t>(cell)].size(), static_cast<std::size_t>(1)));
    for (int factor = 0; factor < data.n_factors; ++factor) {
      data.factor_fractions[static_cast<std::size_t>(cell * data.n_factors + factor)] =
          static_cast<double>(data.factor_counts[static_cast<std::size_t>(cell * data.n_factors + factor)]) /
          denom;
    }
  }
  if (cells.centroid_x.size() == table.num_cells() && cells.centroid_y.size() == table.num_cells()) {
    data.centroid_x = cells.centroid_x;
    data.centroid_y = cells.centroid_y;
    data.centroid_z = cells.centroid_z.size() == table.num_cells()
        ? cells.centroid_z
        : std::vector<double>(table.num_cells(), 0.0);
  } else {
    data.centroid_x = infer_cell_centroids_x(table);
    data.centroid_y = infer_cell_centroids_y(table);
    data.centroid_z = infer_cell_centroids_z(table);
  }
  return data;
}

struct MemoryStream {
  const unsigned char* data = nullptr;
  OPJ_SIZE_T size = 0;
  OPJ_SIZE_T offset = 0;
};

OPJ_SIZE_T memory_stream_read(void* buffer, OPJ_SIZE_T n_bytes, void* user_data) {
  auto* stream = static_cast<MemoryStream*>(user_data);
  if (stream->offset >= stream->size) {
    return static_cast<OPJ_SIZE_T>(-1);
  }
  const OPJ_SIZE_T remaining = stream->size - stream->offset;
  const OPJ_SIZE_T to_read = std::min(n_bytes, remaining);
  std::memcpy(buffer, stream->data + stream->offset, static_cast<std::size_t>(to_read));
  stream->offset += to_read;
  return to_read;
}

OPJ_OFF_T memory_stream_skip(OPJ_OFF_T n_bytes, void* user_data) {
  auto* stream = static_cast<MemoryStream*>(user_data);
  if (n_bytes < 0) {
    return -1;
  }
  const OPJ_SIZE_T skip = std::min(static_cast<OPJ_SIZE_T>(n_bytes), stream->size - stream->offset);
  stream->offset += skip;
  return static_cast<OPJ_OFF_T>(skip);
}

OPJ_BOOL memory_stream_seek(OPJ_OFF_T position, void* user_data) {
  auto* stream = static_cast<MemoryStream*>(user_data);
  if (position < 0 || static_cast<OPJ_SIZE_T>(position) > stream->size) {
    return OPJ_FALSE;
  }
  stream->offset = static_cast<OPJ_SIZE_T>(position);
  return OPJ_TRUE;
}

void openjpeg_noop_callback(const char*, void*) {}

struct TilePixels {
  int width = 0;
  int height = 0;
  std::vector<float> values;
};

TilePixels decode_jpeg2000_tile_once(
    const std::vector<unsigned char>& raw,
    OPJ_CODEC_FORMAT format) {
  MemoryStream mem{raw.data(), static_cast<OPJ_SIZE_T>(raw.size()), 0};
  opj_stream_t* stream = opj_stream_create(4096, OPJ_STREAM_READ);
  if (stream == nullptr) {
    throw std::runtime_error("Failed to create OpenJPEG stream");
  }
  opj_stream_set_user_data(stream, &mem, nullptr);
  opj_stream_set_user_data_length(stream, static_cast<OPJ_UINT64>(raw.size()));
  opj_stream_set_read_function(stream, memory_stream_read);
  opj_stream_set_skip_function(stream, memory_stream_skip);
  opj_stream_set_seek_function(stream, memory_stream_seek);

  opj_dparameters_t params;
  opj_set_default_decoder_parameters(&params);
  opj_codec_t* codec = opj_create_decompress(format);
  if (codec == nullptr) {
    opj_stream_destroy(stream);
    throw std::runtime_error("Failed to create OpenJPEG decoder");
  }
  opj_set_error_handler(codec, openjpeg_noop_callback, nullptr);
  opj_set_warning_handler(codec, openjpeg_noop_callback, nullptr);
  opj_set_info_handler(codec, openjpeg_noop_callback, nullptr);

  opj_image_t* image = nullptr;
  bool ok = opj_setup_decoder(codec, &params) == OPJ_TRUE &&
      opj_read_header(stream, codec, &image) == OPJ_TRUE &&
      opj_decode(codec, stream, image) == OPJ_TRUE &&
      opj_end_decompress(codec, stream) == OPJ_TRUE;
  TilePixels out;
  if (ok && image != nullptr && image->numcomps > 0 && image->comps[0].data != nullptr) {
    const auto& comp = image->comps[0];
    out.width = static_cast<int>(comp.w);
    out.height = static_cast<int>(comp.h);
    out.values.resize(static_cast<std::size_t>(out.width) * static_cast<std::size_t>(out.height));
    const int max_value = comp.prec >= 31U ? std::numeric_limits<int>::max() : ((1 << comp.prec) - 1);
    for (std::size_t i = 0; i < out.values.size(); ++i) {
      const int value = comp.data[i];
      out.values[i] = static_cast<float>(std::max(0, std::min(max_value, value)));
    }
  }
  if (image != nullptr) {
    opj_image_destroy(image);
  }
  opj_destroy_codec(codec);
  opj_stream_destroy(stream);
  if (!ok || out.values.empty()) {
    throw std::runtime_error("OpenJPEG failed to decode tile");
  }
  return out;
}

TilePixels decode_jpeg2000_tile(const std::vector<unsigned char>& raw) {
  const bool looks_jp2 = raw.size() >= 12U &&
      raw[4] == 'j' && raw[5] == 'P' && raw[6] == ' ' && raw[7] == ' ';
  try {
    return decode_jpeg2000_tile_once(raw, looks_jp2 ? OPJ_CODEC_JP2 : OPJ_CODEC_J2K);
  } catch (const std::exception&) {
    return decode_jpeg2000_tile_once(raw, looks_jp2 ? OPJ_CODEC_J2K : OPJ_CODEC_JP2);
  }
}

class TiffStainImage {
 public:
  explicit TiffStainImage(MembraneImageOptions options) : options_(std::move(options)) {
    if (options_.image_path.empty()) {
      throw std::runtime_error("Membrane image path is empty");
    }
    TIFFSetWarningHandler(nullptr);
    tif_ = TIFFOpen(options_.image_path.c_str(), "r");
    if (tif_ == nullptr) {
      throw std::runtime_error("Could not open membrane TIFF image: " + options_.image_path);
    }
    TIFFGetField(tif_, TIFFTAG_IMAGEWIDTH, &width_);
    TIFFGetField(tif_, TIFFTAG_IMAGELENGTH, &height_);
    TIFFGetFieldDefaulted(tif_, TIFFTAG_BITSPERSAMPLE, &bits_per_sample_);
    TIFFGetFieldDefaulted(tif_, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel_);
    TIFFGetFieldDefaulted(tif_, TIFFTAG_COMPRESSION, &compression_);
    TIFFGetFieldDefaulted(tif_, TIFFTAG_SAMPLEFORMAT, &sample_format_);
    if (samples_per_pixel_ != 1) {
      throw std::runtime_error("Membrane scoring expects a single-channel TIFF image");
    }
    tiled_ = TIFFIsTiled(tif_) != 0;
    if (tiled_) {
      TIFFGetField(tif_, TIFFTAG_TILEWIDTH, &tile_width_);
      TIFFGetField(tif_, TIFFTAG_TILELENGTH, &tile_height_);
      TIFFGetField(tif_, TIFFTAG_TILEBYTECOUNTS, &tile_byte_counts_);
      if (tile_width_ == 0 || tile_height_ == 0 || tile_byte_counts_ == nullptr) {
        throw std::runtime_error("Invalid tiled TIFF metadata in membrane image");
      }
    } else {
      tile_width_ = width_;
      tile_height_ = height_;
    }
    tiles_across_ = (width_ + tile_width_ - 1U) / tile_width_;
    const std::size_t n_tiles = tiled_
        ? static_cast<std::size_t>(TIFFNumberOfTiles(tif_))
        : 1U;
    tile_store_.resize(n_tiles);
    tile_ready_ = std::vector<std::atomic<TilePixels*>>(n_tiles);
    for (auto& slot : tile_ready_) {
      slot.store(nullptr, std::memory_order_relaxed);
    }
  }

  ~TiffStainImage() {
    if (tif_ != nullptr) {
      TIFFClose(tif_);
    }
  }

  int width() const { return static_cast<int>(width_); }
  int height() const { return static_cast<int>(height_); }

  double sample(double x, double y) const {
    const double px = (x - options_.x_offset) / options_.pixel_size;
    const double py = (y - options_.y_offset) / options_.pixel_size;
    const double clamped_x = std::clamp(px, 0.0, static_cast<double>(width_ - 1U));
    const double clamped_y = std::clamp(py, 0.0, static_cast<double>(height_ - 1U));
    const int x0 = static_cast<int>(std::floor(clamped_x));
    const int y0 = static_cast<int>(std::floor(clamped_y));
    const int x1 = std::min(x0 + 1, static_cast<int>(width_ - 1U));
    const int y1 = std::min(y0 + 1, static_cast<int>(height_ - 1U));
    const double tx = clamped_x - static_cast<double>(x0);
    const double ty = clamped_y - static_cast<double>(y0);
    const double v00 = value_at_pixel(x0, y0);
    const double v10 = value_at_pixel(x1, y0);
    const double v01 = value_at_pixel(x0, y1);
    const double v11 = value_at_pixel(x1, y1);
    const double v0 = (1.0 - tx) * v00 + tx * v10;
    const double v1 = (1.0 - tx) * v01 + tx * v11;
    return (1.0 - ty) * v0 + ty * v1;
  }

 private:
  TilePixels decode_tile(std::uint32_t tile) const {
    if (!tiled_) {
      return decode_scanline_image();
    }
    if (compression_ == COMPRESSION_JP2000) {
      const auto n_bytes = static_cast<tmsize_t>(tile_byte_counts_[tile]);
      if (n_bytes <= 0) {
        throw std::runtime_error("Invalid TIFF raw tile byte count");
      }
      std::vector<unsigned char> raw(static_cast<std::size_t>(n_bytes));
      const auto got = TIFFReadRawTile(tif_, tile, raw.data(), n_bytes);
      if (got <= 0) {
        throw std::runtime_error("Failed to read raw JPEG2000 TIFF tile");
      }
      raw.resize(static_cast<std::size_t>(got));
      return decode_jpeg2000_tile(raw);
    }

    const auto decoded_size = TIFFTileSize(tif_);
    std::vector<unsigned char> buffer(static_cast<std::size_t>(decoded_size));
    const auto got = TIFFReadEncodedTile(tif_, tile, buffer.data(), decoded_size);
    if (got <= 0) {
      throw std::runtime_error("Failed to read encoded TIFF tile");
    }
    return convert_interleaved_buffer(buffer.data(), tile_width_, tile_height_);
  }

  TilePixels decode_scanline_image() const {
    TilePixels out;
    out.width = static_cast<int>(width_);
    out.height = static_cast<int>(height_);
    out.values.resize(static_cast<std::size_t>(out.width) * static_cast<std::size_t>(out.height));
    const tmsize_t scanline_size = TIFFScanlineSize(tif_);
    std::vector<unsigned char> row(static_cast<std::size_t>(scanline_size));
    for (std::uint32_t y = 0; y < height_; ++y) {
      if (TIFFReadScanline(tif_, row.data(), y, 0) < 0) {
        throw std::runtime_error("Failed to read TIFF scanline");
      }
      const auto converted = convert_interleaved_buffer(row.data(), width_, 1);
      std::copy(
          converted.values.begin(),
          converted.values.end(),
          out.values.begin() + static_cast<std::ptrdiff_t>(y * width_));
    }
    return out;
  }

  TilePixels convert_interleaved_buffer(const unsigned char* buffer, std::uint32_t width, std::uint32_t height) const {
    TilePixels out;
    out.width = static_cast<int>(width);
    out.height = static_cast<int>(height);
    out.values.resize(static_cast<std::size_t>(width) * static_cast<std::size_t>(height));
    if (bits_per_sample_ == 8) {
      for (std::size_t i = 0; i < out.values.size(); ++i) {
        out.values[i] = static_cast<float>(buffer[i]);
      }
    } else if (bits_per_sample_ == 16) {
      const auto* src = reinterpret_cast<const std::uint16_t*>(buffer);
      for (std::size_t i = 0; i < out.values.size(); ++i) {
        out.values[i] = static_cast<float>(src[i]);
      }
    } else if (bits_per_sample_ == 32 && sample_format_ == SAMPLEFORMAT_IEEEFP) {
      const auto* src = reinterpret_cast<const float*>(buffer);
      for (std::size_t i = 0; i < out.values.size(); ++i) {
        out.values[i] = src[i];
      }
    } else {
      throw std::runtime_error("Unsupported TIFF pixel format for membrane image");
    }
    return out;
  }

  // Warm reads are lock-free, and JPEG2000 tiles are decompressed outside
  // the lock (only the raw-byte read needs the shared TIFF handle), so both
  // pixel sampling and tile decoding scale with the worker count.
  const TilePixels& get_tile(std::uint32_t tile) const {
    TilePixels* ready = tile_ready_[tile].load(std::memory_order_acquire);
    if (ready != nullptr) {
      return *ready;
    }

    if (tiled_ && compression_ == COMPRESSION_JP2000) {
      std::vector<unsigned char> raw;
      {
        std::lock_guard<std::mutex> lock(mutex_);
        ready = tile_ready_[tile].load(std::memory_order_relaxed);
        if (ready != nullptr) {
          return *ready;
        }
        raw = read_raw_tile(tile);
      }
      auto decoded = std::make_unique<TilePixels>(decode_jpeg2000_tile(raw));
      std::lock_guard<std::mutex> lock(mutex_);
      ready = tile_ready_[tile].load(std::memory_order_relaxed);
      if (ready == nullptr) {
        tile_store_[tile] = std::move(decoded);
        ready = tile_store_[tile].get();
        tile_ready_[tile].store(ready, std::memory_order_release);
      }
      return *ready;
    }

    std::lock_guard<std::mutex> lock(mutex_);
    ready = tile_ready_[tile].load(std::memory_order_relaxed);
    if (ready == nullptr) {
      tile_store_[tile] = std::make_unique<TilePixels>(decode_tile(tile));
      ready = tile_store_[tile].get();
      tile_ready_[tile].store(ready, std::memory_order_release);
    }
    return *ready;
  }

  std::vector<unsigned char> read_raw_tile(std::uint32_t tile) const {
    const auto n_bytes = static_cast<tmsize_t>(tile_byte_counts_[tile]);
    if (n_bytes <= 0) {
      throw std::runtime_error("Invalid TIFF raw tile byte count");
    }
    std::vector<unsigned char> raw(static_cast<std::size_t>(n_bytes));
    const auto got = TIFFReadRawTile(tif_, tile, raw.data(), n_bytes);
    if (got <= 0) {
      throw std::runtime_error("Failed to read raw JPEG2000 TIFF tile");
    }
    raw.resize(static_cast<std::size_t>(got));
    return raw;
  }

  double value_at_pixel(int x, int y) const {
    const std::uint32_t px = static_cast<std::uint32_t>(std::clamp(x, 0, static_cast<int>(width_ - 1U)));
    const std::uint32_t py = static_cast<std::uint32_t>(std::clamp(y, 0, static_cast<int>(height_ - 1U)));
    const std::uint32_t tile = tiled_
        ? (py / tile_height_) * tiles_across_ + (px / tile_width_)
        : 0;
    const auto& pixels = get_tile(tile);
    const std::uint32_t tile_x = tiled_ ? (px / tile_width_) * tile_width_ : 0;
    const std::uint32_t tile_y = tiled_ ? (py / tile_height_) * tile_height_ : 0;
    const int local_x = std::clamp(static_cast<int>(px - tile_x), 0, std::max(0, pixels.width - 1));
    const int local_y = std::clamp(static_cast<int>(py - tile_y), 0, std::max(0, pixels.height - 1));
    return pixels.values[static_cast<std::size_t>(local_y) * static_cast<std::size_t>(pixels.width) +
                         static_cast<std::size_t>(local_x)];
  }

  MembraneImageOptions options_;
  mutable TIFF* tif_ = nullptr;
  std::uint32_t width_ = 0;
  std::uint32_t height_ = 0;
  std::uint32_t tile_width_ = 0;
  std::uint32_t tile_height_ = 0;
  std::uint16_t bits_per_sample_ = 0;
  std::uint16_t samples_per_pixel_ = 0;
  std::uint16_t compression_ = 0;
  std::uint16_t sample_format_ = SAMPLEFORMAT_UINT;
  std::uint64_t* tile_byte_counts_ = nullptr;
  bool tiled_ = false;
  std::uint32_t tiles_across_ = 1;
  mutable std::mutex mutex_;
  mutable std::vector<std::unique_ptr<TilePixels>> tile_store_;
  mutable std::vector<std::atomic<TilePixels*>> tile_ready_;
};

std::vector<BridgeCandidatePair> discover_membrane_candidates(
    const TranscriptTable& table,
    const MembranePreparedData& data,
    const MembraneTestOptions& options) {
  if (options.cell_candidate_k <= 0) {
    return {};
  }
  const int n_cells = static_cast<int>(data.members_by_cell.size());
  std::vector<int> active_cells;
  active_cells.reserve(data.members_by_cell.size());
  for (int cell = 0; cell < n_cells; ++cell) {
    if (!data.members_by_cell[static_cast<std::size_t>(cell)].empty()) {
      if (data.cell_type_index[static_cast<std::size_t>(cell)] < 0) {
        continue;
      }
      active_cells.push_back(cell);
    }
  }
  if (active_cells.empty()) {
    return {};
  }
  TranscriptTable centroids;
  centroids.x.reserve(active_cells.size());
  centroids.y.reserve(active_cells.size());
  centroids.z.reserve(active_cells.size());
  for (const int cell : active_cells) {
    centroids.x.push_back(data.centroid_x[static_cast<std::size_t>(cell)]);
    centroids.y.push_back(data.centroid_y[static_cast<std::size_t>(cell)]);
    centroids.z.push_back(data.centroid_z[static_cast<std::size_t>(cell)]);
  }
  const SpatialKnnIndex cell_index(centroids);

  std::vector<double> radius(static_cast<std::size_t>(n_cells), 0.0);
  for (int cell = 0; cell < n_cells; ++cell) {
    radius[static_cast<std::size_t>(cell)] = cell_radius_xy(
        table,
        data.members_by_cell[static_cast<std::size_t>(cell)],
        data.centroid_x[static_cast<std::size_t>(cell)],
        data.centroid_y[static_cast<std::size_t>(cell)]);
  }

  const int n_active = static_cast<int>(active_cells.size());
  const int workers = effective_threads(options.num_threads, n_active);
  std::vector<std::unordered_map<std::uint64_t, std::vector<BridgeCandidatePair>>> local_by_type(
      static_cast<std::size_t>(workers));

  subpar::parallelize_range<true>(workers, n_active, [&](int worker, int start, int length) {
    auto& by_type = local_by_type[static_cast<std::size_t>(worker)];
    for (int active_target = start; active_target < start + length; ++active_target) {
      const int target_cell = active_cells[static_cast<std::size_t>(active_target)];
      const auto hits = cell_index.query(active_target, options.cell_candidate_k, false);
      for (const auto& hit : hits) {
        const int source_cell = active_cells[static_cast<std::size_t>(hit.index)];
        if (source_cell == target_cell || source_cell < 0 ||
            source_cell >= n_cells ||
            data.members_by_cell[static_cast<std::size_t>(source_cell)].empty()) {
          continue;
        }
        const double dx = data.centroid_x[static_cast<std::size_t>(target_cell)] -
            data.centroid_x[static_cast<std::size_t>(source_cell)];
        const double dy = data.centroid_y[static_cast<std::size_t>(target_cell)] -
            data.centroid_y[static_cast<std::size_t>(source_cell)];
        const double distance = std::sqrt(dx * dx + dy * dy);
        const double radius_sum =
            radius[static_cast<std::size_t>(target_cell)] +
            radius[static_cast<std::size_t>(source_cell)];
        if (options.cell_candidate_halo >= 0.0 &&
            distance > radius_sum + options.cell_candidate_halo) {
          continue;
        }

        const double normalized_distance = distance / std::max(radius_sum, 1e-6);
        const int contact_rank = std::max(
            1,
            static_cast<int>(std::round(1000000.0 / (1.0 + normalized_distance))));
        const int target_type = data.cell_type_index[static_cast<std::size_t>(target_cell)];
        const int source_type = data.cell_type_index[static_cast<std::size_t>(source_cell)];
        by_type[pair_key(target_type, source_type)].push_back(
            {target_cell, source_cell, target_type, source_type, contact_rank});
      }
    }
  });

  std::unordered_map<std::uint64_t, std::vector<BridgeCandidatePair>> by_type_pair;
  for (auto& local : local_by_type) {
    for (auto& entry : local) {
      auto& dest = by_type_pair[entry.first];
      dest.insert(dest.end(), entry.second.begin(), entry.second.end());
    }
  }

  const int max_per_type_pair = options.candidate_pairs_per_type_pair > 0
      ? options.candidate_pairs_per_type_pair
      : std::max(options.max_cells_per_type_pair * 5, options.max_cells_per_type_pair);
  std::vector<BridgeCandidatePair> out;
  for (auto& entry : by_type_pair) {
    auto& pairs = entry.second;
    std::sort(pairs.begin(), pairs.end(), [](const auto& lhs, const auto& rhs) {
      if (lhs.contact_count != rhs.contact_count) return lhs.contact_count > rhs.contact_count;
      if (lhs.target_cell != rhs.target_cell) return lhs.target_cell < rhs.target_cell;
      return lhs.source_cell < rhs.source_cell;
    });

    std::unordered_set<int> used_targets;
    std::unordered_set<int> used_sources;
    std::unordered_set<int> used_any_same_type;
    int kept = 0;
    for (const auto& candidate : pairs) {
      if (candidate.target_type == candidate.source_type) {
        if (used_any_same_type.count(candidate.target_cell) ||
            used_any_same_type.count(candidate.source_cell)) {
          continue;
        }
        used_any_same_type.insert(candidate.target_cell);
        used_any_same_type.insert(candidate.source_cell);
      } else {
        if (used_targets.count(candidate.target_cell) ||
            used_sources.count(candidate.source_cell)) {
          continue;
        }
        used_targets.insert(candidate.target_cell);
        used_sources.insert(candidate.source_cell);
      }
      out.push_back(candidate);
      ++kept;
      if (max_per_type_pair > 0 && kept >= max_per_type_pair) {
        break;
      }
    }
  }

  std::sort(out.begin(), out.end(), [](const auto& lhs, const auto& rhs) {
    if (lhs.target_type != rhs.target_type) return lhs.target_type < rhs.target_type;
    if (lhs.source_type != rhs.source_type) return lhs.source_type < rhs.source_type;
    if (lhs.contact_count != rhs.contact_count) return lhs.contact_count > rhs.contact_count;
    if (lhs.target_cell != rhs.target_cell) return lhs.target_cell < rhs.target_cell;
    return lhs.source_cell < rhs.source_cell;
  });
  return out;
}

double line_signal(
    const TiffStainImage& image,
    double x0,
    double y0,
    double x1,
    double y1,
    int n_samples) {
  n_samples = std::max(2, n_samples);
  double best = 0.0;
  for (int i = 0; i < n_samples; ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(n_samples - 1);
    best = std::max(best, image.sample(x0 + (x1 - x0) * t, y0 + (y1 - y0) * t));
  }
  return best;
}

double directional_weight(
    double mx,
    double my,
    double target_x,
    double target_y,
    double source_x,
    double source_y) {
  const double tx = target_x - mx;
  const double ty = target_y - my;
  const double sx = source_x - mx;
  const double sy = source_y - my;
  const double tl = std::sqrt(tx * tx + ty * ty);
  const double sl = std::sqrt(sx * sx + sy * sy);
  if (tl <= 1e-9 || sl <= 1e-9) {
    return 0.0;
  }
  const double ux = tx / tl - sx / sl;
  const double uy = ty / tl - sy / sl;
  return 0.5 * std::sqrt(ux * ux + uy * uy);
}

int choose_control_member(
    const TranscriptTable& table,
    const std::vector<int>& members,
    const std::vector<double>& centroid_distance,
    int member_offset,
    double distance_tolerance,
    unsigned int seed) {
  const int query = members[static_cast<std::size_t>(member_offset)];
  const double qdist = centroid_distance[static_cast<std::size_t>(member_offset)];
  const double qz = table.z[static_cast<std::size_t>(query)];
  std::vector<int> candidates;
  candidates.reserve(members.size());
  double best_abs = std::numeric_limits<double>::infinity();
  for (int i = 0; i < static_cast<int>(members.size()); ++i) {
    if (i == member_offset) {
      continue;
    }
    const int row = members[static_cast<std::size_t>(i)];
    if (std::abs(table.z[static_cast<std::size_t>(row)] - qz) > 0.5) {
      continue;
    }
    const double delta = std::abs(centroid_distance[static_cast<std::size_t>(i)] - qdist);
    if (delta <= distance_tolerance) {
      candidates.push_back(i);
    }
    best_abs = std::min(best_abs, delta);
  }
  if (candidates.empty() && std::isfinite(best_abs)) {
    for (int i = 0; i < static_cast<int>(members.size()); ++i) {
      if (i == member_offset) {
        continue;
      }
      const int row = members[static_cast<std::size_t>(i)];
      if (std::abs(table.z[static_cast<std::size_t>(row)] - qz) > 0.5) {
        continue;
      }
      const double delta = std::abs(centroid_distance[static_cast<std::size_t>(i)] - qdist);
      if (std::abs(delta - best_abs) <= 1e-9) {
        candidates.push_back(i);
      }
    }
  }
  if (candidates.empty()) {
    return -1;
  }
  const std::uint64_t hash =
      static_cast<std::uint64_t>(seed) * 1103515245ULL +
      static_cast<std::uint64_t>(query + 1) * 2654435761ULL;
  return candidates[static_cast<std::size_t>(hash % candidates.size())];
}

std::vector<MembranePairScore> score_membrane_candidates(
    const TranscriptTable& table,
    const std::vector<int>& labels,
    const MembranePreparedData& data,
    const TiffStainImage& image,
    const std::vector<BridgeCandidatePair>& candidates,
    const MembraneImageOptions& image_options,
    const MembraneTestOptions& options) {
  const int workers = effective_threads(options.num_threads, static_cast<int>(candidates.size()));
  std::vector<std::vector<MembranePairScore>> local(static_cast<std::size_t>(workers));

  subpar::parallelize_range<true>(workers, static_cast<int>(candidates.size()), [&](int worker, int start, int length) {
    auto& out = local[static_cast<std::size_t>(worker)];
    for (int ci = start; ci < start + length; ++ci) {
      const auto& candidate = candidates[static_cast<std::size_t>(ci)];
      const int target = candidate.target_cell;
      const int source = candidate.source_cell;
      const auto& members = data.members_by_cell[static_cast<std::size_t>(target)];
      if (members.size() < 2U) {
        continue;
      }

      const double target_x = data.centroid_x[static_cast<std::size_t>(target)];
      const double target_y = data.centroid_y[static_cast<std::size_t>(target)];
      const double source_x = data.centroid_x[static_cast<std::size_t>(source)];
      const double source_y = data.centroid_y[static_cast<std::size_t>(source)];
      std::vector<double> centroid_distance(members.size(), 0.0);
      double max_dist = 0.0;
      for (std::size_t i = 0; i < members.size(); ++i) {
        const int row = members[i];
        const double dx = table.x[static_cast<std::size_t>(row)] - target_x;
        const double dy = table.y[static_cast<std::size_t>(row)] - target_y;
        centroid_distance[i] = std::sqrt(dx * dx + dy * dy);
        max_dist = std::max(max_dist, centroid_distance[i]);
      }
      const double tolerance = std::max(max_dist * options.control_distance_fraction, 1e-6);

      for (int factor = 0; factor < data.n_factors; ++factor) {
        const int factor_count =
            data.factor_counts[static_cast<std::size_t>(target * data.n_factors + factor)];
        if (factor_count < options.min_factor_molecules) {
          continue;
        }
        std::vector<double> scores;
        scores.reserve(static_cast<std::size_t>(factor_count));
        double positive = 0.0;
        double weight_sum = 0.0;
        for (int mi = 0; mi < static_cast<int>(members.size()); ++mi) {
          const int row = members[static_cast<std::size_t>(mi)];
          if (labels[static_cast<std::size_t>(row)] != factor) {
            continue;
          }
          const int control_offset = choose_control_member(
              table,
              members,
              centroid_distance,
              mi,
              tolerance,
              options.seed + static_cast<unsigned int>((factor + 1) * 7919 + (ci + 1) * 104729));
          if (control_offset < 0) {
            continue;
          }
          const int control_row = members[static_cast<std::size_t>(control_offset)];
          const double mx = table.x[static_cast<std::size_t>(row)];
          const double my = table.y[static_cast<std::size_t>(row)];
          const double cx = table.x[static_cast<std::size_t>(control_row)];
          const double cy = table.y[static_cast<std::size_t>(control_row)];
          const double real_signal = line_signal(image, mx, my, target_x, target_y, options.line_samples);
          const double control_signal = line_signal(image, cx, cy, target_x, target_y, options.line_samples);
          const double real_weight = directional_weight(mx, my, target_x, target_y, source_x, source_y);
          const double control_weight = directional_weight(cx, cy, target_x, target_y, source_x, source_y);
          const double real = real_signal * std::max(real_weight, 1e-6) + image_options.epsilon;
          const double control = control_signal * std::max(control_weight, 1e-6) + image_options.epsilon;
          const double score = std::log(real / control);
          if (std::isfinite(score)) {
            scores.push_back(score);
            positive += score > 0.0 ? 1.0 : 0.0;
            weight_sum += real_weight;
          }
        }
        if (scores.empty()) {
          continue;
        }
        MembranePairScore row;
        row.target_cell = target;
        row.source_cell = source;
        row.target_type = candidate.target_type;
        row.source_type = candidate.source_type;
        row.factor = factor;
        row.factor_count = factor_count;
        row.scored_molecules = static_cast<int>(scores.size());
        row.mean_score = mean_value(scores);
        row.fraction_positive = positive / static_cast<double>(scores.size());
        row.mean_directional_weight = weight_sum / static_cast<double>(scores.size());
        out.push_back(row);
      }
    }
  });

  std::size_t total = 0;
  for (const auto& piece : local) total += piece.size();
  std::vector<MembranePairScore> out;
  out.reserve(total);
  for (auto& piece : local) {
    out.insert(out.end(), piece.begin(), piece.end());
  }
  std::sort(out.begin(), out.end(), [](const auto& lhs, const auto& rhs) {
    if (lhs.target_type != rhs.target_type) return lhs.target_type < rhs.target_type;
    if (lhs.source_type != rhs.source_type) return lhs.source_type < rhs.source_type;
    if (lhs.factor != rhs.factor) return lhs.factor < rhs.factor;
    if (lhs.mean_score != rhs.mean_score) return lhs.mean_score > rhs.mean_score;
    if (lhs.target_cell != rhs.target_cell) return lhs.target_cell < rhs.target_cell;
    return lhs.source_cell < rhs.source_cell;
  });
  return out;
}

std::vector<MembraneSummary> summarize_membrane_scores(
    std::vector<MembranePairScore>& pair_scores,
    const MembraneTestOptions& options) {
  std::map<std::tuple<int, int, int>, std::vector<int>> groups;
  for (int i = 0; i < static_cast<int>(pair_scores.size()); ++i) {
    const auto& row = pair_scores[static_cast<std::size_t>(i)];
    groups[std::make_tuple(row.target_type, row.source_type, row.factor)].push_back(i);
  }

  std::vector<MembraneSummary> out;
  for (auto& entry : groups) {
    auto& indices = entry.second;
    std::sort(indices.begin(), indices.end(), [&](int lhs, int rhs) {
      const auto& lrow = pair_scores[static_cast<std::size_t>(lhs)];
      const auto& rrow = pair_scores[static_cast<std::size_t>(rhs)];
      if (lrow.mean_score != rrow.mean_score) return lrow.mean_score > rrow.mean_score;
      if (lrow.target_cell != rrow.target_cell) return lrow.target_cell < rrow.target_cell;
      return lrow.source_cell < rrow.source_cell;
    });
    const int keep = options.max_cells_per_type_pair > 0
        ? std::min(options.max_cells_per_type_pair, static_cast<int>(indices.size()))
        : static_cast<int>(indices.size());
    if (keep < options.min_pairs) {
      continue;
    }
    std::vector<double> scores;
    scores.reserve(static_cast<std::size_t>(keep));
    double positive_pairs = 0.0;
    for (int i = 0; i < keep; ++i) {
      auto& row = pair_scores[static_cast<std::size_t>(indices[static_cast<std::size_t>(i)])];
      row.used_in_summary = true;
      scores.push_back(row.mean_score);
      positive_pairs += row.mean_score > 0.0 ? 1.0 : 0.0;
    }
    MembraneSummary summary;
    summary.target_type = std::get<0>(entry.first);
    summary.source_type = std::get<1>(entry.first);
    summary.factor = std::get<2>(entry.first);
    summary.n_pairs = keep;
    summary.mean_score = mean_value(scores);
    summary.q75_score = quantile_75(scores);
    summary.fraction_positive_pairs = positive_pairs / static_cast<double>(keep);
    summary.p_value = one_sample_wilcoxon_greater(scores);
    summary.neg_log10_p = summary.p_value > 0.0 ? -std::log10(summary.p_value) : 300.0;
    out.push_back(summary);
  }

  std::sort(out.begin(), out.end(), [](const auto& lhs, const auto& rhs) {
    if (lhs.target_type != rhs.target_type) return lhs.target_type < rhs.target_type;
    if (lhs.source_type != rhs.source_type) return lhs.source_type < rhs.source_type;
    return lhs.factor < rhs.factor;
  });
  return out;
}

}  // namespace

Image2D read_stain_image_crop(
    const MembraneImageOptions& image_options,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    int max_pixels_per_axis) {
  if (!(xmax > xmin) || !(ymax > ymin)) {
    throw std::runtime_error("Image crop bbox must satisfy xmax > xmin and ymax > ymin");
  }
  if (max_pixels_per_axis <= 0) {
    throw std::runtime_error("max_pixels_per_axis must be positive");
  }
  TiffStainImage image(image_options);
  const double width_pixels = (xmax - xmin) / image_options.pixel_size;
  const double height_pixels = (ymax - ymin) / image_options.pixel_size;
  const double scale = std::max(1.0, std::max(width_pixels, height_pixels) /
      static_cast<double>(max_pixels_per_axis));
  const int out_width = std::max(1, static_cast<int>(std::ceil(width_pixels / scale)));
  const int out_height = std::max(1, static_cast<int>(std::ceil(height_pixels / scale)));

  Image2D out;
  out.width = out_width;
  out.height = out_height;
  out.pixel_size = image_options.pixel_size * scale;
  out.x_origin = xmin;
  out.y_origin = ymin;
  out.values.resize(static_cast<std::size_t>(out_width) * static_cast<std::size_t>(out_height));

  for (int row = 0; row < out_height; ++row) {
    const double y = ymin + (static_cast<double>(row) + 0.5) *
        (ymax - ymin) / static_cast<double>(out_height);
    for (int col = 0; col < out_width; ++col) {
      const double x = xmin + (static_cast<double>(col) + 0.5) *
          (xmax - xmin) / static_cast<double>(out_width);
      out.values[static_cast<std::size_t>(row) * static_cast<std::size_t>(out_width) +
                 static_cast<std::size_t>(col)] = image.sample(x, y);
    }
  }
  return out;
}

MembraneTestResult run_membrane_test(
    const TranscriptTable& table,
    const CellTable& cells,
    const std::vector<int>& labels,
    const MembraneImageOptions& image_options,
    const MembraneTestOptions& options) {
  if (options.min_factor_molecules < 1) {
    throw std::runtime_error("min_factor_molecules must be positive");
  }
  if (options.line_samples < 2) {
    throw std::runtime_error("line_samples must be at least 2");
  }

  const auto start = std::chrono::steady_clock::now();
  auto stage_start = start;
  emit_membrane_progress(
      options,
      "Starting membrane test: molecules=" + std::to_string(table.size()) +
          ", cells=" + std::to_string(table.num_cells()),
      start);

  TiffStainImage image(image_options);
  emit_membrane_progress(
      options,
      "Opened membrane image: width=" + std::to_string(image.width()) +
          ", height=" + std::to_string(image.height()),
      stage_start);

  stage_start = std::chrono::steady_clock::now();
  auto prepared = prepare_membrane_data(table, cells, labels);
  int active_cells = 0;
  int typed_active_cells = 0;
  for (std::size_t cell = 0; cell < prepared.members_by_cell.size(); ++cell) {
    if (!prepared.members_by_cell[cell].empty()) {
      ++active_cells;
      if (prepared.cell_type_index[cell] >= 0) {
        ++typed_active_cells;
      }
    }
  }
  emit_membrane_progress(
      options,
      "Prepared membrane data: cell_types=" + std::to_string(prepared.type_names.size()) +
          ", factors=" + std::to_string(prepared.n_factors) +
          ", active_cells=" + std::to_string(active_cells) +
          ", typed_active_cells=" + std::to_string(typed_active_cells),
      stage_start);

  stage_start = std::chrono::steady_clock::now();
  const auto candidates = discover_membrane_candidates(table, prepared, options);
  emit_membrane_progress(
      options,
      "Discovered membrane candidates: pairs=" + std::to_string(candidates.size()),
      stage_start);

  stage_start = std::chrono::steady_clock::now();
  MembraneTestResult result;
  result.cell_types = prepared.type_names;
  result.pair_scores = score_membrane_candidates(
      table,
      labels,
      prepared,
      image,
      candidates,
      image_options,
      options);
  emit_membrane_progress(
      options,
      "Scored membrane candidates: score_rows=" + std::to_string(result.pair_scores.size()),
      stage_start);

  stage_start = std::chrono::steady_clock::now();
  result.summaries = summarize_membrane_scores(result.pair_scores, options);
  emit_membrane_progress(
      options,
      "Summarized membrane scores: summaries=" + std::to_string(result.summaries.size()),
      stage_start);
  emit_membrane_progress(
      options,
      "Finished membrane test: summaries=" + std::to_string(result.summaries.size()),
      start);
  return result;
}

}  // namespace celladmix
