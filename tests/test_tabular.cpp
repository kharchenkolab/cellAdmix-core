#include "celladmix/tabular.hpp"
#include "celladmix/input_store.hpp"

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <tiffio.h>

#include "test_framework.hpp"

using namespace celladmix;

namespace {

void write_text(const std::filesystem::path& path, const std::string& text) {
  std::ofstream out(path);
  out << text;
}

void write_label_tiff(
    const std::filesystem::path& path,
    uint32_t width,
    uint32_t height,
    const std::vector<std::uint16_t>& labels) {
  TIFF* tif = TIFFOpen(path.string().c_str(), "w");
  if (tif == nullptr) {
    throw std::runtime_error("failed to open TIFF for writing");
  }

  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
  TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, height);

  for (uint32_t row = 0; row < height; ++row) {
    const auto* scanline = reinterpret_cast<const void*>(
        labels.data() + static_cast<std::size_t>(row) * static_cast<std::size_t>(width));
    if (TIFFWriteScanline(tif, const_cast<void*>(scanline), row, 0) < 0) {
      TIFFClose(tif);
      throw std::runtime_error("failed to write TIFF scanline");
    }
  }

  TIFFClose(tif);
}

}  // namespace

TEST_CASE("Tabular loader accepts explicit molecule-level cell ids") {
  const auto dir = std::filesystem::temp_directory_path() / "celladmix_tabular_test_ids";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);

  write_text(
      dir / "molecules.csv",
      "x,y,gene,cell_id,cell_type,qv\n"
      "1,1,G1,c1,fibro,20\n"
      "2,2,G2,c1,fibro,20\n"
      "25,25,G3,c2,malignant,30\n"
      "3,3,G4,0,,30\n");

  TabularSourceSpec source;
  source.molecules_path = (dir / "molecules.csv").string();
  source.gene_col = "gene";
  source.x_col = "x";
  source.y_col = "y";
  source.cell_id_col = "cell_id";
  source.cell_type_col = "cell_type";
  source.qv_col = "qv";

  TabularLoadOptions options;
  options.keep_unassigned = false;
  options.min_qv = 10.0;
  options.crops.push_back(CropBox{"crop_1", 0.0, 10.0, 0.0, 10.0, 0.0, 0.0, false});

  const auto bundle = load_tabular_bundle(source, options);
  REQUIRE_EQ(bundle.used_parquet, false);
  REQUIRE_EQ(bundle.transcripts.size(), 2U);
  REQUIRE_EQ(bundle.transcripts.orig_cell_id[0], std::string("c1"));
  REQUIRE_EQ(bundle.transcripts.orig_cell_id[1], std::string("c1"));
  REQUIRE_EQ(bundle.transcripts.cell_types[0], std::string("fibro"));
  REQUIRE_EQ(bundle.transcript_crop_ids.size(), 2U);
  REQUIRE_EQ(bundle.transcript_crop_ids[0], std::string("crop_1"));
  REQUIRE_EQ(bundle.transcript_crop_ids[1], std::string("crop_1"));

  std::filesystem::remove_all(dir);
}

TEST_CASE("Tabular loader can assign cells from a labeled TIFF mask") {
  const auto dir = std::filesystem::temp_directory_path() / "celladmix_tabular_test_mask";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);

  write_text(
      dir / "molecules.csv",
      "x,y,gene\n"
      "1,1,G1\n"
      "2,2,G2\n"
      "3,3,G3\n");

  write_label_tiff(
      dir / "mask.tiff",
      3,
      3,
      {
          1, 0, 0,
          0, 2, 0,
          0, 0, 0,
      });

  TabularSourceSpec source;
  source.molecules_path = (dir / "molecules.csv").string();
  source.gene_col = "gene";
  source.x_col = "x";
  source.y_col = "y";
  source.segmentation_mask_path = (dir / "mask.tiff").string();

  TabularLoadOptions options;
  options.keep_unassigned = false;

  const auto bundle = load_tabular_bundle(source, options);
  REQUIRE_EQ(bundle.transcripts.size(), 2U);
  REQUIRE_EQ(bundle.transcripts.orig_cell_id[0], std::string("1"));
  REQUIRE_EQ(bundle.transcripts.orig_cell_id[1], std::string("2"));
  REQUIRE_EQ(bundle.transcripts.genes.size(), 2U);

  std::filesystem::remove_all(dir);
}

TEST_CASE("Tabular input store streams labeled TIFF mask assignment") {
  const auto dir = std::filesystem::temp_directory_path() / "celladmix_tabular_store_test_mask";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);

  write_text(
      dir / "molecules.csv",
      "x,y,gene\n"
      "1,1,G1\n"
      "2,2,G2\n"
      "3,3,G3\n");

  write_label_tiff(
      dir / "mask.tiff",
      3,
      3,
      {
          1, 0, 0,
          0, 2, 0,
          0, 0, 0,
      });

  TabularSourceSpec source;
  source.molecules_path = (dir / "molecules.csv").string();
  source.gene_col = "gene";
  source.x_col = "x";
  source.y_col = "y";
  source.segmentation_mask_path = (dir / "mask.tiff").string();

  TabularLoadOptions load_options;
  load_options.keep_unassigned = false;

  InputStoreBuildOptions store_options;
  store_options.store_dir = (dir / "input_store").string();
  store_options.materialize_molecules = true;
  store_options.force = true;

  const auto manifest = build_tabular_input_store(source, load_options, store_options);
  REQUIRE_EQ(manifest.store_mode, std::string("full"));
  REQUIRE_EQ(manifest.has_molecule_rows, true);
  REQUIRE_EQ(manifest.n_molecules, 2U);
  REQUIRE_EQ(manifest.n_cells, 2U);

  const auto counts = load_input_store_counts(store_options.store_dir);
  REQUIRE_EQ(counts.cells.cell_ids.size(), 2U);
  REQUIRE_EQ(counts.cells.cell_ids[0], std::string("1"));
  REQUIRE_EQ(counts.cells.cell_ids[1], std::string("2"));
  REQUIRE_EQ(counts.transcript_counts[0], 1);
  REQUIRE_EQ(counts.transcript_counts[1], 1);

  const auto block = load_input_store_molecule_range(store_options.store_dir, 0, 2);
  REQUIRE_EQ(block.size(), 2U);
  REQUIRE_EQ(block.cell_idx[0], 0);
  REQUIRE_EQ(block.cell_idx[1], 1);

  std::filesystem::remove_all(dir);
}
