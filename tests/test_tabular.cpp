#include "celladmix/tabular.hpp"
#include "celladmix/input_store.hpp"

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <tiffio.h>

#include "test_framework.hpp"

using namespace celladmix;

namespace {

void require_arrow_status(const arrow::Status& status, const char* context) {
  if (!status.ok()) {
    throw std::runtime_error(std::string(context) + ": " + status.ToString());
  }
}

// Write a synthetic molecules parquet table with num_rows rows spread over
// cells distinct cell ids and genes distinct gene names, split into
// row_group_size-sized row groups. The many row groups keep the streaming
// parquet reader busy reading metadata long after the store source returned.
void write_molecules_parquet(
    const std::filesystem::path& path,
    int num_rows,
    int cells,
    int genes,
    std::int64_t row_group_size) {
  arrow::DoubleBuilder x;
  arrow::DoubleBuilder y;
  arrow::StringBuilder gene;
  arrow::StringBuilder cell;
  require_arrow_status(x.Reserve(num_rows), "reserve x");
  require_arrow_status(y.Reserve(num_rows), "reserve y");
  require_arrow_status(gene.Reserve(num_rows), "reserve gene");
  require_arrow_status(cell.Reserve(num_rows), "reserve cell");
  for (int i = 0; i < num_rows; ++i) {
    require_arrow_status(x.Append(static_cast<double>(i % 512)), "append x");
    require_arrow_status(y.Append(static_cast<double>(i % 383)), "append y");
    require_arrow_status(
        gene.Append("G" + std::to_string(i % genes)), "append gene");
    require_arrow_status(
        cell.Append("c" + std::to_string(i % cells)), "append cell");
  }

  auto table = arrow::Table::Make(
      arrow::schema({
          arrow::field("x", arrow::float64()),
          arrow::field("y", arrow::float64()),
          arrow::field("gene", arrow::utf8()),
          arrow::field("cell", arrow::utf8()),
      }),
      {
          x.Finish().ValueOrDie(),
          y.Finish().ValueOrDie(),
          gene.Finish().ValueOrDie(),
          cell.Finish().ValueOrDie(),
      });

  auto sink = arrow::io::FileOutputStream::Open(path.string()).ValueOrDie();
  require_arrow_status(
      parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), sink, row_group_size),
      "write molecules parquet");
}

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

TEST_CASE("Tabular parquet store build keeps the parquet reader alive while streaming") {
  const auto dir = std::filesystem::temp_directory_path() /
      "celladmix_tabular_store_test_parquet";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);

  constexpr int kRows = 50000;
  constexpr int kCells = 500;
  constexpr int kGenes = 40;
  write_molecules_parquet(dir / "molecules.parquet", kRows, kCells, kGenes, 4096);

  TabularSourceSpec source;
  source.molecules_path = (dir / "molecules.parquet").string();
  source.gene_col = "gene";
  source.x_col = "x";
  source.y_col = "y";
  source.cell_id_col = "cell";

  TabularLoadOptions load_options;

  InputStoreBuildOptions store_options;
  store_options.store_dir = (dir / "input_store").string();
  store_options.materialize_molecules = true;
  store_options.force = true;

  const auto manifest = build_tabular_input_store(source, load_options, store_options);
  REQUIRE_EQ(manifest.store_mode, std::string("full"));
  REQUIRE_EQ(manifest.has_molecule_rows, true);
  REQUIRE_EQ(manifest.n_molecules, static_cast<std::size_t>(kRows));
  REQUIRE_EQ(manifest.n_cells, static_cast<std::size_t>(kCells));
  REQUIRE_EQ(manifest.n_genes, static_cast<std::size_t>(kGenes));

  const auto counts = load_input_store_counts(store_options.store_dir);
  REQUIRE_EQ(counts.transcript_counts.size(), static_cast<std::size_t>(kCells));
  std::size_t total = 0;
  for (const auto value : counts.transcript_counts) {
    total += static_cast<std::size_t>(value);
  }
  REQUIRE_EQ(total, static_cast<std::size_t>(kRows));

  std::filesystem::remove_all(dir);
}
