#include "celladmix/xenium.hpp"

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>

#include <filesystem>
#include <fstream>
#include <string>

#include <zlib.h>

#include "test_framework.hpp"

using namespace celladmix;

namespace {

void require_arrow_status(const arrow::Status& status, const char* context) {
  if (!status.ok()) {
    throw std::runtime_error(std::string(context) + ": " + status.ToString());
  }
}

void write_text(const std::filesystem::path& path, const std::string& text) {
  std::ofstream out(path);
  out << text;
}

void write_gzip_text(const std::filesystem::path& path, const std::string& text) {
  gzFile file = gzopen(path.string().c_str(), "wb");
  if (file == nullptr) {
    throw std::runtime_error("failed to open gzip file for writing");
  }
  gzwrite(file, text.data(), static_cast<unsigned int>(text.size()));
  gzclose(file);
}

void write_transcripts_parquet(const std::filesystem::path& path, int64_t row_group_size = 1024) {
  arrow::StringBuilder transcript_id;
  arrow::Int32Builder cell_id;
  arrow::StringBuilder feature_name;
  arrow::DoubleBuilder x_location;
  arrow::DoubleBuilder y_location;
  arrow::DoubleBuilder z_location;
  arrow::DoubleBuilder qv;
  arrow::StringBuilder fov_name;

  require_arrow_status(transcript_id.Append("t1"), "append transcript_id");
  require_arrow_status(transcript_id.Append("t2"), "append transcript_id");
  require_arrow_status(transcript_id.Append("t3"), "append transcript_id");
  require_arrow_status(cell_id.Append(1), "append cell_id");
  require_arrow_status(cell_id.Append(0), "append cell_id");
  require_arrow_status(cell_id.Append(2), "append cell_id");
  require_arrow_status(feature_name.Append("G1"), "append feature_name");
  require_arrow_status(feature_name.Append("G2"), "append feature_name");
  require_arrow_status(feature_name.Append("G3"), "append feature_name");
  require_arrow_status(x_location.Append(10.0), "append x_location");
  require_arrow_status(x_location.Append(50.0), "append x_location");
  require_arrow_status(x_location.Append(110.0), "append x_location");
  require_arrow_status(y_location.Append(10.0), "append y_location");
  require_arrow_status(y_location.Append(50.0), "append y_location");
  require_arrow_status(y_location.Append(110.0), "append y_location");
  require_arrow_status(z_location.Append(0.0), "append z_location");
  require_arrow_status(z_location.Append(0.0), "append z_location");
  require_arrow_status(z_location.Append(0.0), "append z_location");
  require_arrow_status(qv.Append(20.0), "append qv");
  require_arrow_status(qv.Append(20.0), "append qv");
  require_arrow_status(qv.Append(35.0), "append qv");
  require_arrow_status(fov_name.Append("F1"), "append fov_name");
  require_arrow_status(fov_name.Append("F1"), "append fov_name");
  require_arrow_status(fov_name.Append("F2"), "append fov_name");

  auto table = arrow::Table::Make(
      arrow::schema({
          arrow::field("transcript_id", arrow::utf8()),
          arrow::field("cell_id", arrow::int32()),
          arrow::field("feature_name", arrow::utf8()),
          arrow::field("x_location", arrow::float64()),
          arrow::field("y_location", arrow::float64()),
          arrow::field("z_location", arrow::float64()),
          arrow::field("qv", arrow::float64()),
          arrow::field("fov_name", arrow::utf8()),
      }),
      {
          transcript_id.Finish().ValueOrDie(),
          cell_id.Finish().ValueOrDie(),
          feature_name.Finish().ValueOrDie(),
          x_location.Finish().ValueOrDie(),
          y_location.Finish().ValueOrDie(),
          z_location.Finish().ValueOrDie(),
          qv.Finish().ValueOrDie(),
          fov_name.Finish().ValueOrDie(),
      });

  auto sink = arrow::io::FileOutputStream::Open(path.string()).ValueOrDie();
  require_arrow_status(
      parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), sink, row_group_size),
      "write transcripts parquet");
}

void write_cells_parquet(const std::filesystem::path& path, int64_t row_group_size = 1024) {
  arrow::Int32Builder cell_id;
  arrow::DoubleBuilder x_centroid;
  arrow::DoubleBuilder y_centroid;

  require_arrow_status(cell_id.Append(1), "append cell_id");
  require_arrow_status(cell_id.Append(2), "append cell_id");
  require_arrow_status(x_centroid.Append(10.0), "append x_centroid");
  require_arrow_status(x_centroid.Append(110.0), "append x_centroid");
  require_arrow_status(y_centroid.Append(10.0), "append y_centroid");
  require_arrow_status(y_centroid.Append(110.0), "append y_centroid");

  auto table = arrow::Table::Make(
      arrow::schema({
          arrow::field("cell_id", arrow::int32()),
          arrow::field("x_centroid", arrow::float64()),
          arrow::field("y_centroid", arrow::float64()),
      }),
      {
          cell_id.Finish().ValueOrDie(),
          x_centroid.Finish().ValueOrDie(),
          y_centroid.Finish().ValueOrDie(),
      });

  auto sink = arrow::io::FileOutputStream::Open(path.string()).ValueOrDie();
  require_arrow_status(
      parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), sink, row_group_size),
      "write cells parquet");
}

}  // namespace

TEST_CASE("Xenium bundle loader applies crop filtering at load time") {
  const auto dir = std::filesystem::temp_directory_path() / "celladmix_xenium_test_bundle";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);

  write_text(
      dir / "experiment.xenium",
      "{\n"
      "  \"run_name\": \"Mock Xenium\",\n"
      "  \"pixel_size\": 0.2125,\n"
      "  \"z_step_size\": 3.0\n"
      "}\n");

  write_gzip_text(
      dir / "transcripts.csv.gz",
      "transcript_id,cell_id,feature_name,x_location,y_location,z_location,qv,fov_name\n"
      "t1,c1,G1,10,10,0,20,F1\n"
      "t2,UNASSIGNED,G2,50,50,0,20,F1\n"
      "t3,c2,G3,110,110,0,35,F2\n");

  write_gzip_text(
      dir / "cells.csv.gz",
      "cell_id,x_centroid,y_centroid,z_centroid,cell_type\n"
      "c1,10,10,0,fibroblast\n"
      "c2,110,110,0,malignant\n");

  XeniumLoadOptions options;
  options.min_qv = 10.0;
  options.keep_unassigned = false;
  options.read_cells = true;
  options.prefer_parquet = true;
  options.crops.push_back(CropBox{"crop_1", 0.0, 40.0, 0.0, 40.0, 0.0, 0.0, false});

  const auto bundle = load_xenium_bundle(dir.string(), options);
  REQUIRE_EQ(bundle.used_parquet, false);
  REQUIRE_EQ(bundle.manifest.run_name, std::string("Mock Xenium"));
  REQUIRE_EQ(bundle.transcripts.size(), 1U);
  REQUIRE_EQ(bundle.transcripts.orig_cell_id[0], std::string("c1"));
  REQUIRE_EQ(bundle.transcript_crop_ids.size(), 1U);
  REQUIRE_EQ(bundle.transcript_crop_ids[0], std::string("crop_1"));
  REQUIRE_EQ(bundle.cells.size(), 1U);
  REQUIRE_EQ(bundle.cells.cell_ids[0], std::string("c1"));
  REQUIRE_EQ(bundle.cells.cell_types[0], std::string("fibroblast"));
  REQUIRE_EQ(bundle.cell_crop_ids[0], std::string("crop_1"));

  std::filesystem::remove_all(dir);
}

TEST_CASE("Xenium bundle loader prefers parquet and accepts numeric identifiers") {
  const auto dir = std::filesystem::temp_directory_path() / "celladmix_xenium_test_bundle_parquet";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);

  write_text(
      dir / "experiment.xenium",
      "{\n"
      "  \"run_name\": \"Mock Xenium Parquet\",\n"
      "  \"pixel_size\": 0.2125,\n"
      "  \"z_step_size\": 3.0\n"
      "}\n");

  write_transcripts_parquet(dir / "transcripts.parquet", 1);
  write_cells_parquet(dir / "cells.parquet", 1);

  XeniumLoadOptions options;
  options.min_qv = 10.0;
  options.keep_unassigned = false;
  options.read_cells = true;
  options.prefer_parquet = true;
  options.crops.push_back(CropBox{"crop_1", 0.0, 40.0, 0.0, 40.0, 0.0, 0.0, false});

  const auto bundle = load_xenium_bundle(dir.string(), options);
  REQUIRE_EQ(bundle.used_parquet, true);
  REQUIRE_EQ(bundle.manifest.run_name, std::string("Mock Xenium Parquet"));
  REQUIRE_EQ(bundle.transcripts.size(), 1U);
  REQUIRE_EQ(bundle.transcripts.orig_cell_id[0], std::string("1"));
  REQUIRE_EQ(bundle.transcript_crop_ids.size(), 1U);
  REQUIRE_EQ(bundle.transcript_crop_ids[0], std::string("crop_1"));
  REQUIRE_EQ(bundle.cells.size(), 1U);
  REQUIRE_EQ(bundle.cells.cell_ids[0], std::string("1"));
  REQUIRE_EQ(bundle.cell_crop_ids[0], std::string("crop_1"));

  std::filesystem::remove_all(dir);
}

TEST_CASE("Xenium parquet loader handles crops with no matching row groups") {
  const auto dir =
      std::filesystem::temp_directory_path() / "celladmix_xenium_test_bundle_parquet_empty_crop";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);

  write_text(
      dir / "experiment.xenium",
      "{\n"
      "  \"run_name\": \"Mock Xenium Empty Crop\",\n"
      "  \"pixel_size\": 0.2125,\n"
      "  \"z_step_size\": 3.0\n"
      "}\n");

  write_transcripts_parquet(dir / "transcripts.parquet", 1);
  write_cells_parquet(dir / "cells.parquet", 1);

  XeniumLoadOptions options;
  options.min_qv = 10.0;
  options.keep_unassigned = false;
  options.read_cells = true;
  options.prefer_parquet = true;
  options.crops.push_back(CropBox{"crop_1", 200.0, 240.0, 200.0, 240.0, 0.0, 0.0, false});

  const auto bundle = load_xenium_bundle(dir.string(), options);
  REQUIRE_EQ(bundle.used_parquet, true);
  REQUIRE_EQ(bundle.transcripts.size(), 0U);
  REQUIRE_EQ(bundle.transcript_crop_ids.size(), 0U);
  REQUIRE_EQ(bundle.cells.size(), 0U);
  REQUIRE_EQ(bundle.cell_crop_ids.size(), 0U);

  std::filesystem::remove_all(dir);
}
