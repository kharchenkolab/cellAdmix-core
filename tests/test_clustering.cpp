#include "celladmix/clustering.hpp"

#include <set>
#include <string>
#include <vector>

#include "test_framework.hpp"

using namespace celladmix;

namespace {

TranscriptTable make_cell_clustering_fixture() {
  TranscriptTable table;
  const std::vector<std::string> malignant_genes = {"KRT19", "KRT8", "KRT17"};
  const std::vector<std::string> fibro_genes = {"COL1A1", "COL3A1", "DCN"};

  int tx_id = 0;
  for (int cell = 0; cell < 4; ++cell) {
    const std::string cell_id = "malignant_" + std::to_string(cell + 1);
    const double base_x = 10.0 + cell * 20.0;
    const double base_y = 10.0;
    for (int rep = 0; rep < 12; ++rep) {
      for (std::size_t gene = 0; gene < malignant_genes.size(); ++gene) {
        table.gene_key.push_back(malignant_genes[gene]);
        table.orig_cell_id.push_back(cell_id);
        table.cell_types.push_back("malignant");
        table.x.push_back(base_x + static_cast<double>(rep) * 0.1 + gene * 0.01);
        table.y.push_back(base_y + static_cast<double>(rep) * 0.1);
        table.z.push_back(0.0);
        table.transcript_ids.push_back("tx_" + std::to_string(tx_id++));
      }
      table.gene_key.push_back("COL1A1");
      table.orig_cell_id.push_back(cell_id);
      table.cell_types.push_back("malignant");
      table.x.push_back(base_x + static_cast<double>(rep) * 0.1);
      table.y.push_back(base_y + 2.0 + static_cast<double>(rep) * 0.1);
      table.z.push_back(0.0);
      table.transcript_ids.push_back("tx_" + std::to_string(tx_id++));
    }
  }

  for (int cell = 0; cell < 4; ++cell) {
    const std::string cell_id = "fibroblast_" + std::to_string(cell + 1);
    const double base_x = 10.0 + cell * 20.0;
    const double base_y = 60.0;
    for (int rep = 0; rep < 12; ++rep) {
      for (std::size_t gene = 0; gene < fibro_genes.size(); ++gene) {
        table.gene_key.push_back(fibro_genes[gene]);
        table.orig_cell_id.push_back(cell_id);
        table.cell_types.push_back("fibroblast");
        table.x.push_back(base_x + static_cast<double>(rep) * 0.1 + gene * 0.01);
        table.y.push_back(base_y + static_cast<double>(rep) * 0.1);
        table.z.push_back(0.0);
        table.transcript_ids.push_back("tx_" + std::to_string(tx_id++));
      }
      table.gene_key.push_back("KRT19");
      table.orig_cell_id.push_back(cell_id);
      table.cell_types.push_back("fibroblast");
      table.x.push_back(base_x + static_cast<double>(rep) * 0.1);
      table.y.push_back(base_y + 2.0 + static_cast<double>(rep) * 0.1);
      table.z.push_back(0.0);
      table.transcript_ids.push_back("tx_" + std::to_string(tx_id++));
    }
  }

  table.finalize();
  return table;
}

}  // namespace

TEST_CASE("Cell clustering separates simple malignant and fibroblast programs") {
  const auto table = make_cell_clustering_fixture();
  CellClusteringOptions options;
  options.min_molecules = 10;
  options.min_genes = 1;
  options.cells_max = -1;
  options.n_variable_genes = 6;
  options.pca_dims = 4;
  options.graph_k = 3;
  options.cluster_resolution = 1.0;
  options.umap_neighbors = 4;
  options.umap_epochs = 50;
  options.seed = 13;

  const auto result = cluster_cells(table, options, nullptr);
  REQUIRE_EQ(result.cells.size(), 8U);
  REQUIRE_EQ(result.clusters.size(), 8U);
  REQUIRE_EQ(result.umap.rows(), 8);
  REQUIRE_EQ(result.umap.cols(), 2);
  REQUIRE_EQ(result.pcs.rows(), 8);
  REQUIRE_GE(result.pcs.cols(), 2);

  std::set<int> cluster_ids(result.clusters.begin(), result.clusters.end());
  REQUIRE_GE(cluster_ids.size(), 2U);

  const int malignant_cluster = result.clusters[0];
  const int fibro_cluster = result.clusters[4];
  REQUIRE(malignant_cluster != fibro_cluster);
  for (int i = 0; i < 4; ++i) {
    REQUIRE_EQ(result.clusters[static_cast<std::size_t>(i)], malignant_cluster);
  }
  for (int i = 4; i < 8; ++i) {
    REQUIRE_EQ(result.clusters[static_cast<std::size_t>(i)], fibro_cluster);
  }
}

TEST_CASE("Cell clustering can subsample the clustering cell set") {
  const auto table = make_cell_clustering_fixture();
  CellClusteringOptions options;
  options.min_molecules = 10;
  options.min_genes = 1;
  options.cells_max = 4;
  options.n_variable_genes = 6;
  options.pca_dims = 4;
  options.graph_k = 2;
  options.cluster_resolution = 1.0;
  options.umap_neighbors = 3;
  options.umap_epochs = 30;
  options.seed = 7;

  const auto result = cluster_cells(table, options, nullptr);
  REQUIRE_EQ(result.cells.size(), 4U);
  REQUIRE_EQ(result.clusters.size(), 4U);
  REQUIRE_EQ(result.umap.rows(), 4);
  REQUIRE_EQ(result.umap.cols(), 2);
}
