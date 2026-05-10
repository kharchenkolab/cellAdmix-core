#include "celladmix/ncv.hpp"

#include "test_framework.hpp"

using namespace celladmix;

TEST_CASE("NCV counts nearest transcript genes within cell only") {
  TranscriptTable table;
  table.gene_key = {"A", "B", "B", "C", "C"};
  table.orig_cell_id = {"cell1", "cell1", "cell1", "cell1", "cell2"};
  table.x = {0.0, 1.0, 2.0, 8.0, 0.5};
  table.y = {0.0, 0.0, 0.0, 0.0, 0.0};
  table.z = {0.0, 0.0, 0.0, 0.0, 0.0};
  table.finalize();

  NcvOptions options;
  options.k = 2;
  options.include_self = false;
  options.within_cell = true;
  options.query_indices = {0, 1};

  const auto ncv = build_ncv_matrix(table, options);
  REQUIRE_EQ(ncv.rows(), 2);
  REQUIRE_EQ(ncv.cols(), 3);

  const int gene_a = table.gene_index[0];
  const int gene_b = table.gene_index[1];
  const int gene_c = table.gene_index[3];

  REQUIRE_NEAR(ncv(0, gene_b), 2.0, 1e-9);
  REQUIRE_NEAR(ncv(0, gene_a), 0.0, 1e-9);
  REQUIRE_NEAR(ncv(0, gene_c), 0.0, 1e-9);

  REQUIRE_NEAR(ncv(1, gene_a), 1.0, 1e-9);
  REQUIRE_NEAR(ncv(1, gene_b), 1.0, 1e-9);
  REQUIRE_NEAR(ncv(1, gene_c), 0.0, 1e-9);

  const auto sparse = build_sparse_ncv_matrix(table, options);
  REQUIRE_EQ(sparse.rows(), ncv.rows());
  REQUIRE_EQ(sparse.cols(), ncv.cols());
  const auto dense_from_sparse = sparse.to_dense();
  for (int i = 0; i < ncv.rows(); ++i) {
    for (int j = 0; j < ncv.cols(); ++j) {
      REQUIRE_NEAR(dense_from_sparse(i, j), ncv(i, j), 1e-9);
    }
  }
}

TEST_CASE("Sparse NCV matrix preserves explicit query order across cells") {
  TranscriptTable table;
  table.gene_key = {"A", "B", "C", "D", "E", "F"};
  table.orig_cell_id = {"cell1", "cell2", "cell1", "cell2", "cell1", "cell2"};
  table.x = {0.0, 0.0, 1.0, 1.0, 2.0, 2.0};
  table.y = {0.0, 5.0, 0.0, 5.0, 0.0, 5.0};
  table.z = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  table.finalize();

  NcvOptions options;
  options.k = 2;
  options.include_self = true;
  options.within_cell = true;
  options.query_indices = {1, 0, 3, 2, 5, 4};

  const auto dense = build_ncv_matrix(table, options);
  const auto sparse = build_sparse_ncv_matrix(table, options);
  const auto dense_from_sparse = sparse.to_dense();
  REQUIRE_EQ(dense.rows(), dense_from_sparse.rows());
  REQUIRE_EQ(dense.cols(), dense_from_sparse.cols());
  for (int i = 0; i < dense.rows(); ++i) {
    for (int j = 0; j < dense.cols(); ++j) {
      REQUIRE_NEAR(dense(i, j), dense_from_sparse(i, j), 1e-9);
    }
  }
}

TEST_CASE("NCV honors z coordinates when selecting nearest neighbors") {
  TranscriptTable table;
  table.gene_key = {"A", "B", "C"};
  table.orig_cell_id = {"cell1", "cell1", "cell1"};
  table.x = {0.0, 0.0, 1.0};
  table.y = {0.0, 0.0, 0.0};
  table.z = {0.0, 5.0, 0.0};
  table.finalize();

  NcvOptions options;
  options.k = 1;
  options.include_self = false;
  options.within_cell = true;
  options.query_indices = {0};

  const auto ncv = build_ncv_matrix(table, options);
  REQUIRE_EQ(ncv.rows(), 1);
  REQUIRE_EQ(ncv.cols(), 3);
  REQUIRE_NEAR(ncv(0, table.gene_index[2]), 1.0, 1e-9);
  REQUIRE_NEAR(ncv(0, table.gene_index[1]), 0.0, 1e-9);
}

TEST_CASE("NCV factor projection matches between serial and threaded paths") {
  TranscriptTable table;
  table.gene_key = {"A", "B", "C", "A", "C", "B"};
  table.orig_cell_id = {"cell1", "cell1", "cell1", "cell2", "cell2", "cell2"};
  table.x = {0.0, 1.0, 2.0, 0.0, 1.0, 2.0};
  table.y = {0.0, 0.0, 0.0, 5.0, 5.0, 5.0};
  table.z = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  table.finalize();

  DenseMatrix h(2, static_cast<int>(table.num_genes()), 0.0);
  h(0, table.gene_index[0]) = 1.0;
  h(0, table.gene_index[1]) = 0.5;
  h(1, table.gene_index[1]) = 0.5;
  h(1, table.gene_index[2]) = 1.0;

  NcvOptions options;
  options.k = 2;
  options.include_self = true;
  options.within_cell = true;

  const auto serial = project_ncv_to_factors(table, h, options, 1);
  const auto threaded = project_ncv_to_factors(table, h, options, 4);
  REQUIRE_EQ(serial.rows(), threaded.rows());
  REQUIRE_EQ(serial.cols(), threaded.cols());
  for (int i = 0; i < serial.rows(); ++i) {
    for (int j = 0; j < serial.cols(); ++j) {
      REQUIRE_NEAR(serial(i, j), threaded(i, j), 1e-9);
    }
  }
}

TEST_CASE("KL NCV factor projection matches between serial and threaded paths") {
  TranscriptTable table;
  table.gene_key = {"A", "B", "C", "A", "C", "B"};
  table.orig_cell_id = {"cell1", "cell1", "cell1", "cell2", "cell2", "cell2"};
  table.x = {0.0, 1.0, 2.0, 0.0, 1.0, 2.0};
  table.y = {0.0, 0.0, 0.0, 5.0, 5.0, 5.0};
  table.z = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  table.finalize();

  DenseMatrix h(2, static_cast<int>(table.num_genes()), 0.0);
  h(0, table.gene_index[0]) = 1.0;
  h(0, table.gene_index[1]) = 0.5;
  h(1, table.gene_index[1]) = 0.5;
  h(1, table.gene_index[2]) = 1.0;

  NcvOptions options;
  options.k = 2;
  options.include_self = true;
  options.within_cell = true;

  const auto serial = project_ncv_to_factors_kl(table, h, options, 1);
  const auto threaded = project_ncv_to_factors_kl(table, h, options, 4);
  REQUIRE_EQ(serial.rows(), threaded.rows());
  REQUIRE_EQ(serial.cols(), threaded.cols());
  for (int i = 0; i < serial.rows(); ++i) {
    for (int j = 0; j < serial.cols(); ++j) {
      REQUIRE_NEAR(serial(i, j), threaded(i, j), 1e-9);
    }
  }
}

TEST_CASE("KL sparse-row projection matches table-based projection on the same sampled NCVs") {
  TranscriptTable table;
  constexpr int n_cells = 8;
  constexpr int per_cell = 30;
  constexpr int n_genes = 15;
  for (int cell = 0; cell < n_cells; ++cell) {
    for (int i = 0; i < per_cell; ++i) {
      table.gene_key.push_back("G" + std::to_string((7 * i + 3 * cell) % n_genes));
      table.orig_cell_id.push_back("cell" + std::to_string(cell));
      table.x.push_back(cell * 100.0 + ((13 * i) % 17) * 0.3);
      table.y.push_back(cell * 10.0 + ((5 * i) % 11) * 0.2);
      table.z.push_back(((3 * i + cell) % 7) * 0.1);
    }
  }
  table.finalize();

  DenseMatrix h(5, static_cast<int>(table.num_genes()), 0.0);
  for (int factor = 0; factor < h.rows(); ++factor) {
    for (int gene = 0; gene < h.cols(); ++gene) {
      h(factor, gene) =
          static_cast<double>(((factor + 2) * (gene + 3)) % 17 + 1) / 17.0;
    }
  }

  NcvOptions options;
  options.k = 21;
  options.include_self = true;
  options.within_cell = true;
  for (int i = 0; i < static_cast<int>(table.size()); i += 2) {
    options.query_indices.push_back(i);
  }

  const auto sparse = build_sparse_ncv_matrix(table, options);
  const auto via_sparse = project_rows_to_factors_kl(sparse, h, 1);
  const auto via_table = project_ncv_to_factors_kl(table, h, options, 1);

  REQUIRE_EQ(via_sparse.rows(), via_table.rows());
  REQUIRE_EQ(via_sparse.cols(), via_table.cols());
  for (int i = 0; i < via_sparse.rows(); ++i) {
    for (int j = 0; j < via_sparse.cols(); ++j) {
      REQUIRE_NEAR(via_sparse(i, j), via_table(i, j), 1e-9);
    }
  }
}
