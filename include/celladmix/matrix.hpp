// Minimal dense and sparse matrix containers plus small linear-algebra helpers
// used throughout the core algorithms.

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace celladmix {

class DenseMatrix {
 public:
  DenseMatrix() = default;
  DenseMatrix(int rows, int cols, double value = 0.0)
      : rows_(rows), cols_(cols), data_(static_cast<std::size_t>(rows * cols), value) {}

  int rows() const { return rows_; }
  int cols() const { return cols_; }

  double& operator()(int row, int col) {
    return data_.at(static_cast<std::size_t>(row * cols_ + col));
  }

  double operator()(int row, int col) const {
    return data_.at(static_cast<std::size_t>(row * cols_ + col));
  }

  const std::vector<double>& data() const { return data_; }
  std::vector<double>& data() { return data_; }

  // Sum each row into a dense vector.
  std::vector<double> row_sums() const {
    std::vector<double> out(static_cast<std::size_t>(rows_), 0.0);
    for (int i = 0; i < rows_; ++i) {
      double sum = 0.0;
      for (int j = 0; j < cols_; ++j) {
        sum += (*this)(i, j);
      }
      out[static_cast<std::size_t>(i)] = sum;
    }
    return out;
  }

  // Sum each column into a dense vector.
  std::vector<double> col_sums() const {
    std::vector<double> out(static_cast<std::size_t>(cols_), 0.0);
    for (int j = 0; j < cols_; ++j) {
      double sum = 0.0;
      for (int i = 0; i < rows_; ++i) {
        sum += (*this)(i, j);
      }
      out[static_cast<std::size_t>(j)] = sum;
    }
    return out;
  }

  // Return the index of the largest value in one row.
  int row_argmax(int row) const {
    int best = 0;
    double best_value = -std::numeric_limits<double>::infinity();
    for (int j = 0; j < cols_; ++j) {
      if ((*this)(row, j) > best_value) {
        best_value = (*this)(row, j);
        best = j;
      }
    }
    return best;
  }

  // Normalize each row to sum to one when it has non-trivial mass.
  void normalize_rows(double eps = 1e-12) {
    for (int i = 0; i < rows_; ++i) {
      double sum = 0.0;
      for (int j = 0; j < cols_; ++j) {
        sum += (*this)(i, j);
      }
      if (sum <= eps) {
        continue;
      }
      for (int j = 0; j < cols_; ++j) {
        (*this)(i, j) /= sum;
      }
    }
  }

 private:
  int rows_ = 0;
  int cols_ = 0;
  std::vector<double> data_;
};

class SparseRowMatrix {
 public:
  SparseRowMatrix() = default;
  SparseRowMatrix(
      int rows,
      int cols,
      std::vector<int> indptr,
      std::vector<int> indices,
      std::vector<double> values)
      : rows_(rows),
        cols_(cols),
        indptr_(std::move(indptr)),
        indices_(std::move(indices)),
        values_(std::move(values)) {}

  int rows() const { return rows_; }
  int cols() const { return cols_; }
  int nnz() const { return static_cast<int>(values_.size()); }

  const std::vector<int>& indptr() const { return indptr_; }
  const std::vector<int>& indices() const { return indices_; }
  const std::vector<double>& values() const { return values_; }

  // Sum each sparse column into a dense vector.
  std::vector<double> col_sums() const {
    std::vector<double> out(static_cast<std::size_t>(cols_), 0.0);
    for (std::size_t p = 0; p < values_.size(); ++p) {
      out[static_cast<std::size_t>(indices_[p])] += values_[p];
    }
    return out;
  }

  // Materialize the sparse matrix as a dense matrix.
  DenseMatrix to_dense() const {
    DenseMatrix out(rows_, cols_, 0.0);
    for (int row = 0; row < rows_; ++row) {
      for (int p = indptr_[static_cast<std::size_t>(row)];
           p < indptr_[static_cast<std::size_t>(row + 1)];
           ++p) {
        out(row, indices_[static_cast<std::size_t>(p)]) = values_[static_cast<std::size_t>(p)];
      }
    }
    return out;
  }

 private:
  int rows_ = 0;
  int cols_ = 0;
  std::vector<int> indptr_;
  std::vector<int> indices_;
  std::vector<double> values_;
};

// Multiply two dense matrices in row-major storage.
inline DenseMatrix multiply(const DenseMatrix& a, const DenseMatrix& b) {
  if (a.cols() != b.rows()) {
    throw std::runtime_error("Matrix dimensions do not match for multiply");
  }
  DenseMatrix out(a.rows(), b.cols(), 0.0);
  for (int i = 0; i < a.rows(); ++i) {
    for (int k = 0; k < a.cols(); ++k) {
      const double aik = a(i, k);
      if (std::abs(aik) < 1e-15) {
        continue;
      }
      for (int j = 0; j < b.cols(); ++j) {
        out(i, j) += aik * b(k, j);
      }
    }
  }
  return out;
}

// Compute a * b^T without explicitly transposing b.
inline DenseMatrix transpose_multiply_right(const DenseMatrix& a, const DenseMatrix& b) {
  if (a.cols() != b.cols()) {
    throw std::runtime_error("Matrix dimensions do not match for transpose_multiply_right");
  }
  DenseMatrix out(a.rows(), b.rows(), 0.0);
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < b.rows(); ++j) {
      double sum = 0.0;
      for (int k = 0; k < a.cols(); ++k) {
        sum += a(i, k) * b(j, k);
      }
      out(i, j) = sum;
    }
  }
  return out;
}

// Compute the Frobenius distance between two dense matrices.
inline double frobenius_distance(const DenseMatrix& a, const DenseMatrix& b) {
  if (a.rows() != b.rows() || a.cols() != b.cols()) {
    throw std::runtime_error("Matrix dimensions do not match for frobenius_distance");
  }
  double sum = 0.0;
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < a.cols(); ++j) {
      const double diff = a(i, j) - b(i, j);
      sum += diff * diff;
    }
  }
  return std::sqrt(sum);
}

}  // namespace celladmix
