#include "matrix.h"

namespace linalg {

    Matrix::Matrix() : rows{0}, cols{0} {}

    Matrix::Matrix(size_t dim) : rows{dim}, cols{dim}, elem(dim * dim) {}

    Matrix::Matrix(size_t rows, size_t cols) : rows{rows}, cols{cols}, elem(rows * cols) {}

    Matrix::Matrix(size_t rows, size_t cols, const std::vector<double>& elem) : rows{rows}, cols{cols}, elem{elem} {
        if (elem.size() != rows * cols) {
            throw std::invalid_argument("Matrix dimension incompatible with initializing vector.");
        }
    }

    size_t Matrix::nRows() const {
        return rows;
    }

    size_t Matrix::nCols() const {
        return cols;
    }

    double& Matrix::operator()(int row, int col) {
        return const_cast<double&>(static_cast<const Matrix&>(*this)(row, col));
    }

    const double& Matrix::operator()(int row, int col) const {
        if ( (row >= rows) || (col >= cols) ) {
            throw std::out_of_range("Index out of range.");
        }
        return elem[row * cols + col];
    }

} // linalg
