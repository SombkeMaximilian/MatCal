#ifndef MATCAL_MATRIX_HPP
#define MATCAL_MATRIX_HPP

#include <vector>
#include <stdexcept>

namespace linalg {

    template<typename T>
    class Matrix {
    private:
        std::vector<T> elem;
        size_t rows;
        size_t cols;

    public:
        Matrix() : rows{0}, cols{0}, elem() {}

        Matrix(size_t dim) : rows{dim}, cols{dim}, elem(dim * dim) {}

        Matrix(size_t dim, const std::vector<T>& elem) : rows{dim}, cols{dim}, elem{elem}  {
            if ( elem.size() != dim * dim ) {
                throw std::invalid_argument("Matrix dimension incompatible with initializing vector.");
            }
        }

        Matrix(size_t rows, size_t cols) : rows{rows}, cols{cols}, elem(rows * cols) {}

        Matrix(size_t rows, size_t cols, const std::vector<T>& elem) : rows{rows}, cols{cols}, elem{elem} {
            if ( elem.size() != rows * cols ) {
                throw std::invalid_argument("Matrix dimension incompatible with initializing vector.");
            }
        }

        size_t nRows() const {
            return rows;
        }

        size_t nCols() const {
            return cols;
        }

        T& operator()(int row, int col) {
            return const_cast<T&>(static_cast<const Matrix&>(*this)(row, col));
        }

        const T& operator()(int row, int col) const {
            if ( (row >= rows) || (col >= cols) ) {
                throw std::out_of_range("Index out of range.");
            }

            return elem[row * cols + col];
        }

        Matrix operator-() const {
            Matrix<T> result(rows, cols);

            for ( size_t i = 0; i < elem.size(); ++i ) {
                result.elem[i] = -elem[i];
            }

            return result;
        }

        Matrix& operator+=(const Matrix<T>& m) {
            if ( (m.nRows() != rows) || (m.nCols() != cols) ) {
                throw std::invalid_argument("Matrix dimensions do not match.");
            }

            for ( size_t i = 0; i < elem.size(); ++i ) {
                elem[i] += m.elem[i];
            }

            return *this;
        }

        Matrix& operator-=(const Matrix<T>& m) {
            return *this += -m;
        }

    }; // Matrix

} // linalg

#endif //MATCAL_MATRIX_HPP
