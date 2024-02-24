#ifndef MATCAL_MATRIX_HPP
#define MATCAL_MATRIX_HPP

#include <vector>
#include <stdexcept>

namespace linalg {

    template<typename T>
    class Matrix {
    public:
        Matrix();
        Matrix(const Matrix&) = default;
        Matrix(Matrix&&) noexcept = default;
        Matrix& operator=(const Matrix&) = default;
        Matrix& operator=(Matrix&&) noexcept = default;

        explicit Matrix(size_t dim);
        Matrix(size_t dim, const std::vector<T>& elem);
        Matrix(size_t dim, std::vector<T>&& elem);
        Matrix(size_t rows, size_t cols);
        Matrix(size_t rows, size_t cols, const std::vector<T>& elem);
        Matrix(size_t rows, size_t cols, std::vector<T>&& elem);

        [[nodiscard]] size_t nRows() const;
        [[nodiscard]] size_t nCols() const;

        T& operator()(size_t row, size_t col);
        const T& operator()(size_t row, size_t col) const;

        Matrix& operator+=(const Matrix<T>& other);
        Matrix& operator-=(const Matrix<T>& other);
        Matrix& operator*=(const Matrix<T>& other);

        Matrix operator+(const Matrix<T>& other) const;
        Matrix operator-(const Matrix<T>& other) const;
        Matrix operator-() const;
        Matrix operator*(const Matrix<T>& other) const;

        bool operator==(const Matrix<T>& other) const;
        bool operator!=(const Matrix<T>& other) const;

    private:
        std::vector<T> elem;
        size_t rows;
        size_t cols;

    }; // Matrix

    template<typename T>
    Matrix<T>::Matrix() : rows(0), cols(0), elem() {}

    template<typename T>
    Matrix<T>::Matrix(size_t dim) : rows{dim}, cols{dim}, elem(dim * dim) {}

    template<typename T>
    Matrix<T>::Matrix(size_t dim, const std::vector<T> &elem) : rows{dim}, cols{dim}, elem{elem}  {
        if ( elem.size() != dim * dim ) {
            throw std::invalid_argument("Matrix dimension incompatible with initializing vector.");
        }
    }

    template<typename T>
    Matrix<T>::Matrix(size_t rows, size_t cols) : rows{rows}, cols{cols}, elem(rows * cols) {}

    template<typename T>
    Matrix<T>::Matrix(size_t dim, std::vector<T> &&elem) : rows{dim}, cols{dim}, elem{std::move(elem)}  {
        if ( elem.size() != dim * dim ) {
            throw std::invalid_argument("Matrix dimension incompatible with initializing vector.");
        }
    }

    template<typename T>
    Matrix<T>::Matrix(size_t rows, size_t cols, const std::vector<T> &elem) : rows{rows}, cols{cols}, elem{elem} {
        if ( elem.size() != rows * cols ) {
            throw std::invalid_argument("Matrix dimension incompatible with initializing vector.");
        }
    }

    template<typename T>
    Matrix<T>::Matrix(size_t rows, size_t cols, std::vector<T> &&elem) : rows{rows}, cols{cols}, elem{std::move(elem)} {
        if ( elem.size() != rows * cols ) {
            throw std::invalid_argument("Matrix dimension incompatible with initializing vector.");
        }
    }

    template<typename T>
    size_t Matrix<T>::nRows() const {
        return rows;
    }

    template<typename T>
    size_t Matrix<T>::nCols() const {
        return cols;
    }

    template<typename T>
    T &Matrix<T>::operator()(size_t row, size_t col) {
        return const_cast<T&>(static_cast<const Matrix&>(*this)(row, col));
    }

    template<typename T>
    const T &Matrix<T>::operator()(size_t row, size_t col) const {
        if ( (rows <= row) || (cols <= col) ) {
            throw std::out_of_range("Index out of range.");
        }
        return elem[row * cols + col];
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator+=(const Matrix<T> &other) {
        if ( (rows != other.rows) || (cols != other.cols) ) {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
        for ( size_t i = 0; i < elem.size(); ++i ) {
            elem[i] += other.elem[i];
        }
        return *this;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator-=(const Matrix<T> &other) {
        return *this += -other;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator*=(const Matrix<T> &other) {
        Matrix<T> result{(*this) * other};
        *this = std::move(result);
        return *this;
    }

    template<typename T>
    Matrix<T> Matrix<T>::operator+(const Matrix<T> &other) const {
        Matrix<T> result{*this};
        result += other;
        return result;
    }

    template<typename T>
    Matrix<T> Matrix<T>::operator-(const Matrix<T> &other) const {
        return (*this) + (-other);
    }

    template<typename T>
    Matrix<T> Matrix<T>::operator-() const {
        Matrix<T> result(rows, cols);
        for ( size_t i = 0; i < elem.size(); ++i ) {
            result.elem[i] = -elem[i];
        }
        return result;
    }

    template<typename T>
    Matrix<T> Matrix<T>::operator*(const Matrix<T> &other) const {
        Matrix<T> result(rows, other.cols);
        if ( cols != other.rows ) {
            throw std::invalid_argument("Number of columns in A must match number of rows in B.");
        }
        for ( size_t i = 0; i < rows; ++i ) {
            for ( size_t k = 0; k < cols; ++k ) {
                for ( size_t j = 0; j < other.cols; ++j ) {
                    result.elem[i * result.cols + j] += elem[i * cols + k] * other.elem[k * other.cols + j];
                }
            }
        }
        return result;
    }

    template<typename T>
    bool Matrix<T>::operator==(const Matrix<T> &other) const {
        return ( (rows == other.rows) && (cols == other.cols) && (elem == other.elem) );
    }

    template<typename T>
    bool Matrix<T>::operator!=(const Matrix<T> &other) const {
        return ( (rows != other.rows) || (cols != other.cols) || (elem != other.elem) );
    }

} // linalg

#endif //MATCAL_MATRIX_HPP
