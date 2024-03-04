#ifndef MATCAL_MATRIX_HPP
#define MATCAL_MATRIX_HPP

#include <stdexcept>
#include <vector>
#include "../../src/traits/is_complex.hpp"

namespace matcal {

    template<typename T>
    class Matrix {

    public:
        Matrix();
        Matrix(const Matrix&) = default;
        Matrix(Matrix&&) noexcept = default;
        Matrix& operator=(const Matrix&) = default;
        Matrix& operator=(Matrix&&) noexcept = default;
        ~Matrix() = default;

        explicit Matrix(size_t dim);
        Matrix(size_t dim, const std::vector<T>& elem);
        Matrix(size_t dim, std::vector<T>&& elem);
        Matrix(size_t rows, size_t cols);
        Matrix(size_t rows, size_t cols, const std::vector<T>& elem);
        Matrix(size_t rows, size_t cols, std::vector<T>&& elem);

        [[nodiscard]] size_t getRows() const;
        [[nodiscard]] size_t getCols() const;

        void setDimensions(size_t newRowCount, size_t newColCount);
        void setElements(const T& value);
        void setElements(const std::vector<T>& values);
        void setElements(std::vector<T>&& values);

        T& operator()(size_t row, size_t col);
        const T& operator()(size_t row, size_t col) const;

        bool operator==(const Matrix<T>& other) const;
        bool operator!=(const Matrix<T>& other) const;

        Matrix operator-() const &;
        Matrix operator-() &&;

        Matrix& operator+=(const Matrix<T>& other);
        Matrix& operator-=(const Matrix<T>& other);
        Matrix& operator*=(const Matrix<T>& other);
        Matrix& operator*=(const T& scalar);
        Matrix& operator/=(const T& scalar);

        template<typename U> friend Matrix<U> operator+(const Matrix<U>& lhs, const Matrix<U>& rhs);
        template<typename U> friend Matrix<U> operator+(Matrix<U>&& lhs, const Matrix<U>& rhs);
        template<typename U> friend Matrix<U> operator+(const Matrix<U>& lhs, Matrix<U>&& rhs);
        template<typename U> friend Matrix<U> operator+(Matrix<U>&& lhs, Matrix<U>&& rhs);

        template<typename U> friend Matrix<U> operator-(const Matrix<U>& lhs, const Matrix<U>& rhs);
        template<typename U> friend Matrix<U> operator-(Matrix<U>&& lhs, const Matrix<U>& rhs);
        template<typename U> friend Matrix<U> operator-(const Matrix<U>& lhs, Matrix<U>&& rhs);
        template<typename U> friend Matrix<U> operator-(Matrix<U>&& lhs, Matrix<U>&& rhs);

        template<typename U> friend Matrix<U> operator*(const Matrix<U>& lhs, const Matrix<U>& rhs);
        template<typename U> friend Matrix<U> operator*(Matrix<U>&& lhs, const Matrix<U>& rhs);
        template<typename U> friend Matrix<U> operator*(const Matrix<U>& lhs, Matrix<U>&& rhs);
        template<typename U> friend Matrix<U> operator*(Matrix<U>&& lhs, Matrix<U>&& rhs);

        template<typename U> friend Matrix<U> operator*(const Matrix<U>& lhs, const U& rhs);
        template<typename U> friend Matrix<U> operator*(Matrix<U>&& lhs, const U& rhs);
        template<typename U> friend Matrix<U> operator*(const U& lhs, const Matrix<U>& rhs);
        template<typename U> friend Matrix<U> operator*(const U& lhs, Matrix<U>&& rhs);

        template<typename U> friend Matrix<U> operator/(const Matrix<U>& lhs, const U& rhs);
        template<typename U> friend Matrix<U> operator/(Matrix<U>&& lhs, const U& rhs);

        Matrix& hadamardProduct(const Matrix<T>& other);
        template<typename U> friend Matrix<U> hadamardProduct(const Matrix<U>& lhs, const Matrix<U>& rhs);
        template<typename U> friend Matrix<U> hadamardProduct(Matrix<U>&& lhs, const Matrix<U>& rhs);
        template<typename U> friend Matrix<U> hadamardProduct(const Matrix<U>& lhs, Matrix<U>&& rhs);
        template<typename U> friend Matrix<U> hadamardProduct(Matrix<U>&& lhs, Matrix<U>&& rhs);

        Matrix& kroneckerProduct(const Matrix<T>& other);
        template<typename U> friend Matrix<U> kroneckerProduct(const Matrix<U>& lhs, const Matrix<U>& rhs);
        template<typename U> friend Matrix<U> kroneckerProduct(Matrix<U>&& lhs, const Matrix<U>& rhs);
        template<typename U> friend Matrix<U> kroneckerProduct(const Matrix<U>& lhs, Matrix<U>&& rhs);
        template<typename U> friend Matrix<U> kroneckerProduct(Matrix<U>&& lhs, Matrix<U>&& rhs);

        Matrix transpose() const &;
        Matrix transpose() &&;

        Matrix hermitianConjugate() const &;
        Matrix hermitianConjugate() &&;

    private:
        size_t rows;
        size_t cols;
        std::vector<T> elem;

        void validateMatrixDimensions() const;
        void validateIndex(size_t row, size_t col) const;
        bool compareMatrixDimensions(const Matrix& other) const;

        void negateInPlace();
        void transposeInPlace();
        void complexConjugateInPlace();

    }; // Matrix

    template<typename T>
    Matrix<T>::Matrix() : rows(0), cols(0), elem() {}

    template<typename T>
    Matrix<T>::Matrix(size_t dim) : rows{dim}, cols{dim}, elem(dim * dim) {}

    template<typename T>
    Matrix<T>::Matrix(size_t dim, const std::vector<T>& elem) : rows{dim}, cols{dim}, elem{elem} {
        validateMatrixDimensions();
    }

    template<typename T>
    Matrix<T>::Matrix(size_t dim, std::vector<T>&& elem) : rows{dim}, cols{dim}, elem{std::move(elem)} {
        validateMatrixDimensions();
    }

    template<typename T>
    Matrix<T>::Matrix(size_t rows, size_t cols) : rows{rows}, cols{cols}, elem(rows * cols) {}

    template<typename T>
    Matrix<T>::Matrix(size_t rows, size_t cols, const std::vector<T>& elem) : rows{rows}, cols{cols}, elem{elem} {
        validateMatrixDimensions();
    }

    template<typename T>
    Matrix<T>::Matrix(size_t rows, size_t cols, std::vector<T>&& elem) : rows{rows}, cols{cols}, elem{std::move(elem)} {
        validateMatrixDimensions();
    }

    template<typename T>
    size_t Matrix<T>::getRows() const {
        return rows;
    }

    template<typename T>
    size_t Matrix<T>::getCols() const {
        return cols;
    }

    template<typename T>
    void Matrix<T>::setDimensions(size_t newRowCount, size_t newColCount) {
        rows = newRowCount;
        cols = newColCount;
        validateMatrixDimensions();
    }

    template<typename T>
    void Matrix<T>::setElements(const T& value) {
        std::fill(elem.begin(), elem.end(), value);
    }

    template<typename T>
    void Matrix<T>::setElements(const std::vector<T>& values) {
        elem = values;
        validateMatrixDimensions();
    }

    template<typename T>
    void Matrix<T>::setElements(std::vector<T>&& values) {
        elem = std::move(values);
        validateMatrixDimensions();
    }

    template<typename T>
    T& Matrix<T>::operator()(size_t row, size_t col) {
        return const_cast<T&>(static_cast<const Matrix&>(*this)(row, col));
    }

    template<typename T>
    const T& Matrix<T>::operator()(size_t row, size_t col) const {
        validateIndex(row, col);
        return elem[row * cols + col];
    }

    template<typename T>
    bool Matrix<T>::operator==(const Matrix<T>& other) const {
        return ( (rows == other.rows) && (cols == other.cols) && (elem == other.elem) );
    }

    template<typename T>
    bool Matrix<T>::operator!=(const Matrix<T>& other) const {
        return ( (rows != other.rows) || (cols != other.cols) || (elem != other.elem) );
    }

    template<typename T>
    Matrix<T> Matrix<T>::operator-() const & {
        Matrix<T> result{*this};
        result.negateInPlace();
        return result;
    }

    template<typename T>
    Matrix<T> Matrix<T>::operator-() && {
        negateInPlace();
        return std::move(*this);
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& other) {
        if ( compareMatrixDimensions(other) ) {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
        for ( size_t i = 0; i < elem.size(); ++i ) {
            elem[i] += other.elem[i];
        }
        return *this;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& other) {
        return *this += -other;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& other) {
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
        *this = std::move(result);
        return *this;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator*=(const T& scalar) {
        for ( size_t i = 0; i < elem.size(); ++i ) {
            elem[i] *= scalar;
        }
        return *this;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator/=(const T& scalar) {
        for ( size_t i = 0; i < elem.size(); ++i ) {
            elem[i] /= scalar;
        }
        return *this;
    }

    template<typename T>
    Matrix<T> operator+(const Matrix<T>& lhs, const Matrix<T>& rhs) {
        Matrix<T> result{lhs};
        result += rhs;
        return result;
    }

    template<typename T>
    Matrix<T> operator+(Matrix<T>&& lhs, const Matrix<T>& rhs) {
        lhs += rhs;
        return std::move(lhs);
    }

    template<typename T>
    Matrix<T> operator+(const Matrix<T>& lhs, Matrix<T>&& rhs) {
        return rhs + lhs;
    }

    template<typename T>
    Matrix<T> operator+(Matrix<T>&& lhs, Matrix<T>&& rhs) {
        return std::move(lhs) + rhs;
    }

    template<typename T>
    Matrix<T> operator-(const Matrix<T>& lhs, const Matrix<T>& rhs) {
        return lhs + (-rhs);
    }

    template<typename T>
    Matrix<T> operator-(Matrix<T>&& lhs, const Matrix<T>& rhs) {
        return std::move(lhs) + (-rhs);
    }

    template<typename T>
    Matrix<T> operator-(const Matrix<T>& lhs, Matrix<T>&& rhs) {
        return lhs + (-std::move(rhs));
    }

    template<typename T>
    Matrix<T> operator-(Matrix<T>&& lhs, Matrix<T>&& rhs) {
        return std::move(lhs) + (-std::move(rhs));
    }

    template<typename T>
    Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs) {
        Matrix<T> result{lhs};
        result *= rhs;
        return result;
    }

    template<typename T>
    Matrix<T> operator*(Matrix<T>&& lhs, const Matrix<T>& rhs) {
        lhs *= rhs;
        return std::move(lhs);
    }

    template<typename T>
    Matrix<T> operator*(const Matrix<T>& lhs, Matrix<T>&& rhs) {
        rhs = lhs * rhs;
        return std::move(rhs);
    }

    template<typename T>
    Matrix<T> operator*(Matrix<T>&& lhs, Matrix<T>&& rhs) {
        return std::move(lhs) * rhs;
    }

    template<typename T>
    Matrix<T> operator*(const Matrix<T>& lhs, const T& rhs) {
        Matrix<T> result{lhs};
        result *= rhs;
        return result;
    }

    template<typename T>
    Matrix<T> operator*(Matrix<T>&& lhs, const T& rhs) {
        lhs *= rhs;
        return std::move(lhs);
    }

    template<typename T>
    Matrix<T> operator*(const T& lhs, const Matrix<T>& rhs) {
        return rhs * lhs;
    }

    template<typename T>
    Matrix<T> operator*(const T& lhs, Matrix<T>&& rhs) {
        return std::move(rhs) * lhs;
    }

    template<typename T>
    Matrix<T> operator/(const Matrix<T>& lhs, const T& rhs) {
        Matrix<T> result{lhs};
        result /= rhs;
        return result;
    }

    template<typename T>
    Matrix<T> operator/(Matrix<T>&& lhs, const T& rhs) {
        lhs /= rhs;
        return std::move(lhs);
    }

    template<typename T>
    Matrix<T>& Matrix<T>::hadamardProduct(const Matrix<T>& other) {
        if ( compareMatrixDimensions(other) ) {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
        for ( size_t i = 0; i < elem.size(); ++i ) {
            elem[i] *= other.elem[i];
        }
        return *this;
    }

    template<typename T>
    Matrix<T> hadamardProduct(const Matrix<T>& lhs, const Matrix<T>& rhs) {
        Matrix<T> result{lhs};
        result.hadamardProduct(rhs);
        return result;
    }

    template<typename T>
    Matrix<T> hadamardProduct(Matrix<T>&& lhs, const Matrix<T>& rhs) {
        lhs.hadamardProduct(rhs);
        return std::move(lhs);
    }

    template<typename T>
    Matrix<T> hadamardProduct(const Matrix<T>& lhs, Matrix<T>&& rhs) {
        rhs.hadamardProduct(lhs);
        return std::move(rhs);
    }

    template<typename T>
    Matrix<T> hadamardProduct(Matrix<T>&& lhs, Matrix<T>&& rhs) {
        lhs.hadamardProduct(rhs);
        return std::move(lhs);
    }

    template<typename T>
    Matrix<T>& Matrix<T>::kroneckerProduct(const Matrix<T>& other) {
        Matrix<T> result(rows * other.rows, cols * other.cols);
        size_t resultIndex{0};
        for ( size_t i = 0; i < rows; ++i ) {
            for ( size_t k = 0; k < other.rows; ++k ) {
                for ( size_t j = 0; j < cols; ++ j ) {
                    T this_ij = elem[i * cols + j];
                    for ( size_t l = 0; l < other.cols; ++l ) {
                        result.elem[resultIndex++] = this_ij * other.elem[k * other.cols + l];
                    }
                }
            }
        }
        *this = std::move(result);
        return *this;
    }

    template<typename T>
    Matrix<T> kroneckerProduct(const Matrix<T>& lhs, const Matrix<T>& rhs) {
        Matrix<T> result{lhs};
        result.kroneckerProduct(rhs);
        return result;
    }

    template<typename T>
    Matrix<T> kroneckerProduct(Matrix<T>&& lhs, const Matrix<T>& rhs) {
        lhs.kroneckerProduct(rhs);
        return std::move(lhs);
    }

    template<typename T>
    Matrix<T> kroneckerProduct(const Matrix<T>& lhs, Matrix<T>&& rhs) {
        rhs = kroneckerProduct(lhs, rhs);
        return std::move(rhs);
    }

    template<typename T>
    Matrix<T> kroneckerProduct(Matrix<T>&& lhs, Matrix<T>&& rhs) {
        lhs.kroneckerProduct(rhs);
        return std::move(lhs);
    }

    template<typename T>
    Matrix<T> Matrix<T>::transpose() const & {
        Matrix<T> result(cols, rows);
        for ( size_t i = 0; i < rows; ++i ) {
            for ( size_t j = 0; j < cols; ++j ) {
                result.elem[j * result.cols + i] = elem[i * cols + j];
            }
        }
        return result;
    }

    template<typename T>
    Matrix<T> Matrix<T>::transpose() && {
        transposeInPlace();
        return std::move(*this);
    }

    template<typename T>
    Matrix<T> Matrix<T>::hermitianConjugate() const & {
        Matrix<T> result{*this};
        result.transposeInPlace();
        result.complexConjugateInPlace();
        return result;
    }

    template<typename T>
    Matrix<T> Matrix<T>::hermitianConjugate() && {
        transposeInPlace();
        complexConjugateInPlace();
        return std::move(*this);
    }


    template<typename T>
    void Matrix<T>::validateMatrixDimensions() const {
        if ( elem.size() != rows * cols ) {
            throw std::invalid_argument("Matrix dimension incompatible with initializing vector.");
        }
    }

    template<typename T>
    void Matrix<T>::validateIndex(size_t row, size_t col) const {
        if ( (rows <= row) || (cols <= col) ) {
            throw std::out_of_range("Index out of range.");
        }
    }

    template<typename T>
    bool Matrix<T>::compareMatrixDimensions(const Matrix<T> &other) const {
        return ( (rows != other.rows) || (cols != other.cols) );
    }

    template<typename T>
    void Matrix<T>::negateInPlace() {
        for ( size_t i = 0; i < elem.size(); ++i ) {
            elem[i] = -elem[i];
        }
    }

    template<typename T>
    void Matrix<T>::transposeInPlace() {
        std::vector<bool> visited(elem.size(), false);
        size_t currIndex;
        for ( size_t start = 1; start < elem.size() - 1; ++start ) {
            if ( visited[start] ) {
                continue;
            }
            currIndex = start;
            do {
                visited[currIndex] = true;
                currIndex = rows * currIndex % (rows * cols - 1);
                if ( !visited[currIndex] ) {
                    std::swap(elem[start], elem[currIndex]);
                }
            } while ( currIndex != start );
        }
        std::swap(rows, cols);
    }

    template<typename T>
    void Matrix<T>::complexConjugateInPlace() {
        if constexpr ( is_complex<T>::value ) {
            for ( size_t i = 0; i < elem.size(); ++i ) {
                elem[i] = std::conj(elem[i]);
            }
        }
    }

} // matcal

#endif //MATCAL_MATRIX_HPP
