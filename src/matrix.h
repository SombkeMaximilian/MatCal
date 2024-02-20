#ifndef MATCAL_MATRIX_H
#define MATCAL_MATRIX_H

#include <vector>
#include <stdexcept>

namespace linalg {

    class Matrix {
    private:
        std::vector<double> elem;
        size_t rows;
        size_t cols;
    public:
        Matrix();
        Matrix(size_t dim);
        Matrix(size_t rows, size_t cols);
        Matrix(size_t rows, size_t cols, const std::vector<double>& elem);

        size_t nRows() const;
        size_t nCols() const;

        double& operator()(int row, int col);
        const double& operator()(int row, int col) const;
    }; // Matrix

} // linalg

#endif //MATCAL_MATRIX_H
