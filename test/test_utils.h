#ifndef MATCAL_TEST_UTILS_H
#define MATCAL_TEST_UTILS_H

#include <complex>
#include <vector>
#include "matrix.h"

template<typename T>
void EXPECT_TYPE_EQ(const T& actual, const T& expected) {
    EXPECT_EQ(actual, expected);
}

template<>
inline void EXPECT_TYPE_EQ(const float& actual, const float& expected) {
    EXPECT_FLOAT_EQ(actual, expected);
}

template<>
inline void EXPECT_TYPE_EQ(const double& actual, const double& expected) {
    EXPECT_DOUBLE_EQ(actual, expected);
}

template<>
inline void EXPECT_TYPE_EQ(const std::complex<float>& actual, const std::complex<float>& expected) {
    EXPECT_FLOAT_EQ(actual.real(), expected.real());
    EXPECT_FLOAT_EQ(actual.imag(), expected.imag());
}

template<>
inline void EXPECT_TYPE_EQ(const std::complex<double>& actual, const std::complex<double>& expected) {
    EXPECT_DOUBLE_EQ(actual.real(), expected.real());
    EXPECT_DOUBLE_EQ(actual.imag(), expected.imag());
}

template<typename T>
class MatrixTestBase : public ::testing::Test {

protected:
    void testMatrixDimensions(matcal::Matrix<T>& matrixTest, size_t expectedRows, size_t expectedCols) {
        EXPECT_EQ(matrixTest.getRows(), expectedRows);
        EXPECT_EQ(matrixTest.getCols(), expectedCols);
    }

    void testMatrixElements(matcal::Matrix<T>& matrixTest, std::vector<T>& expectedValues) {
        for (size_t i = 0; i < matrixTest.getRows(); ++i) {
            for (size_t j = 0; j < matrixTest.getCols(); ++j) {
                EXPECT_TYPE_EQ(matrixTest(i, j), expectedValues[i * matrixTest.getCols() + j]);
            }
        }
    }

}; // MatrixTestBase

#endif //MATCAL_TEST_UTILS_H
