#include <gtest/gtest.h>
#include <complex>
#include <vector>
#include "matrix.hpp"
#include "test_utils.hpp"

template<typename T>
class MatrixOperators : public ::testing::Test {

protected:

    void testMatrixDimensions(linalg::Matrix<T>& m, size_t expected_rows, size_t expected_cols) {
        EXPECT_EQ(m.getRows(), expected_rows);
        EXPECT_EQ(m.getCols(), expected_cols);
    }

    void testMatrixElements(linalg::Matrix<T>& m, std::vector<T>& expected_values) {
        for ( size_t i = 0; i < m.getRows(); ++i ) {
            for ( size_t j = 0; j < m.getCols(); ++j ) {
                expect_type_eq(m(i, j), expected_values[i * m.getCols() + j]);
            }
        }
    }

};

TYPED_TEST_SUITE(MatrixOperators, MatrixTypes);

TYPED_TEST(MatrixOperators, PlusOperator) {
    const size_t init_dim{2};
    std::vector<TypeParam> init_vec1{1, 2, 3, 4};
    std::vector<TypeParam> init_vec2{5, 6, 7, 8};
    std::vector<TypeParam> expected_vec{6, 8, 10, 12};
    linalg::Matrix<TypeParam> test_matrix1(init_dim, init_vec1);
    linalg::Matrix<TypeParam> test_matrix2(init_dim, init_vec2);
    linalg::Matrix<TypeParam> result_matrix;
    result_matrix = test_matrix1 + test_matrix2;
    this->testMatrixElements(result_matrix, expected_vec);
}
TYPED_TEST(MatrixOperators, MinusOperator) {
    const size_t init_dim{2};
    std::vector<TypeParam> init_vec1{1, 8, 2, 9};
    std::vector<TypeParam> init_vec2{3, 6, 4, 7};
    std::vector<TypeParam> expected_vec{-2, 2, -2, 2};
    linalg::Matrix<TypeParam> test_matrix1(init_dim, init_vec1);
    linalg::Matrix<TypeParam> test_matrix2(init_dim, init_vec2);
    linalg::Matrix<TypeParam> result_matrix;
    result_matrix = test_matrix1 - test_matrix2;
    this->testMatrixElements(result_matrix, expected_vec);
}
TYPED_TEST(MatrixOperators, UnaryMinusOperator) {
    const size_t init_dim{2};
    std::vector<TypeParam> init_vec{1, 2, 3, 4};
    std::vector<TypeParam> expected_vec{-1, -2, -3, -4};
    linalg::Matrix<TypeParam> test_matrix(init_dim, init_vec);
    linalg::Matrix<TypeParam> result_matrix;
    linalg::Matrix<TypeParam> expected_matrix(init_dim, expected_vec);
    result_matrix = -test_matrix;
    this->testMatrixElements(result_matrix, expected_vec);
}
TYPED_TEST(MatrixOperators, MultiplicationOperator) {
    const size_t init_row1{2}, init_row2{3}, init_col1{3}, init_col2{2}, expected_dim{2};
    std::vector<TypeParam> init_vec1{1, 2, 3, 4, 5, 6};
    std::vector<TypeParam> init_vec2{6, 5, 4, 3, 2, 1};
    std::vector<TypeParam> expected_vec{20, 14, 56, 41};
    linalg::Matrix<TypeParam> test_matrix1(init_row1, init_col1, init_vec1);
    linalg::Matrix<TypeParam> test_matrix2(init_row2, init_col2, init_vec2);
    linalg::Matrix<TypeParam> result_matrix;
    result_matrix = test_matrix1 * test_matrix2;
    this->testMatrixDimensions(result_matrix, expected_dim, expected_dim);
    this->testMatrixElements(result_matrix, expected_vec);
}
