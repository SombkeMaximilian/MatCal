#include <gtest/gtest.h>
#include <complex>
#include <vector>
#include "matrix.hpp"
#include "test_utils.hpp"

template<typename T>
class MatrixConstructor : public ::testing::Test {

protected:
    size_t init_dim{3}, init_row{2}, init_col{3};
    std::vector<T> dim_init_vec;
    std::vector<T> rvalue_dim_init_vec;
    std::vector<T> rvalue_dim_expected;
    std::vector<T> row_col_init_vec;
    std::vector<T> rvalue_row_col_init_vec;
    std::vector<T> rvalue_row_col_expected;
    linalg::Matrix<T> default_matrix;
    linalg::Matrix<T> dim_matrix;
    linalg::Matrix<T> dim_init_vec_matrix;
    linalg::Matrix<T> dim_rvalue_vec_matrix;
    linalg::Matrix<T> row_col_matrix;
    linalg::Matrix<T> row_col_init_vec_matrix;
    linalg::Matrix<T> row_col_rvalue_vec_matrix;

    void SetUp() override {
        dim_init_vec = generateRandomVector<T>(init_dim * init_dim);
        rvalue_dim_init_vec = generateRandomVector<T>(init_dim * init_dim);
        rvalue_dim_expected = rvalue_dim_init_vec;
        row_col_init_vec = generateRandomVector<T>(init_row * init_col);
        rvalue_row_col_init_vec = generateRandomVector<T>(init_row * init_col);
        rvalue_row_col_expected = rvalue_row_col_init_vec;

        default_matrix = linalg::Matrix<T>();
        dim_matrix = linalg::Matrix<T>(init_dim);
        dim_init_vec_matrix = linalg::Matrix<T>(init_dim, dim_init_vec);
        dim_rvalue_vec_matrix = linalg::Matrix<T>(init_dim, std::move(rvalue_dim_init_vec));
        row_col_matrix = linalg::Matrix<T>(init_row, init_col);
        row_col_init_vec_matrix = linalg::Matrix<T>(init_row, init_col, row_col_init_vec);
        row_col_rvalue_vec_matrix = linalg::Matrix<T>(init_row, init_col, std::move(rvalue_row_col_init_vec));
    }

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

TYPED_TEST_SUITE(MatrixConstructor, MatrixTypes);

TYPED_TEST(MatrixConstructor, Default) {
    this->testMatrixDimensions(this->default_matrix, 0, 0);
}
TYPED_TEST(MatrixConstructor, Dim) {
    std::vector<TypeParam> expected_values(this->init_dim * this->init_dim, 0);
    this->testMatrixDimensions(this->dim_matrix, this->init_dim, this->init_dim);
    this->testMatrixElements(this->dim_matrix, expected_values);
}
TYPED_TEST(MatrixConstructor, DimInitializerVector) {
    this->testMatrixDimensions(this->dim_init_vec_matrix, this->init_dim, this->init_dim);
    this->testMatrixElements(this->dim_init_vec_matrix, this->dim_init_vec);
}
TYPED_TEST(MatrixConstructor, DimRValueInitializerVector) {
    this->testMatrixDimensions(this->dim_rvalue_vec_matrix, this->init_dim, this->init_dim);
    this->testMatrixElements(this->dim_rvalue_vec_matrix, this->rvalue_dim_expected);
    EXPECT_TRUE(this->rvalue_dim_init_vec.empty());
}
TYPED_TEST(MatrixConstructor, RowCol) {
    std::vector<TypeParam> expected_values(this->init_row * this->init_col, 0);
    this->testMatrixDimensions(this->row_col_matrix, this->init_row, this->init_col);
    this->testMatrixElements(this->row_col_matrix, expected_values);
}
TYPED_TEST(MatrixConstructor, RowColInitializerVector) {
    this->testMatrixDimensions(this->row_col_init_vec_matrix, this->init_row, this->init_col);
    this->testMatrixElements(this->row_col_init_vec_matrix, this->row_col_init_vec);
}
TYPED_TEST(MatrixConstructor, RowColRValueInitializerVector) {
    this->testMatrixDimensions(this->row_col_rvalue_vec_matrix, this->init_row, this->init_col);
    this->testMatrixElements(this->row_col_rvalue_vec_matrix, this->rvalue_row_col_expected);
    EXPECT_TRUE(this->rvalue_row_col_init_vec.empty());
}
TYPED_TEST(MatrixConstructor, DimensionMismatch) {
    const size_t init_dim{4};
    const std::vector<TypeParam> init_vec(init_dim * init_dim + 1, 1);
    EXPECT_THROW({ linalg::Matrix<TypeParam> test_matrix(init_dim, init_vec); }, std::invalid_argument);
}
TYPED_TEST(MatrixConstructor, RowColMismatch) {
    const size_t init_row{3}, init_col{4};
    const std::vector<TypeParam> init_vec(init_row * init_col + 1, 1);
    EXPECT_THROW({ linalg::Matrix<TypeParam> test_matrix(init_row, init_col, init_vec); }, std::invalid_argument);
}
