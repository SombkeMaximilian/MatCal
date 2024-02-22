#include <gtest/gtest.h>
#include <complex>
#include <vector>
#include "matrix.hpp"
#include "test_utils.hpp"

template<typename T>
class MatrixConstructor : public ::testing::Test {};

TYPED_TEST_SUITE(MatrixConstructor, MatrixTypes);

TYPED_TEST(MatrixConstructor, DefaultConstructor) {
    linalg::Matrix<TypeParam> test_matrix;

    EXPECT_EQ(test_matrix.nRows(), 0);
    EXPECT_EQ(test_matrix.nCols(), 0);
}

TYPED_TEST(MatrixConstructor, DimensionConstructor) {
    const size_t init_dim{3};
    linalg::Matrix<TypeParam> test_matrix(init_dim);

    EXPECT_EQ(test_matrix.nRows(), init_dim);
    EXPECT_EQ(test_matrix.nCols(), init_dim);
    for ( size_t i = 0; i < init_dim; ++i ) {
        for ( size_t j = 0; j < init_dim; ++j ) {
            expect_type_eq(test_matrix(i, j), TypeParam(0));
        }
    }
}


TYPED_TEST(MatrixConstructor, InitializerListDimensionConstructor) {
    const size_t init_dim{2};
    const std::vector<TypeParam> init_vec{1, 2, 3, 4};
    linalg::Matrix<TypeParam> test_matrix(init_dim, init_vec);

    EXPECT_EQ(test_matrix.nRows(), init_dim);
    EXPECT_EQ(test_matrix.nCols(), init_dim);
    for ( size_t i = 0; i < init_dim; ++i ) {
        for ( size_t j = 0; j < init_dim; ++j ) {
            expect_type_eq(test_matrix(i, j), init_vec[i * test_matrix.nCols() + j]);
        }
    }
}

TYPED_TEST(MatrixConstructor, RowColConstructor) {
    const size_t init_row{2}, init_col{3};
    linalg::Matrix<TypeParam> test_matrix(init_row, init_col);

    EXPECT_EQ(test_matrix.nRows(), init_row);
    EXPECT_EQ(test_matrix.nCols(), init_col);
    for ( size_t i = 0; i < init_row; ++i ) {
        for ( size_t j = 0; j < init_col; ++j ) {
            expect_type_eq(test_matrix(i, j), TypeParam(0));
        }
    }
}

TYPED_TEST(MatrixConstructor, InitializerListRowColConstructor) {
    const size_t init_row{2}, init_col{3};
    const std::vector<TypeParam> init_vec{1, 2, 3, 4, 5, 6};
    linalg::Matrix<TypeParam> test_matrix(init_row, init_col, init_vec);

    EXPECT_EQ(test_matrix.nRows(), init_row);
    EXPECT_EQ(test_matrix.nCols(), init_col);
    for ( size_t i = 0; i < init_row; ++i ) {
        for ( size_t j = 0; j < init_col; ++j ) {
            expect_type_eq(test_matrix(i, j), init_vec[i * test_matrix.nCols() + j]);
        }
    }
}

TYPED_TEST(MatrixConstructor, DimensionConstructorMismatch) {
    const size_t init_dim{4};
    const std::vector<TypeParam> init_vec(init_dim * init_dim + 1, 1);

    EXPECT_THROW({ linalg::Matrix<TypeParam> test_matrix(init_dim, init_vec); }, std::invalid_argument);
}

TYPED_TEST(MatrixConstructor, RowColConstructorMismatch) {
    const size_t init_row{3}, init_col{4};
    const std::vector<TypeParam> init_vec(init_row * init_col + 1, 1);

    EXPECT_THROW({ linalg::Matrix<TypeParam> test_matrix(init_row, init_col, init_vec); }, std::invalid_argument);
}
