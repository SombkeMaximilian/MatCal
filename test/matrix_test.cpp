#include <gtest/gtest.h>
#include "matrix.h"

TEST(MatrixConstructor, DefaultConstructor) {
    linalg::Matrix m;

    EXPECT_EQ(m.nRows(), 0);
    EXPECT_EQ(m.nCols(), 0);
}

TEST(MatrixConstructor, DimensionConstructor) {
    size_t init_dim{3};
    linalg::Matrix m(init_dim);

    EXPECT_EQ(m.nRows(), init_dim);
    EXPECT_EQ(m.nCols(), init_dim);
    for ( size_t i = 0; i < init_dim; ++i ) {
        for ( size_t j = 0; j < init_dim; ++j ) {
            EXPECT_DOUBLE_EQ(m(i, j), 0);
        }
    }
}


TEST(MatrixConstructor, InitializerListDimensionConstructor) {
    size_t init_dim{2};
    std::vector<double> init_vec{1.0, 2.0, 3.0, 4.0};
    double expected_value = 0.0;
    linalg::Matrix m(init_dim, init_vec);

    EXPECT_EQ(m.nRows(), init_dim);
    EXPECT_EQ(m.nCols(), init_dim);
    for ( size_t i = 0; i < init_dim; ++i ) {
        for ( size_t j = 0; j < init_dim; ++j ) {
            EXPECT_DOUBLE_EQ(m(i, j), ++expected_value);
        }
    }
}

TEST(MatrixConstructor, RowColConstructor) {
    size_t init_row{2}, init_col{3};
    linalg::Matrix m(init_row, init_col);
    EXPECT_EQ(m.nRows(), init_row);
    EXPECT_EQ(m.nCols(), init_col);
    for ( size_t i = 0; i < init_row; ++i ) {
        for ( size_t j = 0; j < init_col; ++j ) {
            EXPECT_DOUBLE_EQ(m(i, j), 0);
        }
    }
}

TEST(MatrixConstructor, InitializerListRowColConstructor) {
    size_t init_row{2}, init_col{3};
    std::vector<double> init_vec{1, 2, 3, 4, 5, 6};
    double expected_value = 0.0;
    linalg::Matrix m(init_row, init_col, init_vec);

    EXPECT_EQ(m.nRows(), init_row);
    EXPECT_EQ(m.nCols(), init_col);
    for ( size_t i = 0; i < init_row; ++i ) {
        for ( size_t j = 0; j < init_col; ++j ) {
            EXPECT_DOUBLE_EQ(m(i, j), ++expected_value);
        }
    }
}

TEST(MatrixConstructor, DimensionConstructorMismatch) {
    size_t init_dim{4};
    std::vector<double> init_vec(init_dim * init_dim + 1, 1.0);

    EXPECT_THROW({ linalg::Matrix m(init_dim, init_vec); }, std::invalid_argument);
}

TEST(MatrixConstructor, RowColConstructorMismatch) {
    size_t init_row{3}, init_col{4};
    std::vector<double> init_vec(init_row * init_col + 1, 1.0);

    EXPECT_THROW({ linalg::Matrix m(init_row, init_col, init_vec); }, std::invalid_argument);
}

TEST(MatrixAccess, AccessElementOutOfRange) {
    size_t init_dim{3};
    linalg::Matrix m(init_dim);

    EXPECT_THROW({ double val{m(3, 0)}; }, std::out_of_range);
}
