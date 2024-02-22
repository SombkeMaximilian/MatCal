#include <gtest/gtest.h>
#include <complex>
#include "matrix.hpp"
#include "test_utils.hpp"

template<typename T>
class MatrixElementAccess : public ::testing::Test {};

TYPED_TEST_SUITE(MatrixElementAccess, MatrixTypes);

TYPED_TEST(MatrixElementAccess, AccessElementOutOfRangeRow) {
    size_t init_dim{3};
    linalg::Matrix<TypeParam> test_matrix(init_dim);

    EXPECT_THROW({ TypeParam val{test_matrix(3, 0)}; }, std::out_of_range);
}

TYPED_TEST(MatrixElementAccess, AccessElementOutOfRangeCol) {
    size_t init_dim{3};
    linalg::Matrix<TypeParam> test_matrix(init_dim);

    EXPECT_THROW({ TypeParam val{test_matrix(0, 3)}; }, std::out_of_range);
}

TYPED_TEST(MatrixElementAccess, AccessElementOutOfRangeRowCol) {
    size_t init_dim{3};
    linalg::Matrix<TypeParam> test_matrix(init_dim);

    EXPECT_THROW({ TypeParam val{test_matrix(3, 3)}; }, std::out_of_range);
}
