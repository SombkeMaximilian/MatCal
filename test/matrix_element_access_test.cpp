#include <gtest/gtest.h>
#include <complex>
#include "matrix.hpp"
#include "test_utils.hpp"

template<typename T>
class MatrixElementAccess : public ::testing::Test {};

TYPED_TEST_SUITE(MatrixElementAccess, MatrixTypes);

TYPED_TEST(MatrixElementAccess, AccessElementOutOfRange) {
    size_t init_dim{3};
    linalg::Matrix<TypeParam> m(init_dim);

    EXPECT_THROW({ TypeParam val{m(3, 0)}; }, std::out_of_range);
    EXPECT_THROW({ TypeParam val{m(0, 3)}; }, std::out_of_range);
    EXPECT_THROW({ TypeParam val{m(3, 3)}; }, std::out_of_range);
}
