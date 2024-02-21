#include <gtest/gtest.h>
#include <complex>
#include "matrix.hpp"
#include "test_utils.hpp"

template<typename T>
class MatrixElementAccess : public ::testing::Test {};

TYPED_TEST_SUITE(MatrixElementAccess, MatrixTypes);

TEST(MatrixAccess, AccessElementOutOfRange) {
    size_t init_dim{3};
    linalg::Matrix<double> m(init_dim);

    EXPECT_THROW({ double val{m(3, 0)}; }, std::out_of_range);
    EXPECT_THROW({ double val{m(0, 3)}; }, std::out_of_range);
    EXPECT_THROW({ double val{m(3, 3)}; }, std::out_of_range);
}
