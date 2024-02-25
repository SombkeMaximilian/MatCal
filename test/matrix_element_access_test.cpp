#include <gtest/gtest.h>
#include <complex>
#include "matrix.hpp"
#include "test_utils.hpp"

template<typename T>
class MatrixElementAccess : public ::testing::Test {
protected:
    size_t init_dim{3};
    linalg::Matrix<T> test_matrix;
    void SetUp() override {
        test_matrix = linalg::Matrix<T>(init_dim);
    }
};

TYPED_TEST_SUITE(MatrixElementAccess, MatrixTypes);

TYPED_TEST(MatrixElementAccess, AccessElementOutOfRangeRow) {
    EXPECT_THROW({ [[maybe_unused]] TypeParam val{this->test_matrix(3, 0)}; }, std::out_of_range);
}
TYPED_TEST(MatrixElementAccess, AccessElementOutOfRangeCol) {
    EXPECT_THROW({ [[maybe_unused]] TypeParam val{this->test_matrix(0, 3)}; }, std::out_of_range);
}
TYPED_TEST(MatrixElementAccess, AccessElementOutOfRangeRowCol) {
    EXPECT_THROW({ [[maybe_unused]] TypeParam val{this->test_matrix(3, 3)}; }, std::out_of_range);
}
