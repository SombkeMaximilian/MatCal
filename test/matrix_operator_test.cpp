#include <gtest/gtest.h>
#include <complex>
#include <vector>
#include "matrix.hpp"
#include "test_utils.hpp"

template<typename T>
class MatrixOperators : public ::testing::Test {};

TYPED_TEST_SUITE(MatrixOperators, MatrixTypes);

TYPED_TEST(MatrixOperators, PlusOperator) {
    const size_t init_dim{2};
    std::vector<TypeParam> init_vec1{1, 2, 3, 4};
    std::vector<TypeParam> init_vec2{5, 6, 7, 8};
    std::vector<TypeParam> expected_vec{6, 8, 10, 12};
    linalg::Matrix<TypeParam> test_matrix1(init_dim, init_vec1);
    linalg::Matrix<TypeParam> test_matrix2(init_dim, init_vec2);
    linalg::Matrix<TypeParam> result_matrix;
    linalg::Matrix<TypeParam> expected_matrix(init_dim, expected_vec);

    result_matrix = test_matrix1 + test_matrix2;
    for ( size_t i = 0; i < init_dim; ++i ) {
        for ( size_t j = 0; j < init_dim; ++j ) {
            expect_type_eq(result_matrix(i, j), expected_matrix(i, j));
        }
    }
}

TYPED_TEST(MatrixOperators, MinusOperator) {
    const size_t init_dim{2};
    std::vector<TypeParam> init_vec1{1, 8, 2, 9};
    std::vector<TypeParam> init_vec2{3, 6, 4, 7};
    std::vector<TypeParam> expected_vec{-2, 2, -2, 2};
    linalg::Matrix<TypeParam> test_matrix1(init_dim, init_vec1);
    linalg::Matrix<TypeParam> test_matrix2(init_dim, init_vec2);
    linalg::Matrix<TypeParam> result_matrix;
    linalg::Matrix<TypeParam> expected_matrix(init_dim, expected_vec);

    result_matrix = test_matrix1 - test_matrix2;
    for ( size_t i = 0; i < init_dim; ++i ) {
        for ( size_t j = 0; j < init_dim; ++j ) {
            expect_type_eq(result_matrix(i, j), expected_matrix(i, j));
        }
    }
}

TYPED_TEST(MatrixOperators, UnaryMinusOperator) {
    const size_t init_dim{2};
    std::vector<TypeParam> init_vec{1, 2, 3, 4};
    std::vector<TypeParam> expected_vec{-1, -2, -3, -4};
    linalg::Matrix<TypeParam> test_matrix(init_dim, init_vec);
    linalg::Matrix<TypeParam> result_matrix;
    linalg::Matrix<TypeParam> expected_matrix(init_dim, expected_vec);

    result_matrix = -test_matrix;
    for ( size_t i = 0; i < init_dim; ++i ) {
        for ( size_t j = 0; j < init_dim; ++j ) {
            expect_type_eq(result_matrix(i, j), expected_matrix(i, j));
        }
    }
}
