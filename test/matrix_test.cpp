#include <gtest/gtest.h>
#include <complex>
#include <vector>
#include "matrix.hpp"
#include "test_utils.hpp"

using MatrixTypes = ::testing::Types<
        int,
        short,
        long,
        float,
        double,
        std::complex<float>,
        std::complex<double>
>;

template<typename T>
class MatrixConstructor : public MatrixTestBase<T> {
}; // MatrixConstructor

TYPED_TEST_SUITE(MatrixConstructor, MatrixTypes);

TYPED_TEST(MatrixConstructor, Default) {
    linalg::Matrix<TypeParam> testMatrix;
    this->testMatrixDimensions(testMatrix, 0, 0);
}

TYPED_TEST(MatrixConstructor, Square) {
    size_t initDim{3};
    std::vector<TypeParam> expectedValues(initDim * initDim, 0);
    linalg::Matrix<TypeParam> testMatrix(initDim);
    this->testMatrixDimensions(testMatrix, initDim, initDim);
    this->testMatrixElements(testMatrix, expectedValues);
}

TYPED_TEST(MatrixConstructor, SquareInitVecCopy) {
    size_t initDim{3};
    std::vector<TypeParam> initVec{1, 2, 3, 4, 5, 6, 7, 8, 9};
    linalg::Matrix<TypeParam> testMatrix(initDim, initVec);
    this->testMatrixDimensions(testMatrix, initDim, initDim);
    this->testMatrixElements(testMatrix, initVec);
}

TYPED_TEST(MatrixConstructor, SquareInitVecMove) {
    size_t initDim{3};
    std::vector<TypeParam> initVec{1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<TypeParam> expectedValues{initVec};
    linalg::Matrix<TypeParam> testMatrix(initDim, std::move(initVec));
    this->testMatrixDimensions(testMatrix, initDim, initDim);
    this->testMatrixElements(testMatrix, expectedValues);
    EXPECT_TRUE(initVec.empty());
}

TYPED_TEST(MatrixConstructor, ByDims) {
    size_t initRow{2};
    size_t initCol{3};
    std::vector<TypeParam> expectedValues(initRow * initCol, 0);
    linalg::Matrix<TypeParam> testMatrix(initRow, initCol);
    this->testMatrixDimensions(testMatrix, initRow, initCol);
    this->testMatrixElements(testMatrix, expectedValues);
}

TYPED_TEST(MatrixConstructor, ByDimsInitVecCopy) {
    size_t initRow{2};
    size_t initCol{3};
    std::vector<TypeParam> initVec{1, 2, 3, 4, 5, 6};
    linalg::Matrix<TypeParam> testMatrix(initRow, initCol, initVec);
    this->testMatrixDimensions(testMatrix, initRow, initCol);
    this->testMatrixElements(testMatrix, initVec);
}

TYPED_TEST(MatrixConstructor, ByDimsInitVecMove) {
    size_t initRow{2};
    size_t initCol{3};
    std::vector<TypeParam> initVec{1, 2, 3, 4, 5, 6};
    std::vector<TypeParam> expectedValues{initVec};
    linalg::Matrix<TypeParam> testMatrix(initRow, initCol, std::move(initVec));
    this->testMatrixDimensions(testMatrix, initRow, initCol);
    this->testMatrixElements(testMatrix, expectedValues);
    EXPECT_TRUE(initVec.empty());
}

TYPED_TEST(MatrixConstructor, SquareMismatch) {
    size_t initDim{3};
    std::vector<TypeParam> initVec(initDim * initDim + 1, 1);
    EXPECT_THROW({
            linalg::Matrix<TypeParam> testMatrix(initDim, initVec);
        },
        std::invalid_argument);
}

TYPED_TEST(MatrixConstructor, ByDimsMismatch) {
    size_t initRow{3};
    size_t initCol{4};
    std::vector<TypeParam> initVec(initRow * initCol + 1, 1);
    EXPECT_THROW({
            linalg::Matrix<TypeParam> testMatrix(initRow, initCol, initVec);
        },
        std::invalid_argument);
}


template<typename T>
class MatrixElements : public MatrixTestBase<T> {

protected:
    size_t initDim{3};
    linalg::Matrix<T> testMatrix;

    void SetUp() override {
        testMatrix = linalg::Matrix<T>(initDim);
    }

}; // MatrixElements

TYPED_TEST_SUITE(MatrixElements, MatrixTypes);

TYPED_TEST(MatrixElements, AccessElementOutOfRangeRow) {
    EXPECT_THROW({
            [[maybe_unused]] TypeParam val{this->testMatrix(3, 0)};
        },
        std::out_of_range);
}
TYPED_TEST(MatrixElements, AccessElementOutOfRangeCol) {
    EXPECT_THROW({
            [[maybe_unused]] TypeParam val{this->testMatrix(0, 3)};
        },
        std::out_of_range);
}
TYPED_TEST(MatrixElements, AccessElementOutOfRangeRowCol) {
    EXPECT_THROW({
            [[maybe_unused]] TypeParam val{this->testMatrix(3, 3)};
        },
        std::out_of_range);
}


template<typename T>
class MatrixArithmetics : public MatrixTestBase<T>  {

protected:

}; // Matrix Arithmetics
