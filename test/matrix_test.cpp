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
    this->testMatrixDimensions(this->matrixDefault, 0, 0);
}
TYPED_TEST(MatrixConstructor, Square) {
    std::vector<TypeParam> expectedValues(this->initDim * this->initDim, 0);
    this->testMatrixDimensions(this->matrixSquare, this->initDim, this->initDim);
    this->testMatrixElements(this->matrixSquare, expectedValues);
}
TYPED_TEST(MatrixConstructor, SquareInitVecCopy) {
    this->testMatrixDimensions(this->matrixSquareCopy, this->initDim, this->initDim);
    this->testMatrixElements(this->matrixSquareCopy, this->expectedVecSquare);
}
TYPED_TEST(MatrixConstructor, SquareInitVecMove) {
    this->testMatrixDimensions(this->matrixSquareMove, this->initDim, this->initDim);
    this->testMatrixElements(this->matrixSquareMove, this->expectedVecSquare);
    EXPECT_TRUE(this->initVecSquare.empty());
}
TYPED_TEST(MatrixConstructor, ByDims) {
    std::vector<TypeParam> expectedValues(this->initRow * this->initCol, 0);
    this->testMatrixDimensions(this->matrixByDims, this->initRow, this->initCol);
    this->testMatrixElements(this->matrixByDims, expectedValues);
}
TYPED_TEST(MatrixConstructor, ByDimsInitVecCopy) {
    this->testMatrixDimensions(this->matrixByDimsCopy, this->initRow, this->initCol);
    this->testMatrixElements(this->matrixByDimsCopy, this->expectedVecByDims);
}
TYPED_TEST(MatrixConstructor, ByDimsInitVecMove) {
    this->testMatrixDimensions(this->matrixByDimsMove, this->initRow, this->initCol);
    this->testMatrixElements(this->matrixByDimsMove, this->expectedVecByDims);
    EXPECT_TRUE(this->initVecByDims.empty());
}
TYPED_TEST(MatrixConstructor, SquareMismatch) {
    EXPECT_THROW({
            linalg::Matrix<TypeParam> test_matrix(this->initDim, this->initVecMismatch1);
        },
        std::invalid_argument);
}
TYPED_TEST(MatrixConstructor, ByDimsMismatch) {
    EXPECT_THROW({
            linalg::Matrix<TypeParam> test_matrix(this->initRow, this->initCol, this->initVecMismatch2);
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
