#include <gtest/gtest.h>
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
