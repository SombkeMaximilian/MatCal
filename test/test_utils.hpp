#ifndef MATCAL_TEST_UTILS_HPP
#define MATCAL_TEST_UTILS_HPP

#include <complex>

template<typename T>
void EXPECT_TYPE_EQ(const T& actual, const T& expected) {
    EXPECT_EQ(actual, expected);
}

template<>
inline void EXPECT_TYPE_EQ(const float& actual, const float& expected) {
    EXPECT_FLOAT_EQ(actual, expected);
}

template<>
inline void EXPECT_TYPE_EQ(const double& actual, const double& expected) {
    EXPECT_DOUBLE_EQ(actual, expected);
}

template<>
inline void EXPECT_TYPE_EQ(const std::complex<float>& actual, const std::complex<float>& expected) {
    EXPECT_FLOAT_EQ(actual.real(), expected.real());
    EXPECT_FLOAT_EQ(actual.imag(), expected.imag());
}

template<>
inline void EXPECT_TYPE_EQ(const std::complex<double>& actual, const std::complex<double>& expected) {
    EXPECT_DOUBLE_EQ(actual.real(), expected.real());
    EXPECT_DOUBLE_EQ(actual.imag(), expected.imag());
}

#endif //MATCAL_TEST_UTILS_HPP
