#ifndef MATCAL_TEST_UTILS_HPP
#define MATCAL_TEST_UTILS_HPP

#include <vector>
#include <complex>
#include <random>
#include <type_traits>

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
std::vector<T> generateRandomVector(size_t n) {
    std::vector<T> vec(n);
    std::mt19937 rng(42);
    if constexpr (std::is_integral<T>::value) {
        std::uniform_int_distribution<T> dist(0, 100);
        for (size_t i = 0; i < n; ++i) {
            vec[i] = dist(rng);
        }
    } else if constexpr (std::is_floating_point<T>::value) {
        std::uniform_real_distribution<T> dist(0.0, 1.0);
        for (size_t i = 0; i < n; ++i) {
            vec[i] = dist(rng);
        }
    } else if constexpr (std::is_same<T, std::complex<float>>::value || std::is_same<T, std::complex<double>>::value) {
        using value_type = typename T::value_type;
        std::uniform_real_distribution<value_type> dist(0.0, 1.0);
        for (size_t i = 0; i < n; ++i) {
            value_type realPart = dist(rng);
            value_type imagPart = dist(rng);
            vec[i] = T(realPart, imagPart);
        }
    }
    return vec;
}

template<typename T>
void expect_type_eq(const T& actual, const T& expected) {
    EXPECT_EQ(actual, expected);
}

template<>
inline void expect_type_eq(const float& actual, const float& expected) {
    EXPECT_FLOAT_EQ(actual, expected);
}

template<>
inline void expect_type_eq(const double& actual, const double& expected) {
    EXPECT_DOUBLE_EQ(actual, expected);
}

template<>
inline void expect_type_eq(const std::complex<float>& actual, const std::complex<float>& expected) {
    EXPECT_FLOAT_EQ(actual.real(), expected.real());
    EXPECT_FLOAT_EQ(actual.imag(), expected.imag());
}

template<>
inline void expect_type_eq(const std::complex<double>& actual, const std::complex<double>& expected) {
    EXPECT_DOUBLE_EQ(actual.real(), expected.real());
    EXPECT_DOUBLE_EQ(actual.imag(), expected.imag());
}

#endif //MATCAL_TEST_UTILS_HPP
