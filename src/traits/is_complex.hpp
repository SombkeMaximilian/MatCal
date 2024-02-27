#ifndef MATCAL_IS_COMPLEX_HPP
#define MATCAL_IS_COMPLEX_HPP

#include <complex>
#include <type_traits>

template<typename T>
struct is_complex : std::false_type {};

template<typename T>
struct is_complex<std::complex<T>> : std::true_type {};

#endif //MATCAL_IS_COMPLEX_HPP
