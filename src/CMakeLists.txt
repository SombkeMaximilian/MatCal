add_library(linalg matrix.cpp matrix.h)
target_include_directories(linalg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
