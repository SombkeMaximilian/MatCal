cmake_minimum_required(VERSION 3.27)
project(matrix_test)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(MATRIX_TEST_FILES
        matrix_test.cpp
        test_utils.h
)

add_executable(MatrixTest ${MATRIX_TEST_FILES})
target_link_libraries(MatrixTest PRIVATE gtest_main)
target_include_directories(MatrixTest PRIVATE ../include/matcal)

add_test(NAME MatrixTest COMMAND MatrixTest)
