cmake_minimum_required(VERSION 3.27)
project(matrix_test)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(TEST_FILES
        matrix_constructor_test.cpp
        matrix_element_access_test.cpp
        test_utils.hpp
)

add_executable(matrix_test ${TEST_FILES})
target_link_libraries(matrix_test PRIVATE gtest_main)
target_include_directories(matrix_test PRIVATE ../src)

add_test(NAME matrix_test COMMAND matrix_test)
