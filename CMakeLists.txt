cmake_minimum_required(VERSION 3.27)
project(MatCal)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wextra -Wpedantic")

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

add_subdirectory(src)
add_subdirectory(test)

add_executable(MatCal examples/main.cpp)
target_include_directories(MatCal PRIVATE ${PROJECT_SOURCE_DIR}/src)

enable_testing()
