cmake_minimum_required(VERSION 3.27)
project(MatCal)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(src)

add_executable(MatCal main.cpp)
target_link_libraries(MatCal PRIVATE linalg)
target_include_directories(MatCal PRIVATE ${PROJECT_SOURCE_DIR}/src)
