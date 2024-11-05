cmake_minimum_required(VERSION 3.10)
project(num_methods_util_algorithms)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

add_library(
        num_methods_util_algorithms
        src/algorithms.cpp)
target_include_directories(
        num_methods_util_checks
        PUBLIC
        ./include)
target_link_libraries(
        num_methods_util_algorithms
        PUBLIC
        num_methods_util_checks)
target_link_libraries(
        num_methods_util_algorithms
        PUBLIC
        num_methods_util_vecn)
target_link_libraries(
        num_methods_util_algorithms
        PUBLIC
        num_methods_util_matrixnxn)
target_link_libraries(
        num_methods_util_algorithms
        PUBLIC
        num_methods_util_polynomial)