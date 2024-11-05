cmake_minimum_required(VERSION 3.10)
project(num_methods_util)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

add_subdirectory(checks)

add_subdirectory(vecN)
add_subdirectory(matrixNxN)
add_subdirectory(polynomial)
add_subdirectory(piecewise_polynomial)

add_subdirectory(algorithms)