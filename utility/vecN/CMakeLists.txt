cmake_minimum_required(VERSION 3.10)
project(num_methods_util_vecn)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

add_library(
        num_methods_util_vecn
        src/vecN.cpp)
target_include_directories(
        num_methods_util_vecn
        PUBLIC
        ./include)
target_link_libraries(
        num_methods_util_vecn
        PUBLIC
        num_methods_util_checks)