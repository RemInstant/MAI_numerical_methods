cmake_minimum_required(VERSION 3.10)
project(num_methods_util_checks)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

add_library(
        num_methods_util_checks
        src/checks.cpp)
target_include_directories(
        num_methods_util_checks
        PUBLIC
        ./include)