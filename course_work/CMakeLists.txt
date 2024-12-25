cmake_minimum_required(VERSION 3.10)
project(num_methods_cw)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

configure_file(input.txt input.txt COPYONLY)
configure_file(graphs.py graphs.py COPYONLY)

add_executable(
        num_methods_cw
        main.cpp)
target_link_libraries(
        num_methods_cw
        PUBLIC
        num_methods_util_checks)
target_link_libraries(
        num_methods_cw
        PUBLIC
        num_methods_util_vecn)
target_link_libraries(
        num_methods_cw
        PUBLIC
        num_methods_util_algorithms)