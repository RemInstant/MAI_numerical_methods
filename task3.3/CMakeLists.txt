cmake_minimum_required(VERSION 3.10)
project(num_methods_tsk_3_3)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

configure_file(input.txt input.txt COPYONLY)
add_executable(
        num_methods_tsk_3_3
        main.cpp)
target_link_libraries(
        num_methods_tsk_3_3
        PUBLIC
        num_methods_util_checks)
target_link_libraries(
        num_methods_tsk_3_3
        PUBLIC
        num_methods_util_polynomial)
target_link_libraries(
        num_methods_tsk_3_3
        PUBLIC
        num_methods_util_piecewise_polynomial)
target_link_libraries(
        num_methods_tsk_3_3
        PUBLIC
        num_methods_util_algorithms)