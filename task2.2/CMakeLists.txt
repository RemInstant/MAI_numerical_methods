cmake_minimum_required(VERSION 3.10)
project(num_methods_tsk_2_2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

configure_file(input.txt input.txt COPYONLY)
add_executable(
        num_methods_tsk_2_2
        main.cpp)
target_link_libraries(
        num_methods_tsk_2_2
        PUBLIC
        num_methods_util_checks)
target_link_libraries(
        num_methods_tsk_2_2
        PUBLIC
        num_methods_util_vecn)
target_link_libraries(
        num_methods_tsk_2_2
        PUBLIC
        num_methods_util_matrixnxn)
target_link_libraries(
        num_methods_tsk_2_2
        PUBLIC
        num_methods_util_algorithms)