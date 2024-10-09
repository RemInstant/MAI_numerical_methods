cmake_minimum_required(VERSION 3.10)
project(num_methods_tsk_1_3)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

configure_file(input.txt input.txt COPYONLY)
configure_file(graphics.py graphics.py COPYONLY)

add_executable(
        num_methods_tsk_1_3
        main.cpp)

target_link_libraries(
        num_methods_tsk_1_3
        PUBLIC
        num_methods_util_vecn)
target_link_libraries(
        num_methods_tsk_1_3
        PUBLIC
        num_methods_util_matrixnxn)