cmake_minimum_required(VERSION 3.10)
project(num_methods_tsk_1_2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

configure_file(input.txt input.txt COPYONLY)

add_executable(
        num_methods_tsk_1_2
        main.cpp)