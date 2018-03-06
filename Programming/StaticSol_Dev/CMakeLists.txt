cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project("all" CXX)

include_directories(./header /usr/include/gsl)
include_directories(./header)

add_executable(tutorial ./src/main.cpp ./src/rhs.cpp)

target_link_libraries(tutorial boost_program_options)
target_link_libraries(tutorial gsl)
target_link_libraries(tutorial gslcblas)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS -Wall\ -O2)

