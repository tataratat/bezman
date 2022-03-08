cmake_minimum_required(VERSION 3.0.0)

# Define local project name
project(basic_bezier_operations)

# Prefix Path to install directory
set(CMAKE_PREFIX_PATH ../../install)
message(STATUS "Searching for install directory in ${CMAKE_PREFIX_PATH}!")

set(CMAKE_CXX_STANDARD 17)

find_package(bezierManipulation)

add_executable(basic_bezier_operations basic_bezier_operations.cpp)

target_link_libraries(basic_bezier_operations PRIVATE bezierManipulation::bezierManipulation)
