cmake_minimum_required(VERSION 3.16.0)

# Define local project name
project(simple_composition)

# Prefix Path to install directory
set(CMAKE_PREFIX_PATH ../../install)
message(STATUS "Searching for install directory in ${CMAKE_PREFIX_PATH}!")

set(CMAKE_CXX_STANDARD 17)

find_package(bezierManipulation)

add_executable(simple_composition simple_composition.cpp)

target_link_libraries(simple_composition PRIVATE bezierManipulation::bezierManipulation)
