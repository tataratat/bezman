cmake_minimum_required(VERSION 3.16.0)

# Define local project name
project(logging_example)

# Prefix Path to install directory
set(CMAKE_PREFIX_PATH ../../install)
message(STATUS "Searching for install directory in ${CMAKE_PREFIX_PATH}!")

set(CMAKE_CXX_STANDARD 17)

find_package(bezman)

add_executable(logging_example logging_example.cpp)

target_link_libraries(logging_example PRIVATE bezman::bezman)
