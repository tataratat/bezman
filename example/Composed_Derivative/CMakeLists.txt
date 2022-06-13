cmake_minimum_required(VERSION 3.16.0)

# Define local project name
project(composition_derivative)

# Prefix Path to install directory
set(CMAKE_PREFIX_PATH ../../install)
message(STATUS "Searching for install directory in ${CMAKE_PREFIX_PATH}!")

set(CMAKE_CXX_STANDARD 17)

find_package(bezman)

add_executable(composition_derivative composition_derivative.cpp)

target_link_libraries(composition_derivative PRIVATE bezman::bezman)
