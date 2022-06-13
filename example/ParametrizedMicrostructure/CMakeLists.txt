cmake_minimum_required(VERSION 3.16.0)

# Define local project name
project(microstructure_build)

# Prefix Path to install directory
set(CMAKE_PREFIX_PATH ../../install)
message(STATUS "Searching for install directory in ${CMAKE_PREFIX_PATH}!")

set(CMAKE_CXX_STANDARD 17)

find_package(bezman)

add_executable(microstructure_build microstructure_build.cpp)

target_link_libraries(microstructure_build PRIVATE bezman::bezman)
