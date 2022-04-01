cmake_minimum_required(VERSION 3.16.0)

# Define local project name
project(microstructure_build)

# Check if user wants parallel assembly
option(OPENMP "Build project with OpenMP support." OFF)

# Prefix Path to install directory
set(CMAKE_PREFIX_PATH ../../install)
message(STATUS "Searching for install directory in ${CMAKE_PREFIX_PATH}!")

set(CMAKE_CXX_STANDARD 17)

find_package(bezierManipulation)

add_executable(microstructure_build microstructure_build.cpp)

target_sources(microstructure_build PRIVATE cross_tile_3d.hpp ring_segments_3d.hpp)

set(LINK_LIBRARIES bezierManipulation::bezierManipulation)


if (OPENMP)
  # Add OpenMP to the build system
  find_package(OpenMP REQUIRED QUIET)
  message("Found OpenMP Version ${OPENMP_CXX_VERSION}")
  list(APPEND LINK_LIBRARIES OpenMP::OpenMP_CXX)
  # Add Preprocessor flag
  add_compile_definitions(ENABLE_OPEN_MP_PARALLEL)
  add_compile_options(-fopenmp)
else()
  message("......................\n"
          "This example can be constructed in parallel. To enable it set the \n" 
          "CMAKE option OPENMP to on and rerun the exampe:\n"
          "For Example:\n"
          "cmake -DOPENMP=ON ..\n"
          "......................")
endif()

message(${LINK_LIBRARIES})

target_link_libraries(microstructure_build PRIVATE ${LINK_LIBRARIES})

