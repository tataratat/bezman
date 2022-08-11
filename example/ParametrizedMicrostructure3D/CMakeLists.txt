cmake_minimum_required(VERSION 3.16.0)

# Define local project name
project(microstructure_build)

option(LOGGING "Show and automatically save logging " ON)

# Check if user wants parallel assembly
option(OPENMP "Build project with OpenMP support." OFF)

# Prefix Path to install directory
# This is to help find_package(bezman).
#
# You can also set this in commendline with:
#  cmake -DCMAKE_PREFIX_PATH=<path-to-install-dir> .....
#
# <path-to-install-dir> should be the path where you set
#   `-DCMAKE_INSTALL_PREFIX` for `bezman`
#
# Here, we append path that corresponds to default install dir of `bezman`.
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../install)
message(STATUS "Searching for install directory in ${CMAKE_PREFIX_PATH}!")

find_package(bezman REQUIRED)

add_executable(microstructure_build microstructure_build.cpp)
target_include_directories(microstructure_build PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(microstructure_build PRIVATE bezman::bezman)

if(LOGGING)
    # link `bezman::logging` to enable logging.
    # note that having LOGGING option is not necessary.
    target_link_libraries(microstructure_build PRIVATE bezman::logging)
endif()

if (OPENMP)
  # Add OpenMP to the build system
  find_package(OpenMP REQUIRED)
  message("Found OpenMP Version ${OPENMP_CXX_VERSION}")
  # Add Preprocessor flag
  add_compile_definitions(ENABLE_OPEN_MP_PARALLEL)
  target_compile_options(microstructure_build PRIVATE -fopenmp)
  target_link_libraries(microstructure_build PRIVATE OpenMP::OpenMP_CXX)
else()
  message("......................\n"
          "This example can be constructed in parallel. To enable it set the \n" 
          "CMAKE option OPENMP to on and rerun the exampe:\n"
          "For Example:\n"
          "cmake -DOPENMP=ON ..\n"
          "......................")
endif()
