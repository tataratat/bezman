cmake_minimum_required(VERSION 3.16.0)

# Define local project name
project(microstructure_build)

option(LOGGING "Show and automatically save logging " ON)

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
target_link_libraries(microstructure_build PRIVATE bezman::bezman)

if(LOGGING)
    # link `bezman::logging` to enable logging.
    # note that having LOGGING option is not necessary.
    target_link_libraries(microstructure_build PRIVATE bezman::logging)
endif()
