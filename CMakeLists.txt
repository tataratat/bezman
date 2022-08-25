####
# General information
####
# Required CMAKE 
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# Define Project
project("bezman"
        VERSION 0.0.1
        DESCRIPTION
            "A header library to facilitate the manipulation of Bezier-splines"
        LANGUAGES CXX)

# Options
option(BEZMAN_GOOGLETEST "Build with GoogleTest" ON)

# Config
set(namespace "${PROJECT_NAME}::")
set(TARGETS_EXPORT_NAME "${PROJECT_NAME}Targets")
set(MY_FANCY_TARGET_NAME "${PROJECT_NAME}")
set(MY_FANCY_LOGGER_TARGET "logging")
set(lib_destination "lib")
set(config_destination "${lib_destination}/cmake/${PROJECT_NAME}")
set(gen_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")
set(version_config "${gen_dir}/${PROJECT_NAME}ConfigVersion.cmake")
set(project_config "${gen_dir}/${PROJECT_NAME}Config.cmake")
set(include_destination "include")
set(include_destination_headers ${include_destination}/${PROJECT_NAME})

# Set defaults
# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Install prefix
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(default_install_dir ${CMAKE_CURRENT_SOURCE_DIR}/install)
    message(STATUS
                "Setting default install prefix: "
                "${default_install_dir}")
    set(CMAKE_INSTALL_PREFIX
            ${default_install_dir} CACHE PATH "" FORCE)
endif()

# Print out information
message(STATUS "Compile using ${CMAKE_CXX_COMPILER_ID} compiler")
if(CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Build debug!")
else()
    message(STATUS "Build release!")
endif()

# Set compile information
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(OPTIMIZATION_FLAGS $<$<NOT:$<CONFIG:Debug>>:-O3 -march=native>)
    set(RUNTIME_CHECKS_DEBUG -D_GLIBCXX_DEBUG)
    set(WARNING_FLAGS -Wall -Wextra -Wpedantic -Wzero-as-null-pointer-constant
                      -Wno-unused>)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(
        WARNING 
            "CLANG is not sufficiently tested. "
            "Please help us by submitting issues (or PR!) to "
            "[github.com/tataratat/bezman] with any problems."
            "Thank you!"
    )
    set(OPTIMIZATION_FLAGS $<$<NOT:$<CONFIG:Debug>>:-O3 -march=native>)
    set(RUNTIME_CHECKS_DEBUG)
    set(WARNING_FLAGS 
            -Wall -Wextra -Wmost -Wextra -Wpedantic
            -Wunreachable-code -Wshadow -Wfloat-equal -Weffc++
            -Wno-unused-parameter -Wno-unused-variable>
            $<$<BOOL:${BEZMAN_GOOGLETEST}>:
                -Wno-gnu-zero-variadic-macro-arguments>)
else()
    message(ERROR "Unsupported compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()

# Add Files by going through directory structure and add support for IDEs
add_subdirectory(src)
set_source_files_properties(${HEADERS}
    PROPERTIES
        LANGUAGE CXX
        HEADER_FILE_ONLY TRUE)
source_group(TREE ${PROJECT_SOURCE_DIR} FILES ${HEADERS})

# Define Targets
# 1. My fancy target, bezman
add_library(${MY_FANCY_TARGET_NAME} INTERFACE)
add_library(${namespace}${MY_FANCY_TARGET_NAME} ALIAS ${MY_FANCY_TARGET_NAME})
target_compile_definitions(${MY_FANCY_TARGET_NAME}
    INTERFACE
        $<$<NOT:$<CONFIG:Debug>>:NDEBUG>)
target_include_directories(${MY_FANCY_TARGET_NAME}
    INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:${incl_dest}>)
target_compile_options(${MY_FANCY_TARGET_NAME}
    INTERFACE
        $<IF:$<CONFIG:Release>,${OPTIMIZATION_FLAGS},
                               ${RUNTIME_CHECK_DEBUG}>
        $<IF:$<CONFIG:Debug>,${WARNING_FLAGS}>)
target_compile_features(${MY_FANCY_TARGET_NAME} INTERFACE cxx_std_17)

# 2. logging
add_library(${MY_FANCY_LOGGER_TARGET} INTERFACE)
add_library(
    ${namespace}${MY_FANCY_LOGGER_TARGET} ALIAS ${MY_FANCY_LOGGER_TARGET}
)
target_compile_definitions(${MY_FANCY_LOGGER_TARGET} INTERFACE ENABLE_LOGGING)

# Install and provide testing 
# Always download and build googletest
if(BEZMAN_GOOGLETEST)
    message(STATUS "Build tests!")
    message(STATUS "---- downloading and building GoogleTest...")
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/AddGoogleTest.cmake)
    add_subdirectory(tests)
else()
    message(STATUS "No tests!")
endif()

# Load cmake Macros
include(CMakePackageConfigHelpers)

# Set Up Version
write_basic_package_version_file(
    "${version_config}" COMPATIBILITY SameMajorVersion
)

# Automatic configuration file
configure_package_config_file(
    "cmake/config.cmake.in"
    "${project_config}"
    INSTALL_DESTINATION "${config_destination}"
)

## INSTALLATION PROCESS
install(
    TARGETS ${MY_FANCY_TARGET_NAME}
    EXPORT ${TARGETS_EXPORT_NAME}
    INCLUDES DESTINATION ${include_destination}
)
install(
    TARGETS ${MY_FANCY_LOGGER_TARGET}
    EXPORT ${TARGETS_EXPORT_NAME}
)
install(
    FILES ${project_config} ${version_config}
    DESTINATION ${config_destination}
)
install(
    EXPORT ${TARGETS_EXPORT_NAME}
    NAMESPACE ${namespace}
    DESTINATION ${config_destination}
)


