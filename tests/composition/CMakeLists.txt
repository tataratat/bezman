set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/composition_test.cpp
)

if (MSVC_VERSION GREATER_EQUAL "1900")
    message("IS 1900")
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("/std:c++20" _cpp_latest_flag_supported)
    if (_cpp_latest_flag_supported)
        message("supports 20")
        add_compile_options("/std:c++20")
    else()
        message("NO 20")
    endif()
else()
    message("NOT 1900")
endif()

# Search for the library
add_executable(composition_test ${TEST_SOURCES})

# Link bezier manipulation library to the executable target
add_dependencies(composition_test ${MY_FANCY_TARGET_NAME})
target_link_libraries(composition_test 
                  PRIVATE ${MY_FANCY_TARGET_NAME} gtest_main)


target_compile_options(composition_test PRIVATE "/std:c++20")

include(GoogleTest)
add_test(NAME composition_test COMMAND composition_test)
