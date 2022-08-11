set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/rational_splines_test.cpp
)

# Search for the library
add_executable(rational_splines_test ${TEST_SOURCES})

# Link bezier manipulation library to the executable target
add_dependencies(rational_splines_test ${MY_FANCY_TARGET_NAME})
target_link_libraries(rational_splines_test 
                  PRIVATE ${MY_FANCY_TARGET_NAME}
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME rational_splines_test COMMAND rational_splines_test)
