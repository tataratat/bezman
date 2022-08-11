set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/spline_splitting_test.cpp
)

# Search for the library
add_executable(spline_splitting_test ${TEST_SOURCES})

# Link bezier manipulation library to the executable target
add_dependencies(spline_splitting_test ${MY_FANCY_TARGET_NAME})
target_link_libraries(spline_splitting_test 
                  PRIVATE ${MY_FANCY_TARGET_NAME}
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME spline_splitting_test COMMAND spline_splitting_test)
