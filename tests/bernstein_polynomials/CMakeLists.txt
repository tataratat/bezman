set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/bernstein_polynomial_test.cpp
)

# Search for the library
add_executable(bernstein_polynomial_test ${TEST_SOURCES})

# Link bezier manipulation library to the executable target
add_dependencies(bernstein_polynomial_test ${MY_FANCY_TARGET_NAME})
target_link_libraries(bernstein_polynomial_test
                  PRIVATE ${MY_FANCY_TARGET_NAME}
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME bernstein_polynomial_test COMMAND bernstein_polynomial_test)
