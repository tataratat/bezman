set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/bezier_group_test.cpp
)

# Search for the library
add_executable(bezier_group_test ${TEST_SOURCES})

# Link bezier manipulation library to the executable target
add_dependencies(bezier_group_test ${MY_FANCY_TARGET_NAME})
target_link_libraries(bezier_group_test
                  PRIVATE ${MY_FANCY_TARGET_NAME}
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME bezier_group_test COMMAND bezier_group_test)
