set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/point_test.cpp
)

# Search for the library
add_executable(point_test ${TEST_SOURCES})

# Link bezier manipulation library to the executable target
add_dependencies(point_test ${MY_FANCY_TARGET_NAME})
target_link_libraries(point_test 
                  PRIVATE ${MY_FANCY_TARGET_NAME}
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME point_test COMMAND point_test)