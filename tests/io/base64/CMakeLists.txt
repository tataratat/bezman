set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/base64_test.cpp
)

# Search for the library
add_executable(base64_test ${TEST_SOURCES})

# Link bezier manipulation library to the executable target
add_dependencies(base64_test ${MY_FANCY_TARGET_NAME})
target_link_libraries(base64_test 
                  PRIVATE ${MY_FANCY_TARGET_NAME}
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME base64_test COMMAND base64_test)