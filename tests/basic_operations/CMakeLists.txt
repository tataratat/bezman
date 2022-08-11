set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/operations_test.cpp
)

# Search for the library
add_executable(operations_test ${TEST_SOURCES})

# Link bezier manipulation library to the executable target
add_dependencies(operations_test ${MY_FANCY_TARGET_NAME})
target_link_libraries(operations_test 
                  PRIVATE ${MY_FANCY_TARGET_NAME}
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME operations_test COMMAND operations_test)
