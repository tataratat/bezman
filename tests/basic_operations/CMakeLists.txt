set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/operations_test.cpp
)

# Search for the library
set(CMAKE_CXX_STANDARD 17)
add_executable(test_operations ${TEST_SOURCES})
find_package(${MY_FANCY_TARGET_NAME})

# Link bezier manipulation library to the executable target
target_link_libraries(test_operations 
                  PRIVATE bezierManipulation::bezierManipulation
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME test_operations COMMAND test_operations)