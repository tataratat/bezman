set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/operations_test.cpp
)

# Search for the library
set(CMAKE_CXX_STANDARD 17)
add_executable(test_operations ${TEST_SOURCES})

# Link bezier manipulation library to the executable target
include_directories(test_operations ${PROJECT_SOURCE_DIR}/..)
add_dependencies(test_operations ${MY_FANCY_TARGET_NAME})
target_link_libraries(test_operations 
                  PRIVATE ${MY_FANCY_TARGET_NAME}
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME test_operations COMMAND test_operations)
