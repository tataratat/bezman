set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/algo_diff_type_test.cpp
)

# Search for the library
set(CMAKE_CXX_STANDARD 17)
add_executable(algo_diff_type_test ${TEST_SOURCES})

# Link bezier manipulation library to the executable target
include_directories(algo_diff_type_test ${PROJECT_SOURCE_DIR}/..)
add_dependencies(algo_diff_type_test ${MY_FANCY_TARGET_NAME})
target_link_libraries(algo_diff_type_test 
                  PRIVATE ${MY_FANCY_TARGET_NAME}
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME algo_diff_type_test COMMAND algo_diff_type_test)
