set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/connectivity_test.cpp
)

# Search for the library
set(CMAKE_CXX_STANDARD 17)
add_executable(connectivity_test ${TEST_SOURCES})

# Link bezier manipulation library to the executable target
include_directories(connectivity_test ${PROJECT_SOURCE_DIR}/..)
add_dependencies(connectivity_test ${MY_FANCY_TARGET_NAME})
target_link_libraries(connectivity_test 
                  PRIVATE ${MY_FANCY_TARGET_NAME}
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME connectivity_test COMMAND connectivity_test)