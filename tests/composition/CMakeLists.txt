set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/composition_test.cpp
)

# Search for the library
set(CMAKE_CXX_STANDARD 17)
add_executable(composition_test ${TEST_SOURCES})
find_package(${MY_FANCY_TARGET_NAME})

# Link bezier manipulation library to the executable target
target_link_libraries(composition_test 
                  PRIVATE bezierManipulation::bezierManipulation
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME composition_test COMMAND composition_test)