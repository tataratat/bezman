set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/composition_test.cpp
)

# Search for the library
set(CMAKE_CXX_STANDARD 17)
add_executable(composition_test ${TEST_SOURCES})

# Link bezier manipulation library to the executable target
include_directories(composition_test ${PROJECT_SOURCE_DIR}/..)
add_dependencies(composition_test ${MY_FANCY_TARGET_NAME})
target_link_libraries(composition_test 
                  PRIVATE ${MY_FANCY_TARGET_NAME}
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME composition_test COMMAND composition_test)