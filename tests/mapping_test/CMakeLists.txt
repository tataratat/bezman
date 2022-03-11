set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/mapping_test.cpp
)

# Search for the library
set(CMAKE_CXX_STANDARD 17)
add_executable(mapping_test ${TEST_SOURCES})

# Link bezier manipulation library to the executable target
include_directories(mapping_test ${PROJECT_SOURCE_DIR}/..)
add_dependencies(mapping_test ${MY_FANCY_TARGET_NAME})
target_link_libraries(mapping_test 
                  PRIVATE ${MY_FANCY_TARGET_NAME}
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME mapping_test COMMAND mapping_test)