set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/order_and_derivation_test.cpp
)

# Search for the library
add_executable(order_and_derivation_test ${TEST_SOURCES})

# Link bezier manipulation library to the executable target
add_dependencies(order_and_derivation_test ${MY_FANCY_TARGET_NAME})
target_link_libraries(order_and_derivation_test 
                  PRIVATE ${MY_FANCY_TARGET_NAME}
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME order_and_derivation_test COMMAND order_and_derivation_test)
