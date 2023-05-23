set(
    TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/composition_sensitivities.cpp
)

# Search for the library
add_executable(composition_sensitivities ${TEST_SOURCES})

# Link bezier manipulation library to the executable target
add_dependencies(composition_sensitivities ${MY_FANCY_TARGET_NAME})
target_link_libraries(composition_sensitivities 
                  PRIVATE ${MY_FANCY_TARGET_NAME}
                  PRIVATE gtest_main)

include(GoogleTest)
add_test(NAME composition_sensitivities COMMAND composition_sensitivities)