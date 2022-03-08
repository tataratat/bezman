# Add files
set(LOCAL_HEADERS 
      ${CMAKE_CURRENT_SOURCE_DIR}/bezier.hpp
      ${CMAKE_CURRENT_SOURCE_DIR}/point.hpp
    )

# Add All Subdirectories
add_subdirectory(utils) # utilities

# Install files in according file structure
install(FILES ${LOCAL_HEADERS} DESTINATION ${include_destination_headers}/src)

# Update Headers
set(HEADERS ${HEADERS} ${LOCAL_HEADERS} PARENT_SCOPE)