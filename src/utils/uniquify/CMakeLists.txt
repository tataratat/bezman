# Add Files
set(LOCAL_HEADERS 
  ${CMAKE_CURRENT_SOURCE_DIR}/point_uniquifier.hpp
)

# Add All Subdirectories

# Install files in according file structure
install(FILES ${LOCAL_HEADERS} DESTINATION ${include_destination_headers}/src/utils/uniquify)

# Update Headers
set(HEADERS ${HEADERS} ${LOCAL_HEADERS} PARENT_SCOPE)