set(SOURCES
  collision_candidate.cpp
  collision_candidate.hpp
  hash_grid.cpp
  hash_grid.hpp
  spatial_hash.cpp
  spatial_hash.hpp
)

prepend_current_path(SOURCES)
ipc_toolkit_copy_headers(${SOURCES})
ipc_toolkit_set_source_group(${SOURCES})
target_sources(${PROJECT_NAME} PRIVATE ${SOURCES})
