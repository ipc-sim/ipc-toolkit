set(SOURCES
  eigen_ext.cpp
  eigen_ext.hpp
  eigen_ext.tpp
  intersection.cpp
  intersection.hpp
  local_to_global.hpp
  logger.cpp
  logger.hpp
  rational.hpp
  unordered_map_and_set.hpp
  world_bbox_diagonal_length.hpp
)

prepend_current_path(SOURCES)
ipc_toolkit_copy_headers(${SOURCES})
ipc_toolkit_set_source_group(${SOURCES})
target_sources(${PROJECT_NAME} PRIVATE ${SOURCES})
