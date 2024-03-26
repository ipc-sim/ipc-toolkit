set(SOURCES
  area_gradient.cpp
  area_gradient.hpp
  eigen_ext.hpp
  eigen_ext.tpp
  intersection.cpp
  intersection.hpp
  interval.cpp
  interval.hpp
  local_to_global.hpp
  logger.cpp
  logger.hpp
  merge_thread_local.hpp
  save_obj.cpp
  save_obj.hpp
  unordered_map_and_set.cpp
  unordered_map_and_set.hpp
  vertex_to_min_edge.cpp
  vertex_to_min_edge.hpp
  world_bbox_diagonal_length.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})