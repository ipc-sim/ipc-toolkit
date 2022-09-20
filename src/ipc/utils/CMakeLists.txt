set(SOURCES
  area_gradient.cpp
  area_gradient.hpp
  eigen_ext.hpp
  eigen_ext.tpp
  intersection.cpp
  intersection.hpp
  local_to_global.hpp
  logger.cpp
  logger.hpp
  merge_thread_local_vectors.hpp
  rational.hpp
  save_obj.cpp
  save_obj.hpp
  unordered_map_and_set.hpp
  vertex_to_min_edge.cpp
  vertex_to_min_edge.hpp
  world_bbox_diagonal_length.hpp
)

ipc_toolkit_prepend_current_path(SOURCES)
ipc_toolkit_set_source_group(${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
