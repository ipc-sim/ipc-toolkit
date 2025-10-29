set(SOURCES
  eigen_ext.hpp
  eigen_ext.tpp
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
  autodiff_types.hpp
  maybe_parallel_for.hpp
  maybe_parallel_for.tpp
  matrix_cache.hpp
  matrix_cache.cpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})