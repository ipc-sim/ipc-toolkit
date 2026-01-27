set(SOURCES
  autodiff_types.hpp
  autodiff_types.hpp
  default_init_allocator.hpp
  eigen_ext.hpp
  eigen_ext.tpp
  local_to_global.hpp
  logger.cpp
  logger.hpp
  matrix_cache.cpp
  matrix_cache.hpp
  maybe_parallel_for.hpp
  maybe_parallel_for.tpp
  merge_thread_local.hpp
  profiler.cpp
  profiler.hpp
  unordered_map_and_set.cpp
  unordered_map_and_set.hpp
  vertex_to_min_edge.cpp
  vertex_to_min_edge.hpp
  world_bbox_diagonal_length.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})
