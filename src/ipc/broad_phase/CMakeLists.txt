set(SOURCES
  aabb.cpp
  aabb.hpp
  broad_phase.cpp
  broad_phase.hpp
  brute_force.cpp
  brute_force.hpp
  bvh.cpp
  bvh.hpp
  default_broad_phase.hpp
  hash_grid.cpp
  hash_grid.hpp
  spatial_hash.cpp
  spatial_hash.hpp
  sweep_and_prune.cpp
  sweep_and_prune.hpp
  sweep_and_tiniest_queue.cpp
  sweep_and_tiniest_queue.hpp
  voxel_size_heuristic.cpp
  voxel_size_heuristic.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})