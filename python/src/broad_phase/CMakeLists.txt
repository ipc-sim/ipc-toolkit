set(SOURCES
  aabb.cpp
  broad_phase.cpp
  brute_force.cpp
  bvh.cpp
  hash_grid.cpp
  spatial_hash.cpp
  sweep_and_prune.cpp
  sweep_and_tiniest_queue.cpp
  voxel_size_heuristic.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipctk PRIVATE ${SOURCES})