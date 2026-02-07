set(SOURCES
  aabb.cpp
  broad_phase.cpp
  brute_force.cpp
  hash_grid.cpp
  lbvh.cpp
  spatial_hash.cpp
  sweep_and_prune.cpp
  sweep_and_tiniest_queue.cpp
  voxel_size_heuristic.cpp
)

target_sources(ipctk PRIVATE ${SOURCES})
