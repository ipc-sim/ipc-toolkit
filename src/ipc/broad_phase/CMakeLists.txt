set(SOURCES
  aabb.cpp
  aabb.hpp
  broad_phase.cpp
  broad_phase.hpp
  brute_force.cpp
  brute_force.hpp
  hash_grid.cpp
  hash_grid.hpp
  spatial_hash.cpp
  spatial_hash.hpp
  sweep_and_tiniest_queue.cpp
  sweep_and_tiniest_queue.hpp
  voxel_size_heuristic.cpp
  voxel_size_heuristic.hpp
  broadmark_util.hpp
  broadmark_util.cpp
  broadmark_util.tpp
  broadmark.hpp
  broadmark.cpp
)

ipc_toolkit_prepend_current_path(SOURCES)
ipc_toolkit_set_source_group(${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
