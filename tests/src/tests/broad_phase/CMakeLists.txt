set(SOURCES
  # Tests
  test_aabb.cpp
  test_broad_phase.cpp
  test_spatial_hash.cpp
  test_stq.cpp
  test_voxel_size_heuristic.cpp

  # Benchmarks
  benchmark_broad_phase.cpp

  # Utilities
  brute_force_comparison.cpp
  brute_force_comparison.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################