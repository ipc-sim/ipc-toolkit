set(SOURCES
  # Tests
  test_ccd.cpp
  test_gpu_ccd.cpp
  test_ccd_benchmark.cpp
  test_edge_edge_ccd.cpp
  test_nonlinear_ccd.cpp
  test_point_edge_ccd.cpp
  test_point_point_ccd.cpp
  test_point_triangle_ccd.cpp

  # Benchmarks
  benchmark_ccd.cpp

  # Utilities
  collision_generator.cpp
  collision_generator.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################