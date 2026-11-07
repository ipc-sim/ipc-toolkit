set(SOURCES
  # Tests
  test_cuda_is_available.cpp
  test_local_to_global.cpp
  test_matrixcache.cpp
  test_vertex_matrix_view.cpp
  test_profiler.cpp
  test_simd_utils.cpp
  test_sinc.cpp
  test_utils.cpp

  # Benchmarks

  # Utilities
)

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################