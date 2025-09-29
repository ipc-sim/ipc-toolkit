set(SOURCES
  # Tests
  test_adaptive_stiffness.cpp
  test_barrier.cpp

  # Benchmarks

  # Utilities
  test_eigen.cpp
)

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################