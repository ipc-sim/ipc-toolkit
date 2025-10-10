set(SOURCES
  # Tests
  test_adaptive_stiffness.cpp
  test_barrier.cpp

  # Benchmarks

  # Utilities
)

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################