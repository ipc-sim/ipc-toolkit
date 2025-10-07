set(SOURCES
  # Tests
  test_feasible_region.cpp
  test_trust_region.cpp

  # Benchmarks

  # Utilities
)

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################