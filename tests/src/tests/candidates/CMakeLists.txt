set(SOURCES
  # Tests
  test_candidates.cpp
  test_coefficients.cpp
  test_normals.cpp

  # Benchmarks

  # Utilities
)

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################