set(SOURCES
  # Tests
  test_adhesion.cpp
  test_simd_adhesion.cpp

  # Benchmarks

  # Utilities
)

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################