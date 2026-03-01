set(SOURCES
  # Tests
  test_utils.cpp
  test_matrixcache.cpp
  test_local_to_global.cpp

  # Benchmarks

  # Utilities
)

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################