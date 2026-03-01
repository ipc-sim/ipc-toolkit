set(SOURCES
  # Tests
  test_local_to_global.cpp
  test_matrixcache.cpp
  test_sinc.cpp
  test_utils.cpp

  # Benchmarks

  # Utilities
)

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################