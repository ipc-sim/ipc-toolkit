set(SOURCES
  # Tests
  test_closest_point.cpp
  test_relative_velocity.cpp
  test_tangent_basis.cpp

  # Benchmarks

  # Utilities
)

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################