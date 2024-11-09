set(SOURCES
  # Tests
  test_closest_point.cpp
  test_relative_velocity.cpp
  test_tangent_basis.cpp

  # Benchmarks

  # Utilities
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################