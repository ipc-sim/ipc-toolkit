set(SOURCES
  # Tests
  test_normal_collisions.cpp
  test_plane_vertex_collisions.cpp

  # Benchmarks

  # Utilities
)

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################