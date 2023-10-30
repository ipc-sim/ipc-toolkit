set(SOURCES
  # Tests
  test_distance_type.cpp
  test_edge_edge_mollifier.cpp
  test_edge_edge.cpp
  test_line_line.cpp
  test_point_edge.cpp
  test_point_line.cpp
  test_point_plane.cpp
  test_point_point.cpp
  test_point_triangle.cpp

  # Benchmarks

  # Utilities
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################