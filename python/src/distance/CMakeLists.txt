set(SOURCES
  distance_type.cpp
  edge_edge_mollifier.cpp
  edge_edge.cpp
  line_line.cpp
  point_edge.cpp
  point_line.cpp
  point_plane.cpp
  point_point.cpp
  point_triangle.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipctk PRIVATE ${SOURCES})