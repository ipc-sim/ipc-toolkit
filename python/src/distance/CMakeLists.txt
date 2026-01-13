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
  signed_distance.cpp
)

target_sources(ipctk PRIVATE ${SOURCES})