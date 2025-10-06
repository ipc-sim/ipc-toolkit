set(SOURCES
  edge_edge.cpp
  edge_vertex.cpp
  face_vertex.cpp
  normal_collision.cpp
  normal_collisions.cpp
  plane_vertex.cpp
  vertex_vertex.cpp
)

target_sources(ipctk PRIVATE ${SOURCES})