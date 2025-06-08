set(SOURCES
  edge_edge.cpp
  edge_vertex.cpp
  face_vertex.cpp
  tangential_collision.cpp
  tangential_collisions.cpp
  vertex_vertex.cpp
)

target_sources(ipctk PRIVATE ${SOURCES})