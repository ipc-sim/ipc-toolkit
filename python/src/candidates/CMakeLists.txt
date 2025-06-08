set(SOURCES
  collision_stencil.cpp
  edge_edge.cpp
  edge_face.cpp
  edge_vertex.cpp
  face_face.cpp
  face_vertex.cpp
  vertex_vertex.cpp
  candidates.cpp
)

target_sources(ipctk PRIVATE ${SOURCES})