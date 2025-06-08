set(SOURCES
  candidates.cpp
  candidates.hpp
  collision_stencil.cpp
  collision_stencil.hpp
  edge_edge.cpp
  edge_edge.hpp
  edge_face.cpp
  edge_face.hpp
  edge_vertex.cpp
  edge_vertex.hpp
  face_face.cpp
  face_face.hpp
  face_vertex.cpp
  face_vertex.hpp
  # plane_vertex.cpp
  # plane_vertex.hpp
  vertex_vertex.cpp
  vertex_vertex.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})