set(SOURCES
  edge_edge.cpp
  edge_edge.hpp
  edge_vertex.hpp
  face_vertex.hpp
  normal_collision.cpp
  normal_collision.hpp
  normal_collisions_builder.cpp
  normal_collisions_builder.hpp
  normal_collisions.cpp
  normal_collisions.hpp
  plane_vertex.cpp
  plane_vertex.hpp
  vertex_vertex.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})