set(SOURCES
  edge_edge.cpp
  edge_edge.hpp
  edge_vertex.cpp
  edge_vertex.hpp
  face_vertex.cpp
  face_vertex.hpp
  tangential_collision.cpp
  tangential_collision.hpp
  tangential_collisions.cpp
  tangential_collisions.hpp
  vertex_vertex.cpp
  vertex_vertex.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})