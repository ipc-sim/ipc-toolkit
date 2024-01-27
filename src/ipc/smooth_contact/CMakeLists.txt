set(SOURCES
  common.hpp
  edge.hpp
  face.hpp
  point.hpp
  smooth_point_point.hpp
  smooth_point_edge.hpp
  smooth_point_face.hpp
  smooth_edge_edge.hpp
  smooth_contact_potential.cpp
  smooth_contact_potential.hpp
  smooth_collisions.cpp
  smooth_collisions.hpp
  smooth_collisions_builder.cpp
  smooth_collisions_builder.hpp
  edge_vertex.cpp
  edge_vertex.hpp
  edge_edge_3d.cpp
  edge_edge_3d.hpp
  edge_vertex_3d.cpp
  edge_vertex_3d.hpp
  face_vertex.cpp
  face_vertex.hpp
  vertex_vertex.cpp
  vertex_vertex.hpp
  vertex_vertex_3d.cpp
  vertex_vertex_3d.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
