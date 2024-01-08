set(SOURCES
  smooth_point_edge.cpp
  smooth_point_edge.hpp
  smooth_point_face.cpp
  smooth_point_face.hpp
  smooth_contact_potential.cpp
  smooth_contact_potential.hpp
  smooth_collisions.cpp
  smooth_collisions.hpp
  smooth_collisions_builder.cpp
  smooth_collisions_builder.hpp
  edge_vertex.cpp
  edge_vertex.hpp
  edge_edge.cpp
  edge_edge.hpp
  face_vertex.cpp
  face_vertex.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
