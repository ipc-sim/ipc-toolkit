set(SOURCES
  smooth_point_edge.cpp
  smooth_point_edge.hpp
  smooth_point_face.cpp
  smooth_point_face.hpp
  smooth_edge_edge.cpp
  smooth_edge_edge.hpp
  smooth_contact_potential.cpp
  smooth_contact_potential.hpp
  smooth_collisions.cpp
  smooth_collisions.hpp
  smooth_collisions_builder.cpp
  smooth_collisions_builder.hpp
  edge_edge.cpp
  edge_edge.hpp
  edge_edge_face.cpp
  edge_edge_face.hpp
  face_face.cpp
  face_face.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
