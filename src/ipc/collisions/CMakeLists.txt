set(SOURCES
  collision_constraint.cpp
  collision_constraint.hpp
  collision_constraints_builder.cpp
  collision_constraints_builder.hpp
  collision_constraints.cpp
  collision_constraints.hpp
  edge_edge.cpp
  edge_edge.hpp
  # edge_vertex.cpp
  edge_vertex.hpp
  # face_vertex.cpp
  face_vertex.hpp
  plane_vertex.cpp
  plane_vertex.hpp
  # vertex_vertex.cpp
  vertex_vertex.hpp
)

ipc_toolkit_prepend_current_path(SOURCES)
ipc_toolkit_set_source_group(${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
