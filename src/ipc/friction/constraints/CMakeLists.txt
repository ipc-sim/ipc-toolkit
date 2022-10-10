set(SOURCES
  edge_edge.cpp
  edge_edge.hpp
  edge_vertex.cpp
  edge_vertex.hpp
  face_vertex.cpp
  face_vertex.hpp
  friction_constraint.cpp
  friction_constraint.hpp
  vertex_vertex.cpp
  vertex_vertex.hpp
)

ipc_toolkit_prepend_current_path(SOURCES)
ipc_toolkit_set_source_group(${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
