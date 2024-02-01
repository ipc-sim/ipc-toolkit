set(SOURCES
  edge_vertex.cpp
  edge_vertex.hpp
  # edge_edge_3d.cpp
  # edge_edge_3d.hpp
  # edge_vertex_3d.cpp
  # edge_vertex_3d.hpp
  # face_vertex.cpp
  # face_vertex.hpp
  vertex_vertex.cpp
  vertex_vertex.hpp
  # vertex_vertex_3d.cpp
  # vertex_vertex_3d.hpp
  smooth_collision.cpp
  smooth_collision.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
