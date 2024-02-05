set(SOURCES
  edge_vertex.cpp
  edge_vertex.hpp
  vertex_vertex.cpp
  vertex_vertex.hpp
  smooth_collision.cpp
  smooth_collision.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
