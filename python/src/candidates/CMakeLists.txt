set(SOURCES
  collision_stencil.cpp
  continuous_collision_candidate.cpp
  edge_edge.cpp
  edge_face.cpp
  edge_vertex.cpp
  face_vertex.cpp
  vertex_vertex.cpp
  candidates.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipctk PRIVATE ${SOURCES})