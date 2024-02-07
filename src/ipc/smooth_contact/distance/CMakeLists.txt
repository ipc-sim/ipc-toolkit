set(SOURCES
  mollifier.cpp
  mollifier.hpp
  mollifier.tpp
  point_edge.cpp
  point_edge.hpp
  point_edge.tpp
  point_face.hpp
  point_face.tpp
  edge_edge.cpp
  edge_edge.hpp
  edge_edge.tpp
  primitive_distance.cpp
  primitive_distance.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
