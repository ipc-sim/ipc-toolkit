set(SOURCES
  autogen.hpp
  autogen.cpp
  mollifier.cpp
  mollifier.hpp
  mollifier.tpp
  point_edge.cpp
  point_edge.hpp
  point_face.cpp
  point_face.hpp
  edge_edge.cpp
  edge_edge.hpp
  primitive_distance.cpp
  primitive_distance.hpp
  primitive_distance.tpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
