set(SOURCES
  edge.cpp
  edge.hpp
  face.cpp
  face.hpp
  primitive.hpp
  point.cpp
  point.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
