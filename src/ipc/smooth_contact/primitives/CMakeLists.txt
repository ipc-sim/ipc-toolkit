set(SOURCES
  autogen.cpp
  autogen.hpp
  edge2.cpp
  edge2.hpp
  edge3.cpp
  edge3.hpp
  face.cpp
  face.hpp
  primitive.hpp
  point3.cpp
  point3.hpp
  point2.cpp
  point2.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
