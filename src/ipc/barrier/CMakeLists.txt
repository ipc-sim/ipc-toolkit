set(SOURCES
  adaptive_stiffness.cpp
  adaptive_stiffness.hpp
  barrier.cpp
  barrier.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})