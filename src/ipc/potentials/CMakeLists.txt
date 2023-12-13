set(SOURCES
  barrier_potential.cpp
  barrier_potential.hpp
  distance_based_potential.cpp
  distance_based_potential.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})