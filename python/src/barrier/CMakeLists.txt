set(SOURCES
  adaptive_stiffness.cpp
  barrier_force_magnitude.cpp
  barrier.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipctk PRIVATE ${SOURCES})