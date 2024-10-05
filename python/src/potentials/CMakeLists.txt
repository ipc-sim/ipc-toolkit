set(SOURCES
  barrier_potential.cpp
  friction_potential.cpp
  normal_adhesion_potential.cpp
  normal_potential.cpp
  tangential_adhesion_potential.cpp
  tangential_potential.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipctk PRIVATE ${SOURCES})