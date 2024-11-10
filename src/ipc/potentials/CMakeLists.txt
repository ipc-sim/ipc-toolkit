set(SOURCES
  barrier_potential.cpp
  barrier_potential.hpp
  friction_potential.cpp
  friction_potential.hpp
  normal_adhesion_potential.cpp
  normal_adhesion_potential.hpp
  normal_potential.cpp
  normal_potential.hpp
  potential.hpp
  potential.tpp
  tangential_adhesion_potential.cpp
  tangential_adhesion_potential.hpp
  tangential_potential.cpp
  tangential_potential.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})