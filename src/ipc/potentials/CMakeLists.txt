set(SOURCES
  barrier_potential.cpp
  barrier_potential.hpp
  friction_potential.cpp
  friction_potential.hpp
  normal_adhesion_potential.cpp
  normal_adhesion_potential.hpp
  normal_potential.cpp
  normal_potential.hpp
  potential.cpp
  potential.hpp
  tangential_adhesion_potential.cpp
  tangential_adhesion_potential.hpp
  tangential_potential.cpp
  tangential_potential.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})