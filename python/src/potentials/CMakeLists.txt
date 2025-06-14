set(SOURCES
  barrier_potential.cpp
  friction_potential.cpp
  normal_adhesion_potential.cpp
  normal_potential.cpp
  tangential_adhesion_potential.cpp
  tangential_potential.cpp
)

target_sources(ipctk PRIVATE ${SOURCES})