set(SOURCES
  inertial_term.cpp
  inertial_term.hpp
  mass.cpp
  mass.hpp
  pose.cpp
  pose.hpp
  rigid_bodies.cpp
  rigid_bodies.hpp
  rigid_body.cpp
  rigid_body.hpp
  simulator.cpp
  time_integrator.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################