set(SOURCES
  affine_dof.cpp
  affine_dof.hpp
  augmented_lagrangian.cpp
  augmented_lagrangian.hpp
  body_forces.cpp
  body_forces.hpp
  inertial_term.cpp
  inertial_term.hpp
  joints.cpp
  joints.hpp
  orthogonality_potential.cpp
  orthogonality_potential.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################