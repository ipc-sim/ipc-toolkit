set(SOURCES
  inertial_term.cpp
  inertial_term.hpp
  orthogonality_potential.cpp
  orthogonality_potential.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################