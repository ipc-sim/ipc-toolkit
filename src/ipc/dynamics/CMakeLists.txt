set(SOURCES
  body_potentials.cpp
  body_potentials.hpp
  to_affine.cpp
  to_affine.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(affine)
add_subdirectory(rigid)
add_subdirectory(time_integration)
