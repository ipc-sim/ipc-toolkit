set(SOURCES
  common.hpp
  gcp_potential.cpp
  gcp_potential.hpp
  gcp_collisions.cpp
  gcp_collisions.hpp
  gcp_collisions_builder.cpp
  gcp_collisions_builder.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(collisions)
add_subdirectory(distance)
add_subdirectory(primitives)
