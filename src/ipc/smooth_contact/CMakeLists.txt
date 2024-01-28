set(SOURCES
  common.hpp
  smooth_contact_potential.cpp
  smooth_contact_potential.hpp
  smooth_collisions.cpp
  smooth_collisions.hpp
  smooth_collisions_builder.cpp
  smooth_collisions_builder.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(collisions)
add_subdirectory(distance)
add_subdirectory(pairs)
add_subdirectory(primitives)
