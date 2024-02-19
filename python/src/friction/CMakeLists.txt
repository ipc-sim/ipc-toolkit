set(SOURCES
  closest_point.cpp
  friction_collisions.cpp
  normal_force_magnitude.cpp
  relative_velocity.cpp
  smooth_friction_mollifier.cpp
  tangent_basis.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipctk PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(collisions)