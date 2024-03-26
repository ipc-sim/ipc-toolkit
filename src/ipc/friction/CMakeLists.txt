set(SOURCES
  closest_point.cpp
  closest_point.hpp
  friction_collisions.cpp
  friction_collisions.hpp
  normal_force_magnitude.cpp
  normal_force_magnitude.hpp
  relative_velocity.cpp
  relative_velocity.hpp
  smooth_friction_mollifier.cpp
  smooth_friction_mollifier.hpp
  tangent_basis.cpp
  tangent_basis.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(collisions)