set(SOURCES
  closest_point.cpp
  closest_point.hpp
  relative_velocity.cpp
  relative_velocity.hpp
  tangent_basis.cpp
  tangent_basis.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################