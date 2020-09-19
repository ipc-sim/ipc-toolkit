set(SOURCES
  collision_constraint.cpp
  collision_constraint.hpp
  ipc.hpp
  ipc.cpp
)

prepend_current_path(SOURCES)
ipc_toolkit_copy_headers(${SOURCES})
ipc_toolkit_set_source_group(${SOURCES})
target_sources(${PROJECT_NAME} PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(barrier)
add_subdirectory(ccd)
add_subdirectory(distance)
add_subdirectory(friction)
add_subdirectory(spatial_hash)
add_subdirectory(utils)
