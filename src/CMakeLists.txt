set(SOURCES
  collision_constraint.cpp
  collision_constraint.hpp
  ipc.hpp
  ipc.cpp
)

ipc_toolkit_prepend_current_path(SOURCES)
ipc_toolkit_copy_headers(${SOURCES})
ipc_toolkit_set_source_group(${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(barrier)
add_subdirectory(broad_phase)
add_subdirectory(ccd)
add_subdirectory(distance)
add_subdirectory(friction)
add_subdirectory(implicits)
add_subdirectory(utils)
