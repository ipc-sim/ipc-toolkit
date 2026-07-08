set(SOURCES
  test_mass.cpp
  test_pose.cpp
  test_rigid_bodies.cpp
  test_rigid_body.cpp
  test_rigid_candidates.cpp
  test_rigid_trajectory.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################
