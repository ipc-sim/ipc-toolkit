set(SOURCES
  test_body_forces.cpp
  test_ground_contact.cpp
  test_inertial_term.cpp
  test_mass.cpp
  test_pose.cpp
  test_rigid_body.cpp
  test_simulator.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################
