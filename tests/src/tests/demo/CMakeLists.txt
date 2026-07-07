# Tests for the demo simulator (demo/); only built when the demo is.
if(NOT IPC_TOOLKIT_BUILD_DEMO)
  return()
endif()

set(SOURCES
  test_affine_simulator.cpp
  test_joints.cpp
  test_rigid_simulator.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})
