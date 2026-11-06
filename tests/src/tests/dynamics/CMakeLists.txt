set(SOURCES
  test_bdf.cpp
  test_body_potentials.cpp
  test_to_affine.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(affine)
add_subdirectory(rigid)