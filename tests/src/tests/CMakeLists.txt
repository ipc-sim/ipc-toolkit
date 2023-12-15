set(SOURCES
  main.cpp
  config.hpp

  # Tests
  test_cfl.cpp
  test_collision_mesh.cpp
  test_has_intersections.cpp
  test_ipc.cpp

  # Benchmarks
  benchmark_eigen.cpp

  # Utilities
  utils.cpp
  utils.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(barrier)
add_subdirectory(broad_phase)
add_subdirectory(candidates)
add_subdirectory(ccd)
add_subdirectory(collisions)
add_subdirectory(distance)
add_subdirectory(friction)
add_subdirectory(potential)
add_subdirectory(utils)