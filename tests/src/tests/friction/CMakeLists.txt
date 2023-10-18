set(SOURCES
  # Tests
  test_closest_point.cpp
  test_force_jacobian.cpp
  test_friction.cpp
  test_normal_force_magnitude.cpp
  test_relative_velocity.cpp
  test_smooth_friction_mollifier.cpp
  test_tangent_basis.cpp

  # Benchmarks

  # Utilities
  friction_data_generator.cpp
  friction_data_generator.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################