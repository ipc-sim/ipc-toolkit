set(SOURCES
  # Tests
  test_adhesion_potentials.cpp
  test_barrier_potential.cpp
  test_smooth_potential.cpp
  test_friction_potential.cpp
  test_distance_vector_methods.cpp
  ipc_test_ee_near_parallel_shape_derivative.cpp

  # Benchmarks

  # Utilities
)

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################