set(SOURCES
  # Tests
  test_adhesion_potentials.cpp
  test_barrier_potential.cpp
  test_smooth_potential.cpp
  test_friction_potential.cpp
  test_distance_vector_methods.cpp

  # Benchmarks
  benchmark_assembly.cpp

  # Utilities
  assembly_scene.cpp
  assembly_scene.hpp
)

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################