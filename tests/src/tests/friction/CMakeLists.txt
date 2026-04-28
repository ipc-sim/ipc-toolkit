set(SOURCES
  # Tests
  test_anisotropic_friction.cpp
  test_force_jacobian.cpp
  test_friction.cpp
  test_smooth_friction_mollifier.cpp
  test_smooth_mu.cpp

  # Benchmarks

  # Utilities
  friction_data_generator.cpp
  friction_data_generator.hpp
)

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################