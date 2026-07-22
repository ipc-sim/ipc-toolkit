set(SOURCES
  # Tests
  test_adhesion_potentials.cpp
  test_barrier_potential.cpp
  test_smooth_potential.cpp
  test_friction_potential.cpp
  test_distance_vector_methods.cpp

  # Benchmarks

  # Utilities
)

if(IPC_TOOLKIT_WITH_CUDA)
  list(APPEND SOURCES
    test_gpu_barrier_potential.cu
    test_gpu_psd_projection.cu
  )
endif()

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################