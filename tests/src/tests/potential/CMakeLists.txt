set(SOURCES
  # Tests
  test_adhesion_potentials.cpp
  test_barrier_potential.cpp
  test_smooth_potential.cpp
  test_friction_potential.cpp
  test_distance_vector_methods.cpp
  test_full_dof_assembly.cpp
  test_gradient_assembly.cpp
  test_meshfem_assembly.cpp

  # Benchmarks
  benchmark_assembly.cpp
  benchmark_gradient_assembly.cpp
  benchmark_simd_barrier_potential.cpp

  # Utilities
  assembly_scene.cpp
  assembly_scene.hpp
  barrier_potential_bench.hpp
  synthetic_contact_scene.cpp
  synthetic_contact_scene.hpp
)

if(IPC_TOOLKIT_WITH_CUDA)
  list(APPEND SOURCES
    test_gpu_barrier_potential.cu
    test_gpu_psd_projection.cu
    benchmark_cuda_barrier_potential.cu
  )
endif()

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################