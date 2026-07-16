set(SOURCES
  # Tests
  test_closest_point.cpp
  test_relative_velocity.cpp
  test_tangent_basis.cpp

  # Benchmarks

  # Utilities
)

if(IPC_TOOLKIT_WITH_CUDA)
  list(APPEND SOURCES
    test_gpu_tangent.cu
  )
endif()

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################