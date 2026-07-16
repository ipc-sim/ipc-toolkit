set(SOURCES
  # Tests
  test_adaptive_stiffness.cpp
  test_barrier.cpp

  # Benchmarks

  # Utilities
)

if(IPC_TOOLKIT_WITH_CUDA)
  list(APPEND SOURCES
    test_gpu_barrier.cu
  )
endif()

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################