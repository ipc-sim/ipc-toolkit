set(SOURCES
  test_angle.cpp
)

if(IPC_TOOLKIT_WITH_CUDA)
  list(APPEND SOURCES
    test_gpu_geometry.cu
  )
endif()

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################