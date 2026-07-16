set(SOURCES
  angle.hpp
  area.hpp
  intersection.cpp
  intersection.hpp
  normal.hpp
)

# With CUDA enabled, these implementations are textually included into their
# headers (see the IPC_TOOLKIT_WITH_CUDA footer in each header) so they can be
# compiled as __host__ __device__ code; compiling them separately would produce
# duplicate symbols.
if(NOT IPC_TOOLKIT_WITH_CUDA)
  list(APPEND SOURCES
    angle.cpp
    area.cpp
    normal.cpp
  )
endif()

target_sources(ipc_toolkit PRIVATE ${SOURCES})

# Definitions shared between host C++ and CUDA device code; under CUDA the
# helper generates a .cu wrapper per file (see ipc_toolkit_shared_device_sources).
ipc_toolkit_target_shared_device_sources(ipc_toolkit
  angle.cpp
  area.cpp
  normal.cpp
)

################################################################################
# Subfolders
################################################################################