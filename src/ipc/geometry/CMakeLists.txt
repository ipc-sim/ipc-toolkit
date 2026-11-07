set(SOURCES
  angle.hpp
  area.hpp
  intersection.cpp
  intersection.hpp
  normal.hpp
)

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