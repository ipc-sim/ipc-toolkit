set(SOURCES
  closest_point.hpp
  relative_velocity.hpp
  tangent_basis.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})

# Definitions shared between host C++ and CUDA device code; under CUDA the
# helper generates a .cu wrapper per file (see ipc_toolkit_shared_device_sources).
ipc_toolkit_target_shared_device_sources(ipc_toolkit
  closest_point.cpp
  relative_velocity.cpp
  tangent_basis.cpp
)

################################################################################
# Subfolders
################################################################################