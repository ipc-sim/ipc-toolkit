set(SOURCES
  line_line.hpp
  point_line.hpp
  point_plane.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})

# Definitions shared between host C++ and CUDA device code; under CUDA the
# helper generates a .cu wrapper per file (see ipc_toolkit_shared_device_sources).
ipc_toolkit_target_shared_device_sources(ipc_toolkit
  line_line.cpp
  point_line.cpp
  point_plane.cpp
)