set(SOURCES
  line_line.hpp
  point_line.hpp
  point_plane.hpp
)

# With CUDA enabled, these implementations are textually included into their
# headers (see the IPC_TOOLKIT_WITH_CUDA footer in each header) so they can be
# compiled as __host__ __device__ code; compiling them separately would produce
# duplicate symbols.
if(NOT IPC_TOOLKIT_WITH_CUDA)
  list(APPEND SOURCES
    line_line.cpp
    point_line.cpp
    point_plane.cpp
  )
endif()

target_sources(ipc_toolkit PRIVATE ${SOURCES})