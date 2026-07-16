set(SOURCES
  distance_type.hpp
  edge_edge.hpp
  edge_edge_mollifier.hpp
  line_line.hpp
  point_edge.hpp
  point_line.hpp
  point_plane.hpp
  point_point.hpp
  point_triangle.hpp
)

# With CUDA enabled, these implementations are textually included into their
# headers (see the IPC_TOOLKIT_WITH_CUDA footer in each header) so they can be
# compiled as __host__ __device__ code; compiling them separately would produce
# duplicate symbols.
if(NOT IPC_TOOLKIT_WITH_CUDA)
  list(APPEND SOURCES
    distance_type.cpp
    edge_edge.cpp
    edge_edge_mollifier.cpp
    line_line.cpp
    point_edge.cpp
    point_line.cpp
    point_plane.cpp
    point_point.cpp
    point_triangle.cpp
  )
endif()

target_sources(ipc_toolkit PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(signed)