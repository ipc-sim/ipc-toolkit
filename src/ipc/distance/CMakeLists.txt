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

target_sources(ipc_toolkit PRIVATE ${SOURCES})

# These definitions are shared between host C++ and CUDA device code. The
# headers only declare them (IPC_TOOLKIT_HOST_DEVICE), so editing a definition
# rebuilds one TU plus a link step rather than every TU that includes the
# header. Under CUDA the helper generates a .cu wrapper per file so nvcc emits
# relocatable device code linkable from other TUs.
ipc_toolkit_target_shared_device_sources(ipc_toolkit
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

################################################################################
# Subfolders
################################################################################

add_subdirectory(signed)