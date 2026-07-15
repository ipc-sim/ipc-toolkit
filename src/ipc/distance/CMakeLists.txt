set(SOURCES
  distance_type.cpp
  distance_type.hpp
  edge_edge.cpp
  edge_edge.hpp
  edge_edge_mollifier.cpp
  edge_edge_mollifier.hpp
  line_line.cpp
  line_line.hpp
  point_edge.cpp
  point_edge.hpp
  point_line.cpp
  point_line.hpp
  point_plane.cpp
  point_plane.hpp
  point_point.hpp
  point_triangle.cpp
  point_triangle.hpp
)

if(NOT IPC_TOOLKIT_WITH_CUDA)
  list(APPEND SOURCES
    point_point.cpp
  )
endif()

target_sources(ipc_toolkit PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(signed)