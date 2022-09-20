set(SOURCES
  distance_type.hpp
  distance_type.tpp
  edge_edge.hpp
  edge_edge_mollifier.cpp
  edge_edge_mollifier.hpp
  line_line.cpp
  line_line.hpp
  point_edge.hpp
  point_line.cpp
  point_line.hpp
  point_plane.cpp
  point_plane.hpp
  point_point.hpp
  point_triangle.hpp
)

ipc_toolkit_prepend_current_path(SOURCES)
ipc_toolkit_set_source_group(${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
