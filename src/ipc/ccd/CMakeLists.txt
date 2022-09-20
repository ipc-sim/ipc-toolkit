set(SOURCES
  aabb.cpp
  aabb.hpp
  ccd.cpp
  ccd.hpp
  inexact_point_edge.cpp
  inexact_point_edge.hpp
  point_static_plane.cpp
  point_static_plane.hpp
)

ipc_toolkit_prepend_current_path(SOURCES)
ipc_toolkit_set_source_group(${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
