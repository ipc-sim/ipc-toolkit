set(SOURCES
  aabb.cpp
  aabb.hpp
  additive_ccd.cpp
  additive_ccd.hpp
  ccd.cpp
  ccd.hpp
  inexact_point_edge.cpp
  inexact_point_edge.hpp
  nonlinear_ccd.cpp
  nonlinear_ccd.hpp
  point_static_plane.cpp
  point_static_plane.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})