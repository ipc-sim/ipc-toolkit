set(SOURCES
  ccd.cpp
  aabb.cpp
  additive_ccd.cpp
  inexact_point_edge.cpp
  nonlinear_ccd.cpp
  point_static_plane.cpp
  tight_inclusion_ccd.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipctk PRIVATE ${SOURCES})