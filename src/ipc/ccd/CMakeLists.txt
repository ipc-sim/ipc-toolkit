set(SOURCES
  aabb.cpp
  aabb.hpp
  additive_ccd.cpp
  additive_ccd.hpp
  check_initial_distance.hpp
  default_narrow_phase_ccd.cpp
  default_narrow_phase_ccd.hpp
  inexact_ccd.cpp
  inexact_ccd.hpp
  inexact_point_edge.cpp
  inexact_point_edge.hpp
  nonlinear_ccd.cpp
  nonlinear_ccd.hpp
  point_static_plane.cpp
  point_static_plane.hpp
  tight_inclusion_ccd.cpp
  tight_inclusion_ccd.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})