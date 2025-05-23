set(SOURCES
  aabb.cpp
  additive_ccd.cpp
  check_initial_distance.cpp
  inexact_ccd.cpp
  inexact_point_edge.cpp
  narrow_phase_ccd.cpp
  nonlinear_ccd.cpp
  point_static_plane.cpp
  tight_inclusion_ccd.cpp
)

target_sources(ipctk PRIVATE ${SOURCES})