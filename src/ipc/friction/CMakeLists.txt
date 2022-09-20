set(SOURCES
  closest_point.cpp
  closest_point.hpp
  friction.cpp
  friction.hpp
  friction.tpp
  friction_constraint.cpp
  friction_constraint.hpp
  normal_force_magnitude.cpp
  normal_force_magnitude.hpp
  relative_displacement.hpp
  smooth_friction_mollifier.hpp
  smooth_friction_mollifier.tpp
  tangent_basis.cpp
  tangent_basis.hpp
)

ipc_toolkit_prepend_current_path(SOURCES)
ipc_toolkit_set_source_group(${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
