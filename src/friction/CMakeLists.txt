set(SOURCES
  closest_point.hpp
  friction.cpp
  friction.hpp
  friction.tpp
  friction_constraint.cpp
  friction_constraint.hpp
  relative_displacement.hpp
  tangent_basis.hpp
)

ipc_toolkit_prepend_current_path(SOURCES)
ipc_toolkit_copy_headers(${SOURCES})
ipc_toolkit_set_source_group(${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
