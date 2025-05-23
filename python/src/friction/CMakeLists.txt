set(SOURCES
  bindings.hpp
  smooth_friction_mollifier.cpp
  material_friction.cpp
)

target_sources(ipctk PRIVATE ${SOURCES})
