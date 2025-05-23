set(SOURCES
  smooth_friction_mollifier.hpp
  smooth_friction_mollifier.cpp
  material_friction.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})