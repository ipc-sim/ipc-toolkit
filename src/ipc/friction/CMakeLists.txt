set(SOURCES
  smooth_friction_mollifier.cpp
  smooth_friction_mollifier.hpp
  smooth_mu.cpp
  smooth_mu.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})