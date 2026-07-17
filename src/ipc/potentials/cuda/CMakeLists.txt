set(SOURCES
  barrier_potential.cu
  barrier_potential.hpp
  local_hessian_assembly.cpp
  local_hessian_assembly.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})
