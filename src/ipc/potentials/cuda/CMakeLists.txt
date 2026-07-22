set(SOURCES
  barrier_potential.cu
  barrier_potential.hpp
  hessian_assembly.cu
  hessian_assembly.cuh
  psd_projection.cu
  psd_projection.cuh
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})
