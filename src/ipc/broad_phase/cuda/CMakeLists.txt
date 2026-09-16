set(SOURCES
  lbvh.cu
  lbvh.hpp
  lbvh_impl.cuh
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})
