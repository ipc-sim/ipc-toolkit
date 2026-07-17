set(SOURCES
  device_normal_collisions.cu
  device_normal_collisions.hpp
  device_normal_collisions_impl.cuh
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})
