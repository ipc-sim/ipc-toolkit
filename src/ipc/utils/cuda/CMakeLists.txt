set(SOURCES
  device_buffer.cuh
  device_utils.cuh
  is_available.cu
  is_available.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})
