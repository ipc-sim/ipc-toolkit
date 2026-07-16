set(SOURCES
  adaptive_stiffness.cpp
  adaptive_stiffness.hpp
  barrier_classes.cpp
  barrier_force_magnitude.cpp
  barrier_force_magnitude.hpp
  barrier.hpp
)

# With CUDA enabled, barrier.cpp (the free functions only) is textually
# included into barrier.hpp (see the IPC_TOOLKIT_WITH_CUDA footer) so it can be
# compiled as __host__ __device__ code; compiling it separately would produce
# duplicate symbols.
if(NOT IPC_TOOLKIT_WITH_CUDA)
  list(APPEND SOURCES
    barrier.cpp
  )
endif()

target_sources(ipc_toolkit PRIVATE ${SOURCES})