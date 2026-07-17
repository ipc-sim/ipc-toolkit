set(SOURCES
  adaptive_stiffness.cpp
  adaptive_stiffness.hpp
  barrier_classes.cpp
  barrier_force_magnitude.cpp
  barrier_force_magnitude.hpp
  barrier_type.hpp
  barrier.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})

# barrier.cpp (the free functions) is shared between host C++ and CUDA device
# code; under CUDA the helper generates a .cu wrapper (see
# ipc_toolkit_shared_device_sources). The other barrier .cpp files are ordinary
# host-only sources listed in SOURCES above.
ipc_toolkit_target_shared_device_sources(ipc_toolkit
  barrier.cpp
)