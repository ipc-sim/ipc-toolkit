set(SOURCES
  adaptive_stiffness.cpp
  adaptive_stiffness.hpp
  barrier.cpp
  barrier.hpp
  barrier.tpp
)

ipc_toolkit_prepend_current_path(SOURCES)
ipc_toolkit_set_source_group(${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
