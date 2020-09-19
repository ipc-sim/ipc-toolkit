set(SOURCES
  adaptive_stiffness.cpp
  adaptive_stiffness.hpp
  barrier.cpp
  barrier.hpp
  barrier.tpp
)

prepend_current_path(SOURCES)
ipc_toolkit_copy_headers(${SOURCES})
ipc_toolkit_set_source_group(${SOURCES})
target_sources(${PROJECT_NAME} PRIVATE ${SOURCES})
