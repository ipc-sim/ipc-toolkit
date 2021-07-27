set(SOURCES
  aabb.cpp
  aabb.hpp
  ccd.cpp
  ccd.hpp
)

ipc_toolkit_prepend_current_path(SOURCES)
ipc_toolkit_copy_headers(${SOURCES})
ipc_toolkit_set_source_group(${SOURCES})
target_sources(${PROJECT_NAME} PRIVATE ${SOURCES})
