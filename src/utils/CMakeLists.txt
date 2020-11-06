set(SOURCES
  eigen_ext.cpp
  eigen_ext.hpp
  local_hessian_to_global_triplets.hpp
  logger.cpp
  logger.hpp
)

prepend_current_path(SOURCES)
ipc_toolkit_copy_headers(${SOURCES})
ipc_toolkit_set_source_group(${SOURCES})
target_sources(${PROJECT_NAME} PRIVATE ${SOURCES})
