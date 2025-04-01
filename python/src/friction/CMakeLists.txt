set(SOURCES
  bindings.hpp
  smooth_friction_mollifier.cpp
  material_friction.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit_python PRIVATE ${SOURCES})