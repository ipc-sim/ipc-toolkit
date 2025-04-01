set(SOURCES
  smooth_friction_mollifier.hpp
  smooth_friction_mollifier.cpp
  material_friction.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})