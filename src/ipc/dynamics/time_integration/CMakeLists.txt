set(SOURCES
  bdf.cpp
  bdf.hpp
  time_integrator.cpp
  time_integrator.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
