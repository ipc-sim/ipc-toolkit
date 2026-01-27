set(SOURCES
  feasible_region.cpp
  feasible_region.hpp
  trust_region.cpp
  trust_region.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################