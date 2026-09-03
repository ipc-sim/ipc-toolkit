set(SOURCES
  adaptive_support.cpp
  adaptive_support.hpp
  arbitrary_point_bvh.cpp
  arbitrary_point_bvh.hpp
  arbitrary_point_potential.cpp
  arbitrary_point_potential.hpp
  esp_collisions.cpp
  esp_collisions.hpp
  esp_collisions_builder.cpp
  esp_collisions_builder.hpp
  esp_potential.hpp
  esp_potential.cpp
  quadrature_potential.cpp
  quadrature_potential.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(collisions)
