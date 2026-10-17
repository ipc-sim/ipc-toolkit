set(SOURCES
  arbitrary_point_bvh.cpp
  arbitrary_point_bvh.hpp
  arbitrary_point_esp.cpp
  arbitrary_point_esp.hpp
  esp_collision_maps.hpp
  esp_collisions.cpp
  esp_collisions.hpp
  esp_collisions_builder.cpp
  esp_collisions_builder.hpp
  esp_distance.hpp
  esp_parameters.hpp
  esp_potential.cpp
  esp_potential.hpp
  esp_tangential_collisions.cpp
  quadrature_potential.cpp
  quadrature_potential.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(collisions)
