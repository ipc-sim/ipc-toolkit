set(SOURCES
  closest_point.cpp
  relative_velocity.cpp
  tangent_basis.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipctk PRIVATE ${SOURCES})