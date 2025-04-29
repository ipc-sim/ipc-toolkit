set(SOURCES
  area_gradient.cpp
  eigen_ext.cpp
  intersection.cpp
  interval.cpp
  logger.cpp
  thread_limiter.cpp
  vertex_to_min_edge.cpp
  world_bbox_diagonal_length.cpp
)

target_sources(ipctk PRIVATE ${SOURCES})