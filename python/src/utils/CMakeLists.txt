set(SOURCES
  eigen_ext.cpp
  logger.cpp
  thread_limiter.cpp
  vertex_to_min_edge.cpp
  world_bbox_diagonal_length.cpp
)

target_sources(ipctk PRIVATE ${SOURCES})