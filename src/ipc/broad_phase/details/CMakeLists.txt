set(SOURCES
  connectivity_filters.hpp
  lbvh_build.hpp
  lbvh_traverse.hpp
  spatial_hash_impl.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})
