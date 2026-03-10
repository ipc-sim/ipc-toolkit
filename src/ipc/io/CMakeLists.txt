set(SOURCES
  read_gltf.cpp
  read_gltf.hpp
  write_candidates_obj.cpp
  write_candidates_obj.hpp
  write_gltf.cpp
  write_gltf.hpp
)

target_sources(ipc_toolkit PRIVATE ${SOURCES})