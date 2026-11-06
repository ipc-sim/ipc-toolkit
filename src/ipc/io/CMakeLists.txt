set(SOURCES
  write_candidates_obj.cpp
  write_candidates_obj.hpp
)

if(IPC_TOOLKIT_WITH_GLTF)
  list(APPEND SOURCES
    read_gltf.cpp
    read_gltf.hpp
    write_gltf.cpp
    write_gltf.hpp
  )
endif()

target_sources(ipc_toolkit PRIVATE ${SOURCES})