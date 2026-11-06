set(SOURCES)

if(IPC_TOOLKIT_WITH_GLTF)
  list(APPEND SOURCES
    test_read_gltf.cpp
    test_write_gltf.cpp
  )
endif()

target_sources(ipc_toolkit_tests PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################