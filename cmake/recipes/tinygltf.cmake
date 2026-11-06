# TinyGLTF (https://github.com/syoyo/tinygltf)
# License: MIT
if(TARGET tinygltf::tinygltf)
    return()
endif()

message(STATUS "Third-party: creating target 'tinygltf::tinygltf'")

include(CPM)
CPMAddPackage(
    URI "gh:syoyo/tinygltf@2.9.6"
    DOWNLOAD_ONLY TRUE
)

add_library(tinygltf)
add_library(tinygltf::tinygltf ALIAS tinygltf)
target_sources(tinygltf PRIVATE ${tinygltf_SOURCE_DIR}/tiny_gltf.cc)
target_include_directories(tinygltf INTERFACE
    $<BUILD_INTERFACE:${tinygltf_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)