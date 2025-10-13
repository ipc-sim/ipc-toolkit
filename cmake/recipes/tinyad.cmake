# TinyAD (https://github.com/patr-schm/tinyad)
# License: MIT
if(TARGET TinyAD::TinyAD)
    return()
endif()

message(STATUS "Third-party: creating target 'TinyAD::TinyAD'")

find_package(Patch REQUIRED)
set(PATCH_COMMAND_ARGS "-rnN")

file(GLOB_RECURSE patches_for_tinyad CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/patches/tinyad_*.patch"
)

set(PATCH_COMMAND_FOR_CPM_BASE "${Patch_EXECUTABLE}" ${PATCH_COMMAND_ARGS} -p1 < )

set(PATCH_COMMAND_FOR_CPM "")
foreach(patch_filename IN LISTS patches_for_tinyad)
    list(APPEND PATCH_COMMAND_FOR_CPM ${PATCH_COMMAND_FOR_CPM_BASE})
    list(APPEND PATCH_COMMAND_FOR_CPM ${patch_filename})
    list(APPEND PATCH_COMMAND_FOR_CPM &&)
endforeach()
list(POP_BACK PATCH_COMMAND_FOR_CPM)

message("Patch command: ${PATCH_COMMAND_FOR_CPM}")

include(CPM)
CPMAddPackage(
        URI "gh:patr-schm/tinyad#4b48d1a1a588874556a692a3abbdecd0db4c23e1"
        PATCH_COMMAND ${PATCH_COMMAND_FOR_CPM})

add_library(TinyAD::TinyAD ALIAS TinyAD)

# Folder name for IDE
set_target_properties(TinyAD PROPERTIES FOLDER "ThirdParty")