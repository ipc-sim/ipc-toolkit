# libigl (https://github.com/libigl/libigl)
# License: MPL-2.0
if(TARGET igl::core)
    return()
endif()

message(STATUS "Third-party: creating target 'igl::core'")

include(eigen)

include(CPM)
CPMAddPackage(
    URI "gh:libigl/libigl@2.6.0"
    OPTIONS "LIBIGL_PREDICATES ON"
)

# Folder name for IDE
foreach(target_name IN ITEMS core predicates)
    set_target_properties(igl_${target_name} PROPERTIES FOLDER "ThirdParty/libigl")
endforeach()
