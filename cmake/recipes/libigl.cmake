# libigl (https://github.com/libigl/libigl)
# License: MPL-2.0
if(TARGET igl::core)
    return()
endif()

message(STATUS "Third-party: creating target 'igl::core'")

set(LIBIGL_PREDICATES ON CACHE BOOL "Use exact predicates" FORCE)
if(IPC_TOOLKIT_WITH_CGAL)
    set(LIBIGL_COPYLEFT_CGAL ON CACHE BOOL "Use CGAL" FORCE)
    set(LIBIGL_DEFAULT_CGAL ON CACHE BOOL "Use CGAL" FORCE)
endif()

# set(LIBIGL_COPYLEFT_CGAL ON CACHE BOOL "Use CGAL" FORCE)
include(eigen)

include(CPM)
CPMAddPackage("gh:libigl/libigl@2.5.0")