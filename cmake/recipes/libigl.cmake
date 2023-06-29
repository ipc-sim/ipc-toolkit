if(TARGET igl::core)
    return()
endif()

message(STATUS "Third-party: creating target 'igl::core'")

set(LIBIGL_PREDICATES ON CACHE BOOL "Use exact predicates" FORCE)

include(eigen)

include(CPM)
CPMAddPackage("gh:libigl/libigl@2.4.0")
