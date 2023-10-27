# BVH
# License: MIT

if(TARGET simple_bvh::simple_bvh)
    return()
endif()

message(STATUS "Third-party: creating target 'simple_bvh::simple_bvh'")

include(CPM)
CPMAddPackage("gh:ipc-sim/SimpleBVH#2117898eb366647d6aacdb82860b9315fb42d6ad")