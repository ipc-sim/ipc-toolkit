# BVH (https://github.com/geometryprocessing/SimpleBVH)
# License: MIT

if(TARGET simple_bvh::simple_bvh)
    return()
endif()

message(STATUS "Third-party: creating target 'simple_bvh::simple_bvh'")

include(CPM)
CPMAddPackage("gh:geometryprocessing/SimpleBVH#e1a931337a9e07e8bd2d2e8bbdfd7e54bc850df5")

# Folder name for IDE
set_target_properties(simple_bvh PROPERTIES FOLDER "ThirdParty")