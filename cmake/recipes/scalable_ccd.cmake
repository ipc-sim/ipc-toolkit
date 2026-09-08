# Scalable-CCD (https://github.com/continuous-collision-detection/scalable-ccd)
# License: Apache 2.0
if(TARGET scalable_ccd::scalable_ccd)
    return()
endif()

message(STATUS "Third-party: creating target 'scalable_ccd::scalable_ccd'")

include(CPM)
CPMAddPackage(
    URI "gh:continuous-collision-detection/scalable-ccd#8f9347c1afc36f2dda17424c15ff5b68087fe8dc"
    OPTIONS "SCALABLE_CCD_WITH_CUDA ${IPC_TOOLKIT_WITH_CUDA}"
)

# Folder name for IDE
set_target_properties(scalable_ccd PROPERTIES FOLDER "ThirdParty")