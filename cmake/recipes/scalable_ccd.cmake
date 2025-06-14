# Scalable-CCD (https://github.com/continuous-collision-detection/scalable-ccd)
# License: Apache 2.0
if(TARGET scalable_ccd::scalable_ccd)
    return()
endif()

message(STATUS "Third-party: creating target 'scalable_ccd::scalable_ccd'")

set(SCALABLE_CCD_WITH_CUDA ${IPC_TOOLKIT_WITH_CUDA} CACHE BOOL "Enable CUDA CCD" FORCE)

include(CPM)
CPMAddPackage("gh:continuous-collision-detection/scalable-ccd#2c82b9ca43fba30b85f7e9aa83283464b1bb7843")

# Folder name for IDE
set_target_properties(scalable_ccd PROPERTIES FOLDER "ThirdParty")