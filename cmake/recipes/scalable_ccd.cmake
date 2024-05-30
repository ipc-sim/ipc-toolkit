# Scalable-CCD (https://github.com/continuous-collision-detection/scalable-ccd)
# License: Apache 2.0
if(TARGET scalable_ccd::scalable_ccd)
    return()
endif()

message(STATUS "Third-party: creating target 'scalable_ccd::scalable_ccd'")

set(SCALABLE_CCD_WITH_CUDA ${IPC_TOOLKIT_WITH_CUDA} CACHE BOOL "Enable CUDA CCD" FORCE)

include(CPM)
CPMAddPackage("gh:continuous-collision-detection/scalable-ccd#60692d7cbc03390a3a198792a3ff2b869c5da05b")