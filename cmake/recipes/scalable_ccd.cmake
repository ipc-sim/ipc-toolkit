# Abseil (https://github.com/continuous-collision-detection/scalable-ccd)
# License: Apache 2.0
if(TARGET scalable_ccd::scalable_ccd)
    return()
endif()

message(STATUS "Third-party: creating target 'scalable_ccd::scalable_ccd'")

set(SCALABLE_CCD_WITH_CUDA ${IPC_TOOLKIT_WITH_CUDA} CACHE BOOL "Enable CUDA CCD" FORCE)
set(SCALABLE_CCD_KEEP_CPU_OVERLAPS ON CACHE BOOL "Keep CPU overlaps after CUDA method" FORCE)

include(CPM)
CPMAddPackage("gh:continuous-collision-detection/scalable-ccd#bc21edc28c1fc2d6321d2cc1e4b7b521c376285e")