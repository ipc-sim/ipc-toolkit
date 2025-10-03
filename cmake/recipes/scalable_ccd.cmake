# Scalable-CCD (https://github.com/continuous-collision-detection/scalable-ccd)
# License: Apache 2.0
if(TARGET scalable_ccd::scalable_ccd)
    return()
endif()

message(STATUS "Third-party: creating target 'scalable_ccd::scalable_ccd'")

include(CPM)
CPMAddPackage(
    URI "gh:continuous-collision-detection/scalable-ccd#c80af01cab083b3eeb8dac80312ec9cfe479a5cf"
    OPTIONS "SCALABLE_CCD_WITH_CUDA ${IPC_TOOLKIT_WITH_CUDA}"
)

# Folder name for IDE
set_target_properties(scalable_ccd PROPERTIES FOLDER "ThirdParty")