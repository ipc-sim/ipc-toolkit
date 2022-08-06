if(TARGET gpu_ccd::gpu_ccd)
    return()
endif()

message(STATUS "Third-party: creating target 'gpu_ccd::gpu_ccd'")

option(STQ_WITH_CPU  "Enable CPU Implementation"  ON)
set(STQ_WITH_CUDA ON CACHE BOOL "Enable CUDA Implementation" FORCE)

include(FetchContent)
FetchContent_Declare(
    gpu_ccd
    GIT_REPOSITORY https://github.com/dbelgrod/CCD-GPU.git
    GIT_TAG 2fbf59063918e644452dcfa3c65e72c428bfa9b5
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(gpu_ccd)

add_library(gpu_ccd::gpu_ccd ALIAS CCDGPU)

set_target_properties(STQ_CPU PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(CCDGPU PROPERTIES POSITION_INDEPENDENT_CODE ON)