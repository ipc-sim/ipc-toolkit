if(TARGET gpu_ccd::gpu_ccd)
    return()
endif()

message(STATUS "Third-party: creating target 'gpu_ccd::gpu_ccd'")

if(EXISTS "${IPC_TOOLKIT_GPU_CCD_PATH}")
    message(STATUS "Using GPU CCD found at: ${IPC_TOOLKIT_GPU_CCD_PATH}")
    add_subdirectory("${IPC_TOOLKIT_GPU_CCD_PATH}" "${PROJECT_BINARY_DIR}/gpu_ccd")
else()
    include(FetchContent)
    FetchContent_Declare(
        gpu_ccd
        GIT_REPOSITORY git@github.com:dbelgrod/CCD-GPU.git
        GIT_TAG c0526d70abbe6d3ee1c2c51b8f9e36f15c208187
        GIT_SHALLOW FALSE
    )
    FetchContent_MakeAvailable(gpu_ccd)
endif()

add_library(gpu_ccd::gpu_ccd ALIAS CCDGPU)

set_target_properties(GPUBF PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(CCDGPU PROPERTIES POSITION_INDEPENDENT_CODE ON)