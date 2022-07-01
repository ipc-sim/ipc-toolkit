if(TARGET gpu_ccd::gpu_ccd)
    return()
endif()

message(STATUS "Third-party: creating target 'gpu_ccd::gpu_ccd'")

option(STQ_WITH_CPU  "Enable CPU Implementation"  ON)
set(STQ_WITH_CUDA ON CACHE BOOL "Enable CUDA Implementation" FORCE)

if(EXISTS "${IPC_TOOLKIT_GPU_CCD_PATH}")
    message(STATUS "Using GPU CCD found at: ${IPC_TOOLKIT_GPU_CCD_PATH}")
    add_subdirectory("${IPC_TOOLKIT_GPU_CCD_PATH}" "${PROJECT_BINARY_DIR}/gpu_ccd")
else()
    include(FetchContent)
    FetchContent_Declare(
        gpu_ccd
        GIT_REPOSITORY https://github.com/dbelgrod/CCD-GPU.git
        GIT_TAG deaf52698f1e8bb0c56d4e04619bb9414a0be746
        GIT_SHALLOW FALSE
    )
    FetchContent_MakeAvailable(gpu_ccd)
endif()

add_library(gpu_ccd::gpu_ccd ALIAS CCDGPU)

set_target_properties(STQ_CPU PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(CCDGPU PROPERTIES POSITION_INDEPENDENT_CODE ON)