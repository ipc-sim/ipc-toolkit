if(TARGET broadmark::broadmark)
    return()
endif()

message(STATUS "Third-party: creating target 'broadmark'")

include(CPM)
CPMAddPackage("gh:dbelgrod/Broadmark#e380a7fed66369164c9c22166bcb97f11e3e601b")

add_library(broadmark::broadmark ALIAS broadmark)

target_include_directories(broadmark PUBLIC
    "${Broadmark_SOURCE_DIR}/Broadmark/Algorithms"
    "${Broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies" # Add Bullet3 include directory
    "${Broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet3/Bullet3Common"
    "${Broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet3"
    "${Broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet3/Bullet3Collision/BroadPhaseCollision"
    "${Broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet2"
)

# ------------------------------------------------------------------------------
# SIMD

# Figure out AVX level support
message(STATUS "Searching for AVX...")
find_package(AVX)
# Add AVX flags to compiler flags
string(REPLACE " " ";" AVX_FLAGS "${AVX_FLAGS}")
target_compile_options(broadmark PRIVATE ${AVX_FLAGS})

# ------------------------------------------------------------------------------
# CUDA

include(cuda)
enable_cuda(broadmark)

# ------------------------------------------------------------------------------
# OpenCL

find_package(OpenCL REQUIRED)
include_directories(${OpenCL_INCLUDE_DIRS})
target_link_libraries(broadmark PRIVATE OpenCL::OpenCL)