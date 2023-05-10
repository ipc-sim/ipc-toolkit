if(TARGET broadmark::broadmark)
    return()
endif()

message(STATUS "Third-party: creating target 'broadmark'")

include(FetchContent)
FetchContent_Declare(
    broadmark
    GIT_REPOSITORY https://github.com/dbelgrod/Broadmark.git
    GIT_TAG master
    GIT_SHALLOW FALSE
)
FetchContent_GetProperties(broadmark)

if(NOT broadmark_POPULATED)
    FetchContent_Populate(broadmark)
    add_subdirectory(${broadmark_SOURCE_DIR} ${broadmark_BINARY_DIR})
endif()
add_library(broadmark::broadmark ALIAS broadmark)

target_include_directories(broadmark PUBLIC
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies" # Add Bullet3 include directory
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet3/Bullet3Common"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet3"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet3/Bullet3Collision/BroadPhaseCollision"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet2"
)

# ------------------------------------------------------------------------------
# SIMD

# Figure out SSE level support
message(STATUS "Seaching for SSE...")
find_package(SSE)
# Figure out AVX level support
message(STATUS "Searching for AVX...")
find_package(AVX)
# Figure out FMA level support
message(STATUS "Searching for FMA...")
find_package(FMA)
# Add SSE, AVX, and FMA flags to compiler flags
string(REPLACE " " ";" SIMD_FLAGS "${SSE_FLAGS} ${AVX_FLAGS} ${FMA_FLAGS}")
target_compile_options(broadmark PRIVATE ${SIMD_FLAGS})

# ------------------------------------------------------------------------------
# CUDA

include(cuda)
enable_cuda(broadmark)

# ------------------------------------------------------------------------------
# OpenCL

find_package(OpenCL REQUIRED)
include_directories(${OpenCL_INCLUDE_DIRS})
target_link_libraries(broadmark PRIVATE OpenCL::OpenCL)