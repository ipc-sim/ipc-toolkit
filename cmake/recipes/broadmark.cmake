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

include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
else()
    message(FATAL_ERROR "Broadmark: No CUDA support found!")
endif()

# We need to explicitly state that we need all CUDA files in the particle
# library to be built with -dc as the member functions could be called by
# other libraries and executables.
set_target_properties(broadmark PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.24.0")
    set(CMAKE_CUDA_ARCHITECTURES "native")
    set_target_properties(broadmark PROPERTIES CUDA_ARCHITECTURES "native")
else()
    include(FindCUDA/select_compute_arch)
    CUDA_DETECT_INSTALLED_GPUS(CUDA_ARCH_LIST)
    string(STRIP "${CUDA_ARCH_LIST}" CUDA_ARCH_LIST)
    string(REPLACE " " ";" CUDA_ARCH_LIST "${CUDA_ARCH_LIST}")
    string(REPLACE "." "" CUDA_ARCH_LIST "${CUDA_ARCH_LIST}")
    set(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_LIST})
    message(STATUS "Compiling for CUDA architectures: ${CUDA_ARCH_LIST}")
    set_target_properties(broadmark PROPERTIES CUDA_ARCHITECTURES "${CUDA_ARCH_LIST}")
endif()

if(APPLE)
    # We need to add the path to the driver (libcuda.dylib) as an rpath,
    # so that the static cuda runtime can find it at runtime.
    set_property(TARGET broadmark
                 PROPERTY
                 BUILD_RPATH ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
endif()

find_package(CUDAToolkit)
target_link_libraries(broadmark PRIVATE CUDA::cudart)

# ------------------------------------------------------------------------------
# OpenCL

find_package(OpenCL REQUIRED)
include_directories(${OpenCL_INCLUDE_DIRS})
target_link_libraries(broadmark PRIVATE OpenCL::OpenCL)