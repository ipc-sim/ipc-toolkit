cmake_minimum_required(VERSION 3.8)
project(IPCToolkit)

################################################################################

# project-options
option(IPC_TOOLKIT_BUILD_UNIT_TESTS  "Build unit-tests"      ON)

################################################################################

### Configuration
set(IPC_TOOLKIT_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")

# Include cmake functions
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Color output
include(UseColors)

# Extra warnings
include(Warnings)

# Prepend function
include(PrependCurrentPath)

# IPC Toolkit utils
include(IPCToolkitCopyHeaders)
include(IPCToolkitSetSourceGroup)

# Download external dependencies
include(IPCToolkitDependencies)

################################################################################
# IPC Toolkit Library
################################################################################

# add_library() can only be called without any source since CMake 3.11 ...
add_library(${PROJECT_NAME} src/ipc.cpp)
# NOTE: The actual list of sources will be added by src/CMakeLists.txt
add_subdirectory(src)

# Public include directory for IPC Toolkit
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR}/include)

target_compile_definitions(${PROJECT_NAME} PUBLIC NOMINMAX)

################################################################################
# Optional Definitions
################################################################################

################################################################################
# Required Libraries
################################################################################

# libigl
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core)

# Extra warnings
target_link_libraries(${PROJECT_NAME} PRIVATE warnings::all)

# TBB
target_link_libraries(${PROJECT_NAME} PUBLIC TBB::tbb)

# Etienne Vouga's CTCD Library for the floating point root finding algorithm
target_link_libraries(${PROJECT_NAME} PUBLIC EVCTCD)

################################################################################
# Compiler options
################################################################################

# ## SIMD support
# # Figure out SSE level support
# message(STATUS "Seaching for SSE...")
# find_package(SSE)
# # Figure out AVX level support
# message(STATUS "Searching for AVX...")
# find_package(AVX)
# # Figure out FMA level support
# message(STATUS "Searching for FMA...")
# find_package(FMA)
# # Add SSE, AVX, and FMA flags to compiler flags
# string(REPLACE " " ";" SIMD_FLAGS "${SSE_FLAGS} ${AVX_FLAGS} ${FMA_FLAGS}")
# target_compile_options(${PROJECT_NAME} PUBLIC ${SIMD_FLAGS})

# Use C++14
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)

################################################################################
# Tests
################################################################################

if(IPC_TOOLKIT_BUILD_UNIT_TESTS)
    # Catch2
    if(NOT TARGET Catch2::Catch2)
        ipc_toolkit_download_catch2()
        add_subdirectory(${IPC_TOOLKIT_EXTERNAL}/Catch2 catch2)
        list(APPEND CMAKE_MODULE_PATH ${IPC_TOOLKIT_EXTERNAL}/Catch2/contrib)
    endif()

    include(CTest)
    enable_testing()
    add_subdirectory(tests)
endif()
