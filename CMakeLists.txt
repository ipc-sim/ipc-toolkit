cmake_minimum_required(VERSION 3.8)
project(IPCToolkit)

################################################################################

# project-options
option(IPC_TOOLKIT_BUILD_UNIT_TESTS  "Build unit-tests"      ON)

################################################################################

### Configuration
set(IPC_TOOLKIT_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(Warnings)
include(UseColors)
include(${PROJECT_NAME}Dependencies)

################################################################################
# IPC Toolkit Library
################################################################################

add_library(${PROJECT_NAME}
  # src/autodiff/autodiff.cpp

  # src/ccd/impact.cpp
  # src/ccd/collision_candidate.cpp
  # src/ccd/hash_grid.cpp
  # src/ccd/collision_detection_narrow_phase.cpp
  # src/ccd/collision_detection_broad_phase.cpp
  #
  # src/geometry/intersection.cpp

  src/ipc.cpp

  src/barrier/barrier.cpp

  src/distance/line_line.cpp
  src/distance/point_line.cpp
  src/distance/point_plane.cpp
  src/distance/edge_edge_mollifier.cpp

  src/spatial_hash/collision_candidate.cpp
  src/spatial_hash/hash_grid.cpp

  src/utils/eigen_ext.cpp
)
target_compile_definitions(${PROJECT_NAME} PUBLIC NOMINMAX)

target_include_directories(${PROJECT_NAME} PUBLIC src)

################################################################################
# Optional Definitions
################################################################################

################################################################################
# Required Libraries
################################################################################

# libigl
# TODO: Replace this with Eigen::Eigen3
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core)

# Extra warnings
target_link_libraries(${PROJECT_NAME} PRIVATE warnings::all)

# TBB
target_link_libraries(${PROJECT_NAME} PUBLIC TBB::tbb)

# Etienne Vouga's CTCD Library for the floating point root finding algorithm
target_link_libraries(${PROJECT_NAME} PUBLIC EVCTCD)


# Teseo's Exact CCD library for exact segment-triangle intersection check
if(TARGET RationalCCD)
  target_link_libraries(${PROJECT_NAME} PUBLIC RationalCCD)
endif()

################################################################################
# Compiler options
################################################################################

## SIMD support
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
target_compile_options(${PROJECT_NAME} PUBLIC ${SIMD_FLAGS})

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
