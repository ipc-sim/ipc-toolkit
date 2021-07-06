cmake_minimum_required(VERSION 3.8)
project(IPCToolkit
        DESCRIPTION "A set of reusable functions to integrate IPC into an existing simulation."
        LANGUAGES CXX)

################################################################################

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(IPC_TOOLKIT_TOPLEVEL_PROJECT ON)
else()
  set(IPC_TOOLKIT_TOPLEVEL_PROJECT OFF)
endif()

# project-options
option(IPC_TOOLKIT_BUILD_UNIT_TESTS           "Build unit-tests"  ${IPC_TOOLKIT_TOPLEVEL_PROJECT})
option(IPC_TOOLKIT_WITH_LOGGER                "Enable logger"                                  ON)
option(IPC_TOOLKIT_WITH_CORRECT_CCD           "Use the TightInclusion CCD"                     ON)
option(IPC_TOOLKIT_WITH_PYTHON                "Build Python bindings"                         OFF)
option(IPC_TOOLKIT_WITH_SIMD                  "Enable SIMD"                                   OFF)
option(IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION "Use rational edge-triangle intersection check" OFF)

################################################################################

### Configuration
set(IPC_TOOLKIT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
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
# Libraries
################################################################################

# libigl
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core igl::predicates)

# Extra warnings
target_link_libraries(${PROJECT_NAME} PRIVATE warnings::all)

# TBB
target_link_libraries(${PROJECT_NAME} PUBLIC TBB::tbb)

# CCD
if(IPC_TOOLKIT_WITH_CORRECT_CCD)
  # Provably conservative CCD of [Wang and Ferguson et al. 2020]
  target_link_libraries(${PROJECT_NAME} PUBLIC TightInclusion)
  target_compile_definitions(${PROJECT_NAME} PUBLIC IPC_TOOLKIT_WITH_CORRECT_CCD)
else()
  # Etienne Vouga's CTCD Library for the floating point root finding algorithm
  target_link_libraries(${PROJECT_NAME} PUBLIC EVCTCD)
endif()

# Logger
if(IPC_TOOLKIT_WITH_LOGGER)
  target_link_libraries(${PROJECT_NAME} PRIVATE spdlog::spdlog)
  target_compile_definitions(${PROJECT_NAME} PUBLIC IPC_TOOLKIT_WITH_LOGGER)
endif()

# Faster unoredered map
target_link_libraries(${PROJECT_NAME} PUBLIC tsl::robin_map)

if(IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION)
  find_package(GMP)
  IF(NOT ${GMP_FOUND})
    MESSAGE(FATAL_ERROR "Unable to find GMP necessary for rational intersection")
  ENDIF()
  target_include_directories(${PROJECT_NAME} PUBLIC ${GMP_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} PUBLIC ${GMP_LIBRARIES})
  target_compile_definitions(${PROJECT_NAME} PUBLIC IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION)
endif()

################################################################################
# Compiler options
################################################################################

## SIMD support
if(IPC_TOOLKIT_WITH_SIMD)
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
  target_compile_options(${PROJECT_NAME} PRIVATE ${SIMD_FLAGS})
endif()

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

  # finite-diff
  if(NOT TARGET FiniteDiff::FiniteDiff)
    ipc_toolkit_download_finite_diff()
    add_subdirectory(${IPC_TOOLKIT_EXTERNAL}/finite-diff EXCLUDE_FROM_ALL)
    add_library(FiniteDiff::FiniteDiff ALIAS FiniteDiff)
  endif()

  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()

################################################################################
# Python bindings
################################################################################

if(IPC_TOOLKIT_WITH_PYTHON)
  add_subdirectory(python)
endif()
