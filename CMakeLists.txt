# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
    set(IPC_TOOLKIT_TOPLEVEL_PROJECT OFF)
else()
    set(IPC_TOOLKIT_TOPLEVEL_PROJECT ON)
endif()

# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.14.0")
if(IPC_TOOLKIT_TOPLEVEL_PROJECT)
    cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
else()
    # Don't use cmake_minimum_required here to avoid implicitly overriding parent policies
    if(${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
        message(FATAL_ERROR "CMake required version to build IPC Toolkit is ${REQUIRED_CMAKE_VERSION}")
    endif()
endif()

# Include user-provided default options if available. We do that before the main
# `project()` so that we can define the C/C++ compilers from the option file.
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/IPCToolkitOptions.cmake)
    message(STATUS "Using local options file: ${CMAKE_CURRENT_SOURCE_DIR}/IPCToolkitOptions.cmake")
    include(${CMAKE_CURRENT_SOURCE_DIR}/IPCToolkitOptions.cmake)
endif()

################################################################################

project(IPCToolkit
        DESCRIPTION "A set of reusable functions to integrate IPC into an existing simulation."
        LANGUAGES CXX
        VERSION "0.0.1")

option(IPC_TOOLKIT_BUILD_TESTS                "Build unit-tests"  ${IPC_TOOLKIT_TOPLEVEL_PROJECT})
option(IPC_TOOLKIT_BUILD_PYTHON               "Build Python bindings"                         OFF)
option(IPC_TOOLKIT_WITH_CORRECT_CCD           "Use the TightInclusion CCD"                     ON)
option(IPC_TOOLKIT_WITH_SIMD                  "Enable SIMD"                                   OFF)
option(IPC_TOOLKIT_WITH_CUDA                  "Enable CUDA CCD"                               OFF)
option(IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION "Use rational edge-triangle intersection check" OFF)

# Set default minimum C++ standard
if(IPC_TOOLKIT_TOPLEVEL_PROJECT)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

### Configuration
set(IPC_TOOLKIT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src/ipc")
set(IPC_TOOLKIT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/ipc_toolkit/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/find/")

# Color output
include(ipc_toolkit_use_colors)

# IPC Toolkit utils
include(ipc_toolkit_prepend_current_path)
include(ipc_toolkit_set_source_group)
include(ipc_toolkit_target_link_system_libraries)

# Generate position independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################
# IPC Toolkit Library
################################################################################

# Add an empty library and fill in the list of sources in `src/CMakeLists.txt`.
add_library(ipc_toolkit)
add_library(ipc::toolkit ALIAS ipc_toolkit)

# Fill in configuration options
configure_file(
  "${IPC_TOOLKIT_SOURCE_DIR}/config.hpp.in"
  "${IPC_TOOLKIT_SOURCE_DIR}/config.hpp")

# Add source and header files to ipc_toolkit
add_subdirectory("${IPC_TOOLKIT_SOURCE_DIR}")

# Public include directory for IPC Toolkit
target_include_directories(ipc_toolkit PUBLIC "${IPC_TOOLKIT_INCLUDE_DIR}")

################################################################################
# Optional Definitions
################################################################################

# For MSVC, do not use the min and max macros.
target_compile_definitions(ipc_toolkit PUBLIC NOMINMAX)

################################################################################
# Dependencies
################################################################################

# libigl
include(eigen)
include(libigl)
ipc_toolkit_target_link_system_libraries(ipc_toolkit PUBLIC
  Eigen3::Eigen
  igl::core
  igl::predicates
)

# TBB
include(onetbb)
ipc_toolkit_target_link_system_libraries(ipc_toolkit PUBLIC TBB::tbb)

# CCD
if(IPC_TOOLKIT_WITH_CORRECT_CCD)
  # Provably conservative CCD of [Wang and Ferguson et al. 2021]
  include(tight_inclusion)
  ipc_toolkit_target_link_system_libraries(ipc_toolkit PUBLIC tight_inclusion::tight_inclusion)
else()
  # Etienne Vouga's CTCD Library for the floating point root finding algorithm
  include(evouga_ccd)
  ipc_toolkit_target_link_system_libraries(ipc_toolkit PUBLIC evouga::ccd)
endif()

# Logger
include(spdlog)
ipc_toolkit_target_link_system_libraries(ipc_toolkit PUBLIC spdlog::spdlog)

if(IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION)
  include(gmp)
  ipc_toolkit_target_link_system_libraries(ipc_toolkit PUBLIC GMP::GMP)
endif()

# Sweep and Tiniest Queue and CCD
if(IPC_TOOLKIT_WITH_CUDA)
  include(gpu_ccd)
  ipc_toolkit_target_link_system_libraries(ipc_toolkit PUBLIC gpu_ccd::gpu_ccd)
  ipc_toolkit_target_link_system_libraries(ipc_toolkit PUBLIC STQ::CPU)
else()
  set(STQ_WITH_CUDA OFF CACHE BOOL "Enable CUDA Implementation" FORCE)
  include(sweep_and_tiniest_queue)
  ipc_toolkit_target_link_system_libraries(ipc_toolkit PUBLIC STQ::CPU)
endif()

# Faster unordered map
include(robin_map)
ipc_toolkit_target_link_system_libraries(ipc_toolkit PUBLIC tsl::robin_map)

# Hashes
include(abseil)
ipc_toolkit_target_link_system_libraries(ipc_toolkit PUBLIC absl::hash)

# Extra warnings (link last for highest priority)
include(ipc_toolkit_warnings)
target_link_libraries(ipc_toolkit PRIVATE ipc::toolkit::warnings)

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
  target_compile_options(ipc_toolkit PRIVATE ${SIMD_FLAGS})
endif()

# Use C++17
target_compile_features(ipc_toolkit PUBLIC cxx_std_17)

################################################################################
# CUDA
################################################################################

if(IPC_TOOLKIT_WITH_CUDA)
  include(CheckLanguage)
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
  else()
    message(FATAL_ERROR "No CUDA support found!")
  endif()

  # We need to explicitly state that we need all CUDA files in the particle
  # library to be built with -dc as the member functions could be called by
  # other libraries and executables.
  set_target_properties(ipc_toolkit PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

  # Nvidia RTX8000 -> compute_75
  # Nvidia V100 -> compute_70
  # Nvidia 1080/1080Ti -> compute_61
  # Nvidia 3080Ti -> compute_86
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 70 75 86)
  endif()
  set_target_properties(ipc_toolkit PROPERTIES CUDA_ARCHITECTURES "70;75;86")

  if(APPLE)
    # We need to add the path to the driver (libcuda.dylib) as an rpath,
    # so that the static cuda runtime can find it at runtime.
    set_property(TARGET ipc_toolkit
                 PROPERTY
                 BUILD_RPATH ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
  endif()

  find_package(CUDAToolkit)
  ipc_toolkit_target_link_system_libraries(ipc_toolkit PRIVATE CUDA::cudart)
endif()

################################################################################
# Tests
################################################################################

# Enable unit testing at the root level
if(IPC_TOOLKIT_TOPLEVEL_PROJECT AND IPC_TOOLKIT_BUILD_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()

################################################################################
# Python bindings
################################################################################

if(IPC_TOOLKIT_TOPLEVEL_PROJECT AND IPC_TOOLKIT_BUILD_PYTHON)
  add_subdirectory(python)
endif()
