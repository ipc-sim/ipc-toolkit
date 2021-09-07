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
        LANGUAGES CXX)

option(IPC_TOOLKIT_BUILD_UNIT_TESTS           "Build unit-tests"  ${IPC_TOOLKIT_TOPLEVEL_PROJECT})
option(IPC_TOOLKIT_WITH_LOGGER                "Enable logger"                                  ON)
option(IPC_TOOLKIT_WITH_CORRECT_CCD           "Use the TightInclusion CCD"                     ON)
option(IPC_TOOLKIT_WITH_PYTHON                "Build Python bindings"                         OFF)
option(IPC_TOOLKIT_WITH_SIMD                  "Enable SIMD"                                   OFF)
option(IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION "Use rational edge-triangle intersection check" OFF)

# Set default minimum C++ standard
if(IPC_TOOLKIT_TOPLEVEL_PROJECT)
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

### Configuration
set(IPC_TOOLKIT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
# set(IPC_TOOLKIT_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/ipc_toolkit/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/find/")

# Color output
include(ipc_toolkit_use_colors)

# IPC Toolkit utils
include(ipc_toolkit_prepend_current_path)
include(ipc_toolkit_copy_headers)
include(ipc_toolkit_set_source_group)

################################################################################
# IPC Toolkit Library
################################################################################

# Add an empty library and fill in the list of sources in `src/CMakeLists.txt`.
add_library(ipc_toolkit)
add_library(ipc::toolkit ALIAS ipc_toolkit)

add_subdirectory(src)

# Public include directory for IPC Toolkit
target_include_directories(ipc_toolkit PUBLIC ${PROJECT_BINARY_DIR}/include)

################################################################################
# Optional Definitions
################################################################################

# For MSVC, do not use the min and max macros.
target_compile_definitions(ipc_toolkit PUBLIC NOMINMAX)

################################################################################
# Dependencies
################################################################################

# Extra warnings
include(ipc_toolkit_warnings)
target_link_libraries(ipc_toolkit PRIVATE ipc::toolkit::warnings)

# libigl
include(eigen)
include(libigl)
target_link_libraries(ipc_toolkit PUBLIC
  Eigen3::Eigen
  igl::core
  igl::predicates
)

# TBB
include(tbb)
target_link_libraries(ipc_toolkit PUBLIC TBB::tbb)

# CCD
if(IPC_TOOLKIT_WITH_CORRECT_CCD)
  # Provably conservative CCD of [Wang and Ferguson et al. 2021]
  include(tight_inclusion)
  target_link_libraries(ipc_toolkit PUBLIC TightInclusion)
  target_compile_definitions(ipc_toolkit PUBLIC IPC_TOOLKIT_WITH_CORRECT_CCD)
else()
  # Etienne Vouga's CTCD Library for the floating point root finding algorithm
  include(evouga_ccd)
  target_link_libraries(ipc_toolkit PUBLIC evouga::ccd)
endif()

# Logger
if(IPC_TOOLKIT_WITH_LOGGER)
  include(spdlog)
  target_link_libraries(ipc_toolkit PRIVATE spdlog::spdlog)
  target_compile_definitions(ipc_toolkit PUBLIC IPC_TOOLKIT_WITH_LOGGER)
else()
  include(fmt)
  target_link_libraries(ipc_toolkit PRIVATE fmt::fmt)
endif()

# Faster unoredered map
include(robin_map)
target_link_libraries(ipc_toolkit PUBLIC tsl::robin_map)

if(IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION)
  find_package(GMP REQUIRED)
  IF(NOT ${GMP_FOUND})
    MESSAGE(FATAL_ERROR "Unable to find GMP necessary for rational intersection")
  ENDIF()
  target_include_directories(ipc_toolkit PUBLIC ${GMP_INCLUDE_DIRS})
  target_link_libraries(ipc_toolkit PUBLIC ${GMP_LIBRARIES})
  target_compile_definitions(ipc_toolkit PUBLIC IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION)
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
  target_compile_options(ipc_toolkit PRIVATE ${SIMD_FLAGS})
endif()

# Use C++14
target_compile_features(ipc_toolkit PUBLIC cxx_std_14)

################################################################################
# Tests
################################################################################

# Enable unit testing at the root level
if(IPC_TOOLKIT_BUILD_UNIT_TESTS)
  include(CTest)
  enable_testing()

  # Include Catch2 and provide function `catch_discover_tests` to register tests.
  include(catch2)
  FetchContent_GetProperties(catch2)
  include("${catch2_SOURCE_DIR}/contrib/Catch.cmake")

  add_subdirectory(tests)
endif()

################################################################################
# Python bindings
################################################################################

if(IPC_TOOLKIT_WITH_PYTHON)
  add_subdirectory(python)
endif()
