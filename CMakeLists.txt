# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
  set(IPC_TOOLKIT_TOPLEVEL_PROJECT OFF)
else()
  set(IPC_TOOLKIT_TOPLEVEL_PROJECT ON)
endif()

# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.18.0")
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

# Enable ccache if available
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  option(IPC_TOOLKIT_WITH_CCACHE "Enable ccache when building IPC Toolkit" ${IPC_TOOLKIT_TOPLEVEL_PROJECT})
else()
  option(IPC_TOOLKIT_WITH_CCACHE "Enable ccache when building IPC Toolkit" OFF)
endif()
if(IPC_TOOLKIT_WITH_CCACHE AND CCACHE_PROGRAM)
  message(STATUS "Enabling Ccache support (${CCACHE_PROGRAM})")
  set(ccacheEnv
    CCACHE_BASEDIR=${CMAKE_BINARY_DIR}
    CCACHE_SLOPPINESS=clang_index_store,include_file_ctime,include_file_mtime,locale,pch_defines,time_macros
  )
  foreach(lang IN ITEMS C CXX)
    set(CMAKE_${lang}_COMPILER_LAUNCHER
      ${CMAKE_COMMAND} -E env ${ccacheEnv} ${CCACHE_PROGRAM}
    )
  endforeach()
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

################################################################################
# CMake Policies
################################################################################

cmake_policy(SET CMP0054 NEW) # Only interpret if() arguments as variables or keywords when unquoted.
cmake_policy(SET CMP0076 NEW) # target_sources() command converts relative paths to absolute.
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24")
  cmake_policy(SET CMP0135 NEW) # Set the timestamps of all extracted contents to the time of the extraction.
endif()
cmake_policy(SET CMP0114 NEW) # Support the Xcode "new build system"

################################################################################

project(IPCToolkit
        DESCRIPTION "A set of reusable functions to integrate IPC into an existing simulation."
        LANGUAGES CXX
        VERSION "1.5.0")

include(CheckLanguage)
check_language(CUDA)

if(IPC_TOOLKIT_TOPLEVEL_PROJECT)
  option(IPC_TOOLKIT_BUILD_TESTS  "Build unit-tests"       ON)
  option(IPC_TOOLKIT_BUILD_PYTHON "Build Python bindings" OFF)
else()
  # If this is not the top-level project, we don't want to build tests or Python
  # bindings. This is useful for projects that use IPC Toolkit as a submodule.
  set(IPC_TOOLKIT_BUILD_TESTS  OFF CACHE BOOL "Build unit-tests"      FORCE)
  set(IPC_TOOLKIT_BUILD_PYTHON OFF CACHE BOOL "Build Python bindings" FORCE)
endif()

if(CMAKE_CUDA_COMPILER)
  option(IPC_TOOLKIT_WITH_CUDA "Enable CUDA CCD" ON)
else()
  option(IPC_TOOLKIT_WITH_CUDA "Enable CUDA CCD" OFF)
endif()
option(IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION "Use rational edge-triangle intersection check" OFF)
option(IPC_TOOLKIT_WITH_ROBIN_MAP             "Use Tessil's robin-map rather than std maps"    ON)
option(IPC_TOOLKIT_WITH_ABSEIL                "Use Abseil's hash functions"                    ON)
option(IPC_TOOLKIT_WITH_FILIB                 "Use filib for interval arithmetic"              ON)
option(IPC_TOOLKIT_WITH_INEXACT_CCD           "Use the original inexact CCD method of IPC"    OFF)

# Advanced options
option(IPC_TOOLKIT_WITH_SIMD                  "Enable SIMD"                                   OFF)
option(IPC_TOOLKIT_WITH_CODE_COVERAGE         "Enable coverage reporting"                     OFF)

mark_as_advanced(IPC_TOOLKIT_WITH_SIMD)          # This does not work reliably
mark_as_advanced(IPC_TOOLKIT_WITH_CODE_COVERAGE) # This is used in GitHub Actions

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

# General CMake utils
include(ipc_toolkit_cpm_cache)
include(ipc_toolkit_use_colors)

# Generate position-independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################
# IPC Toolkit Library
################################################################################

# Add an empty library and fill in the list of sources in `src/ipc/CMakeLists.txt`.
add_library(ipc_toolkit)
add_library(ipc::toolkit ALIAS ipc_toolkit)

# Fill in configuration options
configure_file(
  "${IPC_TOOLKIT_SOURCE_DIR}/config.hpp.in"
  "${PROJECT_BINARY_DIR}/include/ipc/config.hpp")

# Add source and header files to ipc_toolkit
add_subdirectory("${IPC_TOOLKIT_SOURCE_DIR}")

# Public include directory for IPC Toolkit
target_include_directories(ipc_toolkit PUBLIC
  "${IPC_TOOLKIT_INCLUDE_DIR}"    # public headers
  "${PROJECT_BINARY_DIR}/include" # generated config.hpp
)

# Folder name for IDE
set_target_properties(ipc_toolkit PROPERTIES FOLDER "SRC")
get_target_property(IPC_TOOLKIT_SOURCES ipc_toolkit SOURCES)
source_group(TREE "${PROJECT_SOURCE_DIR}" FILES ${IPC_TOOLKIT_SOURCES})

################################################################################
# Dependencies
################################################################################

# Eigen
include(eigen)
target_link_libraries(ipc_toolkit PUBLIC Eigen3::Eigen)

# libigl
include(libigl)
target_link_libraries(ipc_toolkit PRIVATE igl::core igl::predicates)

# TBB
include(onetbb)
target_link_libraries(ipc_toolkit PRIVATE TBB::tbb)

# Provably conservative CCD of [Wang and Ferguson et al. 2021]
include(tight_inclusion)
target_link_libraries(ipc_toolkit PRIVATE tight_inclusion::tight_inclusion)

# Scalable CCD (STQ broad phase and GPU Tight Inclusion)
include(scalable_ccd)
target_link_libraries(ipc_toolkit PRIVATE scalable_ccd::scalable_ccd)

# CCD
if(IPC_TOOLKIT_WITH_INEXACT_CCD)
  # Etienne Vouga's CTCD Library for the floating point root finding algorithm
  include(evouga_ccd)
  target_link_libraries(ipc_toolkit PRIVATE evouga::ccd)
endif()

# SimpleBVH
include(simple_bvh)
target_link_libraries(ipc_toolkit PRIVATE simple_bvh::simple_bvh)

# Logger
include(spdlog)
target_link_libraries(ipc_toolkit PUBLIC spdlog::spdlog)

# rational-cpp (requires GMP)
if(IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION)
  include(rational_cpp)
  target_link_libraries(ipc_toolkit PRIVATE rational::rational)
endif()

# Faster unordered map
if(IPC_TOOLKIT_WITH_ROBIN_MAP)
  include(robin_map)
  target_link_libraries(ipc_toolkit PUBLIC tsl::robin_map)
endif()

# Hashes
if(IPC_TOOLKIT_WITH_ABSEIL)
  include(abseil)
  target_link_libraries(ipc_toolkit PUBLIC absl::hash)
endif()

# Intervals
if(IPC_TOOLKIT_WITH_FILIB)
  include(filib)
  target_link_libraries(ipc_toolkit PUBLIC filib::filib)
endif()

# Extra warnings (link last for highest priority)
include(ipc_toolkit_warnings)
target_link_libraries(ipc_toolkit PRIVATE ipc::toolkit::warnings)

################################################################################
# Compiler options
################################################################################

## SIMD support
if(IPC_TOOLKIT_WITH_SIMD)
  # Figure out SIMD support
  message(STATUS "Testing SIMD capabilities...")
  find_package(SIMD)
  # Add SIMD flags to compiler flags
  message(STATUS "Using SIMD flags: ${SIMD_FLAGS}")
  target_compile_options(ipc_toolkit PRIVATE ${SIMD_FLAGS})
else()
  message(STATUS "SIMD support disabled")
endif()

# For MSVC, do not use the min and max macros.
if(MSVC)
  target_compile_definitions(ipc_toolkit PRIVATE NOMINMAX)
endif()

# Use C++17
target_compile_features(ipc_toolkit PUBLIC cxx_std_17)

################################################################################
# CUDA
################################################################################

if(IPC_TOOLKIT_WITH_CUDA)
  if(NOT CMAKE_CUDA_COMPILER)
    message(FATAL_ERROR "CUDA support requested but no CUDA compiler found!")
  endif()

  # Enable CUDA support
  enable_language(CUDA)

  # We need to explicitly state that we need all CUDA files in the particle
  # library to be built with -dc as the member functions could be called by
  # other libraries and executables.
  set_target_properties(ipc_toolkit PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

  # Use the same CUDA architectures Scalable CCD
  get_target_property(CMAKE_CUDA_ARCHITECTURES scalable_ccd CUDA_ARCHITECTURES)
  set_target_properties(ipc_toolkit PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
endif()

################################################################################
# Tests
################################################################################

# Enable unit testing at the root level
if(IPC_TOOLKIT_BUILD_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()

################################################################################
# Code Coverage
################################################################################

if(IPC_TOOLKIT_WITH_CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(ipc_toolkit PRIVATE
    -g         # generate debug info
    --coverage # sets all required flags
    -fprofile-update=atomic
  )
  target_link_options(ipc_toolkit PUBLIC
    --coverage
    -fprofile-update=atomic
  )
endif()

################################################################################
# Python bindings
################################################################################

if(IPC_TOOLKIT_BUILD_PYTHON)
  add_subdirectory(python)
endif()

################################################################################
# Xcode
################################################################################

if (CMAKE_GENERATOR STREQUAL "Xcode")
  set(CMAKE_XCODE_SCHEME_LAUNCH_CONFIGURATION "${CMAKE_BUILD_TYPE}")
  set_target_properties(ipc_toolkit PROPERTIES XCODE_GENERATE_SCHEME ON)
  if(IPC_TOOLKIT_BUILD_TESTS)
    set_target_properties(ipc_toolkit_tests PROPERTIES XCODE_GENERATE_SCHEME ON)
  endif()
  if(IPC_TOOLKIT_BUILD_PYTHON)
    set_target_properties(ipctk PROPERTIES XCODE_GENERATE_SCHEME ON)
  endif()
endif()