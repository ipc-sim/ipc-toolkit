################################################################################
# Tests
################################################################################

include(CMakeDependentOption)
cmake_dependent_option(IPC_TOOLKIT_TESTS_CCD_BENCHMARK "Enable CCD benchmark test" OFF "IPC_TOOLKIT_BUILD_TESTS" OFF)
if(IPC_TOOLKIT_TESTS_CCD_BENCHMARK)
  mark_as_advanced(CLEAR IPC_TOOLKIT_TESTS_CCD_BENCHMARK_DIR)
  mark_as_advanced(CLEAR IPC_TOOLKIT_TESTS_NEW_CCD_BENCHMARK_DIR)
  set(IPC_TOOLKIT_TESTS_CCD_BENCHMARK_DIR "" CACHE PATH "Path to the CCD benchmark directory")
  set(IPC_TOOLKIT_TESTS_NEW_CCD_BENCHMARK_DIR "" CACHE PATH "Path to the new CCD benchmark directory")
else()
  mark_as_advanced(FORCE IPC_TOOLKIT_TESTS_CCD_BENCHMARK_DIR)
  mark_as_advanced(FORCE IPC_TOOLKIT_TESTS_NEW_CCD_BENCHMARK_DIR)
endif()

### Configuration
set(IPC_TOOLKIT_TESTS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/tests")
set(IPC_TOOLKIT_TESTS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

### Data directory
cmake_path(APPEND IPC_TOOLKIT_TESTS_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}" "data")

################################################################################
# IPC Toolkit Tests Executable
################################################################################

# Add an empty library and fill in the list of sources in `src/tests/CMakeLists.txt`.
add_executable(ipc_toolkit_tests)

# Fill in configuration options
configure_file(
  "${IPC_TOOLKIT_TESTS_SOURCE_DIR}/config.hpp.in"
  "${IPC_TOOLKIT_TESTS_SOURCE_DIR}/config.hpp")

# Add source and header files to ipc_toolkit_tests
add_subdirectory("${IPC_TOOLKIT_TESTS_SOURCE_DIR}")

# Public include directory for IPC Toolkit tests
target_include_directories(ipc_toolkit_tests PUBLIC "${IPC_TOOLKIT_TESTS_INCLUDE_DIR}")

################################################################################
# Required Libraries
################################################################################

target_link_libraries(ipc_toolkit_tests PUBLIC ipc::toolkit)

include(catch2)
target_link_libraries(ipc_toolkit_tests PUBLIC Catch2::Catch2)

include(finite_diff)
target_link_libraries(ipc_toolkit_tests PUBLIC finitediff::finitediff)

include(json)
target_link_libraries(ipc_toolkit_tests PUBLIC nlohmann_json::nlohmann_json)

if (IPC_TOOLKIT_TESTS_CCD_BENCHMARK)
  include(ccd_query_io)
  target_link_libraries(ipc_toolkit_tests PUBLIC ccd_io::ccd_io)
  if(NOT (IPC_TOOLKIT_CCD_BENCHMARK_DIR STREQUAL ""))
    message(STATUS "Using CCD benchmark directory: ${IPC_TOOLKIT_CCD_BENCHMARK_DIR}")
  endif()
  if(NOT (IPC_TOOLKIT_CCD_NEW_BENCHMARK_DIR STREQUAL ""))
    message(STATUS "Using new CCD benchmark directory: ${IPC_TOOLKIT_CCD_NEW_BENCHMARK_DIR}")
  endif()
endif()

# Extra warnings (link last for highest priority)
include(ipc_toolkit_warnings)
target_link_libraries(ipc_toolkit_tests PRIVATE ipc::toolkit::warnings)

################################################################################
# Compiler options
################################################################################

target_compile_definitions(ipc_toolkit_tests PUBLIC CATCH_CONFIG_ENABLE_BENCHMARKING)

################################################################################
# Register tests
################################################################################

FetchContent_GetProperties(catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
catch_discover_tests(ipc_toolkit_tests)