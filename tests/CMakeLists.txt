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

################################################################################
# IPC Toolkit Tests Executable
################################################################################

# Add an empty library and fill in the list of sources in `src/tests/CMakeLists.txt`.
add_executable(ipc_toolkit_tests)

# Fill in configuration options
include(ipc_toolkit_tests_data)
configure_file(
  "${IPC_TOOLKIT_TESTS_SOURCE_DIR}/config.hpp.in"
  "${IPC_TOOLKIT_TESTS_SOURCE_DIR}/config.hpp")

# Add source and header files to ipc_toolkit_tests
add_subdirectory("${IPC_TOOLKIT_TESTS_SOURCE_DIR}")

# Public include directory for IPC Toolkit tests
target_include_directories(ipc_toolkit_tests PRIVATE "${IPC_TOOLKIT_TESTS_INCLUDE_DIR}")

# Folder name for IDE
set_target_properties(ipc_toolkit_tests PROPERTIES FOLDER "SRC")
get_target_property(IPC_TOOLKIT_TESTS_SOURCES ipc_toolkit_tests SOURCES)
source_group(TREE "${PROJECT_SOURCE_DIR}/tests" FILES ${IPC_TOOLKIT_TESTS_SOURCES})

################################################################################
# Required Libraries
################################################################################

target_link_libraries(ipc_toolkit_tests PRIVATE ipc::toolkit)

# libigl (use igl::edges)
include(libigl)
target_link_libraries(ipc_toolkit_tests PRIVATE igl::core)

# OneTBB (use tbb::global_control and tbb::parallel_sort)
include(onetbb)
target_link_libraries(ipc_toolkit_tests PRIVATE TBB::tbb)

include(catch2)
target_link_libraries(ipc_toolkit_tests PRIVATE Catch2::Catch2)

include(finite_diff)
target_link_libraries(ipc_toolkit_tests PRIVATE finitediff::finitediff)

include(json)
target_link_libraries(ipc_toolkit_tests PRIVATE nlohmann_json::nlohmann_json)

if(IPC_TOOLKIT_WITH_CUDA)
  # We need to explicitly state that we need all CUDA files in the particle
  # library to be built with -dc as the member functions could be called by
  # other libraries and executables.
  set_target_properties(ipc_toolkit_tests PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

  # Use the same CUDA architectures Scalable CCD
  # get_target_property(CMAKE_CUDA_ARCHITECTURES scalable_ccd CUDA_ARCHITECTURES)
  set_target_properties(ipc_toolkit_tests PROPERTIES CUDA_ARCHITECTURES "native")
endif()

if (IPC_TOOLKIT_TESTS_CCD_BENCHMARK)
  include(ccd_query_io)
  target_link_libraries(ipc_toolkit_tests PRIVATE ccd_io::ccd_io)
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

target_compile_definitions(ipc_toolkit_tests PRIVATE CATCH_CONFIG_ENABLE_BENCHMARKING)

if (FILIB_BUILD_SHARED_LIB AND WIN32)
  # Copy DLLs to the output directory
  add_custom_command(
    TARGET ipc_toolkit_tests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:ipc_toolkit_tests> $<TARGET_FILE_DIR:ipc_toolkit_tests>
    COMMAND_EXPAND_LISTS
  )
endif()

################################################################################
# Register tests
################################################################################

FetchContent_GetProperties(catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
if(CMAKE_GENERATOR STREQUAL "Xcode")
  # https://github.com/catchorg/Catch2/issues/2411
  # Requires adding `-C Release` to `ctest` command
  catch_discover_tests(ipc_toolkit_tests DISCOVERY_MODE PRE_TEST)
else()
  catch_discover_tests(ipc_toolkit_tests)
endif()