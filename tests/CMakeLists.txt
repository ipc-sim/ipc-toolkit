################################################################################
# Tests
################################################################################

add_executable(ipc_toolkit_tests
  main.cpp

  # Test barrier and dertivatives
  barrier/test_barrier.cpp

  # Test CCD
  broad_phase/benchmark_broad_phase.cpp
  broad_phase/brute_force_comparison.cpp
  broad_phase/test_aabb.cpp
  broad_phase/test_broad_phase.cpp
  broad_phase/test_spatial_hash.cpp
  ccd/benchmark_ccd.cpp
  ccd/collision_generator.cpp
  ccd/test_ccd.cpp
  ccd/test_ccd_benchmark.cpp
  ccd/test_point_edge_ccd.cpp
  ccd/test_point_triangle_ccd.cpp
  ccd/test_edge_edge_ccd.cpp
  ccd/test_ccd.cpp
  test_cfl.cpp

  # Test distances and dertivatives
  distance/test_distance_type.cpp
  distance/test_edge_edge.cpp
  distance/test_line_line.cpp
  distance/test_point_edge.cpp
  distance/test_point_line.cpp
  distance/test_point_plane.cpp
  distance/test_point_point.cpp
  distance/test_point_triangle.cpp

  # Test friction
  friction/friction_data_generator.cpp
  friction/test_closest_point.cpp
  friction/test_force_jacobian.cpp
  friction/test_friction.cpp
  friction/test_normal_force_magnitude.cpp
  friction/test_relative_velocity.cpp
  friction/test_smooth_friction_mollifier.cpp
  friction/test_tangent_basis.cpp

  # Test general interface
  test_ipc.cpp

  # Test intersection checks
  test_has_intersections.cpp

  # Utilities for tests
  test_utils.cpp

  benchmark_eigen.cpp
)

target_include_directories(ipc_toolkit_tests PUBLIC ".")

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

if (IPC_TOOLKIT_TEST_CCD_BENCHMARK)
  include(ccd_query_io)
  target_link_libraries(ipc_toolkit_tests PUBLIC ccd_io::ccd_io)
  target_compile_definitions(ipc_toolkit_tests PRIVATE IPC_TOOLKIT_TEST_CCD_BENCHMARK)
  if(NOT (IPC_TOOLKIT_CCD_BENCHMARK_DIR STREQUAL ""))
    message(STATUS "Using CCD benchmark directory: ${IPC_TOOLKIT_CCD_BENCHMARK_DIR}")
    target_compile_definitions(ipc_toolkit_tests PRIVATE IPC_TOOLKIT_CCD_BENCHMARK_DIR="${IPC_TOOLKIT_CCD_BENCHMARK_DIR}")
  endif()
  if(NOT (IPC_TOOLKIT_CCD_NEW_BENCHMARK_DIR STREQUAL ""))
    message(STATUS "Using new CCD benchmark directory: ${IPC_TOOLKIT_CCD_NEW_BENCHMARK_DIR}")
    target_compile_definitions(ipc_toolkit_tests PRIVATE IPC_TOOLKIT_CCD_NEW_BENCHMARK_DIR="${IPC_TOOLKIT_CCD_NEW_BENCHMARK_DIR}")
  endif()
endif()

if(IPC_TOOLKIT_WITH_CUDA)
  find_package(CUDAToolkit)
  target_link_libraries(ipc_toolkit_tests PRIVATE CUDA::cudart)
endif()

# Extra warnings (link last for highest priority)
include(ipc_toolkit_warnings)
target_link_libraries(ipc_toolkit_tests PRIVATE ipc::toolkit::warnings)

################################################################################
# Compiler options
################################################################################

target_compile_definitions(ipc_toolkit_tests PUBLIC CATCH_CONFIG_ENABLE_BENCHMARKING)

# Test data directory definition
if(WIN32)
  target_compile_definitions(ipc_toolkit_tests PUBLIC TEST_DATA_DIR_CSTR="${CMAKE_CURRENT_SOURCE_DIR}\\\\data\\\\")
else()
  target_compile_definitions(ipc_toolkit_tests PUBLIC TEST_DATA_DIR_CSTR="${CMAKE_CURRENT_SOURCE_DIR}/data/")
endif()

################################################################################
# Register tests
################################################################################

FetchContent_GetProperties(catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
catch_discover_tests(ipc_toolkit_tests)
