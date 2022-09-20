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
  friction/test_relative_displacement.cpp
  friction/test_smooth_friction_mollifier.cpp
  friction/test_tangent_basis.cpp

  # Test general interface
  test_ipc.cpp

  # Test intersection checks
  test_has_intersections.cpp

  # Utilities for tests
  test_utils.cpp

  static_vs_dynamic.cpp
)

target_include_directories(ipc_toolkit_tests PUBLIC ".")

################################################################################
# Required Libraries
################################################################################

target_link_libraries(ipc_toolkit_tests PUBLIC ipc::toolkit)

include(ipc_toolkit_warnings)
target_link_libraries(ipc_toolkit_tests PRIVATE ipc::toolkit::warnings)

include(catch2)
ipc_toolkit_target_link_system_libraries(ipc_toolkit_tests PUBLIC Catch2::Catch2)

include(finite_diff)
ipc_toolkit_target_link_system_libraries(ipc_toolkit_tests PUBLIC finitediff::finitediff)

include(json)
ipc_toolkit_target_link_system_libraries(ipc_toolkit_tests PUBLIC nlohmann_json::nlohmann_json)

if(IPC_TOOLKIT_WITH_CUDA)
  find_package(CUDAToolkit)
  ipc_toolkit_target_link_system_libraries(ipc_toolkit_tests PRIVATE CUDA::cudart)
endif()

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

set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)

# Include Catch2 and provide function `catch_discover_tests` to register tests.
FetchContent_GetProperties(catch2)
include("${catch2_SOURCE_DIR}/contrib/Catch.cmake")

catch_discover_tests(ipc_toolkit_tests)
