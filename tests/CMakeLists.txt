cmake_minimum_required(VERSION 3.1)
################################################################################

add_executable(${PROJECT_NAME}_tests
  main.cpp

  # Test barrier and dertivatives
  barrier/test_barrier.cpp

  # Test CCD
  broad_phase/benchmark_spatial_hash.cpp
  broad_phase/test_broad_phase.cpp
  broad_phase/test_hash_grid.cpp
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
  friction/test_friction.cpp
  friction/test_tangent_basis.cpp
  friction/test_closest_point.cpp
  friction/test_relative_displacement.cpp

  # Test general interface
  test_ipc.cpp

  # Test intersection checks
  test_has_intersections.cpp

  # Utilities for tests
  test_utils.cpp
)

target_include_directories(${PROJECT_NAME}_tests PUBLIC ".")

################################################################################
# Required Libraries
################################################################################

target_link_libraries(${PROJECT_NAME}_tests PUBLIC ${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME}_tests PUBLIC Catch2::Catch2 warnings::all)

target_link_libraries(${PROJECT_NAME}_tests PUBLIC FiniteDiff::FiniteDiff)

if(IPC_TOOLKIT_WITH_LOGGER)
  target_link_libraries(${PROJECT_NAME}_tests PRIVATE spdlog::spdlog)
else()
  target_link_libraries(${PROJECT_NAME}_tests PRIVATE fmt::fmt)
endif()

################################################################################
# Compiler options
################################################################################

target_compile_definitions(${PROJECT_NAME}_tests PUBLIC CATCH_CONFIG_ENABLE_BENCHMARKING)

# Test data directory definition
if(WIN32)
  target_compile_definitions(${PROJECT_NAME}_tests PUBLIC TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}\\\\data\\\\")
else()
  target_compile_definitions(${PROJECT_NAME}_tests PUBLIC TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/")
endif()

################################################################################
# Register tests
################################################################################

foreach(source IN ITEMS ${test_sources})
   source_group("tests" FILES "${source}")
endforeach()

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(${PROJECT_NAME}_tests)
