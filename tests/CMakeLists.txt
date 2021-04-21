cmake_minimum_required(VERSION 3.1)
################################################################################

add_executable(${PROJECT_NAME}_tests
  main.cpp

  # Test barrier and dertivatives
  barrier/test_barrier.cpp

  # Test CCD
  ccd/test_ccd.cpp
  spatial_hash/test_spatial_hash.cpp
  spatial_hash/test_hash_grid.cpp
  spatial_hash/test_spatial_hash.cpp
  spatial_hash/benchmark_spatial_hash.cpp

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

target_link_libraries(${PROJECT_NAME}_tests PUBLIC spdlog::spdlog)

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
