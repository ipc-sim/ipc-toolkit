cmake_minimum_required(VERSION 3.1)
################################################################################

add_executable(${PROJECT_NAME}_tests
  main.cpp
  ${IPC_TOOLKIT_EXTERNAL}/finite-diff/finitediff.cpp

  barrier/test_barrier.cpp

  # spatial_hash/test_spatial_hash.cpp

  distance/test_distance_type.cpp
  distance/test_edge_edge.cpp
  distance/test_line_line.cpp
  distance/test_point_edge.cpp
  distance/test_point_line.cpp
  distance/test_point_plane.cpp
  distance/test_point_point.cpp
  distance/test_point_triangle.cpp

  test_ipc.cpp

  # Test CCD
  # ccd/test_ccd.cpp
  # ccd/test_time_of_impact.cpp
  # ccd/test_hash_grid.cpp
)

target_link_libraries(${PROJECT_NAME}_tests PUBLIC ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME}_tests PUBLIC Catch2::Catch2 warnings::all)
# target_link_libraries(${PROJECT_NAME}_tests PUBLIC FiniteDiff::FiniteDiff)
target_include_directories(${PROJECT_NAME}_tests PUBLIC ${IPC_TOOLKIT_EXTERNAL}/finite-diff)
target_link_libraries(${PROJECT_NAME}_tests PUBLIC spdlog::spdlog)

target_compile_definitions(${PROJECT_NAME}_tests PUBLIC CATCH_CONFIG_ENABLE_BENCHMARKING)

foreach(source IN ITEMS ${test_sources})
   source_group("tests" FILES "${source}")
endforeach()

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(${PROJECT_NAME}_tests)
