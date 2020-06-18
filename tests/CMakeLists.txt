cmake_minimum_required(VERSION 3.1)
################################################################################

add_executable(${PROJECT_NAME}_tests
  main.cpp

  barrier/test_barrier.cpp

  # Test CCD
  # ccd/test_ccd.cpp
  # ccd/test_time_of_impact.cpp
  # ccd/test_hash_grid.cpp
  #
  # geometry/test_distance.cpp
  # geometry/test_intersection.cpp
)

target_link_libraries(${PROJECT_NAME}_tests PUBLIC ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME}_tests PUBLIC Catch2::Catch2 warnings::all)
target_link_libraries(${PROJECT_NAME}_tests PUBLIC FiniteDiff::FiniteDiff)

# target_compile_definitions(${PROJECT_NAME}_tests PUBLIC CATCH_CONFIG_ENABLE_BENCHMARKING)

foreach(source IN ITEMS ${test_sources})
   source_group("tests" FILES "${source}")
endforeach()

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(${PROJECT_NAME}_tests)
