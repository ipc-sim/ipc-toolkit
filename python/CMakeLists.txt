cmake_minimum_required(VERSION 3.8)

###############################################################################
# Download NumpyEigen
function(ipc_toolkit_download_pybind11)
  ipc_toolkit_download_project(pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.6.2
  )
endfunction()

if(NOT TARGET pybind11::module)
  ipc_toolkit_download_pybind11()
  add_subdirectory(
    ${IPC_TOOLKIT_EXTERNAL}/pybind11
    ${CMAKE_CURRENT_BINARY_DIR}/pybind11)
endif()

###############################################################################
# Create Python module for Rigid IPC

pybind11_add_module(ipctk
  src/bindings.cpp

  src/barrier/barrier.cpp

  src/broad_phase/broad_phase.cpp

  src/broad_phase/collision_candidate.cpp

  src/collision_constraint.cpp

  src/distance/distance_type.cpp
  src/distance/edge_edge.cpp
  src/distance/edge_edge_mollifier.cpp
  src/distance/line_line.cpp
  src/distance/point_edge.cpp
  src/distance/point_line.cpp
  src/distance/point_plane.cpp
  src/distance/point_point.cpp
  src/distance/point_triangle.cpp

  src/ipc.cpp

  src/utils/logger.cpp
  src/utils/thread_limiter.cpp
)
target_link_libraries(ipctk PRIVATE ${PROJECT_NAME})
target_link_libraries(ipctk PRIVATE spdlog::spdlog)

# Extra warnings
# target_link_libraries(ipctk PRIVATE warnings::all)

# Move Python library to build root
set_target_properties(ipctk PROPERTIES LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}")
