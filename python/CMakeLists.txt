###############################################################################
# Create Python module for IPC Toolkit
include(pybind11)

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

if(IPC_TOOLKIT_WITH_LOGGER)
  include(spdlog)
  target_link_libraries(ipctk PRIVATE spdlog::spdlog)
else()
  include(fmt)
  target_link_libraries(ipctk PRIVATE fmt::fmt)
endif()

# Extra warnings
# target_link_libraries(ipctk PRIVATE IPCToolkit::warnings)

# Move Python library to build root
set_target_properties(ipctk PROPERTIES LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}")
