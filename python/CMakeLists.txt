###############################################################################
# Create Python module for IPC Toolkit
include(pybind11)

pybind11_add_module(ipctk
  src/bindings.cpp

  src/barrier/adaptive_stiffness.cpp
  src/barrier/barrier.cpp

  src/broad_phase/aabb.cpp
  src/broad_phase/broad_phase.cpp
  src/broad_phase/brute_force.cpp
  src/broad_phase/collision_candidate.cpp
  src/broad_phase/hash_grid.cpp
  src/broad_phase/spatial_hash.cpp
  src/broad_phase/sweep.cpp
  src/broad_phase/voxel_size_heuristic.cpp

  src/ccd/aabb.cpp
  src/ccd/ccd.cpp
  src/ccd/inexact_point_edge.cpp
  src/ccd/point_static_plane.cpp

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

  src/collision_mesh.cpp

  src/ipc.cpp

  src/utils/logger.cpp
  src/utils/thread_limiter.cpp
)

target_link_libraries(ipctk PRIVATE ipc::toolkit)

if(IPC_TOOLKIT_WITH_LOGGER)
  include(spdlog)
  ipc_toolkit_target_link_system_libraries(ipctk PRIVATE spdlog::spdlog)
else()
  include(fmt)
  ipc_toolkit_target_link_system_libraries(ipctk PRIVATE fmt::fmt)
endif()

# Extra warnings
# target_link_libraries(ipctk PRIVATE IPCToolkit::warnings)

# Move Python library to build root
set_target_properties(ipctk PROPERTIES LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}")
