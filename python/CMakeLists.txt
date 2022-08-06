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
  src/broad_phase/sweep_and_tiniest_queue.cpp
  src/broad_phase/voxel_size_heuristic.cpp

  src/ccd/aabb.cpp
  src/ccd/ccd.cpp
  src/ccd/inexact_point_edge.cpp
  src/ccd/point_static_plane.cpp

  src/distance/distance_type.cpp
  src/distance/edge_edge.cpp
  src/distance/edge_edge_mollifier.cpp
  src/distance/line_line.cpp
  src/distance/point_edge.cpp
  src/distance/point_line.cpp
  src/distance/point_plane.cpp
  src/distance/point_point.cpp
  src/distance/point_triangle.cpp

  src/friction/closest_point.cpp
  src/friction/friction_constraint.cpp
  src/friction/friction.cpp
  src/friction/normal_force_magnitude.cpp
  src/friction/relative_displacement.cpp
  src/friction/smooth_friction_mollifier.cpp
  src/friction/tangent_basis.cpp

  src/utils/eigen_ext.cpp
  src/utils/intersection.cpp
  src/utils/logger.cpp
  src/utils/thread_limiter.cpp
  src/utils/world_bbox_diagonal_length.cpp

  src/collision_constraint.cpp
  src/collision_mesh.cpp
  src/ipc.cpp
)

target_link_libraries(ipctk PRIVATE ipc::toolkit)

# Extra warnings
# target_link_libraries(ipctk PRIVATE IPCToolkit::warnings)

# Move Python library to build root
set_target_properties(ipctk PROPERTIES LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}")
