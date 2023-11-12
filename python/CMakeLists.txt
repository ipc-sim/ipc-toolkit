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
  src/broad_phase/bvh.cpp
  src/broad_phase/hash_grid.cpp
  src/broad_phase/spatial_hash.cpp
  src/broad_phase/sweep_and_tiniest_queue.cpp
  src/broad_phase/voxel_size_heuristic.cpp

  src/candidates/candidates.cpp
  src/candidates/collision_stencil.cpp
  src/candidates/continuous_collision_candidate.cpp
  src/candidates/edge_edge.cpp
  src/candidates/edge_face.cpp
  src/candidates/edge_vertex.cpp
  src/candidates/face_vertex.cpp
  src/candidates/vertex_vertex.cpp

  src/ccd/aabb.cpp
  src/ccd/ccd.cpp
  src/ccd/additive_ccd.cpp
  src/ccd/inexact_point_edge.cpp
  src/ccd/nonlinear_ccd.cpp
  src/ccd/point_static_plane.cpp

  src/collisions/collision_constraint.cpp
  src/collisions/collision_constraints.cpp
  src/collisions/edge_edge.cpp
  src/collisions/edge_vertex.cpp
  src/collisions/face_vertex.cpp
  src/collisions/plane_vertex.cpp
  src/collisions/vertex_vertex.cpp

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
  src/friction/friction_constraints.cpp
  src/friction/normal_force_magnitude.cpp
  src/friction/relative_velocity.cpp
  src/friction/smooth_friction_mollifier.cpp
  src/friction/tangent_basis.cpp

  src/friction/constraints/edge_edge.cpp
  src/friction/constraints/edge_vertex.cpp
  src/friction/constraints/face_vertex.cpp
  src/friction/constraints/friction_constraint.cpp
  src/friction/constraints/vertex_vertex.cpp

  src/implicits/plane.cpp

  src/utils/area_gradient.cpp
  src/utils/eigen_ext.cpp
  src/utils/interval.cpp
  src/utils/intersection.cpp
  src/utils/logger.cpp
  src/utils/thread_limiter.cpp
  src/utils/vertex_to_min_edge.cpp
  src/utils/world_bbox_diagonal_length.cpp

  src/collision_mesh.cpp
  src/ipc.cpp
)

target_link_libraries(ipctk PRIVATE ipc::toolkit)
target_include_directories(ipctk PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Extra warnings
# target_link_libraries(ipctk PRIVATE IPCToolkit::warnings)

# Move Python library to build root
set_target_properties(ipctk PROPERTIES LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}")
