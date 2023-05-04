#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Compute the gradient of an edge's length.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @return The gradient of the edge's length wrt e0, and e1.
VectorMax6d edge_length_gradient(
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1);

/// @brief Compute the gradient of the area of a triangle.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The gradient of the triangle's area t0, t1, and t2.
Vector9d triangle_area_gradient(
    const Eigen::Ref<const Eigen::Vector3d>& t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2);

namespace autogen {

    // dA is (9×1) flattened in column-major order
    void triangle_area_gradient(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double dA[9]);

} // namespace autogen

} // namespace ipc
