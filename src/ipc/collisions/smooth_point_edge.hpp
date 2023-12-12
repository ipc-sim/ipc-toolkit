#pragma once

#include <ipc/distance/distance_type.hpp>

namespace ipc {

    /// @brief Compute pointwise potential for a point p and a point specified by uv on edge [e0, e1]
    /// @param p One point outside of edge
    /// @param e0 One end point of the edge
    /// @param e1 One end point of the edge
    /// @param uv Barycentric coordinate
    /// @param dhat The effective distance of barrier
    /// @param alpha The effective angle of barrier
    double smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3d>& p,
        const Eigen::Ref<const VectorMax3d>& e0,
        const Eigen::Ref<const VectorMax3d>& e1,
        const double &uv,
        const double &dhat,
        const double &alpha);

    /// @brief Compute potential for a point p and an edge [e0, e1], integrated over the edge
    double smooth_point_edge_potential(
        const Eigen::Ref<const VectorMax3d>& p,
        const Eigen::Ref<const VectorMax3d>& e0,
        const Eigen::Ref<const VectorMax3d>& e1,
        const double &dhat,
        const double &alpha);
}