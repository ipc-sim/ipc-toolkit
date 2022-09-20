#pragma once

#include <Eigen/Core>

namespace ipc {

bool is_edge_intersecting_triangle(
    const Eigen::Vector3d& e0,
    const Eigen::Vector3d& e1,
    const Eigen::Vector3d& t0,
    const Eigen::Vector3d& t1,
    const Eigen::Vector3d& t2);

} // namespace ipc
