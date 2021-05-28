#pragma once

#include <Eigen/Core>

namespace ipc {

/// Construct a matrix that maps from the faces' edges to rows in the edges
/// matrix
Eigen::MatrixXi
faces_to_edges(const Eigen::MatrixXi& F, const Eigen::MatrixXi& E);

} // namespace ipc
