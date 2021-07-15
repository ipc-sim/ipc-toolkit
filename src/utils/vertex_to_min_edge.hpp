#pragma once

#include <vector>

#include <Eigen/Core>

namespace ipc {

std::vector<size_t>
vertex_to_min_edge(size_t num_vertices, const Eigen::MatrixXi& E);

} // namespace ipc
