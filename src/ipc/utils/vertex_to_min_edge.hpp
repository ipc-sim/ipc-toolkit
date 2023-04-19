#pragma once

#include <Eigen/Core>

#include <vector>

namespace ipc {

std::vector<size_t>
vertex_to_min_edge(size_t num_vertices, const Eigen::MatrixXi& edges);

} // namespace ipc
