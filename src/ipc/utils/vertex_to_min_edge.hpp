#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <vector>

namespace ipc {

std::vector<size_t>
vertex_to_min_edge(size_t num_vertices, Eigen::ConstRef<Eigen::MatrixXi> edges);

} // namespace ipc
