#include "vertex_to_min_edge.hpp"

#include <algorithm> // std::min/max

namespace ipc {

std::vector<size_t>
vertex_to_min_edge(size_t num_vertices, const Eigen::MatrixXi& E)
{
    std::vector<size_t> V2E(num_vertices, E.rows() + 1);
    // Column first because colmajor
    for (size_t ej = 0; ej < E.cols(); ej++) {
        for (size_t ei = 0; ei < E.rows(); ei++) {
            const size_t& vi = E(ei, ej);
            V2E[vi] = std::min(V2E[vi], ei);
        }
    }
    return V2E;
}

} // namespace ipc
