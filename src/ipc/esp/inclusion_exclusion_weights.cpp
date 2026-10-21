#include "inclusion_exclusion_weights.hpp"

namespace ipc {

void compute_inclusion_exclusion_weights(
    const CollisionMesh& mesh,
    std::vector<int>& edge_weights,
    std::vector<int>& vertex_weights)
{
    edge_weights.assign(mesh.num_edges(), 1);
    vertex_weights.assign(mesh.num_vertices(), 1);
    for (int f = 0; f < mesh.faces().rows(); f++) {
        for (int j = 0; j < 3; j++) {
            --edge_weights[mesh.faces_to_edges()(f, j)];
            ++vertex_weights[mesh.faces()(f, j)];
        }
    }
    for (int e = 0; e < mesh.edges().rows(); e++) {
        --vertex_weights[mesh.edges()(e, 0)];
        --vertex_weights[mesh.edges()(e, 1)];
    }
}

} // namespace ipc
