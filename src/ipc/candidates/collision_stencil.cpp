#include "collision_stencil.hpp"

namespace ipc {

std::array<VectorMax3d, 4> CollisionStencil::vertices(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
    const VectorMax3d nan_vector = VectorMax3d::Constant(
        std::numeric_limits<double>::signaling_NaN(), vertices.cols());

    const std::array<long, 4> vertex_ids = this->vertex_ids(edges, faces);

    std::array<VectorMax3d, 4> stencil_vertices;
    for (int i = 0; i < 4; i++) {
        stencil_vertices[i] =
            vertex_ids[i] < 0 ? nan_vector : vertices.row(vertex_ids[i]);
    }

    return stencil_vertices;
}

} // namespace ipc