#include "esp_collision.hpp"

#include "ipc/gcp/distance/point_edge.hpp"

#include <ipc/config.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_point.hpp>

namespace ipc {

std::vector<index_t> ESPCollision::vertex_ids() const
{
    std::vector<index_t> ids;
    ids.reserve(num_vertices());
    for (int i = 0; i < num_vertices(); ++i) {
        ids.push_back(vertex_id(i));
    }
    return ids;
}

Eigen::VectorXd ESPCollision::dof(Eigen::ConstRef<Eigen::MatrixXd> X) const
{
    const int DIM = X.cols();
    Eigen::VectorXd x(num_vertices() * DIM);
    if (DIM == 2) {
        for (int i = 0; i < num_vertices(); i++) {
            x.segment<2>(i * 2) = X.row(vertex_id(i));
        }
    } else if (DIM == 3) {
        for (int i = 0; i < num_vertices(); i++) {
            x.segment<3>(i * 3) = X.row(vertex_id(i));
        }
    } else {
        throw std::runtime_error("Invalid dimension!");
    }
    return x;
}

Eigen::VectorXd ESPCollision::dof(VertexMatrixView<3> X_extended) const
{
    Eigen::VectorXd x(num_vertices() * 3);
    for (int i = 0; i < num_vertices(); i++) {
        assert(vertex_id(i) < X_extended.rows());
        x.segment<3>(i * 3) = X_extended(vertex_id(i));
    }
    return x;
}

Eigen::VectorXd ESPCollision::dof(VertexMatrixView<2> X_extended) const
{
    Eigen::VectorXd x(num_vertices() * 2);
    for (int i = 0; i < num_vertices(); i++) {
        assert(vertex_id(i) < X_extended.rows());
        x.segment<2>(i * 2) = X_extended(vertex_id(i));
    }
    return x;
}
} // namespace ipc