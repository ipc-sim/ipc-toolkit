#pragma once

#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct PlaneVertexConstraint : CollisionConstraint {
    PlaneVertexConstraint(
        const VectorMax3d& plane_origin,
        const VectorMax3d& plane_normal,
        const long vertex_id);

    int num_vertices() const override { return 1; };
    std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return { { vertex_id, -1, -1, -1 } };
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override;

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override;

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override;

    VectorMax3d plane_origin;
    VectorMax3d plane_normal;
    long vertex_id;
};

} // namespace ipc
