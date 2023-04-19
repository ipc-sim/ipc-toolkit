#pragma once

#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class PlaneVertexConstraint : public CollisionConstraint {
public:
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

    VectorMax3d plane_origin;
    VectorMax3d plane_normal;
    long vertex_id;

protected:
    double compute_distance(const VectorMax12d& point) const override;

    VectorMax12d
    compute_distance_gradient(const VectorMax12d& point) const override;

    MatrixMax12d
    compute_distance_hessian(const VectorMax12d& point) const override;
};

} // namespace ipc
