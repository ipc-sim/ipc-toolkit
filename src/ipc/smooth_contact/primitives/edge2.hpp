#pragma once

#include "primitive.hpp"

#include <ipc/smooth_contact/distance/mollifier.hpp>

namespace ipc {
class Edge2 : public Primitive {
public:
    static constexpr int N_CORE_POINTS = 2;
    static constexpr int DIM = 2;
    // d is a vector from closest point on the edge to the point outside of the
    // edge
    Edge2(
        const index_t id,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const VectorMax3d& d,
        const SmoothContactParameters& params);

    int n_vertices() const override;
    int n_dofs() const override { return n_vertices() * DIM; }

    double potential(
        Eigen::ConstRef<Eigen::Vector2d> d,
        Eigen::ConstRef<Eigen::Vector4d> x) const;
    Vector6d grad(
        Eigen::ConstRef<Eigen::Vector2d> d,
        Eigen::ConstRef<Eigen::Vector4d> x) const;
    Matrix6d hessian(
        Eigen::ConstRef<Eigen::Vector2d> d,
        Eigen::ConstRef<Eigen::Vector4d> x) const;
};
} // namespace ipc
