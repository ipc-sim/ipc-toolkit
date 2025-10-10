#pragma once

#include "primitive.hpp"

#include <ipc/smooth_contact/distance/mollifier.hpp>

namespace ipc {
class Edge2 : public Primitive {
public:
    constexpr static int N_CORE_POINTS = 2;
    constexpr static int DIM = 2;
    // d is a vector from closest point on the edge to the point outside of the
    // edge
    Edge2(
        const long& id,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const VectorMax3d& d,
        const ParameterType& param);

    int n_vertices() const override;
    int n_dofs() const override { return n_vertices() * DIM; }

    double potential(const Vector2d& d, const Vector4d& x) const;
    Vector6d grad(const Vector2d& d, const Vector4d& x) const;
    Matrix6d hessian(const Vector2d& d, const Vector4d& x) const;
};
} // namespace ipc
