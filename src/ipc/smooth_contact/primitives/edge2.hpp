#pragma once

#include <ipc/smooth_contact/distance/mollifier.hpp>
#include "primitive.hpp"

namespace ipc {
    class Edge2 : public Primitive {
    public:
        constexpr static int n_core_points = 2;
        constexpr static int dim = 2;
        // d is a vector from closest point on the edge to the point outside of the edge
        Edge2(
            const long& id,
            const CollisionMesh& mesh,
            const Eigen::MatrixXd& vertices,
            const VectorMax3d& d,
            const double& alpha,
            const double& beta);

        int n_vertices() const override;
        int n_dofs() const override { return n_vertices() * dim; }

        double potential(const Vector2d& d, const Vector4d& x) const;
        Vector6d grad(const Vector2d& d, const Vector4d& x) const;
        Matrix6d hessian(const Vector2d& d, const Vector4d& x) const;
    };
}
