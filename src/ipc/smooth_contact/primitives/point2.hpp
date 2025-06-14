#pragma once

#include "primitive.hpp"

#include <ipc/smooth_contact/distance/mollifier.hpp>

namespace ipc {
class Point2 : public Primitive {
public:
    constexpr static int n_core_points = 1;
    constexpr static int dim = 2;
    constexpr static int max_size = n_vert_neighbors_2d * dim;
    // d is a vector from this point to the other primitive
    Point2(
        const long& id,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const VectorMax3d& d,
        const ParameterType& param);

    Point2(
        const long& id,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices);

    int n_vertices() const override;
    int n_dofs() const override { return n_vertices() * dim; }

    // assume the following functions are only called if active
    double potential(
        const Vector<double, dim>& d,
        const Vector<double, -1, max_size>& x) const;
    // derivatives including wrt. d (the closest direction) in front
    Vector<double, -1, max_size + dim> grad(
        const Vector<double, dim>& d,
        const Vector<double, -1, max_size>& x) const;
    MatrixMax<double, max_size + dim, max_size + dim> hessian(
        const Vector<double, dim>& d,
        const Vector<double, -1, max_size>& x) const;

private:
    bool has_neighbor_1, has_neighbor_2;
    bool orientable;
};
} // namespace ipc
