#pragma once

#include "primitive.hpp"

#include <ipc/smooth_contact/distance/mollifier.hpp>

namespace ipc {
class Point2 : public Primitive {
public:
    static constexpr int N_CORE_POINTS = 1;
    static constexpr int DIM = 2;
    static constexpr int MAX_SIZE = N_VERT_NEIGHBORS_2D * DIM;
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
    int n_dofs() const override { return n_vertices() * DIM; }

    // assume the following functions are only called if active
    double potential(
        const Vector<double, DIM>& d,
        const Vector<double, -1, MAX_SIZE>& x) const;
    // derivatives including wrt. d (the closest direction) in front
    Vector<double, -1, MAX_SIZE + DIM> grad(
        const Vector<double, DIM>& d,
        const Vector<double, -1, MAX_SIZE>& x) const;
    MatrixMax<double, MAX_SIZE + DIM, MAX_SIZE + DIM> hessian(
        const Vector<double, DIM>& d,
        const Vector<double, -1, MAX_SIZE>& x) const;

private:
    bool has_neighbor_1, has_neighbor_2;
    bool orientable;
};
} // namespace ipc
