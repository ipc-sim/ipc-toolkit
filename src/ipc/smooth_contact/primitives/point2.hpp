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
        const double& alpha,
        const double& beta);

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
};

/// @brief
/// @param v
/// @param direc points from v
/// @param e0
/// @param e1
/// @param alpha
/// @param beta
/// @return
template <class scalar>
scalar smooth_point2_term(
    const Eigen::Ref<const Vector2<scalar>>& v,
    const Eigen::Ref<const Vector2<scalar>>& direc,
    const Eigen::Ref<const Vector2<scalar>>& e0,
    const Eigen::Ref<const Vector2<scalar>>& e1,
    const double& alpha,
    const double& beta);

/// @brief
/// @param v
/// @param direc points from v
/// @param e0
/// @param e1
/// @param alpha
/// @param beta
/// @return
bool smooth_point2_term_type(
    const Eigen::Ref<const Vector2<double>>& v,
    const Eigen::Ref<const Vector2<double>>& direc,
    const Eigen::Ref<const Vector2<double>>& e0,
    const Eigen::Ref<const Vector2<double>>& e1,
    const double& alpha,
    const double& beta);

} // namespace ipc
