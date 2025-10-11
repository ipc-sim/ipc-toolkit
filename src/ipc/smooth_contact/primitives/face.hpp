#pragma once

#include "primitive.hpp"

#include <ipc/smooth_contact/distance/mollifier.hpp>

namespace ipc {
class Face : public Primitive {
public:
    static constexpr int N_CORE_POINTS = 3;
    static constexpr int DIM = 3;
    // d is a vector from closest point on the face to the point outside of the
    // face
    Face(
        const long& id,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const VectorMax3d& d,
        const ParameterType& param);

    int n_vertices() const override;
    int n_dofs() const override { return n_vertices() * DIM; }

    double potential(const Vector3d& d, const Vector9d& x) const;
    Vector12d grad(const Vector3d& d, const Vector9d& x) const;
    Matrix12d hessian(const Vector3d& d, const Vector9d& x) const;
};

/// @brief d points from triangle to the point
template <typename scalar>
scalar smooth_face_term(
    Eigen::ConstRef<Vector3<scalar>> v0,
    Eigen::ConstRef<Vector3<scalar>> v1,
    Eigen::ConstRef<Vector3<scalar>> v2)
{
    return 0.5 * (v1 - v0).cross(v2 - v0).norm(); // area of triangle
}
} // namespace ipc