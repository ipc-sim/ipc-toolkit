#include "ipc/geometry/normal.hpp"

namespace ipc {

std::tuple<VectorMax3d, MatrixMax3d>
normalization_and_jacobian(Eigen::ConstRef<VectorMax3d> x)
{
    const double norm = x.norm();
    const VectorMax3d xhat = x / norm;
    const auto I = MatrixMax3d::Identity(x.size(), x.size());
    const MatrixMax3d J = (I - (xhat * xhat.transpose())) / norm;
    return std::make_tuple(xhat, J);
}

std::tuple<VectorMax3d, MatrixMax3d, std::array<MatrixMax3d, 3>>
normalization_and_jacobian_and_hessian(Eigen::ConstRef<VectorMax3d> x)
{
    // const auto [xhat, J] = normalization_and_jacobian(x);
    const double norm = x.norm();
    const VectorMax3d xhat = x / norm;
    const auto I = MatrixMax3d::Identity(x.size(), x.size());
    const MatrixMax3d J = (I - (xhat * xhat.transpose())) / norm;

    std::array<MatrixMax3d, 3> H;
    for (int i = 0; i < x.size(); i++) {
        H[i] = (xhat(i) * J + xhat * J.row(i) + J.col(i) * xhat.transpose())
            / (-norm);
    }

    return std::make_tuple(xhat, J, H);
}

// --- edge-vertex normal functions -------------------------------------------

VectorMax3d edge_vertex_unnormalized_normal(
    Eigen::ConstRef<VectorMax3d> v,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1)
{
    assert(v.size() == e0.size() && v.size() == e1.size());
    assert(v.size() == 2 || v.size() == 3);
    if (v.size() == 2) {
        // In 2D, the normal is simply the perpendicular vector to the edge
        const Eigen::Vector2d e = e1.tail<2>() - e0.tail<2>();
        return Eigen::Vector2d(-e.y(), e.x());
    } else {
        // Use triple product expansion of the cross product -e × (e × d)
        // (https://en.wikipedia.org/wiki/Cross_product#Triple_product_expansion)
        // NOTE: This would work in 2D as well, but we handle that case above.
        const Eigen::Vector3d e = e1 - e0;
        const Eigen::Vector3d d = v - e0;
        return d * e.dot(e) - e * e.dot(d);
    }
}

MatrixMax<double, 3, 9> edge_vertex_unnormalized_normal_jacobian(
    Eigen::ConstRef<VectorMax3d> v,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1)
{
    assert(v.size() == e0.size() && v.size() == e1.size());
    assert(v.size() == 2 || v.size() == 3);

    MatrixMax<double, 3, 9> dn(v.size(), 3 * v.size());

    if (v.size() == 2) {
        // In 2D, the normal is simply the perpendicular vector to the edge
        dn.leftCols<2>().setZero();
        dn.middleCols<2>(2) << 0, 1, -1, 0;
        dn.rightCols<2>() << 0, -1, 1, 0;
        return dn;
    } else {
        const Eigen::Vector3d e = e1 - e0;
        const Eigen::Vector3d d = v - e0;

        const auto I = Eigen::Matrix3d::Identity();

        // ∂n/∂v
        dn.leftCols<3>() = e.dot(e) * I - e * e.transpose();
        // ∂n/∂e1
        dn.rightCols<3>() =
            -e.dot(d) * I - e * d.transpose() + (2 * d) * e.transpose();
        // ∂n/∂e0
        dn.middleCols<3>(3) = -dn.leftCols<3>() - dn.rightCols<3>();
    }

    return dn;
}

} // namespace ipc