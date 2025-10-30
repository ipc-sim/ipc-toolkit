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

// --- triangle normal functions ----------------------------------------------

namespace {
    void set_cross_product_matrix_jacobian(
        Eigen::Ref<Eigen::Matrix<double, 9, 3>> Jx, double chain_rule = 1.0)
    {
        Jx(2, 1) = Jx(3, 2) = Jx(7, 0) = -chain_rule;
        Jx(1, 2) = Jx(5, 0) = Jx(6, 1) = chain_rule;
    }
} // namespace

Eigen::Matrix<double, 3, 9> cross_product_matrix_jacobian()
{
    Eigen::Matrix<double, 9, 3> J = Eigen::Matrix<double, 9, 3>::Zero();
    set_cross_product_matrix_jacobian(J);
    return J.reshaped(3, 9);
}

Eigen::Matrix<double, 3, 81> triangle_unnormalized_normal_hessian(
    Eigen::ConstRef<Eigen::Vector3d> a,
    Eigen::ConstRef<Eigen::Vector3d> b,
    Eigen::ConstRef<Eigen::Vector3d> c)
{
    Eigen::Matrix<double, 27, 9> H = Eigen::Matrix<double, 27, 9>::Zero();

    // ∂²n/∂a² = 0
    set_cross_product_matrix_jacobian(H.block<9, 3>(0, 3), -1.0); // ∂²n/∂a∂b
    set_cross_product_matrix_jacobian(H.block<9, 3>(0, 6), 1.0);  // ∂²n/∂a∂c
    /**/
    set_cross_product_matrix_jacobian(H.block<9, 3>(9, 0), 1.0); // ∂²n/∂b∂a
    // ∂²n/∂b² = 0
    set_cross_product_matrix_jacobian(H.block<9, 3>(9, 6), -1.0); // ∂²n/∂b∂c
    /**/
    set_cross_product_matrix_jacobian(H.block<9, 3>(18, 0), -1.0); // ∂²n/∂c∂a
    set_cross_product_matrix_jacobian(H.block<9, 3>(18, 3), 1.0);  // ∂²n/∂c∂b
    // ∂²n/∂c² = 0

    return H.reshaped(3, 81);
}

Eigen::Matrix<double, 3, 81> triangle_normal_hessian(
    Eigen::ConstRef<Eigen::Vector3d> a,
    Eigen::ConstRef<Eigen::Vector3d> b,
    Eigen::ConstRef<Eigen::Vector3d> c)
{
    const Eigen::Vector3d z = triangle_unnormalized_normal(a, b, c);
    const double z_norm2 = z.squaredNorm();
    const double z_norm = std::sqrt(z_norm2);
    const double z_norm3 = z_norm2 * z_norm;

    const auto dz_dx = triangle_unnormalized_normal_jacobian(a, b, c);
    const auto d2z_dx2 = triangle_unnormalized_normal_hessian(a, b, c);

    Eigen::Matrix<double, 3, 81> d2n_dx2;
    for (int k = 0; k < 81; ++k) {
        const int i = k / 9, j = k % 9;

        const double alpha = z.dot(dz_dx.col(i)) / z_norm3;

        const double dalpha_dxj =
            (dz_dx.col(j).dot(dz_dx.col(i)) + z.dot(d2z_dx2.col(k))
             - 3 * z.dot(dz_dx.col(i)) * dz_dx.col(j).dot(z) / z_norm2)
            / z_norm3;

        d2n_dx2.col(k) = -dz_dx.col(i) * z.transpose() * dz_dx.col(j) / z_norm3
            + d2z_dx2.col(k) / z_norm - alpha * dz_dx.col(j) - dalpha_dxj * z;
    }

    return d2n_dx2;
}

// --- edge-edge normal functions ---------------------------------------------

Eigen::Matrix<double, 3, 144> edge_edge_unnormalized_normal_hessian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    Eigen::Matrix<double, 36, 12> H = Eigen::Matrix<double, 36, 12>::Zero();

    // ∂²n/∂a² = 0
    // ∂²n/∂a∂b = 0
    set_cross_product_matrix_jacobian(H.block<9, 3>(0, 6), -1.0); // ∂²n/∂a∂c
    set_cross_product_matrix_jacobian(H.block<9, 3>(0, 9), 1.0);  // ∂²n/∂a∂d
    /**/
    // ∂²n/∂b∂a = 0
    // ∂²n/∂b² = 0
    set_cross_product_matrix_jacobian(H.block<9, 3>(9, 6), 1.0);  // ∂²n/∂b∂c
    set_cross_product_matrix_jacobian(H.block<9, 3>(9, 9), -1.0); // ∂²n/∂b∂d
    /**/
    set_cross_product_matrix_jacobian(H.block<9, 3>(18, 0), 1.0);  // ∂²n/∂c∂a
    set_cross_product_matrix_jacobian(H.block<9, 3>(18, 3), -1.0); // ∂²n/∂c∂b
    // ∂²n/∂c² = 0
    // ∂²n/∂d² = 0
    /**/
    set_cross_product_matrix_jacobian(H.block<9, 3>(27, 0), -1.0); // ∂²n/∂d∂a
    set_cross_product_matrix_jacobian(H.block<9, 3>(27, 3), 1.0);  // ∂²n/∂d∂b
    // ∂²n/∂d∂c = 0
    // ∂²n/∂d² = 0

    return H.reshaped(3, 144);
}

Eigen::Matrix<double, 3, 144> edge_edge_normal_hessian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    const Eigen::Vector3d z = edge_edge_unnormalized_normal(ea0, ea1, eb0, eb1);
    const double z_norm2 = z.squaredNorm();
    const double z_norm = std::sqrt(z_norm2);
    const double z_norm3 = z_norm2 * z_norm;

    const auto dz_dx =
        edge_edge_unnormalized_normal_jacobian(ea0, ea1, eb0, eb1);
    const auto d2z_dx2 =
        edge_edge_unnormalized_normal_hessian(ea0, ea1, eb0, eb1);

    Eigen::Matrix<double, 3, 144> d2n_dx2;
    for (int k = 0; k < 144; ++k) {
        const int i = k / 12, j = k % 12;

        const double alpha = z.dot(dz_dx.col(i)) / z_norm3;

        const double dalpha_dxj =
            (dz_dx.col(j).dot(dz_dx.col(i)) + z.dot(d2z_dx2.col(k))
             - 3 * z.dot(dz_dx.col(i)) * dz_dx.col(j).dot(z) / z_norm2)
            / z_norm3;

        d2n_dx2.col(k) = -dz_dx.col(i) * z.transpose() * dz_dx.col(j) / z_norm3
            + d2z_dx2.col(k) / z_norm - alpha * dz_dx.col(j) - dalpha_dxj * z;
    }

    return d2n_dx2;
}

} // namespace ipc