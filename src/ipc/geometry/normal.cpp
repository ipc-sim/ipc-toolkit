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

// --- point-line normal functions -------------------------------------------

VectorMax3d point_line_unnormalized_normal(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1)
{
    assert(p.size() == e0.size() && p.size() == e1.size());
    assert(p.size() == 2 || p.size() == 3);
    if (p.size() == 2) {
        // In 2D, the normal is simply the perpendicular vector to the line
        const Eigen::Vector2d e = e1.head<2>() - e0.head<2>();
        return Eigen::Vector2d(-e.y(), e.x());
    } else {
        // Use triple product expansion of the cross product -e × (e × d)
        // (https://en.wikipedia.org/wiki/Cross_product#Triple_product_expansion)
        // NOTE: This would work in 2D as well, but we handle that case above.
        const Eigen::Vector3d e = e1 - e0;
        const Eigen::Vector3d d = p - e0;
        return d * e.dot(e) - e * e.dot(d);
    }
}

MatrixMax<double, 3, 9> point_line_unnormalized_normal_jacobian(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1)
{
    assert(p.size() == e0.size() && p.size() == e1.size());
    assert(p.size() == 2 || p.size() == 3);

    MatrixMax<double, 3, 9> dn(p.size(), 3 * p.size());

    if (p.size() == 2) {
        // In 2D, the normal is simply the perpendicular vector to the line
        dn.leftCols<2>().setZero();
        dn.middleCols<2>(2) << 0, 1, -1, 0;
        dn.rightCols<2>() << 0, -1, 1, 0;
        return dn;
    } else {
        const Eigen::Vector3d e = e1 - e0;
        const Eigen::Vector3d d = p - e0;

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

MatrixMax<double, 27, 9> point_line_unnormalized_normal_hessian(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1)
{
    assert(p.size() == e0.size() && p.size() == e1.size());
    assert(p.size() == 2 || p.size() == 3);
    if (p.size() == 2) {
        // In 2D, the normal is simply the perpendicular vector to the line
        return MatrixMax<double, 27, 9>::Zero(12, 6);
    } else {
        const Eigen::Vector3d e = e1 - e0;
        const Eigen::Vector3d d = p - e0;
        const auto I = Eigen::Matrix3d::Identity();

        // The full Hessian is 27x9 (3 components of n, each with 9 variables)
        MatrixMax<double, 27, 9> hess = MatrixMax<double, 27, 9>::Zero(27, 9);

        // Compute the Hessian for each component k of the normal vector n
        for (int k = 0; k < 3; ++k) {
            // Standard basis vector for component k
            Eigen::Vector3d ek = Eigen::Vector3d::Zero();
            ek(k) = 1.0;

            // The full 9x9 Hessian for component n[k]
            // Order of variables: p (0-2), e0 (3-5), e1 (6-8)
            Matrix9d Hk = Matrix9d::Zero();

            // -------------------------------------------------------------
            // Block (p, e1): Derivative of ∂n/∂p w.r.t e1
            // -------------------------------------------------------------
            // ∂n/∂p = ||e||^2 * I - e * e^T
            // Differentiating row k w.r.t e:
            // H_pe1 = 2 * ek * e^T - e * ek^T - e[k] * I
            Eigen::Matrix3d H_pe1 =
                2.0 * ek * e.transpose() - e * ek.transpose() - e(k) * I;

            // -------------------------------------------------------------
            // Block (e1, e1): Derivative of ∂n/∂e1 w.r.t e1
            // -------------------------------------------------------------
            // ∂n/∂e1 = -(e.d)I - e*d^T + 2d*e^T
            // Differentiating row k w.r.t e (keeping d constant):
            // H_e1e1 = -ek * d^T - d * ek^T + 2 * d[k] * I
            Eigen::Matrix3d H_e1e1 =
                -ek * d.transpose() - d * ek.transpose() + 2.0 * d(k) * I;

            // -------------------------------------------------------------
            // Block (e1, e0): Derivative of ∂n/∂e1 w.r.t e0
            // -------------------------------------------------------------
            // Chain rule: ∂/∂e0 = ∂/∂e * (-1) + ∂/∂d * (-1)

            // Part 1: Contribution from e (-H_e1e1)

            // Part 2: Contribution from d
            // Differentiating row k w.r.t d (keeping e constant):
            // M = -ek * e^T - e[k] * I + 2 * e * ek^T
            Eigen::Matrix3d term_d =
                -ek * e.transpose() - e(k) * I + 2.0 * e * ek.transpose();

            // Combine: -1 * H_e1e1 - 1 * term_d
            Eigen::Matrix3d H_e1e0 = -H_e1e1 - term_d;

            // -------------------------------------------------------------
            // Assemble 9x9 Matrix (Symmetry & Translation Invariance)
            // -------------------------------------------------------------

            // (p, p) is Zero because ∂n/∂p is linear in p

            // (p, e1) and (e1, p)
            Hk.block<3, 3>(0, 6) = H_pe1;
            Hk.block<3, 3>(6, 0) = H_pe1.transpose();

            // (p, e0) = - (p, e1)
            // Because ∂e/∂e0 = -∂e/∂e1, and ∂n/∂p depends only on e
            Hk.block<3, 3>(0, 3) = -H_pe1;
            Hk.block<3, 3>(3, 0) = -H_pe1.transpose();

            // (e1, e1)
            Hk.block<3, 3>(6, 6) = H_e1e1;

            // (e1, e0) and (e0, e1)
            Hk.block<3, 3>(6, 3) = H_e1e0;
            Hk.block<3, 3>(3, 6) = H_e1e0.transpose();

            // (e0, e0)
            // Translation invariance: Sum of rows (and cols) must be zero.
            // H_e0e0 = -H_pe0 - H_e1e0
            //        = -(-H_pe1) - H_e1e0
            //        = H_pe1 - H_e1e0
            Hk.block<3, 3>(3, 3) = H_pe1 - H_e1e0;

            // Assign Hₖ to the strided rows of the full Hessian
            for (int i = 0; i < 9; ++i) {
                hess.row(3 * i + k) = Hk.row(i);
            }
        }

        return hess;
    }
}

MatrixMax<double, 27, 9> point_line_normal_hessian(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1)
{
    assert(p.size() == e0.size() && p.size() == e1.size());
    assert(p.size() == 2 || p.size() == 3);

    const VectorMax3d z = point_line_unnormalized_normal(p, e0, e1);
    const double z_norm2 = z.squaredNorm();
    const double z_norm = std::sqrt(z_norm2);
    const double z_norm3 = z_norm2 * z_norm;

    const int DIM = z.size(); // dimension (2 or 3)
    const int DOF = 3 * DIM;  // total dof (6 or 9)

    const auto dz_dx = point_line_unnormalized_normal_jacobian(p, e0, e1);
    const auto d2z_dx2 = point_line_unnormalized_normal_hessian(p, e0, e1);

    MatrixMax<double, 27, 9> d2n_dx2(DOF * DIM, DOF);
    for (int k = 0; k < DOF * DOF; ++k) {
        const int i = k / DOF, j = k % DOF;

        const double alpha = z.dot(dz_dx.col(i)) / z_norm3;

        const VectorMax3d d2z_dxidxj = d2z_dx2.block(DIM * i, j, DIM, 1);

        const double dalpha_dxj =
            (dz_dx.col(j).dot(dz_dx.col(i)) + z.dot(d2z_dxidxj)
             - 3 * z.dot(dz_dx.col(i)) * dz_dx.col(j).dot(z) / z_norm2)
            / z_norm3;

        d2n_dx2.block(DIM * i, j, DIM, 1) =
            -dz_dx.col(i) * z.transpose() * dz_dx.col(j) / z_norm3
            + d2z_dxidxj / z_norm - alpha * dz_dx.col(j) - dalpha_dxj * z;
    }

    return d2n_dx2;
}

// --- triangle normal functions ----------------------------------------------

namespace {
    void set_cross_product_matrix_jacobian(
        Eigen::Ref<Eigen::Matrix<double, 9, 3>> Jx, double chain_rule = 1.0)
    {
        Jx(1, 2) = Jx(5, 0) = Jx(6, 1) = -chain_rule;
        Jx(2, 1) = Jx(3, 2) = Jx(7, 0) = chain_rule;
    }
} // namespace

Eigen::Matrix<double, 9, 3> cross_product_matrix_jacobian()
{
    Eigen::Matrix<double, 9, 3> J = Eigen::Matrix<double, 9, 3>::Zero();
    set_cross_product_matrix_jacobian(J);
    return J;
}

Eigen::Matrix<double, 27, 9> triangle_unnormalized_normal_hessian(
    Eigen::ConstRef<Eigen::Vector3d> a,
    Eigen::ConstRef<Eigen::Vector3d> b,
    Eigen::ConstRef<Eigen::Vector3d> c)
{
    Eigen::Matrix<double, 27, 9> H = Eigen::Matrix<double, 27, 9>::Zero();

    // ∂²n/∂a² = 0
    set_cross_product_matrix_jacobian(H.block<9, 3>(0, 3), 1.0);  // ∂²n/∂a∂b
    set_cross_product_matrix_jacobian(H.block<9, 3>(0, 6), -1.0); // ∂²n/∂a∂c
    /**/
    set_cross_product_matrix_jacobian(H.block<9, 3>(9, 0), -1.0); // ∂²n/∂b∂a
    // ∂²n/∂b² = 0
    set_cross_product_matrix_jacobian(H.block<9, 3>(9, 6), 1.0); // ∂²n/∂b∂c
    /**/
    set_cross_product_matrix_jacobian(H.block<9, 3>(18, 0), 1.0);  // ∂²n/∂c∂a
    set_cross_product_matrix_jacobian(H.block<9, 3>(18, 3), -1.0); // ∂²n/∂c∂b
    // ∂²n/∂c² = 0

    return H;
}

Eigen::Matrix<double, 27, 9> triangle_normal_hessian(
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

    Eigen::Matrix<double, 27, 9> d2n_dx2;
    for (int k = 0; k < 81; ++k) {
        const int i = k / 9, j = k % 9;

        const double alpha = z.dot(dz_dx.col(i)) / z_norm3;

        const auto d2z_dxidxj = d2z_dx2.block<3, 1>(3 * i, j);

        const double dalpha_dxj =
            (dz_dx.col(j).dot(dz_dx.col(i)) + z.dot(d2z_dxidxj)
             - 3 * z.dot(dz_dx.col(i)) * dz_dx.col(j).dot(z) / z_norm2)
            / z_norm3;

        d2n_dx2.block<3, 1>(3 * i, j) =
            -dz_dx.col(i) * z.transpose() * dz_dx.col(j) / z_norm3
            + d2z_dxidxj / z_norm - alpha * dz_dx.col(j) - dalpha_dxj * z;
    }

    return d2n_dx2;
}

// --- line-line normal functions ---------------------------------------------

Eigen::Matrix<double, 36, 12> line_line_unnormalized_normal_hessian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    Eigen::Matrix<double, 36, 12> H = Eigen::Matrix<double, 36, 12>::Zero();

    // ∂²n/∂a² = 0
    // ∂²n/∂a∂b = 0
    set_cross_product_matrix_jacobian(H.block<9, 3>(0, 6), 1.0);  // ∂²n/∂a∂c
    set_cross_product_matrix_jacobian(H.block<9, 3>(0, 9), -1.0); // ∂²n/∂a∂d
    /**/
    // ∂²n/∂b∂a = 0
    // ∂²n/∂b² = 0
    set_cross_product_matrix_jacobian(H.block<9, 3>(9, 6), -1.0); // ∂²n/∂b∂c
    set_cross_product_matrix_jacobian(H.block<9, 3>(9, 9), 1.0);  // ∂²n/∂b∂d
    /**/
    set_cross_product_matrix_jacobian(H.block<9, 3>(18, 0), -1.0); // ∂²n/∂c∂a
    set_cross_product_matrix_jacobian(H.block<9, 3>(18, 3), 1.0);  // ∂²n/∂c∂b
    // ∂²n/∂c² = 0
    // ∂²n/∂d² = 0
    /**/
    set_cross_product_matrix_jacobian(H.block<9, 3>(27, 0), 1.0);  // ∂²n/∂d∂a
    set_cross_product_matrix_jacobian(H.block<9, 3>(27, 3), -1.0); // ∂²n/∂d∂b
    // ∂²n/∂d∂c = 0
    // ∂²n/∂d² = 0

    return H;
}

Eigen::Matrix<double, 36, 12> line_line_normal_hessian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    const Eigen::Vector3d z = line_line_unnormalized_normal(ea0, ea1, eb0, eb1);
    const double z_norm2 = z.squaredNorm();
    const double z_norm = std::sqrt(z_norm2);
    const double z_norm3 = z_norm2 * z_norm;

    const Eigen::Matrix<double, 3, 12> dz_dx =
        line_line_unnormalized_normal_jacobian(ea0, ea1, eb0, eb1);
    const Eigen::Matrix<double, 36, 12> d2z_dx2 =
        line_line_unnormalized_normal_hessian(ea0, ea1, eb0, eb1);

    Eigen::Matrix<double, 36, 12> d2n_dx2;
    for (int k = 0; k < 144; ++k) {
        const int i = k / 12, j = k % 12;

        const double alpha = z.dot(dz_dx.col(i)) / z_norm3;

        const auto d2z_dxidxj = d2z_dx2.block<3, 1>(3 * i, j);

        const double dalpha_dxj =
            (dz_dx.col(j).dot(dz_dx.col(i)) + z.dot(d2z_dxidxj)
             - 3 * z.dot(dz_dx.col(i)) * dz_dx.col(j).dot(z) / z_norm2)
            / z_norm3;

        d2n_dx2.block<3, 1>(3 * i, j) =
            -dz_dx.col(i) * z.transpose() * dz_dx.col(j) / z_norm3
            + d2z_dxidxj / z_norm - alpha * dz_dx.col(j) - dalpha_dxj * z;
    }

    return d2n_dx2;
}

} // namespace ipc