#include "angle.hpp"

#include <ipc/geometry/normal.hpp>

namespace ipc {

double dihedral_angle(
    Eigen::ConstRef<Eigen::Vector3d> x0,
    Eigen::ConstRef<Eigen::Vector3d> x1,
    Eigen::ConstRef<Eigen::Vector3d> x2,
    Eigen::ConstRef<Eigen::Vector3d> x3)
{
    const Eigen::Vector3d n0 = triangle_normal(x0, x1, x2);
    const Eigen::Vector3d n1 = triangle_normal(x1, x0, x3);
    const Eigen::Vector3d e = (x1 - x0).normalized();

    const double sin_theta = n0.cross(n1).dot(e);
    const double cos_theta = n0.dot(n1);

    return std::atan2(sin_theta, cos_theta);
}

Eigen::Vector<double, 12> dihedral_angle_gradient(
    Eigen::ConstRef<Eigen::Vector3d> x0,
    Eigen::ConstRef<Eigen::Vector3d> x1,
    Eigen::ConstRef<Eigen::Vector3d> x2,
    Eigen::ConstRef<Eigen::Vector3d> x3)
{
    const Eigen::Vector3d n0 = triangle_normal(x0, x1, x2);
    const Eigen::Vector3d n1 = triangle_normal(x1, x0, x3);
    const Eigen::Vector3d e = (x1 - x0).normalized();

    // --- Normal gradients ---

    Eigen::Matrix<double, 3, 12> dn0_dx;
    dn0_dx.leftCols<9>() = triangle_normal_jacobian(x0, x1, x2);
    dn0_dx.rightCols<3>().setZero();

    Eigen::Matrix<double, 3, 12> dn1_dx;
    const std::array<int, 9> idx = { { 3, 4, 5, 0, 1, 2, 9, 10, 11 } };
    dn1_dx(Eigen::all, idx) = triangle_normal_jacobian(x1, x0, x3);
    dn1_dx.middleCols<3>(6).setZero();

    // --- Edge gradient ---

    Eigen::Matrix<double, 3, 12> de_dx;
    de_dx.middleCols<3>(3) = normalization_jacobian(x1 - x0);
    de_dx.leftCols<3>() = -de_dx.middleCols<3>(3);
    de_dx.rightCols<6>().setZero();

    // --- Angle gradient ---

    const Eigen::Vector<double, 12> dcos_dx =
        dn0_dx.transpose() * n1 + dn1_dx.transpose() * n0;

    const Eigen::Vector<double, 12> dsin_dx =
        (cross_product_matrix(n0) * dn1_dx - cross_product_matrix(n1) * dn0_dx)
            .transpose()
        * e;

    // --- Product rule ---

    const double sin_theta = n0.cross(n1).dot(e);
    const double cos_theta = n0.dot(n1);

    return dsin_dx * cos_theta - dcos_dx * sin_theta;
}

} // namespace ipc