#include "angle.hpp"

#include <ipc/geometry/normal.hpp>

namespace ipc {

double dihedral_angle(
    Eigen::ConstRef<Eigen::Vector3d> x0,
    Eigen::ConstRef<Eigen::Vector3d> x1,
    Eigen::ConstRef<Eigen::Vector3d> x2,
    Eigen::ConstRef<Eigen::Vector3d> x3)
{
    double cos_theta =
        triangle_normal(x0, x1, x2).dot(triangle_normal(x1, x0, x3));
    assert(std::abs(cos_theta) <= 1.0 + 1e-6); // Numerical tolerance
    return std::acos(std::clamp(cos_theta, -1.0, 1.0));
}

Eigen::Vector<double, 12> dihedral_angle_gradient(
    Eigen::ConstRef<Eigen::Vector3d> x0,
    Eigen::ConstRef<Eigen::Vector3d> x1,
    Eigen::ConstRef<Eigen::Vector3d> x2,
    Eigen::ConstRef<Eigen::Vector3d> x3)
{
    const Eigen::Vector3d n0 = triangle_normal(x0, x1, x2);
    const Eigen::Vector3d n1 = triangle_normal(x1, x0, x3);

    const Eigen::Matrix<double, 3, 9> dn0_dx =
        triangle_normal_jacobian(x0, x1, x2);
    const Eigen::Matrix<double, 3, 9> dn1_dx =
        triangle_normal_jacobian(x1, x0, x3);

    const Eigen::Matrix3d dn0_dx0 = dn0_dx.block<3, 3>(0, 0);
    const Eigen::Matrix3d dn0_dx1 = dn0_dx.block<3, 3>(0, 3);
    const Eigen::Matrix3d dn0_dx2 = dn0_dx.block<3, 3>(0, 6);
    const Eigen::Matrix3d dn1_dx0 = dn1_dx.block<3, 3>(0, 3);
    const Eigen::Matrix3d dn1_dx1 = dn1_dx.block<3, 3>(0, 0);
    const Eigen::Matrix3d dn1_dx3 = dn1_dx.block<3, 3>(0, 6);

    Eigen::Vector<double, 12> gradient;
    // Product rule
    gradient.segment<3>(0) =
        dn0_dx0.transpose() * n1 + dn1_dx0.transpose() * n0;
    gradient.segment<3>(3) =
        dn0_dx1.transpose() * n1 + dn1_dx1.transpose() * n0;
    gradient.segment<3>(6) = dn0_dx2.transpose() * n1;
    gradient.segment<3>(9) = dn1_dx3.transpose() * n0;

    // Chain rule through acos
    const double cos_theta = n0.dot(n1);
    gradient *= -1 / sqrt(1.0 - cos_theta * cos_theta);

    return gradient;
}

} // namespace ipc