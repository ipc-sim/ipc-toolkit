#include "point_triangle.hpp"

#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_plane.hpp>
#include <ipc/distance/point_point.hpp>

#include <stdexcept> // std::invalid_argument

namespace ipc {

double point_triangle_distance(
    const Eigen::Ref<const Eigen::Vector3d>& p,
    const Eigen::Ref<const Eigen::Vector3d>& t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2,
    PointTriangleDistanceType dtype)
{
    if (dtype == PointTriangleDistanceType::AUTO) {
        dtype = point_triangle_distance_type(p, t0, t1, t2);
    }

    switch (dtype) {
    case PointTriangleDistanceType::P_T0:
        return point_point_distance(p, t0);

    case PointTriangleDistanceType::P_T1:
        return point_point_distance(p, t1);

    case PointTriangleDistanceType::P_T2:
        return point_point_distance(p, t2);

    case PointTriangleDistanceType::P_E0:
        return point_line_distance(p, t0, t1);

    case PointTriangleDistanceType::P_E1:
        return point_line_distance(p, t1, t2);

    case PointTriangleDistanceType::P_E2:
        return point_line_distance(p, t2, t0);

    case PointTriangleDistanceType::P_T:
        return point_plane_distance(p, t0, t1, t2);

    default:
        throw std::invalid_argument(
            "Invalid distance type for point-triangle distance!");
    }
}

Vector12d point_triangle_distance_gradient(
    const Eigen::Ref<const Eigen::Vector3d>& p,
    const Eigen::Ref<const Eigen::Vector3d>& t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2,
    PointTriangleDistanceType dtype)
{
    if (dtype == PointTriangleDistanceType::AUTO) {
        dtype = point_triangle_distance_type(p, t0, t1, t2);
    }

    Vector12d grad = Vector12d::Zero();

    switch (dtype) {
    case PointTriangleDistanceType::P_T0:
        grad.head<6>() = point_point_distance_gradient(p, t0);
        break;

    case PointTriangleDistanceType::P_T1: {
        const Vector6d local_grad = point_point_distance_gradient(p, t1);
        grad.head<3>() = local_grad.head<3>();
        grad.segment<3>(6) = local_grad.tail<3>();
        break;
    }

    case PointTriangleDistanceType::P_T2: {
        const Vector6d local_grad = point_point_distance_gradient(p, t2);
        grad.head<3>() = local_grad.head<3>();
        grad.tail<3>() = local_grad.tail<3>();
        break;
    }

    case PointTriangleDistanceType::P_E0:
        grad.head<9>() = point_line_distance_gradient(p, t0, t1);
        break;

    case PointTriangleDistanceType::P_E1: {
        const Vector9d local_grad = point_line_distance_gradient(p, t1, t2);
        grad.head<3>() = local_grad.head<3>();
        grad.tail<6>() = local_grad.tail<6>();
        break;
    }

    case PointTriangleDistanceType::P_E2: {
        const Vector9d local_grad = point_line_distance_gradient(p, t2, t0);
        grad.head<3>() = local_grad.head<3>();     // ∇_p
        grad.segment<3>(3) = local_grad.tail<3>(); // ∇_{t0}
        grad.tail<3>() = local_grad.segment<3>(3); // ∇_{t2}
        break;
    }

    case PointTriangleDistanceType::P_T:
        grad = point_plane_distance_gradient(p, t0, t1, t2);
        break;

    default:
        throw std::invalid_argument(
            "Invalid distance type for point-triangle distance gradient!");
    }

    return grad;
}

Matrix12d point_triangle_distance_hessian(
    const Eigen::Ref<const Eigen::Vector3d>& p,
    const Eigen::Ref<const Eigen::Vector3d>& t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2,
    PointTriangleDistanceType dtype)
{
    if (dtype == PointTriangleDistanceType::AUTO) {
        dtype = point_triangle_distance_type(p, t0, t1, t2);
    }

    Matrix12d hess = Matrix12d::Zero();

    switch (dtype) {
    case PointTriangleDistanceType::P_T0:
        hess.topLeftCorner<6, 6>() = point_point_distance_hessian(p, t0);
        break;

    case PointTriangleDistanceType::P_T1: {
        const Matrix6d local_hess = point_point_distance_hessian(p, t1);
        hess.topLeftCorner<3, 3>() = local_hess.topLeftCorner<3, 3>();
        hess.block<3, 3>(0, 6) = local_hess.topRightCorner<3, 3>();
        hess.block<3, 3>(6, 0) = local_hess.bottomLeftCorner<3, 3>();
        hess.block<3, 3>(6, 6) = local_hess.bottomRightCorner<3, 3>();
        break;
    }

    case PointTriangleDistanceType::P_T2: {
        const Matrix6d local_hess = point_point_distance_hessian(p, t2);
        hess.topLeftCorner<3, 3>() = local_hess.topLeftCorner<3, 3>();
        hess.topRightCorner<3, 3>() = local_hess.topRightCorner<3, 3>();
        hess.bottomLeftCorner<3, 3>() = local_hess.bottomLeftCorner<3, 3>();
        hess.bottomRightCorner<3, 3>() = local_hess.bottomRightCorner<3, 3>();
        break;
    }

    case PointTriangleDistanceType::P_E0:
        hess.topLeftCorner<9, 9>() = point_line_distance_hessian(p, t0, t1);
        break;

    case PointTriangleDistanceType::P_E1: {
        const Matrix9d local_hess = point_line_distance_hessian(p, t1, t2);
        hess.topLeftCorner<3, 3>() = local_hess.topLeftCorner<3, 3>();
        hess.topRightCorner<3, 6>() = local_hess.topRightCorner<3, 6>();
        hess.bottomLeftCorner<6, 3>() = local_hess.bottomLeftCorner<6, 3>();
        hess.bottomRightCorner<6, 6>() = local_hess.bottomRightCorner<6, 6>();
        break;
    }

    case PointTriangleDistanceType::P_E2: {
        const Matrix9d local_hess = point_line_distance_hessian(p, t2, t0);
        hess.topLeftCorner<3, 3>() = local_hess.topLeftCorner<3, 3>();
        hess.block<3, 3>(0, 3) = local_hess.topRightCorner<3, 3>();
        hess.topRightCorner<3, 3>() = local_hess.block<3, 3>(0, 3);
        hess.block<3, 3>(3, 0) = local_hess.bottomLeftCorner<3, 3>();
        hess.block<3, 3>(3, 3) = local_hess.bottomRightCorner<3, 3>();
        hess.block<3, 3>(3, 9) = local_hess.block<3, 3>(6, 3);
        hess.bottomLeftCorner<3, 3>() = local_hess.block<3, 3>(3, 0);
        hess.block<3, 3>(9, 3) = local_hess.block<3, 3>(3, 6);
        hess.bottomRightCorner<3, 3>() = local_hess.block<3, 3>(3, 3);
        break;
    }

    case PointTriangleDistanceType::P_T:
        hess = point_plane_distance_hessian(p, t0, t1, t2);
        break;

    default:
        throw std::invalid_argument(
            "Invalid distance type for point-triangle distance hessian!");
    }

    return hess;
}

} // namespace ipc
