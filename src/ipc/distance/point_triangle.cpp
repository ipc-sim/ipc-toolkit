#include "point_triangle.hpp"

#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_plane.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/autodiff_types.hpp>

#include <stdexcept> // std::invalid_argument

namespace ipc {

template <typename T>
T point_triangle_distance(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2,
    PointTriangleDistanceType dtype)
{
    if constexpr (std::is_same_v<T, double> || std::is_same_v<T, float>) {
        if (dtype == PointTriangleDistanceType::AUTO) {
            dtype = point_triangle_distance_type(p, t0, t1, t2);
        }
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

template <typename T>
Eigen::Vector<T, 12> point_triangle_distance_gradient(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2,
    PointTriangleDistanceType dtype)
{
    if (dtype == PointTriangleDistanceType::AUTO) {
        dtype = point_triangle_distance_type(p, t0, t1, t2);
    }

    Eigen::Vector<T, 12> grad = Eigen::Vector<T, 12>::Zero();

    switch (dtype) {
    case PointTriangleDistanceType::P_T0:
        grad.template head<6>() = point_point_distance_gradient(p, t0);
        break;

    case PointTriangleDistanceType::P_T1: {
        const Eigen::Vector<T, 6> local_grad =
            point_point_distance_gradient(p, t1);
        grad.template head<3>() = local_grad.template head<3>();
        grad.template segment<3>(6) = local_grad.template tail<3>();
        break;
    }

    case PointTriangleDistanceType::P_T2: {
        const Eigen::Vector<T, 6> local_grad =
            point_point_distance_gradient(p, t2);
        grad.template head<3>() = local_grad.template head<3>();
        grad.template tail<3>() = local_grad.template tail<3>();
        break;
    }

    case PointTriangleDistanceType::P_E0:
        grad.template head<9>() = point_line_distance_gradient(p, t0, t1);
        break;

    case PointTriangleDistanceType::P_E1: {
        const Eigen::Vector<T, 9> local_grad =
            point_line_distance_gradient(p, t1, t2);
        grad.template head<3>() = local_grad.template head<3>();
        grad.template tail<6>() = local_grad.template tail<6>();
        break;
    }

    case PointTriangleDistanceType::P_E2: {
        const Eigen::Vector<T, 9> local_grad =
            point_line_distance_gradient(p, t2, t0);
        grad.template head<3>() = local_grad.template head<3>();     // ∇_p
        grad.template segment<3>(3) = local_grad.template tail<3>(); // ∇_{t0}
        grad.template tail<3>() = local_grad.template segment<3>(3); // ∇_{t2}
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

template <typename T>
Eigen::Matrix<T, 12, 12> point_triangle_distance_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2,
    PointTriangleDistanceType dtype)
{
    if (dtype == PointTriangleDistanceType::AUTO) {
        dtype = point_triangle_distance_type(p, t0, t1, t2);
    }

    Eigen::Matrix<T, 12, 12> hess = Eigen::Matrix<T, 12, 12>::Zero();

    switch (dtype) {
    case PointTriangleDistanceType::P_T0:
        hess.template topLeftCorner<6, 6>() =
            point_point_distance_hessian(p, t0);
        break;

    case PointTriangleDistanceType::P_T1: {
        const Eigen::Matrix<T, 6, 6> local_hess =
            point_point_distance_hessian(p, t1);
        hess.template topLeftCorner<3, 3>() =
            local_hess.template topLeftCorner<3, 3>();
        hess.template block<3, 3>(0, 6) =
            local_hess.template topRightCorner<3, 3>();
        hess.template block<3, 3>(6, 0) =
            local_hess.template bottomLeftCorner<3, 3>();
        hess.template block<3, 3>(6, 6) =
            local_hess.template bottomRightCorner<3, 3>();
        break;
    }

    case PointTriangleDistanceType::P_T2: {
        const Eigen::Matrix<T, 6, 6> local_hess =
            point_point_distance_hessian(p, t2);
        hess.template topLeftCorner<3, 3>() =
            local_hess.template topLeftCorner<3, 3>();
        hess.template topRightCorner<3, 3>() =
            local_hess.template topRightCorner<3, 3>();
        hess.template bottomLeftCorner<3, 3>() =
            local_hess.template bottomLeftCorner<3, 3>();
        hess.template bottomRightCorner<3, 3>() =
            local_hess.template bottomRightCorner<3, 3>();
        break;
    }

    case PointTriangleDistanceType::P_E0:
        hess.template topLeftCorner<9, 9>() =
            point_line_distance_hessian(p, t0, t1);
        break;

    case PointTriangleDistanceType::P_E1: {
        const Eigen::Matrix<T, 9, 9> local_hess =
            point_line_distance_hessian(p, t1, t2);
        hess.template topLeftCorner<3, 3>() =
            local_hess.template topLeftCorner<3, 3>();
        hess.template topRightCorner<3, 6>() =
            local_hess.template topRightCorner<3, 6>();
        hess.template bottomLeftCorner<6, 3>() =
            local_hess.template bottomLeftCorner<6, 3>();
        hess.template bottomRightCorner<6, 6>() =
            local_hess.template bottomRightCorner<6, 6>();
        break;
    }

    case PointTriangleDistanceType::P_E2: {
        const Eigen::Matrix<T, 9, 9> local_hess =
            point_line_distance_hessian(p, t2, t0);
        hess.template topLeftCorner<3, 3>() =
            local_hess.template topLeftCorner<3, 3>();
        hess.template block<3, 3>(0, 3) =
            local_hess.template topRightCorner<3, 3>();
        hess.template topRightCorner<3, 3>() =
            local_hess.template block<3, 3>(0, 3);
        hess.template block<3, 3>(3, 0) =
            local_hess.template bottomLeftCorner<3, 3>();
        hess.template block<3, 3>(3, 3) =
            local_hess.template bottomRightCorner<3, 3>();
        hess.template block<3, 3>(3, 9) = local_hess.template block<3, 3>(6, 3);
        hess.template bottomLeftCorner<3, 3>() =
            local_hess.template block<3, 3>(3, 0);
        hess.template block<3, 3>(9, 3) = local_hess.template block<3, 3>(3, 6);
        hess.template bottomRightCorner<3, 3>() =
            local_hess.template block<3, 3>(3, 3);
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

// clang-format off
template float point_triangle_distance<float>(Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, PointTriangleDistanceType);
template double point_triangle_distance<double>(Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, PointTriangleDistanceType);
template ADGrad<12> point_triangle_distance<ADGrad<12>>(Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>>, Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>>, Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>>, Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>>, PointTriangleDistanceType);
template ADHessian<12> point_triangle_distance<ADHessian<12>>(Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>>, Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>>, Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>>, Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>>, PointTriangleDistanceType);
template ADGrad<13> point_triangle_distance<ADGrad<13>>(Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>>, Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>>, Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>>, Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>>, PointTriangleDistanceType);
template ADHessian<13> point_triangle_distance<ADHessian<13>>(Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>>, Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>>, Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>>, Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>>, PointTriangleDistanceType);
template Vector12f point_triangle_distance_gradient<float>(Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, PointTriangleDistanceType);
template Vector12d point_triangle_distance_gradient<double>(Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, PointTriangleDistanceType);
template Matrix12f point_triangle_distance_hessian<float>(Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, PointTriangleDistanceType);
template Matrix12d point_triangle_distance_hessian<double>(Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, PointTriangleDistanceType);
// clang-format on

} // namespace ipc
