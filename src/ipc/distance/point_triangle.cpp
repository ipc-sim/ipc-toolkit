#include "point_triangle.hpp"

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_plane.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/autodiff_types.hpp>

namespace ipc::detail {

template <typename T>
T point_triangle_distance(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2,
    PointTriangleDistanceType dtype)
{
    if constexpr (std::is_floating_point_v<T>) {
        if (dtype == PointTriangleDistanceType::AUTO) {
            dtype = point_triangle_distance_type(p, t0, t1, t2);
        }
    } else if (dtype == PointTriangleDistanceType::AUTO) {
        throw_auto_requires_explicit_dtype("point_triangle_distance");
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
        throw_invalid_distance_type("point_triangle_distance");
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
    if constexpr (std::is_floating_point_v<T>) {
        if (dtype == PointTriangleDistanceType::AUTO) {
            dtype = point_triangle_distance_type(p, t0, t1, t2);
        }
    } else if (dtype == PointTriangleDistanceType::AUTO) {
        throw_auto_requires_explicit_dtype("point_triangle_distance_gradient");
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
        throw_invalid_distance_type("point_triangle_distance_gradient");
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
    if constexpr (std::is_floating_point_v<T>) {
        if (dtype == PointTriangleDistanceType::AUTO) {
            dtype = point_triangle_distance_type(p, t0, t1, t2);
        }
    } else if (dtype == PointTriangleDistanceType::AUTO) {
        throw_auto_requires_explicit_dtype("point_triangle_distance_hessian");
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
        throw_invalid_distance_type("point_triangle_distance_hessian");
    }

    return hess;
}

#define IPC_INSTANTIATE_POINT_TRIANGLE(T)                                      \
    template T point_triangle_distance<T>(                                     \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>, PointTriangleDistanceType);        \
    template Eigen::Vector<T, 12> point_triangle_distance_gradient<T>(         \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>, PointTriangleDistanceType);        \
    template Eigen::Matrix<T, 12, 12> point_triangle_distance_hessian<T>(      \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>, PointTriangleDistanceType)

// Autodiff scalars: value only. Gradients/Hessians are intentionally not
// supported for them.
#define IPC_INSTANTIATE_POINT_TRIANGLE_VALUE(T)                                \
    template T point_triangle_distance<T>(                                     \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>, PointTriangleDistanceType)

IPC_INSTANTIATE_POINT_TRIANGLE(float);
IPC_INSTANTIATE_POINT_TRIANGLE(double);
IPC_INSTANTIATE_POINT_TRIANGLE_VALUE(ADGrad<12>);
IPC_INSTANTIATE_POINT_TRIANGLE_VALUE(ADHessian<12>);
IPC_INSTANTIATE_POINT_TRIANGLE_VALUE(ADGrad<13>);
IPC_INSTANTIATE_POINT_TRIANGLE_VALUE(ADHessian<13>);

#undef IPC_INSTANTIATE_POINT_TRIANGLE
#undef IPC_INSTANTIATE_POINT_TRIANGLE_VALUE

} // namespace ipc::detail
