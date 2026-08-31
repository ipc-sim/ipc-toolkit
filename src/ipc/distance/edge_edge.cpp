#include "edge_edge.hpp"

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/line_line.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/autodiff_types.hpp>

namespace ipc::detail {

template <typename T>
T edge_edge_distance(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1,
    EdgeEdgeDistanceType dtype)
{
    if constexpr (std::is_same_v<T, double> || std::is_same_v<T, float>) {
        if (dtype == EdgeEdgeDistanceType::AUTO) {
            dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
        }
    }

    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0:
        return point_point_distance(ea0, eb0);

    case EdgeEdgeDistanceType::EA0_EB1:
        return point_point_distance(ea0, eb1);

    case EdgeEdgeDistanceType::EA1_EB0:
        return point_point_distance(ea1, eb0);

    case EdgeEdgeDistanceType::EA1_EB1:
        return point_point_distance(ea1, eb1);

    case EdgeEdgeDistanceType::EA_EB0:
        return point_line_distance(eb0, ea0, ea1);

    case EdgeEdgeDistanceType::EA_EB1:
        return point_line_distance(eb1, ea0, ea1);

    case EdgeEdgeDistanceType::EA0_EB:
        return point_line_distance(ea0, eb0, eb1);

    case EdgeEdgeDistanceType::EA1_EB:
        return point_line_distance(ea1, eb0, eb1);

    case EdgeEdgeDistanceType::EA_EB:
        return line_line_distance(ea0, ea1, eb0, eb1);

    default:
        throw_invalid_distance_type("edge_edge_distance");
    }
}

template <typename T>
Eigen::Vector<T, 12> edge_edge_distance_gradient(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1,
    EdgeEdgeDistanceType dtype)
{
    using Vector6T = Eigen::Vector<T, 6>;
    using Vector9T = Eigen::Vector<T, 9>;
    using Vector12T = Eigen::Vector<T, 12>;

    if (dtype == EdgeEdgeDistanceType::AUTO) {
        dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
    }

    Vector12T grad = Vector12T::Zero();

    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0: {
        const Vector6T local_grad = point_point_distance_gradient(ea0, eb0);
        grad.template head<3>() = local_grad.template head<3>();
        grad.template segment<3>(6) = local_grad.template tail<3>();
        break;
    }

    case EdgeEdgeDistanceType::EA0_EB1: {
        const Vector6T local_grad = point_point_distance_gradient(ea0, eb1);
        grad.template head<3>() = local_grad.template head<3>();
        grad.template tail<3>() = local_grad.template tail<3>();
        break;
    }

    case EdgeEdgeDistanceType::EA1_EB0:
        grad.template segment<6>(3) = point_point_distance_gradient(ea1, eb0);
        break;

    case EdgeEdgeDistanceType::EA1_EB1: {
        const Vector6T local_grad = point_point_distance_gradient(ea1, eb1);
        grad.template segment<3>(3) = local_grad.template head<3>();
        grad.template tail<3>() = local_grad.template tail<3>();
        break;
    }

    case EdgeEdgeDistanceType::EA_EB0: {
        const Vector9T local_grad = point_line_distance_gradient(eb0, ea0, ea1);
        grad.template head<6>() = local_grad.template tail<6>();
        grad.template segment<3>(6) = local_grad.template head<3>();
        break;
    }

    case EdgeEdgeDistanceType::EA_EB1: {
        const Vector9T local_grad = point_line_distance_gradient(eb1, ea0, ea1);
        grad.template head<6>() = local_grad.template tail<6>();
        grad.template tail<3>() = local_grad.template head<3>();
        break;
    }

    case EdgeEdgeDistanceType::EA0_EB: {
        const Vector9T local_grad = point_line_distance_gradient(ea0, eb0, eb1);
        grad.template head<3>() = local_grad.template head<3>();
        grad.template tail<6>() = local_grad.template tail<6>();
        break;
    }

    case EdgeEdgeDistanceType::EA1_EB:
        grad.template tail<9>() = point_line_distance_gradient(ea1, eb0, eb1);
        break;

    case EdgeEdgeDistanceType::EA_EB:
        grad = line_line_distance_gradient(ea0, ea1, eb0, eb1);
        break;

    default:
        throw_invalid_distance_type("edge_edge_distance_gradient");
    }

    return grad;
}

template <typename T>
Eigen::Matrix<T, 12, 12> edge_edge_distance_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1,
    EdgeEdgeDistanceType dtype)
{
    using Matrix6T = Eigen::Matrix<T, 6, 6>;
    using Matrix9T = Eigen::Matrix<T, 9, 9>;
    using Matrix12T = Eigen::Matrix<T, 12, 12>;

    if (dtype == EdgeEdgeDistanceType::AUTO) {
        dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
    }

    Matrix12T hess = Matrix12T::Zero();

    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0: {
        const Matrix6T local_hess = point_point_distance_hessian(ea0, eb0);
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

    case EdgeEdgeDistanceType::EA0_EB1: {
        const Matrix6T local_hess = point_point_distance_hessian(ea0, eb1);
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

    case EdgeEdgeDistanceType::EA1_EB0:
        hess.template block<6, 6>(3, 3) =
            point_point_distance_hessian(ea1, eb0);
        break;

    case EdgeEdgeDistanceType::EA1_EB1: {
        const Matrix6T local_hess = point_point_distance_hessian(ea1, eb1);
        hess.template block<3, 3>(3, 3) =
            local_hess.template topLeftCorner<3, 3>();
        hess.template block<3, 3>(3, 9) =
            local_hess.template topRightCorner<3, 3>();
        hess.template block<3, 3>(9, 3) =
            local_hess.template bottomLeftCorner<3, 3>();
        hess.template bottomRightCorner<3, 3>() =
            local_hess.template bottomRightCorner<3, 3>();
        break;
    }

    case EdgeEdgeDistanceType::EA_EB0: {
        const Matrix9T local_hess = point_line_distance_hessian(eb0, ea0, ea1);
        hess.template topLeftCorner<6, 6>() =
            local_hess.template bottomRightCorner<6, 6>();
        hess.template block<3, 6>(6, 0) =
            local_hess.template topRightCorner<3, 6>();
        hess.template block<6, 3>(0, 6) =
            local_hess.template bottomLeftCorner<6, 3>();
        hess.template block<3, 3>(6, 6) =
            local_hess.template topLeftCorner<3, 3>();
        break;
    }

    case EdgeEdgeDistanceType::EA_EB1: {
        const Matrix9T local_hess = point_line_distance_hessian(eb1, ea0, ea1);
        hess.template topLeftCorner<6, 6>() =
            local_hess.template bottomRightCorner<6, 6>();
        hess.template topRightCorner<6, 3>() =
            local_hess.template bottomLeftCorner<6, 3>();
        hess.template bottomLeftCorner<3, 6>() =
            local_hess.template topRightCorner<3, 6>();
        hess.template bottomRightCorner<3, 3>() =
            local_hess.template topLeftCorner<3, 3>();
        break;
    }

    case EdgeEdgeDistanceType::EA0_EB: {
        const Matrix9T local_hess = point_line_distance_hessian(ea0, eb0, eb1);
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

    case EdgeEdgeDistanceType::EA1_EB:
        hess.template bottomRightCorner<9, 9>() =
            point_line_distance_hessian(ea1, eb0, eb1);
        break;

    case EdgeEdgeDistanceType::EA_EB:
        hess = line_line_distance_hessian(ea0, ea1, eb0, eb1);
        break;

    default:
        throw_invalid_distance_type("edge_edge_distance_hessian");
    }

    return hess;
}

#define IPC_INSTANTIATE_EDGE_EDGE(T)                                           \
    template T edge_edge_distance<T>(                                          \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>, EdgeEdgeDistanceType);             \
    template Eigen::Vector<T, 12> edge_edge_distance_gradient<T>(              \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>, EdgeEdgeDistanceType);             \
    template Eigen::Matrix<T, 12, 12> edge_edge_distance_hessian<T>(           \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>, EdgeEdgeDistanceType)

// Autodiff scalars: value only. Gradients/Hessians are intentionally not
// supported for them.
#define IPC_INSTANTIATE_EDGE_EDGE_VALUE(T)                                     \
    template T edge_edge_distance<T>(                                          \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>, EdgeEdgeDistanceType)

IPC_INSTANTIATE_EDGE_EDGE(float);
IPC_INSTANTIATE_EDGE_EDGE(double);
IPC_INSTANTIATE_EDGE_EDGE_VALUE(ADGrad<9>);
IPC_INSTANTIATE_EDGE_EDGE_VALUE(ADHessian<9>);
IPC_INSTANTIATE_EDGE_EDGE_VALUE(ADGrad<12>);
IPC_INSTANTIATE_EDGE_EDGE_VALUE(ADHessian<12>);
IPC_INSTANTIATE_EDGE_EDGE_VALUE(ADGrad<13>);
IPC_INSTANTIATE_EDGE_EDGE_VALUE(ADHessian<13>);

#undef IPC_INSTANTIATE_EDGE_EDGE
#undef IPC_INSTANTIATE_EDGE_EDGE_VALUE

} // namespace ipc::detail
