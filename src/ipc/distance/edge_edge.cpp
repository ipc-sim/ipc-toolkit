#include "edge_edge.hpp"

#include <ipc/distance/line_line.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/autodiff_types.hpp>

#include <stdexcept> // std::invalid_argument

namespace ipc {

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
        throw std::invalid_argument(
            "Invalid distance type for edge-edge distance!");
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
        throw std::invalid_argument(
            "Invalid distance type for edge-edge distance gradient!");
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
        throw std::invalid_argument(
            "Invalid distance type for edge-edge distance hessian!");
    }

    return hess;
}

// clang-format off
template float edge_edge_distance<float>(Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, EdgeEdgeDistanceType);
template double edge_edge_distance<double>(Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, EdgeEdgeDistanceType);
template ADGrad<9> edge_edge_distance<ADGrad<9>>(Eigen::ConstRef<Eigen::Vector3<ADGrad<9>>>, Eigen::ConstRef<Eigen::Vector3<ADGrad<9>>>, Eigen::ConstRef<Eigen::Vector3<ADGrad<9>>>, Eigen::ConstRef<Eigen::Vector3<ADGrad<9>>>, EdgeEdgeDistanceType);
template ADHessian<9> edge_edge_distance<ADHessian<9>>(Eigen::ConstRef<Eigen::Vector3<ADHessian<9>>>, Eigen::ConstRef<Eigen::Vector3<ADHessian<9>>>, Eigen::ConstRef<Eigen::Vector3<ADHessian<9>>>, Eigen::ConstRef<Eigen::Vector3<ADHessian<9>>>, EdgeEdgeDistanceType);
template ADGrad<12> edge_edge_distance<ADGrad<12>>(Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>>, Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>>, Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>>, Eigen::ConstRef<Eigen::Vector3<ADGrad<12>>>, EdgeEdgeDistanceType);
template ADHessian<12> edge_edge_distance<ADHessian<12>>(Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>>, Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>>, Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>>, Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>>, EdgeEdgeDistanceType);
template ADGrad<13> edge_edge_distance<ADGrad<13>>(Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>>, Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>>, Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>>, Eigen::ConstRef<Eigen::Vector3<ADGrad<13>>>, EdgeEdgeDistanceType);
template ADHessian<13> edge_edge_distance<ADHessian<13>>(Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>>, Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>>, Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>>, Eigen::ConstRef<Eigen::Vector3<ADHessian<13>>>, EdgeEdgeDistanceType);
template Vector12f edge_edge_distance_gradient<float>(Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, EdgeEdgeDistanceType);
template Vector12d edge_edge_distance_gradient<double>(Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, EdgeEdgeDistanceType);
template Matrix12f edge_edge_distance_hessian<float>(Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, EdgeEdgeDistanceType);
template Matrix12d edge_edge_distance_hessian<double>(Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, EdgeEdgeDistanceType);
// clang-format on

} // namespace ipc
