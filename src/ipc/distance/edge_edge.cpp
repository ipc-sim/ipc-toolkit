

#include "edge_edge.hpp"

#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/line_line.hpp>

#include <stdexcept> // std::invalid_argument

namespace ipc {

double edge_edge_distance(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1,
    EdgeEdgeDistanceType dtype)
{
    if (dtype == EdgeEdgeDistanceType::AUTO) {
        dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
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

Vector12d edge_edge_distance_gradient(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1,
    EdgeEdgeDistanceType dtype)
{
    if (dtype == EdgeEdgeDistanceType::AUTO) {
        dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
    }

    Vector12d grad = Vector12d::Zero();

    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0: {
        const Vector6d local_grad = point_point_distance_gradient(ea0, eb0);
        grad.head<3>() = local_grad.head<3>();
        grad.segment<3>(6) = local_grad.tail<3>();
        break;
    }

    case EdgeEdgeDistanceType::EA0_EB1: {
        const Vector6d local_grad = point_point_distance_gradient(ea0, eb1);
        grad.head<3>() = local_grad.head<3>();
        grad.tail<3>() = local_grad.tail<3>();
        break;
    }

    case EdgeEdgeDistanceType::EA1_EB0:
        grad.segment<6>(3) = point_point_distance_gradient(ea1, eb0);
        break;

    case EdgeEdgeDistanceType::EA1_EB1: {
        const Vector6d local_grad = point_point_distance_gradient(ea1, eb1);
        grad.segment<3>(3) = local_grad.head<3>();
        grad.tail<3>() = local_grad.tail<3>();
        break;
    }

    case EdgeEdgeDistanceType::EA_EB0: {
        const Vector9d local_grad = point_line_distance_gradient(eb0, ea0, ea1);
        grad.head<6>() = local_grad.tail<6>();
        grad.segment<3>(6) = local_grad.head<3>();
        break;
    }

    case EdgeEdgeDistanceType::EA_EB1: {
        const Vector9d local_grad = point_line_distance_gradient(eb1, ea0, ea1);
        grad.head<6>() = local_grad.tail<6>();
        grad.tail<3>() = local_grad.head<3>();
        break;
    }

    case EdgeEdgeDistanceType::EA0_EB: {
        const Vector9d local_grad = point_line_distance_gradient(ea0, eb0, eb1);
        grad.head<3>() = local_grad.head<3>();
        grad.tail<6>() = local_grad.tail<6>();
        break;
    }

    case EdgeEdgeDistanceType::EA1_EB:
        grad.tail<9>() = point_line_distance_gradient(ea1, eb0, eb1);
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

Matrix12d edge_edge_distance_hessian(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1,
    EdgeEdgeDistanceType dtype)
{
    if (dtype == EdgeEdgeDistanceType::AUTO) {
        dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
    }

    Matrix12d hess = Matrix12d::Zero();

    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0: {
        const Matrix6d local_hess = point_point_distance_hessian(ea0, eb0);
        hess.topLeftCorner<3, 3>() = local_hess.topLeftCorner<3, 3>();
        hess.block<3, 3>(0, 6) = local_hess.topRightCorner<3, 3>();
        hess.block<3, 3>(6, 0) = local_hess.bottomLeftCorner<3, 3>();
        hess.block<3, 3>(6, 6) = local_hess.bottomRightCorner<3, 3>();
        break;
    }

    case EdgeEdgeDistanceType::EA0_EB1: {
        const Matrix6d local_hess = point_point_distance_hessian(ea0, eb1);
        hess.topLeftCorner<3, 3>() = local_hess.topLeftCorner<3, 3>();
        hess.topRightCorner<3, 3>() = local_hess.topRightCorner<3, 3>();
        hess.bottomLeftCorner<3, 3>() = local_hess.bottomLeftCorner<3, 3>();
        hess.bottomRightCorner<3, 3>() = local_hess.bottomRightCorner<3, 3>();
        break;
    }

    case EdgeEdgeDistanceType::EA1_EB0:
        hess.block<6, 6>(3, 3) = point_point_distance_hessian(ea1, eb0);
        break;

    case EdgeEdgeDistanceType::EA1_EB1: {
        const Matrix6d local_hess = point_point_distance_hessian(ea1, eb1);
        hess.block<3, 3>(3, 3) = local_hess.topLeftCorner<3, 3>();
        hess.block<3, 3>(3, 9) = local_hess.topRightCorner<3, 3>();
        hess.block<3, 3>(9, 3) = local_hess.bottomLeftCorner<3, 3>();
        hess.bottomRightCorner<3, 3>() = local_hess.bottomRightCorner<3, 3>();
        break;
    }

    case EdgeEdgeDistanceType::EA_EB0: {
        const Matrix9d local_hess = point_line_distance_hessian(eb0, ea0, ea1);
        hess.topLeftCorner<6, 6>() = local_hess.bottomRightCorner<6, 6>();
        hess.block<3, 6>(6, 0) = local_hess.topRightCorner<3, 6>();
        hess.block<6, 3>(0, 6) = local_hess.bottomLeftCorner<6, 3>();
        hess.block<3, 3>(6, 6) = local_hess.topLeftCorner<3, 3>();
        break;
    }

    case EdgeEdgeDistanceType::EA_EB1: {
        const Matrix9d local_hess = point_line_distance_hessian(eb1, ea0, ea1);
        hess.topLeftCorner<6, 6>() = local_hess.bottomRightCorner<6, 6>();
        hess.topRightCorner<6, 3>() = local_hess.bottomLeftCorner<6, 3>();
        hess.bottomLeftCorner<3, 6>() = local_hess.topRightCorner<3, 6>();
        hess.bottomRightCorner<3, 3>() = local_hess.topLeftCorner<3, 3>();
        break;
    }

    case EdgeEdgeDistanceType::EA0_EB: {
        const Matrix9d local_hess = point_line_distance_hessian(ea0, eb0, eb1);
        hess.topLeftCorner<3, 3>() = local_hess.topLeftCorner<3, 3>();
        hess.topRightCorner<3, 6>() = local_hess.topRightCorner<3, 6>();
        hess.bottomLeftCorner<6, 3>() = local_hess.bottomLeftCorner<6, 3>();
        hess.bottomRightCorner<6, 6>() = local_hess.bottomRightCorner<6, 6>();
        break;
    }

    case EdgeEdgeDistanceType::EA1_EB:
        hess.bottomRightCorner<9, 9>() =
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

} // namespace ipc
