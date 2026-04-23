#include "point_edge.hpp"

#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_point.hpp>

#include <stdexcept> // std::invalid_argument

namespace ipc {

template <typename T>
T point_edge_distance(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1,
    PointEdgeDistanceType dtype)
{
    assert(p.size() == 2 || p.size() == 3);
    assert(e0.size() == 2 || e0.size() == 3);
    assert(e1.size() == 2 || e1.size() == 3);

    if constexpr (std::is_same_v<T, double> || std::is_same_v<T, float>) {
        if (dtype == PointEdgeDistanceType::AUTO) {
            dtype = point_edge_distance_type(p, e0, e1);
        }
    }

    switch (dtype) {
    case PointEdgeDistanceType::P_E0:
        return point_point_distance(p, e0);

    case PointEdgeDistanceType::P_E1:
        return point_point_distance(p, e1);

    case PointEdgeDistanceType::P_E:
        return point_line_distance(p, e0, e1);

    default:
        throw std::invalid_argument(
            "Invalid distance type for point-edge distance!");
    }
}

template <typename T>
VectorMax9<T> point_edge_distance_gradient(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1,
    PointEdgeDistanceType dtype)
{
    const int dim = p.size();
    assert(e0.size() == dim);
    assert(e1.size() == dim);

    if (dtype == PointEdgeDistanceType::AUTO) {
        dtype = point_edge_distance_type(p, e0, e1);
    }

    VectorMax9<T> grad = VectorMax9<T>::Zero(3 * dim);

    switch (dtype) {
    case PointEdgeDistanceType::P_E0:
        grad.head(2 * dim) = point_point_distance_gradient(p, e0);
        break;

    case PointEdgeDistanceType::P_E1: {
        const VectorMax6<T> local_grad = point_point_distance_gradient(p, e1);
        grad.head(dim) = local_grad.head(dim);
        grad.tail(dim) = local_grad.tail(dim);
        break;
    }

    case PointEdgeDistanceType::P_E:
        grad = point_line_distance_gradient(p, e0, e1);
        break;

    default:
        throw std::invalid_argument(
            "Invalid distance type for point-edge distance gradient!");
    }

    return grad;
}

template <typename T>
MatrixMax9<T> point_edge_distance_hessian(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1,
    PointEdgeDistanceType dtype)
{
    const int dim = p.size();
    assert(e0.size() == dim);
    assert(e1.size() == dim);

    if (dtype == PointEdgeDistanceType::AUTO) {
        dtype = point_edge_distance_type(p, e0, e1);
    }

    MatrixMax9<T> hess = MatrixMax9<T>::Zero(3 * dim, 3 * dim);

    switch (dtype) {
    case PointEdgeDistanceType::P_E0:
        hess.topLeftCorner(2 * dim, 2 * dim) =
            point_point_distance_hessian(p, e0);
        break;

    case PointEdgeDistanceType::P_E1: {
        const MatrixMax6<T> local_hess = point_point_distance_hessian(p, e1);
        hess.topLeftCorner(dim, dim) = local_hess.topLeftCorner(dim, dim);
        hess.topRightCorner(dim, dim) = local_hess.topRightCorner(dim, dim);
        hess.bottomLeftCorner(dim, dim) = local_hess.bottomLeftCorner(dim, dim);
        hess.bottomRightCorner(dim, dim) =
            local_hess.bottomRightCorner(dim, dim);
        break;
    }

    case PointEdgeDistanceType::P_E:
        hess = point_line_distance_hessian(p, e0, e1);
        break;

    default:
        throw std::invalid_argument(
            "Invalid distance type for point-edge distance hessian!");
    }

    return hess;
}

// clang-format off
template float point_edge_distance<float>(Eigen::ConstRef<VectorMax3f>, Eigen::ConstRef<VectorMax3f>, Eigen::ConstRef<VectorMax3f>, PointEdgeDistanceType);
template double point_edge_distance<double>(Eigen::ConstRef<VectorMax3d>, Eigen::ConstRef<VectorMax3d>, Eigen::ConstRef<VectorMax3d>, PointEdgeDistanceType);
template VectorMax9<float> point_edge_distance_gradient<float>(Eigen::ConstRef<VectorMax3f>, Eigen::ConstRef<VectorMax3f>, Eigen::ConstRef<VectorMax3f>, PointEdgeDistanceType);
template VectorMax9<double> point_edge_distance_gradient<double>(Eigen::ConstRef<VectorMax3d>, Eigen::ConstRef<VectorMax3d>, Eigen::ConstRef<VectorMax3d>, PointEdgeDistanceType);
template MatrixMax9<float> point_edge_distance_hessian<float>(Eigen::ConstRef<VectorMax3f>, Eigen::ConstRef<VectorMax3f>, Eigen::ConstRef<VectorMax3f>, PointEdgeDistanceType);
template MatrixMax9<double> point_edge_distance_hessian<double>(Eigen::ConstRef<VectorMax3d>, Eigen::ConstRef<VectorMax3d>, Eigen::ConstRef<VectorMax3d>, PointEdgeDistanceType);
// clang-format on

} // namespace ipc
