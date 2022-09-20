#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_plane.hpp>
#include <ipc/distance/point_point.hpp>

#include <stdexcept> // std::invalid_argument

namespace ipc {

/// @brief Compute the distance between a points and a triangle.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The distance between the point and triangle.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
auto point_triangle_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    assert(p.size() == 3);
    assert(t0.size() == 3);
    assert(t1.size() == 3);
    assert(t2.size() == 3);

    return point_triangle_distance(
        p, t0, t1, t2, point_triangle_distance_type(p, t0, t1, t2));
}

/// @brief Compute the distance between a points and a triangle.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @param dtype The point-triangle distance type to compute.
/// @return The distance between the point and triangle.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
auto point_triangle_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2,
    const PointTriangleDistanceType dtype)
{
    assert(p.size() == 3);
    assert(t0.size() == 3);
    assert(t1.size() == 3);
    assert(t2.size() == 3);

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

/// @brief Compute the gradient of the distance between a points and a triangle.
/// @note The distance is actually squared distance.
/// @param[in] p The point.
/// @param[in] t0 The first vertex of the triangle.
/// @param[in] t1 The second vertex of the triangle.
/// @param[in] t2 The third vertex of the triangle.
/// @param[out] grad The gradient of the distance wrt p, t0, t1, and t2.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2,
    typename DerivedGrad>
void point_triangle_distance_gradient(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2,
    Eigen::PlainObjectBase<DerivedGrad>& grad)
{
    assert(p.size() == 3);
    assert(t0.size() == 3);
    assert(t1.size() == 3);
    assert(t2.size() == 3);

    return point_triangle_distance_gradient(
        p, t0, t1, t2, point_triangle_distance_type(p, t0, t1, t2), grad);
}

/// @brief Compute the gradient of the distance between a points and a triangle.
/// @note The distance is actually squared distance.
/// @param[in] p The point.
/// @param[in] t0 The first vertex of the triangle.
/// @param[in] t1 The second vertex of the triangle.
/// @param[in] t2 The third vertex of the triangle.
/// @param dtype The point-triangle distance type to compute.
/// @param[out] grad The gradient of the distance wrt p, t0, t1, and t2.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2,
    typename DerivedGrad>
void point_triangle_distance_gradient(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2,
    const PointTriangleDistanceType dtype,
    Eigen::PlainObjectBase<DerivedGrad>& grad)
{
    int dim = p.size();
    assert(t0.size() == dim);
    assert(t1.size() == dim);
    assert(t2.size() == dim);

    grad.resize(4 * dim);
    grad.setZero();

    VectorMax9<typename DerivedGrad::Scalar> local_grad;
    switch (dtype) {
    case PointTriangleDistanceType::P_T0:
        point_point_distance_gradient(p, t0, local_grad);
        grad.head(2 * dim) = local_grad;
        break;

    case PointTriangleDistanceType::P_T1:
        point_point_distance_gradient(p, t1, local_grad);
        grad.head(dim) = local_grad.head(dim);
        grad.segment(2 * dim, dim) = local_grad.tail(dim);
        break;

    case PointTriangleDistanceType::P_T2:
        point_point_distance_gradient(p, t2, local_grad);
        grad.head(dim) = local_grad.head(dim);
        grad.tail(dim) = local_grad.tail(dim);
        break;

    case PointTriangleDistanceType::P_E0:
        point_line_distance_gradient(p, t0, t1, local_grad);
        grad.head(3 * dim) = local_grad;
        break;

    case PointTriangleDistanceType::P_E1:
        point_line_distance_gradient(p, t1, t2, local_grad);
        grad.head(dim) = local_grad.head(dim);
        grad.tail(2 * dim) = local_grad.tail(2 * dim);
        break;

    case PointTriangleDistanceType::P_E2:
        point_line_distance_gradient(p, t2, t0, local_grad);
        grad.head(dim) = local_grad.head(dim);         // ∇_p
        grad.segment(dim, dim) = local_grad.tail(dim); // ∇_{t0}
        grad.tail(dim) = local_grad.segment(dim, dim); // ∇_{t2}
        break;

    case PointTriangleDistanceType::P_T:
        point_plane_distance_gradient(p, t0, t1, t2, grad);
        break;

    default:
        throw std::invalid_argument(
            "Invalid distance type for point-triangle distance gradient!");
    }
}

/// @brief Compute the hessian of the distance between a points and a triangle.
/// @note The distance is actually squared distance.
/// @param[in] p The point.
/// @param[in] t0 The first vertex of the triangle.
/// @param[in] t1 The second vertex of the triangle.
/// @param[in] t2 The third vertex of the triangle.
/// @param[out] hess The hessian of the distance wrt p, t0, t1, and t2.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2,
    typename DerivedHess>
void point_triangle_distance_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2,
    Eigen::PlainObjectBase<DerivedHess>& hess)
{
    assert(p.size() == 3);
    assert(t0.size() == 3);
    assert(t1.size() == 3);
    assert(t2.size() == 3);

    return point_triangle_distance_hessian(
        p, t0, t1, t2, point_triangle_distance_type(p, t0, t1, t2), hess);
}

/// @brief Compute the hessian of the distance between a points and a triangle.
/// @note The distance is actually squared distance.
/// @param[in] p The point.
/// @param[in] t0 The first vertex of the triangle.
/// @param[in] t1 The second vertex of the triangle.
/// @param[in] t2 The third vertex of the triangle.
/// @param dtype The point-triangle distance type to compute.
/// @param[out] hess The hessian of the distance wrt p, t0, t1, and t2.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2,
    typename DerivedHess>
void point_triangle_distance_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2,
    const PointTriangleDistanceType dtype,
    Eigen::PlainObjectBase<DerivedHess>& hess)
{
    int dim = p.size();
    assert(t0.size() == dim);
    assert(t1.size() == dim);
    assert(t2.size() == dim);

    hess.resize(4 * dim, 4 * dim);
    hess.setZero();

    MatrixMax9<typename DerivedHess::Scalar> local_hess;
    switch (dtype) {
    case PointTriangleDistanceType::P_T0:
        point_point_distance_hessian(p, t0, local_hess);
        hess.topLeftCorner(2 * dim, 2 * dim) = local_hess;
        break;

    case PointTriangleDistanceType::P_T1:
        point_point_distance_hessian(p, t1, local_hess);
        hess.topLeftCorner(dim, dim) = local_hess.topLeftCorner(dim, dim);
        hess.block(0, 2 * dim, dim, dim) = local_hess.topRightCorner(dim, dim);
        hess.block(2 * dim, 0, dim, dim) =
            local_hess.bottomLeftCorner(dim, dim);
        hess.block(2 * dim, 2 * dim, dim, dim) =
            local_hess.bottomRightCorner(dim, dim);
        break;

    case PointTriangleDistanceType::P_T2:
        point_point_distance_hessian(p, t2, local_hess);
        hess.topLeftCorner(dim, dim) = local_hess.topLeftCorner(dim, dim);
        hess.topRightCorner(dim, dim) = local_hess.topRightCorner(dim, dim);
        hess.bottomLeftCorner(dim, dim) = local_hess.bottomLeftCorner(dim, dim);
        hess.bottomRightCorner(dim, dim) =
            local_hess.bottomRightCorner(dim, dim);
        break;

    case PointTriangleDistanceType::P_E0:
        point_line_distance_hessian(p, t0, t1, local_hess);
        hess.topLeftCorner(3 * dim, 3 * dim) = local_hess;
        break;

    case PointTriangleDistanceType::P_E1:
        point_line_distance_hessian(p, t1, t2, local_hess);
        hess.topLeftCorner(dim, dim) = local_hess.topLeftCorner(dim, dim);
        hess.topRightCorner(dim, 2 * dim) =
            local_hess.topRightCorner(dim, 2 * dim);
        hess.bottomLeftCorner(2 * dim, dim) =
            local_hess.bottomLeftCorner(2 * dim, dim);
        hess.bottomRightCorner(2 * dim, 2 * dim) =
            local_hess.bottomRightCorner(2 * dim, 2 * dim);
        break;

    case PointTriangleDistanceType::P_E2:
        point_line_distance_hessian(p, t2, t0, local_hess);
        hess.topLeftCorner(dim, dim) = local_hess.topLeftCorner(dim, dim);
        hess.block(0, dim, dim, dim) = local_hess.topRightCorner(dim, dim);
        hess.topRightCorner(dim, dim) = local_hess.block(0, dim, dim, dim);
        hess.block(dim, 0, dim, dim) = local_hess.bottomLeftCorner(dim, dim);
        hess.block(dim, dim, dim, dim) = local_hess.bottomRightCorner(dim, dim);
        hess.block(dim, 3 * dim, dim, dim) =
            local_hess.block(2 * dim, dim, dim, dim);
        hess.bottomLeftCorner(dim, dim) = local_hess.block(dim, 0, dim, dim);
        hess.block(3 * dim, dim, dim, dim) =
            local_hess.block(dim, 2 * dim, dim, dim);
        hess.bottomRightCorner(dim, dim) = local_hess.block(dim, dim, dim, dim);
        break;

    case PointTriangleDistanceType::P_T:
        point_plane_distance_hessian(p, t0, t1, t2, hess);
        break;

    default:
        throw std::invalid_argument(
            "Invalid distance type for point-triangle distance hessian!");
    }
}

} // namespace ipc
