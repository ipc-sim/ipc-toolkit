#pragma once

#include <ipc/utils/eigen_ext.hpp>
#include <Eigen/Geometry>

namespace ipc {

/// @brief Compute the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param origin The origin of the plane.
/// @param normal The normal of the plane.
/// @return The distance between the point and plane.
template <typename DerivedP, typename DerivedOrigin, typename DerivedNormal>
auto point_plane_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedOrigin>& origin,
    const Eigen::MatrixBase<DerivedNormal>& normal)
{
    auto point_to_plane = (p - origin).dot(normal);
    return point_to_plane * point_to_plane / normal.squaredNorm();
}

/// @brief Compute the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param t0 The first vertex of the triangle.
/// @param t1 The second vertex of the triangle.
/// @param t2 The third vertex of the triangle.
/// @return The distance between the point and plane.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
auto point_plane_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    assert(p.size() == 3);
    assert(t0.size() == 3);
    assert(t1.size() == 3);
    assert(t2.size() == 3);

    auto normal = cross(t1 - t0, t2 - t0);
    return point_plane_distance(p, t0, normal);
}

// Symbolically generated derivatives;
namespace autogen {
    void point_plane_distance_gradient(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double v31,
        double v32,
        double v33,
        double g[12]);

    void point_plane_distance_hessian(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double v31,
        double v32,
        double v33,
        double H[144]);
} // namespace autogen

/// @brief Compute the gradient of the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param[in] p The point.
/// @param[in] origin The origin of the plane.
/// @param[in] normal The normal of the plane.
/// @param[out] grad The gradient of the distance wrt p.
template <
    typename DerivedP,
    typename DerivedOrigin,
    typename DerivedNormal,
    typename DerivedGrad>
void point_plane_distance_gradient(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedOrigin>& origin,
    const Eigen::MatrixBase<DerivedNormal>& normal,
    Eigen::PlainObjectBase<DerivedGrad>& grad)
{
    grad = (2 * (p - origin).dot(normal)) / normal.squaredNorm() * normal;
}

/// @brief Compute the gradient of the distance between a point and a plane.
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
void point_plane_distance_gradient(
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

    grad.resize(p.size() + t0.size() + t1.size() + t2.size());
    autogen::point_plane_distance_gradient(
        p[0], p[1], p[2], t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0],
        t2[1], t2[2], grad.data());
}

/// @brief Compute the hessian of the distance between a point and a plane.
/// @note The distance is actually squared distance.
/// @param[in] p The point.
/// @param[in] origin The origin of the plane.
/// @param[in] normal The normal of the plane.
/// @param[out] hess The hessian of the distance wrt p.
template <
    typename DerivedP,
    typename DerivedOrigin,
    typename DerivedNormal,
    typename DerivedHess>
void point_plane_distance_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedOrigin>& origin,
    const Eigen::MatrixBase<DerivedNormal>& normal,
    Eigen::PlainObjectBase<DerivedHess>& hess)
{
    if (normal.cols() == 1) {
        // (n×1)(n×1)ᵀ = (n×n)
        hess = 2 / normal.squaredNorm() * normal * normal.transpose();
    } else {
        assert(normal.rows() == 1);
        // (1×n)ᵀ(1×n) = (n×n)
        hess = 2 / normal.squaredNorm() * normal.transpose() * normal;
    }
}

/// @brief Compute the hessian of the distance between a point and a plane.
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
void point_plane_distance_hessian(
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

    hess.resize(
        p.size() + t0.size() + t1.size() + t2.size(),
        p.size() + t0.size() + t1.size() + t2.size());
    autogen::point_plane_distance_hessian(
        p[0], p[1], p[2], t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0],
        t2[1], t2[2], hess.data());
}

} // namespace ipc
