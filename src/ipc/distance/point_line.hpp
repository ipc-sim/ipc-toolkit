#pragma once

#include <ipc/distance/point_point.hpp>

#include <Eigen/Core>

namespace ipc {

/// @brief Compute the distance between a point and line in 2D or 3D.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge defining the line.
/// @param e1 The second vertex of the edge defining the line.
/// @return The distance between the point and line.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
auto point_line_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    assert(p.size() == 2 || p.size() == 3);
    assert(e0.size() == 2 || e0.size() == 3);
    assert(e1.size() == 2 || e1.size() == 3);

    if (p.size() == 2) {
        auto e = e1 - e0;
        auto numerator =
            (e[1] * p[0] - e[0] * p[1] + e1[0] * e0[1] - e1[1] * e0[0]);
        return numerator * numerator / e.squaredNorm();
    } else {
        return cross(e0 - p, e1 - p).squaredNorm() / (e1 - e0).squaredNorm();
    }
}

// Symbolically generated derivatives;
namespace autogen {
    void point_line_distance_gradient_2D(
        double v01,
        double v02,
        double v11,
        double v12,
        double v21,
        double v22,
        double g[6]);

    void point_line_distance_gradient_3D(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double g[9]);

    void point_line_distance_hessian_2D(
        double v01,
        double v02,
        double v11,
        double v12,
        double v21,
        double v22,
        double H[36]);

    void point_line_distance_hessian_3D(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double H[81]);
} // namespace autogen

/// @brief Compute the gradient of the distance between a point and line.
/// @note The distance is actually squared distance.
/// @param[in] p The point.
/// @param[in] e0 The first vertex of the edge defining the line.
/// @param[in] e1 The second vertex of the edge defining the line.
/// @param[out] grad The gradient of the distance wrt p, e0, and e1.
template <
    typename DerivedP,
    typename DerivedE0,
    typename DerivedE1,
    typename DerivedGrad>
void point_line_distance_gradient(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1,
    Eigen::PlainObjectBase<DerivedGrad>& grad)
{
    assert(p.size() == 2 || p.size() == 3);
    assert(e0.size() == 2 || e0.size() == 3);
    assert(e1.size() == 2 || e1.size() == 3);

    grad.resize(p.size() + e0.size() + e1.size());
    if (p.size() == 2) {
        autogen::point_line_distance_gradient_2D(
            p[0], p[1], e0[0], e0[1], e1[0], e1[1], grad.data());
    } else {
        autogen::point_line_distance_gradient_3D(
            p[0], p[1], p[2], e0[0], e0[1], e0[2], e1[0], e1[1], e1[2],
            grad.data());
    }
}

/// @brief Compute the hessian of the distance between a point and line.
/// @note The distance is actually squared distance.
/// @param[in] p The point.
/// @param[in] e0 The first vertex of the edge defining the line.
/// @param[in] e1 The second vertex of the edge defining the line.
/// @param[out] hess The hessian of the distance wrt p, e0, and e1.
template <
    typename DerivedP,
    typename DerivedE0,
    typename DerivedE1,
    typename DerivedHess>
void point_line_distance_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1,
    Eigen::PlainObjectBase<DerivedHess>& hess)
{
    assert(p.size() == 2 || p.size() == 3);
    assert(e0.size() == 2 || e0.size() == 3);
    assert(e1.size() == 2 || e1.size() == 3);

    hess.resize(
        p.size() + e0.size() + e1.size(), p.size() + e0.size() + e1.size());
    if (p.size() == 2) {
        autogen::point_line_distance_hessian_2D(
            p[0], p[1], e0[0], e0[1], e1[0], e1[1], hess.data());
    } else {
        autogen::point_line_distance_hessian_3D(
            p[0], p[1], p[2], e0[0], e0[1], e0[2], e1[0], e1[1], e1[2],
            hess.data());
    }
}

} // namespace ipc
