#pragma once

#include <Eigen/Core>

#include <ipc/distance/point_point.hpp>

namespace ipc {

/// @brief Compute the distance between a point and a line (defined by an edge).
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0,e1 The points of the edge defining the line.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
auto point_line_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    if (p.size() == 2) {
        auto e = e1 - e0;
        auto numerator =
            (e[1] * p[0] - e[0] * p[1] + e1[0] * e0[1] - e1[1] * e0[0]);
        return numerator * numerator / e.squaredNorm();
    } else {
        return Eigen::cross(e0 - p, e1 - p).squaredNorm()
            / (e1 - e0).squaredNorm();
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

template <
    typename DerivedP,
    typename DerivedE0,
    typename DerivedE1,
    typename DerivedHess>
void point_line_distance_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1,
    Eigen::PlainObjectBase<DerivedHess>& hess,
    bool project_to_psd = false)
{
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
    if (project_to_psd) {
        Eigen::project_to_psd(hess);
    }
}

} // namespace ipc
