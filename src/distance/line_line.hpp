#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <utils/eigen_ext.hpp>

namespace ipc {

/// @brief Compute the distance between a two infinite lines in 3D.
/// @note The distance is actually squared distance.
/// @warning If the lines are parallel this function returns a distance of zero.
/// @param ea0,ea1 The points of the edge defining the first line.
/// @param eb0,eb1 The points of the edge defining the second line.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
auto line_line_distance(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    const auto normal = (ea1 - ea0).cross(eb1 - eb0);
    const auto line_to_line = (eb0 - ea0).dot(normal);
    return line_to_line * line_to_line / normal.squaredNorm();
}

// Symbolically generated derivatives;
namespace autogen {
    void line_line_distance_gradient(
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

    void line_line_distance_hessian(
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

/// @brief Compute the gradient of the distance between a two lines in 3D.
/// @note The distance is actually squared distance.
/// @warning If the lines are parallel this function returns a distance of zero.
/// @param ea0,ea1 The points of the edge defining the first line.
/// @param eb0,eb1 The points of the edge defining the second line.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    typename DerivedGrad>
void line_line_distance_gradient(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    Eigen::PlainObjectBase<DerivedGrad>& grad)
{
    grad.resize(ea0.size() + ea1.size() + eb0.size() + eb1.size());
    autogen::line_line_distance_gradient(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], grad.data());
}

/// @brief Compute the hessian of the distance between a two lines in 3D.
/// @note The distance is actually squared distance.
/// @warning If the lines are parallel this function returns a distance of zero.
/// @param[in] ea0,ea1 The points of the edge defining the first line.
/// @param[in] eb0,eb1 The points of the edge defining the second line.
/// @param[out] hess The computed hessian.
/// @param[in] project_to_psd True if the hessian should be projected to
///                           positive semi-definite.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    typename DerivedHess>
void line_line_distance_hessian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    Eigen::PlainObjectBase<DerivedHess>& hess,
    bool project_to_psd = false)
{
    hess.resize(
        ea0.size() + ea1.size() + eb0.size() + eb1.size(),
        ea0.size() + ea1.size() + eb0.size() + eb1.size());
    autogen::line_line_distance_hessian(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], hess.data());
    if (project_to_psd) {
        Eigen::project_to_psd(hess);
    }
}

} // namespace ipc
